
#devtools::load_all("C:/Users/pk02/Packages/spatind") # spatind package under development

# 1. Preamble ------------------------------------------------------------------

rm(list=ls())

if (!requireNamespace("spatind", quietly = TRUE)) {
  install.packages(
    here::here("Packages", "spatind_0.0.0.9000.tar.gz"),
    repos = NULL,
    type = "source"
  )
}


if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
  install.packages("rnaturalearth")
}

library(spatind)
library(ggplot2)
library(dplyr)
library(here)

figure_dir <- here::here("Outputs", "Figures")

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Load data -----------------------------------------------------------------

scallop_comor <- readRDS(here::here("Data", "comor_2000_2025_density.RDS"))


# 3. Basic checks --------------------------------------------------------------

str(scallop_comor)

ggplot(data = scallop_comor) +
  geom_point(aes(x = lon, y = lat, size = dens_abun_m2)) +
  facet_wrap(~year) +
  theme_bw()


# 4. COMOR Survey --------------------------------------------------------------

## 4.1. Boundary ---------------------------------------------------------------

# Unique sampling locations across the time series
survey_locations <- scallop_comor %>%
  dplyr::distinct(lon, lat) %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  sf::st_transform(3035)

# Concave stock boundary with a 1 km buffer
hull_concave <- survey_locations %>%
  sf::st_geometry() %>%
  sf::st_union() %>%
  sf::st_concave_hull(ratio = 0.3, allow_holes = FALSE) %>%
  sf::st_buffer(1000) %>%
  sf::st_as_sf()

## 4.2. Survey area ------------------------------------------------------------

# French coastline
france <- rnaturalearth::ne_countries(country = "France", scale = "medium", returnclass = "sf") %>%
  sf::st_transform(3035)

# Add 5 km around the buffered boundary for plotting
map_bbox <- hull_concave %>%
  sf::st_buffer(5000) %>%
  sf::st_bbox()

p_survey_area <- ggplot() +
  geom_sf(data = france, fill = "grey90", colour = "grey50", linewidth = 0.3) +
  geom_sf(data = hull_concave, fill = "#56B4E9", colour = "#0072B2", linewidth = 0.7, alpha = 0.25) +
  geom_sf(data = survey_locations, colour = "black", size = 0.7, alpha = 0.45) +
  coord_sf(xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
           ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
           crs = sf::st_crs(3035), datum = sf::st_crs(4326), expand = FALSE) +
  labs(title = "COMOR survey area in the Bay of Seine",
       subtitle = "Sampling locations used during 2000–2025", x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.25),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30"))

p_survey_area

ggsave(filename = file.path(figure_dir, "scallop_COMOR_survey_area.png"),
       plot = p_survey_area, width = 6, height = 6, units = "in", dpi = 400, bg = "white")


# 5. Voronoi weight function ---------------------------------------------------

# debug run through
# df = scallop_comor
# poly_bound = hull_buffed
# crs = 3035

# performs per timestep
get_Voronoi <- function(
    lon,
    lat,
    z,
    poly_bound,
    crs_in = 4326,
    crs_out = 3035
) {
  
  if (!all(c(length(lon), length(lat), length(z)) == length(lon))) {
    stop("lon, lat and z must have equal lengths.")
  }
  
  if (length(lon) < 3) {
    stop("At least three observations are required.")
  }
  
  if (anyNA(lon) || anyNA(lat)) {
    stop("lon and lat cannot contain missing values.")
  }
  
  if (is.na(sf::st_crs(poly_bound))) {
    stop("poly_bound must have a defined CRS.")
  }
  
  # Combine observations at identical locations
  df <- tibble::tibble(
    lon = lon,
    lat = lat,
    z = z
  ) %>%
    dplyr::group_by(lon, lat) %>%
    dplyr::summarise(
      z = if (all(is.na(z))) {
        NA_real_
      } else {
        mean(z, na.rm = TRUE)
      },
      n_observations = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      loc_id = dplyr::row_number()
    )
  
  if (nrow(df) < 3) {
    stop(
      "At least three unique sampling locations are required."
    )
  }
  
  # Create projected sampling points
  loc_sf <- df %>%
    sf::st_as_sf(
      coords = c("lon", "lat"),
      crs = crs_in,
      remove = FALSE
    ) %>%
    sf::st_transform(crs_out)
  
  coordinates <- sf::st_coordinates(loc_sf)
  
  df <- df %>%
    dplyr::mutate(
      x = coordinates[, "X"],
      y = coordinates[, "Y"]
    )
  
  # Prepare clipping boundary
  poly_geom <- poly_bound %>%
    sf::st_transform(crs_out) %>%
    sf::st_make_valid() %>%
    sf::st_geometry() %>%
    sf::st_union()
  
  # Construct initial rectangular Voronoi tessellation
  voronoi_geom <- loc_sf %>%
    sf::st_geometry() %>%
    sf::st_union() %>%
    sf::st_voronoi(
      envelope = sf::st_as_sfc(
        sf::st_bbox(poly_geom)
      )
    ) %>%
    sf::st_collection_extract("POLYGON")
  
  voronoi_sf <- sf::st_sf(
    geometry = voronoi_geom
  )
  
  # Match each polygon to its generating sampling location
  voronoi_sf <- sf::st_join(
    voronoi_sf,
    loc_sf["loc_id"],
    join = sf::st_contains
  )
  
  if (anyNA(voronoi_sf$loc_id)) {
    stop(
      "One or more Voronoi polygons could not be matched ",
      "to a sampling location."
    )
  }
  
  # Clip the rectangular tessellation to poly_bound
  voronoi_sf <- suppressWarnings(
    sf::st_intersection(
      voronoi_sf,
      poly_geom
    )
  )
  
  # Combine fragments belonging to the same location
  voronoi_sf <- voronoi_sf %>%
    dplyr::group_by(loc_id) %>%
    dplyr::summarise(
      geometry = sf::st_union(geometry),
      .groups = "drop"
    )
  
  # Calculate clipped areas and weights
  voronoi_sf$area_km2 <- as.numeric(
    sf::st_area(voronoi_sf)
  ) / 1e6
  
  voronoi_sf <- voronoi_sf %>%
    dplyr::mutate(
      weight = area_km2 / sum(area_km2)
    ) %>%
    dplyr::left_join(
      df,
      by = "loc_id"
    )
  
  return(voronoi_sf)
}


# 6. Get Voronoi weights -------------------------------------------------------

scallop_by_year <- split(scallop_comor, scallop_comor$year)

scallop_voronoi_list <- lapply(scallop_by_year, function(dat) {
  
  polygons <- get_Voronoi(lon = dat$lon, lat = dat$lat, z = dat$dens_abun_m2,
                          poly_bound = hull_concave)
  
  polygons$year <- dat$year[1]
  polygons
})

scallop_voronoi <- do.call(rbind, scallop_voronoi_list)


# 7. Check the weights ----------------------------------------------------------

# Shared styling and reference values
theme_checks <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
        panel.border = element_rect(colour = "grey65", linewidth = 0.4),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey35", margin = margin(b = 12)),
        plot.caption = element_text(colour = "grey40", hjust = 0, margin = margin(t = 10)),
        plot.title.position = "plot", plot.caption.position = "plot")

scallop_voronoi <- scallop_voronoi %>%
  dplyr::mutate(year = factor(year, levels = sort(unique(year))))

weight_mean <- mean(scallop_voronoi$weight) * 100
weight_median <- median(scallop_voronoi$weight) * 100
weight_caption <- "Dashed blue line: pooled mean. Dotted orange line: pooled median."

# 1. Voronoi maps
p_weights_map <- ggplot(scallop_voronoi) +
  geom_sf(aes(fill = weight), colour = "grey50", linewidth = 0.1) +
  facet_wrap(~ year, ncol = 5) +
  viridis::scale_fill_viridis(name = "Voronoi weight", labels = scales::label_percent()) +
  coord_sf(datum = NA) +
  labs(title = "Scallop survey Voronoi weights",
       subtitle = "COMOR survey, Bay of Seine | 2000–2025",
       caption = "Each weight is the proportion of the analysis boundary represented by a sampling location.") +
  theme_void(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9, margin = margin(b = 5)),
        panel.spacing = grid::unit(0.8, "lines"),
        legend.position = "bottom",
        legend.key.width = grid::unit(2, "cm"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey35", margin = margin(b = 12)),
        plot.caption = element_text(colour = "grey40", hjust = 0, margin = margin(t = 10)),
        plot.margin = margin(12, 12, 12, 12))

p_weights_map

ggsave(filename = file.path(figure_dir, "scallop_COMOR_Voronoi_weights_map.png"),
       plot = p_weights_map, width = 9, height = 8, units = "in", dpi = 400, bg = "white")

# 2. Weight distributions by year
p_weights_box <- ggplot(scallop_voronoi, aes(x = year, y = weight * 100)) +
  geom_boxplot(fill = "#56B4E9", alpha = 0.6, colour = "#253746", linewidth = 0.4,
               width = 0.65, outlier.size = 1, outlier.alpha = 0.5) +
  geom_hline(yintercept = weight_mean, linetype = "dashed", colour = "#0072B2", linewidth = 0.6) +
  geom_hline(yintercept = weight_median, linetype = "dotted", colour = "#D55E00", linewidth = 0.6) +
  coord_flip() +
  labs(title = "Scallop survey Voronoi weight distributions",
       subtitle = "COMOR survey, Bay of Seine | 2000–2025",
       x = "Year", y = "Voronoi weight (%)", caption = weight_caption) +
  theme_checks +
  theme(panel.grid.major.y = element_blank())

p_weights_box

ggsave(filename = file.path(figure_dir, "scallop_COMOR_Voronoi_weights_boxplots.png"),
       plot = p_weights_box, width = 7, height = 8, units = "in", dpi = 400, bg = "white")

# 3. Weight density curves
p_weights_density <- ggplot(scallop_voronoi, aes(x = weight * 100, colour = year)) +
  geom_density(linewidth = 0.6, alpha = 0.8, bounds = c(0, Inf), key_glyph = "path") +
  scale_colour_viridis_d(name = "Year") +
  geom_vline(xintercept = weight_mean, linetype = "dashed", colour = "#0072B2", linewidth = 0.6) +
  geom_vline(xintercept = weight_median, linetype = "dotted", colour = "#D55E00", linewidth = 0.6) +
  labs(title = "Scallop survey Voronoi weight densities",
       subtitle = "COMOR survey, Bay of Seine | 2000–2025",
       x = "Voronoi weight (%)", y = "Density", caption = weight_caption) +
  theme_checks +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE))

p_weights_density

ggsave(filename = file.path(figure_dir, "scallop_COMOR_Voronoi_weights_density.png"),
       plot = p_weights_density, width = 8, height = 6, units = "in", dpi = 400, bg = "white")


# 8. Calculate spatial indicators ----------------------------------------------

scallop_spatinds <- scallop_voronoi %>%
  sf::st_drop_geometry() %>%
  group_by(year) %>%
  summarise(sa_km2  = spatind::spreadingarea(z = z, w = area_km2),
            pa_km2  = spatind::positivearea(z = z, w = area_km2),
            ea_km2  = spatind::equivalentarea(z = z, w = area_km2),
            inert   = spatind::inertia(x = x, y = y, z = z*area_km2),
            mean_density = stats::weighted.mean(z, w = area_km2),
            .groups = "drop") %>%
  tidyr::unpack(inert)


# 9. Plot spatial indicators ---------------------------------------------------

ind_order <- c("cog_x", "cog_y", "inert", "isotropy",
               "pa_km2", "sa_km2", "ea_km2", "mean_density")

ind_labels <- c(mean_density = "Mean density (individuals m⁻²)",
                pa_km2 = "Positive area (km²)",
                sa_km2 = "Spreading area (km²)",
                ea_km2 = "Equivalent area (km²)",
                cog_x = "East–west centre of gravity (km)",
                cog_y = "North–south centre of gravity (km)",
                inert = "Spatial inertia (km²)",
                isotropy = "Isotropy")

scallop_spatinds_l <- scallop_spatinds %>%
  tidyr::pivot_longer(cols = -year, names_to = "Indicator", values_to = "Value") %>%
  dplyr::mutate(year = as.numeric(as.character(year))) %>%
  tidyr::complete(Indicator, year = 2000:2025) %>%
  dplyr::mutate(Indicator = factor(Indicator, levels = ind_order)) %>%
  dplyr::group_by(Indicator) %>%
  dplyr::mutate(
    deviation = Value - mean(Value, na.rm = TRUE),
    # Convert CoG from metres to kilometres and inertia from m² to km²
    deviation = dplyr::if_else(Indicator %in% c("cog_x", "cog_y"), deviation / 1000, deviation),
    deviation = dplyr::if_else(Indicator == "inert", deviation / 1e6, deviation)
  ) %>%
  dplyr::ungroup()

p_spatinds <- ggplot(scallop_spatinds_l, aes(x = year, y = deviation)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
  geom_line(colour = "#253746", linewidth = 0.5) +
  geom_point(colour = "#253746", size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "#0072B2", linewidth = 0.7) +
  facet_wrap(~ Indicator, ncol = 2, scales = "free_y",
             labeller = ggplot2::as_labeller(ind_labels)) +
  scale_x_continuous(breaks = seq(2000, 2025, 5), minor_breaks = NULL) +
  labs(title = "King scallop spatial indicators",
       subtitle = expression(italic("Pecten maximus") * " | COMOR survey, Bay of Seine | 2000–2025"),
       x = "Year", y = "Deviation from long-term mean",
       caption = "Blue lines show linear trends. Dashed lines indicate the long-term mean.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
        panel.border = element_rect(colour = "grey65", linewidth = 0.4),
        strip.background = element_rect(fill = "grey95", colour = "grey65", linewidth = 0.4),
        strip.text = element_text(face = "bold", size = 10, margin = margin(6, 4, 6, 4)),
        panel.spacing = grid::unit(1, "lines"),
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(colour = "grey35", margin = margin(b = 12)),
        plot.caption = element_text(colour = "grey40", hjust = 0, margin = margin(t = 10)),
        plot.title.position = "plot",
        plot.caption.position = "plot")

p_spatinds

ggsave(filename = file.path(figure_dir, "scallop_COMOR_spatial_indicators.png"),
       plot = p_spatinds, width = 8, height = 8, units = "in", dpi = 400, bg = "white")


# 10. Checks --------------------------------------------------------------------

scallop_checks <- scallop_voronoi %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(n_locations = dplyr::n(),
                   n_missing = sum(is.na(z)),
                   n_zero = sum(z == 0, na.rm = TRUE),
                   proportion_zero = mean(z == 0, na.rm = TRUE),
                   total_area_km2 = sum(area_km2),
                   maximum_weight = max(weight), .groups = "drop") %>%
  dplyr::mutate(year = as.numeric(as.character(year))) %>%
  tidyr::complete(year = 2000:2025)

print(scallop_checks, n = nrow(scallop_checks))

dev.new()
ggplot(scallop_checks, aes(x = year, y = maximum_weight)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "Year", y = "Largest Voronoi weight") +
  theme_bw()

scallop_zero_checks <- scallop_voronoi %>%
  sf::st_drop_geometry() %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(minimum_density = min(z),
                   lower_5_percent = quantile(z, 0.05),
                   zero_area_km2 = sum(area_km2[z == 0]),
                   zero_area_percent = 100 * sum(weight[z == 0]), .groups = "drop")

print(scallop_zero_checks, n = nrow(scallop_zero_checks))
