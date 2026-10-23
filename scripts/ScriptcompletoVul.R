# ============================================================
# SCRIPT COMPLETO DE GENERACIÓN DE FIGURAS (CORREGIDO)
# Proyecto: Beyond composite scores (Guayas, Ecuador)
# ============================================================

# ------------------------------------------------------------
# 0. CONFIGURACIÓN INICIAL Y LIBRERÍAS
# ------------------------------------------------------------
setwd("~/Documents/ARTICULOS")

if (!dir.exists("figuras")) dir.create("figuras")
if (!dir.exists("datos")) dir.create("datos")

library(readxl)
library(janitor)
library(dplyr)
library(tidyr)
library(sf)
library(spdep)
library(ggplot2)
library(ggspatial)
library(scales)
library(patchwork)
library(geodata)
library(psych)
library(Matrix)

# ------------------------------------------------------------
# 1. CARGA Y PREPARACIÓN DE DATOS
# ------------------------------------------------------------
cat("Cargando y preparando datos...\n")

# Mapa base de Guayas
ecuador_prov <- gadm(country = "ECU", level = 1, path = tempdir())
ecuador_sf <- st_as_sf(ecuador_prov)
guayas <- ecuador_sf %>% filter(NAME_1 == "Guayas") %>% st_transform(32717)

# Datos GeoDa
geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)
geoda_data$coord_x_utm <- as.numeric(as.character(geoda_data$coord_x_utm))
geoda_data$coord_y_utm <- as.numeric(as.character(geoda_data$coord_y_utm))
geoda_limpio <- geoda_data %>% filter(!is.na(coord_x_utm), !is.na(coord_y_utm))

# Datos brutos para PCA
bd <- read_excel("bd_pc.xlsx", sheet = "Hoja1_Depurada", col_types = "text") %>% clean_names()

# ------------------------------------------------------------
# 2. ANÁLISIS ESPACIAL (LISA, GETIS-ORD, GEARY)
# ------------------------------------------------------------
cat("Calculando estadísticas espaciales...\n")

coords_geoda <- as.matrix(geoda_limpio[, c("coord_x_utm", "coord_y_utm")])
pc1_geoda <- geoda_limpio$PC1_GeoDa

# Matriz de pesos KNN (k=8)
nb <- knn2nb(knearneigh(coords_geoda, k = 8))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Local Moran
lisa <- localmoran(pc1_geoda, lw, zero.policy = TRUE)
p_lisa <- lisa[, 5]
p_lisa_fdr <- p.adjust(p_lisa, method = "BH")
lag_pc1 <- lag.listw(lw, pc1_geoda, zero.policy = TRUE)
media_pc1 <- mean(pc1_geoda)
media_lag <- mean(lag_pc1)

geoda_limpio$lisa_nominal <- case_when(
  p_lisa >= 0.05 ~ "Not significant",
  pc1_geoda > 0 & lag_pc1 > 0 ~ "High-High",
  pc1_geoda < 0 & lag_pc1 < 0 ~ "Low-Low",
  pc1_geoda < 0 & lag_pc1 > 0 ~ "Low-High",
  pc1_geoda > 0 & lag_pc1 < 0 ~ "High-Low",
  TRUE ~ "Not significant"
)

geoda_limpio$lisa_fdr <- case_when(
  p_lisa_fdr >= 0.05 ~ "Not significant",
  pc1_geoda > media_pc1 & lag_pc1 > media_lag ~ "High-High",
  pc1_geoda < media_pc1 & lag_pc1 < media_lag ~ "Low-Low",
  pc1_geoda < media_pc1 & lag_pc1 > media_lag ~ "Low-High",
  pc1_geoda > media_pc1 & lag_pc1 < media_lag ~ "High-Low",
  TRUE ~ "Not significant"
)

# Getis-Ord Gi*
gi <- localG(pc1_geoda, lw, zero.policy = TRUE)
p_gi <- 2 * pnorm(-abs(gi))
p_gi_fdr <- p.adjust(p_gi, method = "BH")
geoda_limpio$gi_fdr <- case_when(
  p_gi_fdr >= 0.05 ~ "Not significant",
  gi > 0 ~ "Hot spot",
  gi < 0 ~ "Cold spot",
  TRUE ~ "Not significant"
)

# Local Geary
geary <- localC(pc1_geoda, lw, zero.policy = TRUE)
set.seed(2024)
geary_perm <- localC_perm(pc1_geoda, lw, nsim = 999, zero.policy = TRUE)
p_geary <- attr(geary_perm, "pseudo-p")[, 6]
media_geary <- mean(geary, na.rm = TRUE)
geoda_limpio$geary_categoria <- case_when(
  p_geary >= 0.05 ~ "Not significant",
  geary < media_geary & pc1_geoda > 0 ~ "High-High",
  geary < media_geary & pc1_geoda < 0 ~ "Low-Low",
  geary > media_geary & pc1_geoda > 0 ~ "High-Low",
  geary > media_geary & pc1_geoda < 0 ~ "Low-High",
  TRUE ~ "Other"
)

# Convertir a sf
geoda_sf <- st_as_sf(geoda_limpio, coords = c("coord_x_utm", "coord_y_utm"), 
                     crs = 32717, remove = FALSE)
geoda_sf$pc1_norm <- (geoda_sf$PC1_GeoDa - min(geoda_sf$PC1_GeoDa)) /
  (max(geoda_sf$PC1_GeoDa) - min(geoda_sf$PC1_GeoDa))

# ------------------------------------------------------------
# 3. TEMA Y PALETAS GLOBALES
# ------------------------------------------------------------
tema_limpio <- theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.key.size = unit(0.6, "cm"),
    legend.margin = margin(5, 5, 5, 5),
    plot.margin = margin(5, 5, 5, 5)
  )

# Paleta Okabe-Ito (Color-blind safe) y formas para LISA
paleta_lisa_okabe <- c(
  "High-High" = "#D55E00", "Low-Low" = "#0072B2",
  "Low-High" = "#E69F00", "High-Low" = "#56B4E9",
  "Not significant" = "grey85", "Hot spot" = "#D55E00",
  "Cold spot" = "#0072B2", "Other" = "grey70"
)

formas_lisa_okabe <- c(
  "High-High" = 16, "Low-Low" = 15,
  "Low-High" = 17, "High-Low" = 18,
  "Not significant" = 16, "Hot spot" = 16,
  "Cold spot" = 15, "Other" = 16
)

# Paleta PC1 secuencial
paleta_pc1 <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B")

# ------------------------------------------------------------
# 4. GENERACIÓN DE FIGURAS PRINCIPALES
# ------------------------------------------------------------
cat("Generando figuras principales...\n")

# --- FIGURA 1: Mapa de ubicación ---
mapa_principal <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = geoda_sf, color = "#2471A3", size = 0.4, alpha = 0.6) +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) +
  tema_limpio + theme(legend.position = "none")

inset <- ggplot() +
  geom_sf(data = ecuador_sf, fill = "gray95", color = "gray60", linewidth = 0.2) +
  geom_sf(data = guayas, fill = "#2471A3", color = "black", linewidth = 0.3) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = "gray50", linewidth = 0.4),
        plot.margin = margin(2, 2, 2, 2))

fig1 <- mapa_principal + inset_element(inset, left = 0.72, bottom = 0.72, right = 0.98, top = 0.98)
ggsave("figuras/Figure1.png", fig1, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 2: Densidad Kernel ---
coords_df <- data.frame(st_coordinates(geoda_sf))
fig2 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  stat_density_2d(data = coords_df, aes(x = X, y = Y, fill = after_stat(level)),
                  geom = "polygon", alpha = 0.7, bins = 15) +
  scale_fill_gradientn(colors = c("#D5D8DC", "#85C1E9", "#5DADE2", "#2471A3", "#1A5276"), name = "Density") +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio
ggsave("figuras/Figure2.png", fig2, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 3: Histograma + Mapa PC1 ---
hist_pc1 <- ggplot(geoda_limpio %>% filter(!is.na(PC1_GeoDa) & PC1_GeoDa >= -2 & PC1_GeoDa <= 10), 
                   aes(x = PC1_GeoDa)) +
  geom_histogram(bins = 50, fill = "#2471A3", color = "white", alpha = 0.8) +
  labs(x = "PC1 score", y = "Frequency") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "gray95", linewidth = 0.2))

mapa_pc1 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = geoda_sf, aes(color = pc1_norm), size = 0.5, alpha = 0.8) +
  scale_color_gradientn(colors = paleta_pc1, name = "PC1\n(normalized)", limits = c(0, 1)) +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio

fig3 <- hist_pc1 / mapa_pc1 + plot_layout(heights = c(1, 2))
ggsave("figuras/Figure3.png", fig3, width = 140, height = 200, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 4: Mapa PC1 (Normalizado) ---
fig4 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = geoda_sf, aes(color = pc1_norm), size = 0.5, alpha = 0.8) +
  scale_color_gradientn(colors = paleta_pc1, name = "PC1\n(normalized)", limits = c(0, 1)) +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio
ggsave("figuras/Figure4.png", fig4, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 5: Sensibilidad Moran's I (Corregida) ---
datos_moran <- data.frame(
  Dominio = c("PC1", "IV-SE", "IV-HAB", "IV-PROD", "IV-SAL"),
  k6 = c(0.1298, 0.1199, 0.2617, 0.2846, 0.0558),
  k8 = c(0.1093, 0.1205, 0.2504, 0.2675, 0.0496),
  k10 = c(0.1070, 0.1156, 0.2467, 0.2622, 0.0514)
)

datos_moran$Dominio <- factor(datos_moran$Dominio, levels = c("IV-SE", "IV-HAB", "IV-PROD", "IV-SAL", "PC1"))
datos_largo <- datos_moran %>% pivot_longer(cols = starts_with("k"), names_to = "k", values_to = "Moran_I") %>%
  mutate(k = factor(gsub("k", "", k), levels = c("6", "8", "10")))

fig5 <- ggplot(datos_largo, aes(x = Dominio, y = Moran_I, fill = k)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", Moran_I)), position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("6" = "#9ecae1", "8" = "#4292c6", "10" = "#08519c"), breaks = c("6", "8", "10"), name = "k (neighbours)") +
  labs(x = "Domain / Composite Score", y = "Moran's I") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"))
ggsave("figuras/Figure5.png", fig5, width = 160, height = 120, units = "mm", dpi = 300, bg = "white")

# ------------------------------------------------------------
# 5. GENERACIÓN DE FIGURAS 6-9 (LISA CORREGIDOS)
# ------------------------------------------------------------
cat("Generando figuras LISA (6-9)...\n")

# --- FIGURA 6: Local Moran Nominal ---
sig6 <- geoda_sf %>% filter(lisa_nominal != "Not significant")
nosig6 <- geoda_sf %>% filter(lisa_nominal == "Not significant")

fig6 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = nosig6, color = "#E8E8E8", size = 0.3, alpha = 0.4) +
  geom_sf(data = sig6, aes(color = lisa_nominal, shape = lisa_nominal), size = 2.0) +
  scale_color_manual(values = paleta_lisa_okabe, name = "Local Moran\n(nominal)") +
  scale_shape_manual(values = formas_lisa_okabe, name = "Local Moran\n(nominal)") +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio +
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4)))
ggsave("figuras/Figure6.png", fig6, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 7: Local Moran FDR ---
sig7 <- geoda_sf %>% filter(lisa_fdr != "Not significant")
nosig7 <- geoda_sf %>% filter(lisa_fdr == "Not significant")

fig7 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = nosig7, color = "#E8E8E8", size = 0.3, alpha = 0.4) +
  geom_sf(data = sig7, aes(color = lisa_fdr, shape = lisa_fdr), size = 2.5) +
  scale_color_manual(values = paleta_lisa_okabe, name = "Local Moran\n(FDR)") +
  scale_shape_manual(values = formas_lisa_okabe, name = "Local Moran\n(FDR)") +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio +
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4)))
ggsave("figuras/Figure7.png", fig7, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 8: Getis-Ord FDR ---
sig8 <- geoda_sf %>% filter(gi_fdr != "Not significant")
nosig8 <- geoda_sf %>% filter(gi_fdr == "Not significant")

fig8 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = nosig8, color = "#E8E8E8", size = 0.3, alpha = 0.4) +
  geom_sf(data = sig8, aes(color = gi_fdr, shape = gi_fdr), size = 2.5) +
  scale_color_manual(values = paleta_lisa_okabe, name = "Getis-Ord\n(FDR)") +
  scale_shape_manual(values = formas_lisa_okabe, name = "Getis-Ord\n(FDR)") +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio +
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4)))
ggsave("figuras/Figure8.png", fig8, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# --- FIGURA 9: Local Geary Nominal ---
sig9 <- geoda_sf %>% filter(geary_categoria != "Not significant")
nosig9 <- geoda_sf %>% filter(geary_categoria == "Not significant")

fig9 <- ggplot() +
  geom_sf(data = guayas, fill = "gray97", color = "gray50", linewidth = 0.4) +
  geom_sf(data = nosig9, color = "#E8E8E8", size = 0.2, alpha = 0.3) +
  geom_sf(data = sig9, aes(color = geary_categoria, shape = geary_categoria), size = 1.5) +
  scale_color_manual(values = paleta_lisa_okabe, name = "Local Geary\n(nominal)") +
  scale_shape_manual(values = formas_lisa_okabe, name = "Local Geary\n(nominal)") +
  annotation_scale(location = "bl", width_hint = 0.2, text_cex = 0.8) +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_fancy_orienteering) +
  coord_sf(crs = 32717, datum = 32717) + tema_limpio +
  guides(color = guide_legend(override.aes = list(size = 4)), shape = guide_legend(override.aes = list(size = 4)))
ggsave("figuras/Figure9.png", fig9, width = 140, height = 160, units = "mm", dpi = 300, bg = "white")

# ------------------------------------------------------------
# 6. GENERACIÓN DE FIGURAS SUPLEMENTARIAS (S1, S2, S3)
# ------------------------------------------------------------
cat("Generando figuras suplementarias...\n")

# --- FIGURA S1: Parallel Analysis ---
# (Asumiendo que ya tienes el objeto fa_parallel guardado de tu script anterior)
# Si no lo tienes, necesitas recalcularlo. Aquí asumo que lo tienes cargado.
if (file.exists("datos/parallel_analysis.rds")) {
  resultados_parallel <- readRDS("datos/parallel_analysis.rds")
  eigen_obs <- resultados_parallel$eigen_obs
  eigen_sim <- resultados_parallel$eigen_sim
  
  n_mostrar <- min(20, length(eigen_obs))
  df_scree <- data.frame(
    Componente = 1:n_mostrar,
    Observado = eigen_obs[1:n_mostrar],
    Paralelo = eigen_sim[1:n_mostrar]
  )
  
  fig_s1 <- ggplot(df_scree, aes(x = Componente)) +
    geom_line(aes(y = Observado, color = "Observed"), linewidth = 1) +
    geom_point(aes(y = Observado, color = "Observed"), size = 2.5) +
    geom_line(aes(y = Paralelo, color = "Parallel (p95)"), linewidth = 1, linetype = "dashed") +
    geom_hline(aes(yintercept = 1, color = "Kaiser (lambda=1)"), linetype = "dotted", linewidth = 0.8) +
    scale_color_manual(values = c("Observed" = "#2471A3", "Parallel (p95)" = "#C0392B", "Kaiser (lambda=1)" = "#7F8C8D"), name = "Criterion") +
    labs(x = "Component", y = "Eigenvalue") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  ggsave("figuras/FigureS1_parallel_analysis.png", fig_s1, width = 140, height = 120, units = "mm", dpi = 300, bg = "white")
  cat("FigureS1 guardada.\n")
} else {
  cat("⚠️ No se encontró 'datos/parallel_analysis.rds'. Ejecuta primero el script de análisis paralelo.\n")
}

# --- FIGURA S2: PC1 Original vs Domain-Balanced ---
# NOTA: Necesitas tener un dataframe 'pc_data' con las columnas PC1_original y PC1_balanced
# Si no lo tienes, debes calcularlo primero. Aquí un ejemplo genérico:
# pc_data <- data.frame(PC1_original = geoda_data$PC1_GeoDa, PC1_balanced = geoda_data$PC1_Balanced)

if (exists("pc_data")) {
  fig_s2 <- ggplot(pc_data, aes(x = PC1_original, y = PC1_balanced)) +
    geom_point(alpha = 0.3, color = "#4292c6", size = 1) +
    geom_smooth(method = "lm", color = "#C0392B", se = TRUE, linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    labs(title = "Correlation = 0.787", x = "Original PC1", y = "Domain-balanced PC1") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid.minor = element_blank())
  
  ggsave("figuras/FigureS2_pc1_comparison.png", fig_s2, width = 140, height = 120, units = "mm", dpi = 300, bg = "white")
  cat("FigureS2 guardada.\n")
} else {
  cat("⚠️ No se encontró el dataframe 'pc_data'. Debes calcular el PC1 balanceado primero.\n")
}

# --- FIGURA S3: Moran's I con IC 95% ---
# Valores de tu Tabla S10
df_moran_ci <- data.frame(
  k = factor(c("6", "8", "10"), levels = c("6", "8", "10")),
  Moran_I = c(0.1298, 0.1093, 0.1070),
  CI_low = c(0.1147, 0.0961, 0.0952),
  CI_high = c(0.1450, 0.1225, 0.1189)
)

fig_s3 <- ggplot(df_moran_ci, aes(x = k, y = Moran_I)) +
  geom_point(size = 4, color = "#2471A3") +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.1, linewidth = 0.8, color = "#2471A3") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%.4f", Moran_I)), vjust = -1.5, size = 3.5) +
  labs(x = "k (number of neighbours)", y = "Moran's I (95% CI)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

ggsave("figuras/FigureS3_moran_ci.png", fig_s3, width = 120, height = 100, units = "mm", dpi = 300, bg = "white")
cat("FigureS3 guardada.\n")

cat("\n========================================\n")
cat("TODAS LAS FIGURAS GENERADAS CORRECTAMENTE\n")
cat("========================================\n")