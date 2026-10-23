# ============================================================
# R2-9: SPATIAL BLOCK BOOTSTRAP (VERSION FINAL CORREGIDA)
# Sin ruido en coordenadas, eliminando duplicados
# ============================================================

setwd("~/Documents/ARTICULOS")
library(readxl)
library(dplyr)
library(spdep)
library(ggplot2)

cat("===== CARGANDO DATOS =====\n")

# 1. Cargar datos
geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)
geoda_data$coord_x_utm <- as.numeric(as.character(geoda_data$coord_x_utm))
geoda_data$coord_y_utm <- as.numeric(as.character(geoda_data$coord_y_utm))
geoda_limpio <- geoda_data %>% filter(!is.na(coord_x_utm), !is.na(coord_y_utm))

# Necesitamos el canton
bd <- read_excel("bd_pc.xlsx", sheet = "Hoja1_Depurada", col_types = "text") %>% janitor::clean_names()
bd_4747 <- bd[geoda_limpio$id_fila, ]
geoda_limpio$canton <- bd_4747$canton

coords_geoda <- as.matrix(geoda_limpio[, c("coord_x_utm", "coord_y_utm")])
pc1_geoda <- geoda_limpio$PC1_GeoDa
canton_vec <- geoda_limpio$canton
cantones_validos <- names(table(canton_vec))[table(canton_vec) >= 50]

# Moran Original (Referencia)
nb_orig <- knn2nb(knearneigh(coords_geoda, k = 8))
lw_orig <- nb2listw(nb_orig, style = "W", zero.policy = TRUE)
moran_orig <- moran.test(pc1_geoda, lw_orig, zero.policy = TRUE)

cat("\nMoran I original:", round(moran_orig$estimate[1], 4), "\n")
cat("SE analitico:", round(sqrt(moran_orig$estimate[3]), 4), "\n")

# ------------------------------------------------------------
# FUNCION AUXILIAR: Calcular Moran sin duplicados
# ------------------------------------------------------------
calc_moran <- function(idx) {
  pc1_sub <- pc1_geoda[idx]
  coords_sub <- coords_geoda[idx, ]
  
  # Eliminar puntos duplicados (mismas coordenadas)
  df_sub <- data.frame(x = coords_sub[,1], y = coords_sub[,2], pc1 = pc1_sub)
  df_sub <- df_sub[!duplicated(df_sub[, c("x", "y")]), ]
  
  if (nrow(df_sub) < 20) return(NULL)
  
  # Recalcular pesos
  coords_clean <- as.matrix(df_sub[, c("x", "y")])
  pc1_clean <- df_sub$pc1
  
  nb <- tryCatch(knn2nb(knearneigh(coords_clean, k = 8)), error = function(e) NULL)
  if (is.null(nb)) return(NULL)
  
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  moran <- tryCatch(moran.test(pc1_clean, lw, zero.policy = TRUE), error = function(e) NULL)
  if (is.null(moran)) return(NULL)
  
  return(moran$estimate[1])
}

# ------------------------------------------------------------
# A. IID BOOTSTRAP (100 iteraciones)
# ------------------------------------------------------------
cat("\n===== IID BOOTSTRAP (Rompe la dependencia) =====\n")
set.seed(2024)
n_iter <- 100
moran_iid <- numeric(n_iter)

for (i in 1:n_iter) {
  idx <- sample(1:nrow(coords_geoda), replace = TRUE)
  res <- calc_moran(idx)
  if (!is.null(res)) moran_iid[i] <- res
}
moran_iid <- moran_iid[moran_iid != 0]

# ------------------------------------------------------------
# B. BLOCK BOOTSTRAP (100 iteraciones)
# ------------------------------------------------------------
cat("===== BLOCK BOOTSTRAP (Preserva la dependencia) =====\n")
set.seed(2024)
moran_block <- numeric(n_iter)

for (i in 1:n_iter) {
  cantones_boot <- sample(cantones_validos, length(cantones_validos), replace = TRUE)
  idx <- unlist(lapply(cantones_boot, function(c) which(canton_vec == c)))
  res <- calc_moran(idx)
  if (!is.null(res)) moran_block[i] <- res
}
moran_block <- moran_block[moran_block != 0]

# ------------------------------------------------------------
# C. RESULTADOS
# ------------------------------------------------------------
se_iid <- sd(moran_iid)
se_block <- sd(moran_block)
ci_iid <- quantile(moran_iid, c(0.025, 0.975))
ci_block <- quantile(moran_block, c(0.025, 0.975))

cat("\n========================================\n")
cat("RESULTADOS: IID vs BLOCK BOOTSTRAP\n")
cat("========================================\n")
cat("Observado:", round(moran_orig$estimate[1], 4), "\n")
cat("SE Analitico:", round(sqrt(moran_orig$estimate[3]), 4), "\n\n")

cat("IID Bootstrap (Singh, 1981):\n")
cat("  Media:", round(mean(moran_iid), 4), "\n")
cat("  SE:", round(se_iid, 4), "\n")
cat("  IC 95%:", round(ci_iid[1], 4), "-", round(ci_iid[2], 4), "\n\n")

cat("Block Bootstrap (Lahiri & Zhu, 2006):\n")
cat("  Media:", round(mean(moran_block), 4), "\n")
cat("  SE:", round(se_block, 4), "\n")
cat("  IC 95%:", round(ci_block[1], 4), "-", round(ci_block[2], 4), "\n\n")

cat("Ratio (Block / IID):", round(se_block/se_iid, 2), "x\n")
cat("Conclusion: El IID subestima el error por un factor de", round(se_block/se_iid, 1), "x.\n")

# ------------------------------------------------------------
# D. GRAFICO
# ------------------------------------------------------------
df_plot <- rbind(
  data.frame(Metodo = "IID Bootstrap", Moran_I = moran_iid),
  data.frame(Metodo = "Block Bootstrap", Moran_I = moran_block)
)

fig <- ggplot(df_plot, aes(x = Moran_I, fill = Metodo)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(xintercept = moran_orig$estimate[1], color = "black", linewidth = 1.2, linetype = "dashed") +
  scale_fill_manual(values = c("IID Bootstrap" = "#85C1E9", "Block Bootstrap" = "#C0392B")) +
  labs(
    x = "Moran's I",
    y = "Density",
    title = "IID vs Spatial Block Bootstrap",
    subtitle = paste0("Dashed line = Observed (", round(moran_orig$estimate[1], 4), ")")
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("figuras/FigureS5_bootstrap_comparison.png", fig, width = 140, height = 100, units = "mm", dpi = 300, bg = "white")
cat("\nFigureS5_bootstrap_comparison.png guardada\n")

# ------------------------------------------------------------
# E. TABLA
# ------------------------------------------------------------
tabla <- data.frame(
  Metodo = c("Observado", "IID Bootstrap", "Block Bootstrap"),
  Moran_I = c(moran_orig$estimate[1], mean(moran_iid), mean(moran_block)),
  SE = c(sqrt(moran_orig$estimate[3]), se_iid, se_block),
  CI_Lower = c(NA, ci_iid[1], ci_block[1]),
  CI_Upper = c(NA, ci_iid[2], ci_block[2])
)
write.csv(tabla, "figuras/TableS12_bootstrap_comparison.csv", row.names = FALSE)
cat("TableS12_bootstrap_comparison.csv guardada\n")

cat("\n========================================\n")
cat("COMPLETADO\n")
cat("========================================\n")

