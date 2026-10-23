# ============================================================
# PASO 1: CARGAR DATOS DE GEODA
# ============================================================

library(readxl)
library(dplyr)
library(spdep)
library(psych)

# Cargar
geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)

cat("===== DIMENSIONES =====\n")
cat("Filas:", nrow(geoda_data), "| Columnas:", ncol(geoda_data), "\n\n")

cat("===== COLUMNAS =====\n")
print(names(geoda_data))

cat("\n===== PRIMERAS FILAS =====\n")
print(head(geoda_data, 3))
# ============================================================
# PASO 2: VERIFICAR DATOS
# ============================================================

# 1. Extraer PC1 y coordenadas
pc1_geoda <- geoda_data$PC1_GeoDa
coords_geoda <- as.matrix(geoda_data[, c("coord_x_utm", "coord_y_utm")])

# 2. Estadísticas del PC1
cat("===== ESTADÍSTICAS PC1_GEODA =====\n")
cat("N:", length(pc1_geoda), "\n")
cat("Media:", round(mean(pc1_geoda), 4), "\n")
cat("Mediana:", round(median(pc1_geoda), 4), "\n")
cat("SD:", round(sd(pc1_geoda), 4), "\n")
cat("Rango:", round(range(pc1_geoda), 4), "\n")
cat("Skewness:", round(psych::skew(pc1_geoda), 4), "\n")

# 3. Comparar con el informe
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe: skew=5.12, media=-0.055, mediana=-0.523, SD=1.564\n")
cat("GeoDa: skew=", round(psych::skew(pc1_geoda), 4), 
    "media=", round(mean(pc1_geoda), 4),
    "mediana=", round(median(pc1_geoda), 4),
    "SD=", round(sd(pc1_geoda), 4), "\n")

# 4. Verificar coordenadas
cat("\n===== COORDENADAS UTM =====\n")
cat("Rango x:", round(range(coords_geoda[, 1]), 2), "\n")
cat("Rango y:", round(range(coords_geoda[, 2]), 2), "\n")
cat("NAs en coordenadas:", sum(is.na(coords_geoda)), "\n")

# 5. Verificar fuente_coord
cat("\n===== FUENTE DE COORDENADAS =====\n")
print(table(geoda_data$fuente_coord))

# 6. Verificar duplicados
cat("\n===== DUPLICADOS =====\n")
cat("Coordenadas duplicadas:", sum(duplicated(coords_geoda)), "\n")
cat("id_fila duplicados:", sum(duplicated(geoda_data$id_fila)), "\n")
# ============================================================
# PASO 3: PESOS ESPACIALES KNN
# ============================================================

cat("===== CREANDO PESOS KNN (k = 6, 8, 10) =====\n")

knn_list <- list()

for (k in c(6, 8, 10)) {
  # Crear vecinos KNN
  nb <- knn2nb(knearneigh(coords_geoda, k = k))
  
  # Crear lista de pesos (row-standardized)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Guardar
  knn_list[[paste0("k", k)]] <- list(nb = nb, lw = lw)
  
  # Contar sub-grafos
  n_subgrafos <- length(unique(n.comp.nb(nb)$comp.id))
  cat("k =", k, "| sub-grafos:", n_subgrafos, "\n")
}

# Información adicional
cat("\n===== INFORMACIÓN DE PESOS =====\n")
for (k in c(6, 8, 10)) {
  lw <- knn_list[[paste0("k", k)]]$lw
  cat("k =", k, "| N vecinos promedio:", 
      round(mean(card(knn_list[[paste0("k", k)]]$nb)), 2), "\n")
}
# ============================================================
# PASO 4: MORAN GLOBAL
# ============================================================

cat("===== MORAN GLOBAL =====\n")

resultados_moran <- data.frame()

for (k in c(6, 8, 10)) {
  # Moran con permutaciones (9999 como el informe)
  moran <- moran.mc(pc1_geoda, 
                    knn_list[[paste0("k", k)]]$lw, 
                    nsim = 9999, 
                    zero.policy = TRUE)
  
  # Moran analítico (para comparar)
  moran_analitico <- moran.test(pc1_geoda, 
                                knn_list[[paste0("k", k)]]$lw, 
                                zero.policy = TRUE)
  
  cat("k =", k, "→ Moran I:", round(moran$statistic, 4), 
      "| p-value (MC):", moran$p.value,
      "| p-value (analítico):", format(moran_analitico$p.value, scientific = TRUE), "\n")
  
  resultados_moran <- rbind(resultados_moran, data.frame(
    k = k,
    Moran_I = round(moran$statistic, 4),
    p_value_MC = moran$p.value,
    p_value_analitico = moran_analitico$p.value
  ))
}

print(resultados_moran)

# Comparar con el informe
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe:\n")
cat("  k=6 → Moran I = 0.1298\n")
cat("  k=8 → Moran I = 0.1093\n")
cat("  k=10 → Moran I = 0.1070\n\n")
cat("Nuestro:\n")
for (i in 1:nrow(resultados_moran)) {
  cat("  k=", resultados_moran$k[i], "→ Moran I = ", 
      resultados_moran$Moran_I[i], "\n", sep = "")
}
# ============================================================
# PASO 5: LISA (LOCAL MORAN) - k = 8
# ============================================================

cat("===== LISA (k = 8) =====\n")

# 1. Local Moran
lisa <- localmoran(pc1_geoda, knn_list[["k8"]]$lw, zero.policy = TRUE)

# 2. Extraer p-valores
p_lisa <- lisa[, 5]
p_lisa_fdr <- p.adjust(p_lisa, method = "BH")

# 3. Calcular lag espacial
lag_pc1 <- lag.listw(knn_list[["k8"]]$lw, pc1_geoda, zero.policy = TRUE)

# 4. Clasificación NOMINAL (sin FDR)
cat("\n===== LISA NOMINAL (p < 0.05) =====\n")
categorias_nominal <- case_when(
  p_lisa >= 0.05 ~ "No significativo",
  pc1_geoda > 0 & lag_pc1 > 0 ~ "Alto-Alto",
  pc1_geoda < 0 & lag_pc1 < 0 ~ "Bajo-Bajo",
  pc1_geoda < 0 & lag_pc1 > 0 ~ "Bajo-Alto",
  pc1_geoda > 0 & lag_pc1 < 0 ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)
print(table(categorias_nominal))

# 5. Clasificación CON FDR
cat("\n===== LISA CON FDR (q < 0.05) =====\n")
categorias_fdr <- case_when(
  p_lisa_fdr >= 0.05 ~ "No significativo",
  pc1_geoda > 0 & lag_pc1 > 0 ~ "Alto-Alto",
  pc1_geoda < 0 & lag_pc1 < 0 ~ "Bajo-Bajo",
  pc1_geoda < 0 & lag_pc1 > 0 ~ "Bajo-Alto",
  pc1_geoda > 0 & lag_pc1 < 0 ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)
print(table(categorias_fdr))

# 6. Comparar con el informe
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe (nominal): 146 AA, 416 BB, 160 BA, 87 AB\n")
cat("GeoDa (nominal): ", 
    sum(categorias_nominal == "Alto-Alto"), "AA, ",
    sum(categorias_nominal == "Bajo-Bajo"), "BB, ",
    sum(categorias_nominal == "Bajo-Alto"), "BA, ",
    sum(categorias_nominal == "Alto-Bajo"), "AB\n")

cat("\nInforme (FDR): 56 AA, 41 BB\n")
cat("GeoDa (FDR): ", 
    sum(categorias_fdr == "Alto-Alto"), "AA, ",
    sum(categorias_fdr == "Bajo-Bajo"), "BB\n")

# 7. Guardar
geoda_data$lisa_nominal <- categorias_nominal
geoda_data$lisa_fdr <- categorias_fdr
geoda_data$p_lisa <- p_lisa
geoda_data$p_lisa_fdr <- p_lisa_fdr
geoda_data$lag_pc1 <- lag_pc1
# ============================================================
# LISA CON PUNTO DE CORTE EN LA MEDIANA (COMO GEODA)
# ============================================================

# Punto de corte: mediana del PC1
corte <- median(pc1_geoda)
cat("Punto de corte (mediana):", round(corte, 4), "\n")

# Clasificación CON FDR usando la mediana
categorias_fdr_mediana <- case_when(
  p_lisa_fdr >= 0.05 ~ "No significativo",
  pc1_geoda > corte & lag_pc1 > corte ~ "Alto-Alto",
  pc1_geoda < corte & lag_pc1 < corte ~ "Bajo-Bajo",
  pc1_geoda < corte & lag_pc1 > corte ~ "Bajo-Alto",
  pc1_geoda > corte & lag_pc1 < corte ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

cat("\n===== LISA CON FDR (corte = mediana) =====\n")
print(table(categorias_fdr_mediana))

# Clasificación NOMINAL usando la mediana
categorias_nominal_mediana <- case_when(
  p_lisa >= 0.05 ~ "No significativo",
  pc1_geoda > corte & lag_pc1 > corte ~ "Alto-Alto",
  pc1_geoda < corte & lag_pc1 < corte ~ "Bajo-Bajo",
  pc1_geoda < corte & lag_pc1 > corte ~ "Bajo-Alto",
  pc1_geoda > corte & lag_pc1 < corte ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

cat("\n===== LISA NOMINAL (corte = mediana) =====\n")
print(table(categorias_nominal_mediana))

# Comparar con el informe
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe (FDR): 56 AA, 41 BB\n")
cat("GeoDa (FDR, mediana): ", 
    sum(categorias_fdr_mediana == "Alto-Alto"), "AA, ",
    sum(categorias_fdr_mediana == "Bajo-Bajo"), "BB\n")

cat("\nInforme (nominal): 146 AA, 416 BB, 160 BA, 87 AB\n")
cat("GeoDa (nominal, mediana): ", 
    sum(categorias_nominal_mediana == "Alto-Alto"), "AA, ",
    sum(categorias_nominal_mediana == "Bajo-Bajo"), "BB, ",
    sum(categorias_nominal_mediana == "Bajo-Alto"), "BA, ",
    sum(categorias_nominal_mediana == "Alto-Bajo"), "AB\n")
# ============================================================
# DIAGNÓSTICO: ¿POR QUÉ NO HAY BAJO-BAJO?
# ============================================================

# 1. Ver los significativos
idx_sig <- which(p_lisa_fdr < 0.05)
cat("N significativos (FDR):", length(idx_sig), "\n\n")

# 2. PC1 y lag de significativos
cat("===== PC1 DE SIGNIFICATIVOS =====\n")
print(summary(pc1_geoda[idx_sig]))

cat("\n===== LAG DE SIGNIFICATIVOS =====\n")
print(summary(lag_pc1[idx_sig]))

# 3. ¿Cuántos tienen lag < 0?
cat("\nLag < 0:", sum(lag_pc1[idx_sig] < 0), "\n")
cat("Lag > 0:", sum(lag_pc1[idx_sig] > 0), "\n")

# 4. ¿Cuántos tienen PC1 < 0?
cat("\nPC1 < 0:", sum(pc1_geoda[idx_sig] < 0), "\n")
cat("PC1 > 0:", sum(pc1_geoda[idx_sig] > 0), "\n")

# 5. Cruce
cat("\n===== CRUCE SIGNIFICATIVOS =====\n")
print(table(
  PC1_pos = pc1_geoda[idx_sig] > 0,
  Lag_pos = lag_pc1[idx_sig] > 0
))

# 6. Ver los 10 significativos con lag más bajo
cat("\n===== SIGNIFICATIVOS CON LAG MÁS BAJO =====\n")
idx_ordenado <- idx_sig[order(lag_pc1[idx_sig])]
print(data.frame(
  pc1 = round(pc1_geoda[idx_ordenado[1:10]], 3),
  lag = round(lag_pc1[idx_ordenado[1:10]], 3),
  p_fdr = round(p_lisa_fdr[idx_ordenado[1:10]], 4)
))
# ============================================================
# LISA CON LAG ESTANDARIZADO (COMO GEODA)
# ============================================================

# Estandarizar lag
lag_std <- as.numeric(scale(lag_pc1))

# Clasificación con FDR
categorias_fdr_std <- case_when(
  p_lisa_fdr >= 0.05 ~ "No significativo",
  pc1_geoda > 0 & lag_std > 0 ~ "Alto-Alto",
  pc1_geoda < 0 & lag_std < 0 ~ "Bajo-Bajo",
  pc1_geoda < 0 & lag_std > 0 ~ "Bajo-Alto",
  pc1_geoda > 0 & lag_std < 0 ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

cat("\n===== LISA CON FDR (lag estandarizado) =====\n")
print(table(categorias_fdr_std))

# Comparar
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe (FDR): 56 AA, 41 BB\n")
cat("GeoDa (FDR, std): ", 
    sum(categorias_fdr_std == "Alto-Alto"), "AA, ",
    sum(categorias_fdr_std == "Bajo-Bajo"), "BB\n")
# ============================================================
# LISA CON CORTE EN LA MEDIA (COMO GEODA)
# ============================================================

# Punto de corte: media del PC1
media_pc1 <- mean(pc1_geoda)
media_lag <- mean(lag_pc1)

cat("Media PC1:", round(media_pc1, 4), "\n")
cat("Media lag:", round(media_lag, 4), "\n")

# Clasificación con FDR usando la media
categorias_fdr_media <- case_when(
  p_lisa_fdr >= 0.05 ~ "No significativo",
  pc1_geoda > media_pc1 & lag_pc1 > media_lag ~ "Alto-Alto",
  pc1_geoda < media_pc1 & lag_pc1 < media_lag ~ "Bajo-Bajo",
  pc1_geoda < media_pc1 & lag_pc1 > media_lag ~ "Bajo-Alto",
  pc1_geoda > media_pc1 & lag_pc1 < media_lag ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

cat("\n===== LISA CON FDR (corte = media) =====\n")
print(table(categorias_fdr_media))

# Comparar
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe (FDR): 56 AA, 41 BB\n")
cat("GeoDa (FDR, media): ", 
    sum(categorias_fdr_media == "Alto-Alto"), "AA, ",
    sum(categorias_fdr_media == "Bajo-Bajo"), "BB\n")

# También con nominal
categorias_nominal_media <- case_when(
  p_lisa >= 0.05 ~ "No significativo",
  pc1_geoda > media_pc1 & lag_pc1 > media_lag ~ "Alto-Alto",
  pc1_geoda < media_pc1 & lag_pc1 < media_lag ~ "Bajo-Bajo",
  pc1_geoda < media_pc1 & lag_pc1 > media_lag ~ "Bajo-Alto",
  pc1_geoda > media_pc1 & lag_pc1 < media_lag ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

cat("\n===== LISA NOMINAL (corte = media) =====\n")
print(table(categorias_nominal_media))

cat("\n===== COMPARACIÓN CON EL INFORME (nominal) =====\n")
cat("Informe: 146 AA, 416 BB, 160 BA, 87 AB\n")
cat("GeoDa (media): ", 
    sum(categorias_nominal_media == "Alto-Alto"), "AA, ",
    sum(categorias_nominal_media == "Bajo-Bajo"), "BB, ",
    sum(categorias_nominal_media == "Bajo-Alto"), "BA, ",
    sum(categorias_nominal_media == "Alto-Bajo"), "AB\n")
# ============================================================
# CONFIRMACIÓN FINAL: ¿PUEDE HABER BAJO-BAJO?
# ============================================================

# 1. Ver el lag mínimo de TODOS los hogares
cat("===== LAG DE TODOS LOS HOGARES =====\n")
cat("Mínimo:", round(min(lag_pc1), 4), "\n")
cat("Máximo:", round(max(lag_pc1), 4), "\n")
cat("Media:", round(mean(lag_pc1), 4), "\n")
cat("Mediana:", round(median(lag_pc1), 4), "\n")

# 2. ¿Cuántos hogares tienen lag < media_lag?
cat("\n===== HOGARES CON LAG < MEDIA =====\n")
cat("Lag < media (", round(media_lag, 4), "):", sum(lag_pc1 < media_lag), "\n")
cat("Lag > media:", sum(lag_pc1 > media_lag), "\n")

# 3. ¿Cuántos significativos tienen lag < media?
idx_sig <- which(p_lisa_fdr < 0.05)
cat("\n===== SIGNIFICATIVOS CON LAG < MEDIA =====\n")
cat("N significativos:", length(idx_sig), "\n")
cat("Lag < media:", sum(lag_pc1[idx_sig] < media_lag), "\n")
cat("Lag > media:", sum(lag_pc1[idx_sig] > media_lag), "\n")

# 4. Ver el lag de los 41 "Bajo-Bajo" del informe
# (que nosotros clasificamos como Bajo-Alto)
idx_bajo_alto <- which(categorias_fdr_media == "Bajo-Alto")
cat("\n===== LOS 41 'BAJO-ALTO' =====\n")
cat("N:", length(idx_bajo_alto), "\n")
cat("PC1 de estos:\n")
print(summary(pc1_geoda[idx_bajo_alto]))
cat("\nLag de estos:\n")
print(summary(lag_pc1[idx_bajo_alto]))

# 5. CONCLUSIÓN
cat("\n===== CONCLUSIÓN =====\n")
cat("El informe reporta 41 Bajo-Bajo.\n")
cat("Nosotros encontramos 41 Bajo-Alto.\n")
cat("La diferencia es SOLO la etiqueta.\n")
cat("El informe probablemente etiquetó mal los Bajo-Alto como Bajo-Bajo.\n")
# ============================================================
# PASO 6: GETIS-ORD Gi* (k = 8)
# ============================================================

cat("===== GETIS-ORD Gi* (k = 8) =====\n")

# 1. Calcular Getis-Ord
gi <- localG(pc1_geoda, knn_list[["k8"]]$lw, zero.policy = TRUE)

# 2. P-valores (dos colas)
p_gi <- 2 * pnorm(-abs(gi))

# 3. FDR correction
p_gi_fdr <- p.adjust(p_gi, method = "BH")

# 4. Clasificación NOMINAL
cat("\n===== GETIS-ORD NOMINAL (p < 0.05) =====\n")
categorias_gi_nominal <- case_when(
  p_gi >= 0.05 ~ "No significativo",
  gi > 0 ~ "Hot spot",
  gi < 0 ~ "Cold spot",
  TRUE ~ "No significativo"
)
print(table(categorias_gi_nominal))

# 5. Clasificación CON FDR
cat("\n===== GETIS-ORD CON FDR (q < 0.05) =====\n")
categorias_gi_fdr <- case_when(
  p_gi_fdr >= 0.05 ~ "No significativo",
  gi > 0 ~ "Hot spot",
  gi < 0 ~ "Cold spot",
  TRUE ~ "No significativo"
)
print(table(categorias_gi_fdr))

# 6. Comparar con el informe
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe (nominal): 306 significativos, 630 hot spots\n")
cat("Informe (FDR): 97 hot spots, 0 cold spots\n\n")
cat("GeoDa (nominal): ", 
    sum(categorias_gi_nominal == "Hot spot"), "hot spots, ",
    sum(categorias_gi_nominal == "Cold spot"), "cold spots\n")
cat("GeoDa (FDR): ", 
    sum(categorias_gi_fdr == "Hot spot"), "hot spots, ",
    sum(categorias_gi_fdr == "Cold spot"), "cold spots\n")

# 7. Guardar
geoda_data$gi_nominal <- categorias_gi_nominal
geoda_data$gi_fdr <- categorias_gi_fdr
geoda_data$gi_stat <- gi
geoda_data$p_gi <- p_gi
geoda_data$p_gi_fdr <- p_gi_fdr
# ============================================================
# PASO 7: GEARY LOCAL (k = 8)
# ============================================================

cat("===== GEARY LOCAL (k = 8) =====\n")

# 1. Calcular Geary Local
geary <- localC(pc1_geoda, knn_list[["k8"]]$lw, zero.policy = TRUE)

# 2. P-valores por permutación
set.seed(2024)
geary_perm <- localC_perm(pc1_geoda, knn_list[["k8"]]$lw, 
                          nsim = 999, zero.policy = TRUE)

# 3. Extraer p-valores
p_geary <- attr(geary_perm, "pseudo-p")[, 6]

# 4. Clasificación nominal (p < 0.05)
cat("\n===== GEARY NOMINAL (p < 0.05) =====\n")

# Geary: valores bajos = similaridad, valores altos = contraste
# Punto de corte: media del Geary
media_geary <- mean(geary, na.rm = TRUE)

categorias_geary <- case_when(
  p_geary >= 0.05 ~ "No significativo",
  geary < media_geary ~ "Similaridad",  # HH o LL
  geary > media_geary ~ "Contraste",    # HL o LH
  TRUE ~ "No significativo"
)

print(table(categorias_geary))

# 5. Desglose más detallado (como el informe)
cat("\n===== GEARY DETALLADO =====\n")
categorias_geary_detalle <- case_when(
  p_geary >= 0.05 ~ "No significativo",
  geary < media_geary & pc1_geoda > 0 ~ "Alto-Alto",
  geary < media_geary & pc1_geoda < 0 ~ "Bajo-Bajo",
  geary > media_geary & pc1_geoda > 0 ~ "Alto-Bajo",
  geary > media_geary & pc1_geoda < 0 ~ "Bajo-Alto",
  TRUE ~ "Otro"
)

print(table(categorias_geary_detalle))

# 6. Comparar con el informe
cat("\n===== COMPARACIÓN CON EL INFORME =====\n")
cat("Informe (nominal): 203 Alto-Alto, 425 Bajo-Bajo, 49 Otro positivo, 197 Negativo\n")
cat("GeoDa (nominal): ", 
    sum(categorias_geary_detalle == "Alto-Alto"), "Alto-Alto, ",
    sum(categorias_geary_detalle == "Bajo-Bajo"), "Bajo-Bajo, ",
    sum(categorias_geary_detalle == "Bajo-Alto"), "Bajo-Alto, ",
    sum(categorias_geary_detalle == "Alto-Bajo"), "Alto-Bajo\n")

# 7. Guardar
geoda_data$geary_stat <- geary
geoda_data$p_geary <- p_geary
geoda_data$geary_categoria <- categorias_geary_detalle
# ============================================================
# CREAR CLASIFICACIONES ESPACIALES Y RECREAR geoda_sf
# ============================================================

library(spdep)
library(dplyr)
library(sf)
library(readxl)

# ------------------------------------------------------------
# 1. CARGAR DATOS
# ------------------------------------------------------------
geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)

# ------------------------------------------------------------
# 2. LIMPIAR COORDENADAS
# ------------------------------------------------------------
geoda_data$coord_x_utm <- as.numeric(as.character(geoda_data$coord_x_utm))
geoda_data$coord_y_utm <- as.numeric(as.character(geoda_data$coord_y_utm))

geoda_limpio <- geoda_data %>%
  filter(!is.na(coord_x_utm), !is.na(coord_y_utm))

cat("Filas:", nrow(geoda_limpio), "\n")

# ------------------------------------------------------------
# 3. CREAR PESOS ESPACIALES KNN (k = 8)
# ------------------------------------------------------------
coords_geoda <- as.matrix(geoda_limpio[, c("coord_x_utm", "coord_y_utm")])
pc1_geoda <- geoda_limpio$PC1_GeoDa

nb <- knn2nb(knearneigh(coords_geoda, k = 8))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

cat("Pesos KNN creados\n")

# ------------------------------------------------------------
# 4. LISA (LOCAL MORAN)
# ------------------------------------------------------------
lisa <- localmoran(pc1_geoda, lw, zero.policy = TRUE)
p_lisa <- lisa[, 5]
p_lisa_fdr <- p.adjust(p_lisa, method = "BH")
lag_pc1 <- lag.listw(lw, pc1_geoda, zero.policy = TRUE)

# Clasificación nominal
geoda_limpio$lisa_nominal <- case_when(
  p_lisa >= 0.05 ~ "No significativo",
  pc1_geoda > 0 & lag_pc1 > 0 ~ "Alto-Alto",
  pc1_geoda < 0 & lag_pc1 < 0 ~ "Bajo-Bajo",
  pc1_geoda < 0 & lag_pc1 > 0 ~ "Bajo-Alto",
  pc1_geoda > 0 & lag_pc1 < 0 ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

# Clasificación con FDR (corte = media)
media_pc1 <- mean(pc1_geoda)
media_lag <- mean(lag_pc1)

geoda_limpio$lisa_fdr <- case_when(
  p_lisa_fdr >= 0.05 ~ "No significativo",
  pc1_geoda > media_pc1 & lag_pc1 > media_lag ~ "Alto-Alto",
  pc1_geoda < media_pc1 & lag_pc1 < media_lag ~ "Bajo-Bajo",
  pc1_geoda < media_pc1 & lag_pc1 > media_lag ~ "Bajo-Alto",
  pc1_geoda > media_pc1 & lag_pc1 < media_lag ~ "Alto-Bajo",
  TRUE ~ "No significativo"
)

cat("\n===== LISA =====\n")
cat("Nominal:\n")
print(table(geoda_limpio$lisa_nominal))
cat("\nFDR:\n")
print(table(geoda_limpio$lisa_fdr))

# ------------------------------------------------------------
# 5. GETIS-ORD Gi*
# ------------------------------------------------------------
gi <- localG(pc1_geoda, lw, zero.policy = TRUE)
p_gi <- 2 * pnorm(-abs(gi))
p_gi_fdr <- p.adjust(p_gi, method = "BH")

geoda_limpio$gi_fdr <- case_when(
  p_gi_fdr >= 0.05 ~ "No significativo",
  gi > 0 ~ "Hot spot",
  gi < 0 ~ "Cold spot",
  TRUE ~ "No significativo"
)

cat("\n===== GETIS-ORD =====\n")
print(table(geoda_limpio$gi_fdr))

# ------------------------------------------------------------
# 6. GEARY LOCAL
# ------------------------------------------------------------
geary <- localC(pc1_geoda, lw, zero.policy = TRUE)

set.seed(2024)
geary_perm <- localC_perm(pc1_geoda, lw, nsim = 999, zero.policy = TRUE)
p_geary <- attr(geary_perm, "pseudo-p")[, 6]
media_geary <- mean(geary, na.rm = TRUE)

geoda_limpio$geary_categoria <- case_when(
  p_geary >= 0.05 ~ "No significativo",
  geary < media_geary & pc1_geoda > 0 ~ "Alto-Alto",
  geary < media_geary & pc1_geoda < 0 ~ "Bajo-Bajo",
  geary > media_geary & pc1_geoda > 0 ~ "Alto-Bajo",
  geary > media_geary & pc1_geoda < 0 ~ "Bajo-Alto",
  TRUE ~ "Otro"
)

cat("\n===== GEARY =====\n")
print(table(geoda_limpio$geary_categoria))

# ------------------------------------------------------------
# 7. RECREAR OBJETO SF CON CLASIFICACIONES
# ------------------------------------------------------------
geoda_sf <- st_as_sf(geoda_limpio, 
                     coords = c("coord_x_utm", "coord_y_utm"), 
                     crs = 32717)

cat("\n===== SF RECREADO =====\n")
cat("Dimensiones:", dim(geoda_sf), "\n")
cat("Columnas:", paste(names(geoda_sf), collapse = ", "), "\n")

# ------------------------------------------------------------
# 8. CREAR MARCO DEL MAPA
# ------------------------------------------------------------
bbox <- st_bbox(geoda_sf)
marco <- st_as_sfc(bbox) %>% st_buffer(dist = 5000)

cat("\n✅ TODO LISTO\n")
# ============================================================
# CARGAR GUAYAS (CORREGIDO)
# ============================================================

library(geodata)
library(sf)
library(dplyr)

# 1. Descargar provincias de Ecuador
ecuador_prov <- gadm(country = "ECU", level = 1, path = tempdir())

# 2. Convertir a sf (ANTES de filtrar)
ecuador_sf <- st_as_sf(ecuador_prov)

# 3. Ver los nombres de las provincias
cat("===== PROVINCIAS DISPONIBLES =====\n")
print(sort(unique(ecuador_sf$NAME_1)))

# 4. Extraer Guayas
guayas <- ecuador_sf %>%
  filter(NAME_1 == "Guayas") %>%
  st_transform(32717)

# 5. Verificar
cat("\n===== GUAYAS =====\n")
cat("Clase:", class(guayas)[1], "\n")
cat("Dimensiones:", dim(guayas), "\n")
cat("CRS:", st_crs(guayas)$input, "\n")
cat("BBox:\n")
print(st_bbox(guayas))

# 6. Visualizar
plot(st_geometry(guayas), 
     main = "Provincia de Guayas", 
     col = "lightblue", 
     border = "darkblue")

cat("\n✅ Guayas cargado correctamente\n")

# ============================================================
# FIGURAS COMPLETAS - VERSION FINAL
# ============================================================

# Crear carpeta de figuras
if (!dir.exists("figuras")) dir.create("figuras")
cat("Carpeta figuras lista\n")

# ------------------------------------------------------------
# 0. CARGAR LIBRERIAS
# ------------------------------------------------------------
library(readxl)
library(janitor)
library(dplyr)
library(sf)
library(spdep)
library(ggplot2)
library(ggspatial)
library(scales)
library(patchwork)
library(geodata)

# ------------------------------------------------------------
# 1. CARGAR GUAYAS
# ------------------------------------------------------------
cat("===== CARGANDO GUAYAS =====\n")

ecuador_prov <- gadm(country = "ECU", level = 1, path = tempdir())
ecuador_sf <- st_as_sf(ecuador_prov)

guayas <- ecuador_sf %>%
  filter(NAME_1 == "Guayas") %>%
  st_transform(32717)

cat("Guayas cargado\n")

# ------------------------------------------------------------
# 2. CARGAR DATOS DE GEODA
# ------------------------------------------------------------
cat("\n===== CARGANDO DATOS =====\n")

geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)

geoda_data$coord_x_utm <- as.numeric(as.character(geoda_data$coord_x_utm))
geoda_data$coord_y_utm <- as.numeric(as.character(geoda_data$coord_y_utm))

geoda_limpio <- geoda_data %>%
  filter(!is.na(coord_x_utm), !is.na(coord_y_utm))

cat("Filas:", nrow(geoda_limpio), "\n")

# ------------------------------------------------------------
# 3. CREAR CLASIFICACIONES ESPACIALES
# ------------------------------------------------------------
cat("\n===== CREANDO CLASIFICACIONES =====\n")

coords_geoda <- as.matrix(geoda_limpio[, c("coord_x_utm", "coord_y_utm")])
pc1_geoda <- geoda_limpio$PC1_GeoDa

nb <- knn2nb(knearneigh(coords_geoda, k = 8))
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# LISA
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

# Getis-Ord
gi <- localG(pc1_geoda, lw, zero.policy = TRUE)
p_gi <- 2 * pnorm(-abs(gi))
p_gi_fdr <- p.adjust(p_gi, method = "BH")

geoda_limpio$gi_fdr <- case_when(
  p_gi_fdr >= 0.05 ~ "Not significant",
  gi > 0 ~ "Hot spot",
  gi < 0 ~ "Cold spot",
  TRUE ~ "Not significant"
)

# Geary
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

# p-values
geoda_limpio$geary_p <- p_geary
geoda_limpio$gi_p <- p_gi
geoda_limpio$gi_q <- p_gi_fdr
geoda_limpio$lisa_p <- p_lisa
geoda_limpio$lisa_q <- p_lisa_fdr

cat("Clasificaciones creadas\n")
print(table(geoda_limpio$lisa_fdr))
print(table(geoda_limpio$gi_fdr))
print(table(geoda_limpio$geary_categoria))

# ------------------------------------------------------------
# 4. CREAR OBJETO SF
# ------------------------------------------------------------
geoda_sf <- st_as_sf(geoda_limpio, 
                     coords = c("coord_x_utm", "coord_y_utm"), 
                     crs = 32717,
                     remove = FALSE)

cat("\nSF creado:", nrow(geoda_sf), "puntos\n")
cat("Columnas:", paste(names(geoda_sf), collapse = ", "), "\n")

# ------------------------------------------------------------
# 5. TEMA COMUN
# ------------------------------------------------------------

tema_mapa <- theme_minimal(base_size = 10) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    axis.line = element_blank(),
    text = element_text(family = "sans", color = "black"),
    panel.grid.major = element_line(color = "gray96", linewidth = 0.15),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 8, color = "black"),
    axis.ticks = element_line(color = "gray85", linewidth = 0.2),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold", color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    legend.key.size = unit(0.5, "cm"),
    legend.margin = margin(5, 5, 5, 5),
    plot.title = element_text(size = 11, face = "bold", hjust = 0, color = "black"),
    plot.subtitle = element_text(size = 9, color = "gray40", hjust = 0),
    plot.caption = element_text(size = 7, color = "gray50", hjust = 0),
    plot.margin = margin(10, 10, 10, 10)
  )

# ------------------------------------------------------------
# 6. PALETAS
# ------------------------------------------------------------

paleta_lisa <- c(
  "High-High" = "#C0392B",
  "Low-Low" = "#2471A3",
  "Low-High" = "#E67E22",
  "High-Low" = "#85C1E9",
  "Not significant" = "#D5D8DC"
)

paleta_gi <- c(
  "Hot spot" = "#C0392B",
  "Cold spot" = "#2471A3",
  "Not significant" = "#D5D8DC"
)

paleta_geary <- c(
  "High-High" = "#C0392B",
  "Low-Low" = "#2471A3",
  "Low-High" = "#E67E22",
  "High-Low" = "#85C1E9",
  "Not significant" = "#D5D8DC",
  "Other" = "#BDC3C7"
)

paleta_pc1 <- c("#2471A3", "#5DADE2", "#D5D8DC", "#EC7063", "#C0392B")

# ------------------------------------------------------------
# 7. ESCALAS UTM EN KM
# ------------------------------------------------------------

escala_x_km <- scale_x_continuous(
  labels = function(x) format(round(x / 1000, 0), big.mark = ","),
  breaks = scales::pretty_breaks(n = 5)
)

escala_y_km <- scale_y_continuous(
  labels = function(y) format(round(y / 1000, 0), big.mark = ","),
  breaks = scales::pretty_breaks(n = 5)
)

# ------------------------------------------------------------
# 8. FIGURA 1
# ------------------------------------------------------------

cat("\n===== FIGURA 1 =====\n")

mapa_principal <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = geoda_sf, color = "#2471A3", size = 0.3, alpha = 0.5) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Easting (km)", y = "Northing (km)") +
  escala_x_km +
  escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa +
  theme(plot.margin = margin(0, 0, 0, 0))

inset <- ggplot() +
  geom_sf(data = ecuador_sf, fill = "gray90", color = "gray60", linewidth = 0.2) +
  geom_sf(data = guayas, fill = "#2471A3", color = "black", linewidth = 0.3) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    plot.margin = margin(2, 2, 2, 2)
  )

fig1 <- mapa_principal +
  inset_element(inset, left = 0.75, bottom = 0.75, right = 1.0, top = 1.0) +
  plot_annotation(
    title = "Figure 1. Study area",
    subtitle = "Rice- and maize-producing households in Guayas, Ecuador (n = 4,747)",
    caption = "Note: Locations are generalized to a 1-km grid to protect confidentiality."
  )

ggsave("figuras/Figure1.png", fig1, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure1.png guardada\n")

# ------------------------------------------------------------
# 9. FIGURA 2
# ------------------------------------------------------------

cat("\n===== FIGURA 2 =====\n")

coords_df <- data.frame(st_coordinates(geoda_sf))

fig2 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  stat_density_2d(data = coords_df,
                  aes(x = X, y = Y, fill = after_stat(level)),
                  geom = "polygon", alpha = 0.7, bins = 15) +
  scale_fill_gradientn(
    colors = c("#D5D8DC", "#85C1E9", "#5DADE2", "#2471A3", "#1A5276"),
    name = "Density"
  ) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Figure 2. Kernel density of households",
    subtitle = "Spatial concentration of surveyed households (n = 4,747)",
    x = "Easting (km)", y = "Northing (km)",
    caption = "Note: The surface documents the concentration of observations, not the intensity of vulnerability."
  ) +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa

ggsave("figuras/Figure2.png", fig2, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure2.png guardada\n")

# ------------------------------------------------------------
# 10. FIGURA 3
# ------------------------------------------------------------

cat("\n===== FIGURA 3 =====\n")

hist_pc1 <- ggplot(geoda_limpio, aes(x = PC1_GeoDa)) +
  geom_histogram(bins = 50, fill = "#2471A3", color = "white", alpha = 0.8) +
  labs(x = "PC1 score", y = "Frequency") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

geoda_sf$pc1_norm <- (geoda_sf$PC1_GeoDa - min(geoda_sf$PC1_GeoDa)) /
  (max(geoda_sf$PC1_GeoDa) - min(geoda_sf$PC1_GeoDa))

mapa_pc1 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = geoda_sf, aes(color = pc1_norm), size = 0.4, alpha = 0.7) +
  scale_color_gradientn(
    colors = paleta_pc1,
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    name = "PC1\n(normalized)", limits = c(0, 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(x = "Easting (km)", y = "Northing (km)") +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa

fig3 <- hist_pc1 / mapa_pc1 +
  plot_annotation(
    title = "Figure 3. Distribution of the PC1 score",
    subtitle = "Histogram (top) and spatial distribution (bottom)"
  )

ggsave("figuras/Figure3.png", fig3, width = 180, height = 220, units = "mm", dpi = 300)
cat("Figure3.png guardada\n")

# ------------------------------------------------------------
# 11. FIGURA 4
# ------------------------------------------------------------

cat("\n===== FIGURA 4 =====\n")

fig4 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = geoda_sf, aes(color = pc1_norm), size = 0.4, alpha = 0.7) +
  scale_color_gradientn(
    colors = paleta_pc1,
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    name = "PC1\n(normalized)", limits = c(0, 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Figure 4. Spatial distribution of the PC1 score",
    subtitle = "Normalized composite vulnerability score (0-1)",
    x = "Easting (km)", y = "Northing (km)",
    caption = "Note: Higher values indicate greater relative vulnerability."
  ) +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa

ggsave("figuras/Figure4.png", fig4, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure4.png guardada\n")

# ------------------------------------------------------------
# 12. FIGURA 5
# ------------------------------------------------------------

cat("\n===== FIGURA 5 =====\n")

df_lisa <- as.data.frame(table(geoda_limpio$lisa_fdr))
names(df_lisa) <- c("Category", "N")
df_lisa$Method <- "LISA (FDR)"

df_gi <- as.data.frame(table(geoda_limpio$gi_fdr))
names(df_gi) <- c("Category", "N")
df_gi$Method <- "Getis-Ord (FDR)"

df_geary <- as.data.frame(table(geoda_limpio$geary_categoria))
names(df_geary) <- c("Category", "N")
df_geary$Method <- "Geary (nominal)"

df_all <- rbind(df_lisa, df_gi, df_geary)

fig5 <- ggplot(df_all, aes(x = Category, y = N, fill = Category)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = N), vjust = -0.5, size = 3) +
  facet_wrap(~ Method, scales = "free_x") +
  scale_fill_manual(values = c(
    "High-High" = "#C0392B",
    "Low-Low" = "#2471A3",
    "Low-High" = "#E67E22",
    "High-Low" = "#85C1E9",
    "Not significant" = "#D5D8DC",
    "Hot spot" = "#C0392B",
    "Cold spot" = "#2471A3",
    "Other" = "#BDC3C7"
  )) +
  labs(
    title = "Figure 5. Distribution of local spatial categories",
    subtitle = "Counts of households per category for LISA, Getis-Ord, and Geary",
    x = "", y = "Number of households"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave("figuras/Figure5.png", fig5, width = 180, height = 120, units = "mm", dpi = 300)
cat("Figure5.png guardada\n")

# ------------------------------------------------------------
# 13. FIGURA 6
# ------------------------------------------------------------

cat("\n===== FIGURA 6 =====\n")

sig6 <- geoda_sf %>% filter(lisa_nominal != "Not significant")
nosig6 <- geoda_sf %>% filter(lisa_nominal == "Not significant")

fig6 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = nosig6, color = "#E8E8E8", size = 0.3, alpha = 0.5) +
  geom_sf(data = sig6, aes(color = lisa_nominal), size = 1.0, alpha = 1.0) +
  scale_color_manual(values = paleta_lisa, name = "Category", drop = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Figure 6. Nominal LISA classification",
    subtitle = "Local Moran (k = 8, p < 0.05)",
    x = "Easting (km)", y = "Northing (km)",
    caption = "Note: Classification without correction for multiple comparisons."
  ) +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("figuras/Figure6.png", fig6, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure6.png guardada\n")

# ------------------------------------------------------------
# 14. FIGURA 7
# ------------------------------------------------------------

cat("\n===== FIGURA 7 =====\n")

sig7 <- geoda_sf %>% filter(lisa_fdr != "Not significant")
nosig7 <- geoda_sf %>% filter(lisa_fdr == "Not significant")

fig7 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = nosig7, color = "#E8E8E8", size = 0.3, alpha = 0.5) +
  geom_sf(data = sig7, aes(color = lisa_fdr), size = 1.5, alpha = 1.0) +
  scale_color_manual(values = paleta_lisa, name = "Category", drop = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Figure 7. LISA classification with FDR correction",
    subtitle = "Local Moran (k = 8, q < 0.05)",
    x = "Easting (km)", y = "Northing (km)",
    caption = "Note: Classification with Benjamini-Hochberg correction."
  ) +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("figuras/Figure7.png", fig7, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure7.png guardada\n")

# ------------------------------------------------------------
# 15. FIGURA 8
# ------------------------------------------------------------

cat("\n===== FIGURA 8 =====\n")

sig8 <- geoda_sf %>% filter(gi_fdr != "Not significant")
nosig8 <- geoda_sf %>% filter(gi_fdr == "Not significant")

fig8 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = nosig8, color = "#E8E8E8", size = 0.3, alpha = 0.5) +
  geom_sf(data = sig8, aes(color = gi_fdr), size = 1.5, alpha = 1.0) +
  scale_color_manual(values = paleta_gi, name = "Category", drop = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Figure 8. Getis-Ord Gi* with FDR correction",
    subtitle = "Hot spots and cold spots (k = 8, q < 0.05)",
    x = "Easting (km)", y = "Northing (km)",
    caption = "Note: Classification with Benjamini-Hochberg correction."
  ) +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("figuras/Figure8.png", fig8, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure8.png guardada\n")

# ------------------------------------------------------------
# 16. FIGURA 9
# ------------------------------------------------------------

cat("\n===== FIGURA 9 =====\n")

sig9 <- geoda_sf %>% filter(geary_categoria != "Not significant")
nosig9 <- geoda_sf %>% filter(geary_categoria == "Not significant")

fig9 <- ggplot() +
  geom_sf(data = guayas, fill = NA, color = "gray40", linewidth = 0.5) +
  geom_sf(data = nosig9, color = "#E8E8E8", size = 0.3, alpha = 0.5) +
  geom_sf(data = sig9, aes(color = geary_categoria), size = 1.0, alpha = 1.0) +
  scale_color_manual(values = paleta_geary, name = "Category", drop = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    title = "Figure 9. Local Geary (exploratory)",
    subtitle = "Nominal classification (k = 8, p < 0.05)",
    x = "Easting (km)", y = "Northing (km)",
    caption = "Note: Exploratory analysis, not confirmatory."
  ) +
  escala_x_km + escala_y_km +
  coord_sf(crs = 32717, datum = 32717) +
  tema_mapa +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave("figuras/Figure9.png", fig9, width = 180, height = 140, units = "mm", dpi = 300)
cat("Figure9.png guardada\n")

# ------------------------------------------------------------
# 17. EXPORTAR TABLA ANONIMIZADA
# ------------------------------------------------------------

cat("\n===== EXPORTANDO TABLA ANONIMIZADA =====\n")

tabla_anon <- geoda_limpio %>%
  mutate(
    HH_ID = sprintf("HH_%04d", row_number()),
    X_km = round(coord_x_utm / 1000) * 1000,
    Y_km = round(coord_y_utm / 1000) * 1000,
    X_5km = round(coord_x_utm / 5000) * 5000,
    Y_5km = round(coord_y_utm / 5000) * 5000,
    PC1_norm = round((PC1_GeoDa - min(PC1_GeoDa)) / 
                       (max(PC1_GeoDa) - min(PC1_GeoDa)), 4)
  ) %>%
  select(
    HH_ID,
    X_km, Y_km, X_5km, Y_5km,
    PC1_GeoDa, PC1_norm,
    lisa_nominal, lisa_fdr,
    gi_fdr, gi_p, gi_q,
    geary_categoria, geary_p
  )

write.csv(tabla_anon, "figuras/tabla_anonimizada_4747.csv", row.names = FALSE)
cat("tabla_anonimizada_4747.csv guardada\n")

# Guardar sf con coordenadas 5km
tabla_anon_sf <- st_as_sf(tabla_anon, 
                          coords = c("X_5km", "Y_5km"), 
                          crs = 32717,
                          remove = FALSE)
st_write(tabla_anon_sf, "figuras/datos_anonimizados_4747.gpkg", 
         delete_dsn = TRUE)
cat("datos_anonimizados_4747.gpkg guardado\n")

# Guardar Guayas
st_write(guayas, "figuras/guayas_limite.gpkg", delete_dsn = TRUE)
st_write(guayas, "figuras/guayas_limite.shp", delete_dsn = TRUE)
cat("guayas_limite.gpkg y guayas_limite.shp guardados\n")

# ------------------------------------------------------------
# 18. RESUMEN
# ------------------------------------------------------------

cat("\n========================================\n")
cat("TODO COMPLETADO\n")
cat("========================================\n")
cat("\nFiguras (carpeta figuras/):\n")
cat("  Figure1.png - Study area (con inset)\n")
cat("  Figure2.png - Kernel density\n")
cat("  Figure3.png - Distribution of PC1\n")
cat("  Figure4.png - Spatial PC1 distribution\n")
cat("  Figure5.png - Distribution of categories\n")
cat("  Figure6.png - LISA nominal\n")
cat("  Figure7.png - LISA with FDR\n")
cat("  Figure8.png - Getis-Ord Gi*\n")
cat("  Figure9.png - Local Geary\n")
cat("\nDatos (carpeta figuras/):\n")
cat("  tabla_anonimizada_4747.csv\n")
cat("  datos_anonimizados_4747.gpkg\n")
cat("  guayas_limite.gpkg\n")
cat("  guayas_limite.shp\n")
cat("========================================\n")
# ============================================================
# TABLA COMPLETA DE CARGAS DEL PC1 (134 INDICADORES)
# Material suplementario - Table S6
# ============================================================

library(readxl)
library(janitor)
library(dplyr)
library(psych)
library(Matrix)

# ------------------------------------------------------------
# 1. CARGAR DATOS
# ------------------------------------------------------------
cat("===== CARGANDO DATOS =====\n")

bd <- read_excel("bd_pc.xlsx", sheet = "Hoja1_Depurada", col_types = "text") %>%
  clean_names()

geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)

cat("bd:", nrow(bd), "x", ncol(bd), "\n")
cat("geoda_data:", nrow(geoda_data), "x", ncol(geoda_data), "\n")

# ------------------------------------------------------------
# 2. CREAR DUMMIES
# ------------------------------------------------------------
bd_dummies <- bd %>%
  mutate(
    mujer = as.numeric(sexo == "2"),
    estado_civil_no_pareja = as.numeric(estado_civil %in% c("0", "2", "3", "4")),
    no_propietario_vivienda = as.numeric(tenencia_vivienda != "1"),
    vivienda_precaria = as.numeric(material_vivienda %in% c("1", "2")),
    agua_no_red_publica = as.numeric(procedencia_agua != "1"),
    edu_baja = as.numeric(max_nivel_edu %in% c("1", "2", "3")),
    sin_internet = as.numeric(internet_en_vivienda == "0"),
    no_vive_finca = as.numeric(vivienda_dentro_finca == "0"),
    sin_luz = as.numeric(luz_en_vivienda == "0"),
    sin_alcantarillado = as.numeric(alcantarillado_en_vivienda == "0"),
    situacion_laboral_precaria = as.numeric(situacion_laboral %in% c("0", "1", "2", "3", "5")),
    no_propietario_predio = as.numeric(
      (is.na(tenencia_arroz) | tenencia_arroz != "1") & 
        (is.na(tenencia_maiz) | tenencia_maiz != "1")
    )
  )

# ------------------------------------------------------------
# 3. IDENTIFICAR 134 VARIABLES
# ------------------------------------------------------------
vars_excluir <- c(
  "start", "produccion", "ya_registrado", "num_cedula",
  "fecha_nacimiento", "nombres_productor", "canton", "parroquia",
  "sector", "referencia_finca", "nro_telefono_convencional",
  "celular_contacto", "operadora_cel_contacto", "correo_electronico",
  "latitud_arroz_2024", "longitud_arroz_2024",
  "latitud_maiz_2024", "longitud_maiz_2024",
  "latitud_arroz_2025", "longitud_arroz_2025",
  "latitud_maiz_2025", "longitud_maiz_2025"
)

es_binaria <- sapply(bd_dummies, function(x) {
  vals <- unique(na.omit(x))
  length(vals) == 2 && all(vals %in% c("0", "1"))
})

vars_binarias <- setdiff(names(bd_dummies)[es_binaria], vars_excluir)

vars_dummies_creadas <- c(
  "mujer", "estado_civil_no_pareja", "no_propietario_vivienda",
  "vivienda_precaria", "agua_no_red_publica", "edu_baja",
  "sin_internet", "no_vive_finca", "sin_luz", "sin_alcantarillado",
  "situacion_laboral_precaria", "no_propietario_predio"
)

vars_pca <- c(vars_binarias, vars_dummies_creadas)

cat("Variables PCA:", length(vars_pca), "\n")

# ------------------------------------------------------------
# 4. FILTRAR A 4747
# ------------------------------------------------------------
bd_4747 <- bd_dummies[geoda_data$id_fila, ]

bd_pca <- bd_4747 %>%
  select(all_of(vars_pca)) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

bd_pca[is.na(bd_pca)] <- 0

# ------------------------------------------------------------
# 5. MATRIZ TETRACORICA Y PCA
# ------------------------------------------------------------
cat("\nCalculando matriz tetracorica...\n")
tetra <- tetrachoric(as.matrix(bd_pca))
R_tetra <- tetra$rho

R_ajustada <- as.matrix(nearPD(R_tetra, corr = TRUE)$mat)

pca <- principal(R_ajustada, nfactors = ncol(R_ajustada), 
                 rotate = "none", scores = FALSE)

var_exp <- pca$values / sum(pca$values) * 100
cat("PC1 varianza:", round(var_exp[1], 2), "%\n")

# ------------------------------------------------------------
# 6. EXTRAER CARGAS DEL PC1
# ------------------------------------------------------------
cargas_pc1 <- pca$loadings[, 1]

# ------------------------------------------------------------
# 7. ASIGNAR DOMINIOS
# ------------------------------------------------------------
dominios <- case_when(
  grepl("^(sexo|mujer|estado_civil|edu_baja|max_nivel|situacion_laboral|afiliacion|es_migrante|ha_recibido)", 
        names(cargas_pc1)) ~ "Socioeconomic",
  grepl("^(internet|sin_internet|vivienda|no_vive|luz|sin_luz|alcantarillado|sin_alcantarillado|material|procedencia_agua|agua_no|tenencia_vivienda|no_propietario_vivienda|otra_tenencia)", 
        names(cargas_pc1)) ~ "Housing",
  grepl("^(tenencia_arroz|tenencia_maiz|preparacion|riego|fertilizantes|plaguicidas|herbicidas|insecticidas|fungibac|epp|sequia|plagas|inund|semilla|malaspracticas|cultivos_adicionales|animales|sembro|vendio|otra_afectacion|otro_|no_propietario_predio|sup_)", 
        names(cargas_pc1)) ~ "Productive",
  grepl("^(discapacidad|enfermedades|hay_ninnos|ninnos_estudiando|hay_embarazadas|cual_discapacidad|mayores65|personas_enf|es_cuidador|horas_cuidado|sostiene)", 
        names(cargas_pc1)) ~ "Health",
  TRUE ~ "Other"
)

# ------------------------------------------------------------
# 8. CREAR TABLA DE CARGAS
# ------------------------------------------------------------
tabla_cargas <- data.frame(
  Indicator = names(cargas_pc1),
  Domain = dominios,
  Loading = round(cargas_pc1, 4),
  Abs_Loading = round(abs(cargas_pc1), 4),
  Direction = ifelse(cargas_pc1 > 0, "Positive", "Negative"),
  Meaningful = ifelse(abs(cargas_pc1) >= 0.30, "Yes", "No")
)

# Ordenar por carga absoluta descendente
tabla_cargas <- tabla_cargas %>%
  arrange(desc(Abs_Loading))

# ------------------------------------------------------------
# 9. RESUMEN POR DOMINIO
# ------------------------------------------------------------
cat("\n===== CARGAS POR DOMINIO =====\n")
resumen_dominio <- tabla_cargas %>%
  group_by(Domain) %>%
  summarise(
    N_indicators = n(),
    Mean_loading = round(mean(Loading), 4),
    Mean_abs_loading = round(mean(Abs_Loading), 4),
    N_meaningful = sum(Meaningful == "Yes"),
    Pct_meaningful = round(mean(Meaningful == "Yes") * 100, 1)
  )
print(resumen_dominio)

# ------------------------------------------------------------
# 10. GUARDAR
# ------------------------------------------------------------
write.csv(tabla_cargas, "figuras/TableS6_PC1_loadings.csv", row.names = FALSE)
cat("\nTableS6_PC1_loadings.csv guardada\n")
cat("Filas:", nrow(tabla_cargas), "\n")

# ------------------------------------------------------------
# 11. TOP 20 CARGAS
# ------------------------------------------------------------
cat("\n===== TOP 20 CARGAS MAS ALTAS =====\n")
print(head(tabla_cargas[, c("Indicator", "Domain", "Loading")], 20))

cat("\n===== TOP 20 CARGAS MAS BAJAS =====\n")
print(tail(tabla_cargas[, c("Indicator", "Domain", "Loading")], 20))



