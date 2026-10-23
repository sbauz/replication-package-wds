# ============================================================
# HIGHAM Δ CORREGIDO
# ============================================================

# Asumiendo que tienes bd_pca_filtrado del script anterior

# 1. Matriz tetracórica SIN suavizar
tetra_raw <- tetrachoric(as.matrix(bd_pca_filtrado), smooth = FALSE)
R_raw <- tetra_raw$rho

# 2. Verificar si es PD
eigen_raw <- eigen(R_raw)$values
cat("Matriz raw es PD:", all(eigen_raw > 1e-10), "\n")
cat("Eigenvalue mínimo raw:", round(min(eigen_raw), 6), "\n")

# 3. Aplicar nearPD
R_ajustada <- as.matrix(nearPD(R_raw, corr = TRUE)$mat)

# 4. Calcular Δ
delta_max <- max(abs(R_raw - R_ajustada))
delta_mean <- mean(abs(R_raw - R_ajustada))

cat("Δ máximo (nearPD):", round(delta_max, 6), "\n")
cat("Δ promedio (nearPD):", round(delta_mean, 6), "\n")

# 5. Eigenvalue después del ajuste
eigen_ajustada <- eigen(R_ajustada)$values
cat("Eigenvalue mínimo ajustada:", round(min(eigen_ajustada), 6), "\n")

# 6. Comparar loadings
pca_raw <- principal(R_raw, nfactors = 1, rotate = "none")
pca_ajustada <- principal(R_ajustada, nfactors = 1, rotate = "none")
cor_loadings <- cor(pca_raw$loadings[,1], pca_ajustada$loadings[,1])
cat("Correlación loadings:", round(cor_loadings, 6), "\n")

# 7. Guardar
diagnostico <- data.frame(
  Metric = c("Is raw PD?", "Min eigenvalue (raw)", 
             "Max abs change (nearPD)", "Mean abs change (nearPD)",
             "Min eigenvalue (adjusted)", "Correlation loadings"),
  Value = c(all(eigen_raw > 1e-10), round(min(eigen_raw), 6), 
            round(delta_max, 6), round(delta_mean, 6),
            round(min(eigen_ajustada), 6), round(cor_loadings, 6))
)
print(diagnostico)
write.csv(diagnostico, "outputs_revisores/TableS_Higham_diagnostic_FINAL.csv", 
          row.names = FALSE)
cat("✓ Guardado\n")

# Δ de cor.smooth interno
R_smoothed <- cor.smooth(R_raw)
delta_cor_smooth <- max(abs(R_raw - R_smoothed))
cat("Δ cor.smooth:", round(delta_cor_smooth, 6), "\n")

# ============================================================
# VERIFICACIÓN COMPLETA DEL Δ HIGHAM
# ¿De dónde salió 0.0018?
# ============================================================

setwd("~/Documents/ARTICULOS")
library(readxl)
library(janitor)
library(dplyr)
library(psych)
library(Matrix)

# ------------------------------------------------------------
# 1. CARGAR DATOS Y RECREAR bd_4747
# ------------------------------------------------------------
bd <- read_excel("bd_pc.xlsx", sheet = "Hoja1_Depurada", col_types = "text") %>%
  clean_names()
geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)

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

bd_4747 <- bd_dummies[geoda_data$id_fila, ]

# Variables PCA
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

es_binaria <- sapply(bd_4747, function(x) {
  vals <- unique(na.omit(x))
  length(vals) == 2 && all(vals %in% c("0", "1"))
})

vars_binarias <- setdiff(names(bd_4747)[es_binaria], vars_excluir)
vars_dummies_creadas <- c(
  "mujer", "estado_civil_no_pareja", "no_propietario_vivienda",
  "vivienda_precaria", "agua_no_red_publica", "edu_baja",
  "sin_internet", "no_vive_finca", "sin_luz", "sin_alcantarillado",
  "situacion_laboral_precaria", "no_propietario_predio"
)
vars_pca <- c(vars_binarias, vars_dummies_creadas)

bd_pca <- bd_4747 %>%
  select(all_of(vars_pca)) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))
bd_pca[is.na(bd_pca)] <- 0

# Filtrado (88 variables)
varianza <- sapply(bd_pca, var, na.rm = TRUE)
prevalencia <- sapply(bd_pca, function(x) {
  p <- mean(x == 1, na.rm = TRUE)
  min(p, 1 - p)
})
vars_ok <- names(varianza[varianza > 0 & prevalencia >= 0.005])
bd_pca_filtrado <- bd_pca %>% select(all_of(vars_ok))

cat("Dimensiones bd_pca_filtrado:", dim(bd_pca_filtrado), "\n\n")

# ------------------------------------------------------------
# 2. CALCULAR Δ DE 7 FORMAS DISTINTAS
# ------------------------------------------------------------

# --- Variante 1: tetrachoric() raw + nearPD ---
tetra1 <- tetrachoric(as.matrix(bd_pca_filtrado), smooth = FALSE)
R1 <- tetra1$rho
R1_adj <- nearPD(R1)$mat
delta1 <- max(abs(R1 - R1_adj))

# --- Variante 2: tetrachoric() smooth (default) + nearPD ---
tetra2 <- tetrachoric(as.matrix(bd_pca_filtrado), smooth = TRUE)
R2 <- tetra2$rho
R2_adj <- nearPD(R2)$mat
delta2 <- max(abs(R2 - R2_adj))

# --- Variante 3: Δ de cor.smooth() interno ---
R_raw <- tetra1$rho
R_smooth <- cor.smooth(R_raw)
delta3 <- max(abs(R_raw - R_smooth))

# --- Variante 4: nearPD() sobre la matriz smoothed ---
R4_adj <- nearPD(R_smooth)$mat
delta4 <- max(abs(R_smooth - R4_adj))

# --- Variante 5: polychoric() en lugar de tetrachoric() ---
tetra5 <- polychoric(as.matrix(bd_pca_filtrado))
R5 <- tetra5$rho
R5_adj <- nearPD(R5)$mat
delta5 <- max(abs(R5 - R5_adj))

# --- Variante 6: nearPD() con maxit = 100 (default) ---
tetra6 <- tetrachoric(as.matrix(bd_pca_filtrado), smooth = FALSE)
R6 <- tetra6$rho
R6_adj <- nearPD(R6, maxit = 100)$mat
delta6 <- max(abs(R6 - R6_adj))

# --- Variante 7: nearPD() sobre la matriz con 134 variables (sin filtrar) ---
bd_pca_full <- bd_pca
tetra7 <- tetrachoric(as.matrix(bd_pca_full), smooth = FALSE)
R7 <- tetra7$rho
R7_adj <- nearPD(R7)$mat
delta7 <- max(abs(R7 - R7_adj))

# ------------------------------------------------------------
# 3. MOSTRAR RESULTADOS
# ------------------------------------------------------------
cat("========================================\n")
cat("VERIFICACIÓN DEL Δ HIGHAM\n")
cat("========================================\n\n")

resultados <- data.frame(
  Variante = c(
    "1. tetrachoric(raw) + nearPD (88 vars)",
    "2. tetrachoric(smooth) + nearPD (88 vars)",
    "3. Δ de cor.smooth() interno",
    "4. nearPD sobre smoothed",
    "5. polychoric + nearPD",
    "6. nearPD(maxit=100)",
    "7. tetrachoric(raw) + nearPD (134 vars)"
  ),
  Delta_max = round(c(delta1, delta2, delta3, delta4, delta5, delta6, delta7), 6),
  Delta_mean = round(c(
    mean(abs(R1 - R1_adj)),
    mean(abs(R2 - R2_adj)),
    mean(abs(R_raw - R_smooth)),
    mean(abs(R_smooth - R4_adj)),
    mean(abs(R5 - R5_adj)),
    mean(abs(R6 - R6_adj)),
    mean(abs(R7 - R7_adj))
  ), 6),
  Coincide_0.0018 = c(
    abs(delta1 - 0.0018) < 0.0001,
    abs(delta2 - 0.0018) < 0.0001,
    abs(delta3 - 0.0018) < 0.0001,
    abs(delta4 - 0.0018) < 0.0001,
    abs(delta5 - 0.0018) < 0.0001,
    abs(delta6 - 0.0018) < 0.0001,
    abs(delta7 - 0.0018) < 0.0001
  )
)

print(resultados)

# ------------------------------------------------------------
# 4. CONCLUSIÓN
# ------------------------------------------------------------
cat("\n========================================\n")
cat("CONCLUSIÓN\n")
cat("========================================\n")

if (any(resultados$Coincide_0.0018)) {
  variantes_coinciden <- resultados$Variante[resultados$Coincide_0.0018]
  cat("✅ El valor 0.0018 SÍ corresponde a:\n")
  for (v in variantes_coinciden) cat("   -", v, "\n")
  cat("\nEl manuscrito original estaba en lo correcto para esa variante.\n")
  cat("Debes especificar en el texto CUÁL variante usaste.\n")
} else {
  cat("❌ El valor 0.0018 NO corresponde a NINGUNA variante probada.\n")
  cat("El manuscrito original tenía un error.\n")
  cat("Debes reportar el valor correcto según la variante que uses.\n")
  cat("\nValores reales obtenidos:\n")
  for (i in 1:nrow(resultados)) {
    cat("  -", resultados$Variante[i], ":", resultados$Delta_max[i], "\n")
  }
}

cat("\n========================================\n")
cat("COMPLETADO\n")
cat("========================================\n")
