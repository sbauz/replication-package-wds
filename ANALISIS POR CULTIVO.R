# ============================================================
# ANALISIS POR CULTIVO (R1-4)
# Compara cargas del PC1 entre productores de arroz, maiz y ambos
# ============================================================

# ------------------------------------------------------------
# 0. CONFIGURAR
# ------------------------------------------------------------
setwd("~/Documents/ARTICULOS")
getwd()

# ------------------------------------------------------------
# 1. CARGAR LIBRERIAS
# ------------------------------------------------------------
library(readxl)
library(janitor)
library(dplyr)
library(psych)
library(Matrix)
library(spdep)

# ------------------------------------------------------------
# 2. CARGAR DATOS
# ------------------------------------------------------------
cat("===== CARGANDO DATOS =====\n")

bd <- read_excel("bd_pc.xlsx", sheet = "Hoja1_Depurada", col_types = "text") %>%
  clean_names()

geoda_data <- read_excel("geoda_pca.xlsx", sheet = 1)

cat("bd:", nrow(bd), "x", ncol(bd), "\n")
cat("geoda_data:", nrow(geoda_data), "x", ncol(geoda_data), "\n")

# ------------------------------------------------------------
# 3. CREAR DUMMIES
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
# 4. FILTRAR A 4747
# ------------------------------------------------------------
bd_4747 <- bd_dummies[geoda_data$id_fila, ]

# Identificar cultivo
bd_4747$cultivo <- case_when(
  bd_4747$produccion == "1" ~ "Rice",
  bd_4747$produccion == "2" ~ "Maize",
  bd_4747$produccion == "3" ~ "Both",
  TRUE ~ "None"
)

cat("\n===== DISTRIBUCION POR CULTIVO =====\n")
print(table(bd_4747$cultivo))

# ------------------------------------------------------------
# 5. IDENTIFICAR VARIABLES DEL PCA
# ------------------------------------------------------------
vars_excluir <- c(
  "start", "produccion", "ya_registrado", "num_cedula",
  "fecha_nacimiento", "nombres_productor", "canton", "parroquia",
  "sector", "referencia_finca", "nro_telefono_convencional",
  "celular_contacto", "operadora_cel_contacto", "correo_electronico",
  "latitud_arroz_2024", "longitud_arroz_2024",
  "latitud_maiz_2024", "longitud_maiz_2024",
  "latitud_arroz_2025", "longitud_arroz_2025",
  "latitud_maiz_2025", "longitud_maiz_2025",
  "cultivo"
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

cat("\nVariables PCA:", length(vars_pca), "\n")

# ------------------------------------------------------------
# 6. FUNCION PARA PCA POR CULTIVO
# ------------------------------------------------------------
calcular_pca_cultivo <- function(data, cultivo_nombre, vars) {
  sub <- data %>% filter(cultivo == cultivo_nombre)
  
  if (nrow(sub) < 30) {
    cat("\n", cultivo_nombre, ": muy pocos casos (", nrow(sub), ")\n")
    return(NULL)
  }
  
  bd_pca <- sub %>%
    select(all_of(vars)) %>%
    mutate(across(everything(), ~ as.numeric(as.character(.))))
  bd_pca[is.na(bd_pca)] <- 0
  
  # Filtrar variables problemáticas
  varianza <- sapply(bd_pca, var, na.rm = TRUE)
  prevalencia <- sapply(bd_pca, function(x) {
    p <- mean(x == 1, na.rm = TRUE)
    min(p, 1 - p)
  })
  vars_ok <- names(varianza[varianza > 0 & prevalencia >= 0.005])
  
  bd_pca <- bd_pca %>% select(all_of(vars_ok))
  
  # PCA
  tetra <- tetrachoric(as.matrix(bd_pca))
  R <- as.matrix(nearPD(tetra$rho, corr = TRUE)$mat)
  pca <- principal(R, nfactors = ncol(R), rotate = "none", scores = FALSE)
  
  var_exp <- pca$values / sum(pca$values) * 100
  cargas <- pca$loadings[, 1]
  
  cat("\n===== PCA:", cultivo_nombre, "=====\n")
  cat("N:", nrow(sub), "\n")
  cat("Variables:", ncol(bd_pca), "\n")
  cat("PC1 varianza:", round(var_exp[1], 2), "%\n")
  cat("Componentes lambda>1:", sum(pca$values > 1), "\n")
  
  return(list(
    cultivo = cultivo_nombre,
    n = nrow(sub),
    n_vars = ncol(bd_pca),
    var_exp = var_exp[1],
    n_componentes = sum(pca$values > 1),
    cargas = cargas
  ))
}

# ------------------------------------------------------------
# 7. CALCULAR PCA POR CULTIVO
# ------------------------------------------------------------
pca_rice <- calcular_pca_cultivo(bd_4747, "Rice", vars_pca)
pca_maize <- calcular_pca_cultivo(bd_4747, "Maize", vars_pca)
pca_both <- calcular_pca_cultivo(bd_4747, "Both", vars_pca)

# ------------------------------------------------------------
# 8. CORRELACION DE CARGAS ENTRE CULTIVOS
# ------------------------------------------------------------
cat("\n========================================\n")
cat("CORRELACION DE CARGAS (TUCKER CONGRUENCE)\n")
cat("========================================\n")

# Rice vs Maize
if (!is.null(pca_rice) && !is.null(pca_maize)) {
  vars_comunes <- intersect(names(pca_rice$cargas), names(pca_maize$cargas))
  if (length(vars_comunes) > 10) {
    cor_rice_maize <- cor(pca_rice$cargas[vars_comunes], 
                          pca_maize$cargas[vars_comunes],
                          use = "complete.obs")
    cat("\nRice vs Maize:\n")
    cat("  Variables comunes:", length(vars_comunes), "\n")
    cat("  Correlacion:", round(cor_rice_maize, 4), "\n")
    cat("  Tucker congruence > 0.90:", cor_rice_maize > 0.90, "\n")
  }
}

# Rice vs Both
if (!is.null(pca_rice) && !is.null(pca_both)) {
  vars_comunes <- intersect(names(pca_rice$cargas), names(pca_both$cargas))
  if (length(vars_comunes) > 10) {
    cor_rice_both <- cor(pca_rice$cargas[vars_comunes], 
                         pca_both$cargas[vars_comunes],
                         use = "complete.obs")
    cat("\nRice vs Both:\n")
    cat("  Correlacion:", round(cor_rice_both, 4), "\n")
  }
}

# Maize vs Both
if (!is.null(pca_maize) && !is.null(pca_both)) {
  vars_comunes <- intersect(names(pca_maize$cargas), names(pca_both$cargas))
  if (length(vars_comunes) > 10) {
    cor_maize_both <- cor(pca_maize$cargas[vars_comunes], 
                          pca_both$cargas[vars_comunes],
                          use = "complete.obs")
    cat("\nMaize vs Both:\n")
    cat("  Correlacion:", round(cor_maize_both, 4), "\n")
  }
}

# ------------------------------------------------------------
# 9. MORAN'S I POR CULTIVO
# ------------------------------------------------------------
cat("\n========================================\n")
cat("MORAN'S I POR CULTIVO\n")
cat("========================================\n")

# Cargar coordenadas
coords_geoda <- as.matrix(geoda_data[, c("coord_x_utm", "coord_y_utm")])
pc1_geoda <- geoda_data$PC1_GeoDa

# Identificar cultivo por id_fila
geoda_data$cultivo <- bd_4747$cultivo[match(geoda_data$id_fila, 
                                            1:nrow(bd_4747))]

# Función para Moran por cultivo
moran_por_cultivo <- function(cultivo_nombre) {
  idx <- which(geoda_data$cultivo == cultivo_nombre)
  
  if (length(idx) < 30) {
    cat("\n", cultivo_nombre, ": muy pocos casos (", length(idx), ")\n")
    return(NULL)
  }
  
  coords_sub <- coords_geoda[idx, ]
  pc1_sub <- pc1_geoda[idx]
  
  nb <- knn2nb(knearneigh(coords_sub, k = 8))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  moran <- moran.mc(pc1_sub, lw, nsim = 999, zero.policy = TRUE)
  
  cat("\n", cultivo_nombre, ":\n")
  cat("  N:", length(idx), "\n")
  cat("  Moran I:", round(moran$statistic, 4), "\n")
  cat("  p-value:", moran$p.value, "\n")
  
  return(data.frame(
    Cultivo = cultivo_nombre,
    N = length(idx),
    Moran_I = round(moran$statistic, 4),
    p_value = moran$p.value
  ))
}

resultados_moran <- rbind(
  moran_por_cultivo("Rice"),
  moran_por_cultivo("Maize"),
  moran_por_cultivo("Both")
)

# ------------------------------------------------------------
# 10. GUARDAR RESULTADOS
# ------------------------------------------------------------
cat("\n========================================\n")
cat("GUARDANDO RESULTADOS\n")
cat("========================================\n")

# Tabla S8: PCA por cultivo
tabla_cultivo <- data.frame(
  Cultivo = c(pca_rice$cultivo, pca_maize$cultivo, pca_both$cultivo),
  N = c(pca_rice$n, pca_maize$n, pca_both$n),
  N_vars = c(pca_rice$n_vars, pca_maize$n_vars, pca_both$n_vars),
  PC1_var = round(c(pca_rice$var_exp, pca_maize$var_exp, pca_both$var_exp), 2),
  Componentes = c(pca_rice$n_componentes, pca_maize$n_componentes, pca_both$n_componentes)
)
write.csv(tabla_cultivo, "figuras/TableS8_crop_stratified.csv", row.names = FALSE)
cat("TableS8_crop_stratified.csv guardada\n")

# Tabla S8b: Moran por cultivo
write.csv(resultados_moran, "figuras/TableS8b_moran_crop.csv", row.names = FALSE)
cat("TableS8b_moran_crop.csv guardada\n")

# Guardar objetos
saveRDS(list(
  pca_rice = pca_rice,
  pca_maize = pca_maize,
  pca_both = pca_both,
  moran = resultados_moran
), "datos/crop_stratified.rds")
cat("crop_stratified.rds guardado\n")

cat("\n========================================\n")
cat("COMPLETADO\n")
cat("========================================\n")

# ============================================================
# CORRECCIÓN FINAL DE LA TABLA S6
# Reemplaza los pasos anteriores
# ============================================================

# 1. Cargar el CSV ya corregido (con Labels y renombrado)
tabla_cargas <- read.csv("outputs_revisores/TableS6_PCA_loadings_corrected.csv",
                         stringsAsFactors = FALSE)

# 2. Correcciones finales de dominio
tabla_cargas$Domain[tabla_cargas$Indicator == "hay_discapacidad_hogar"] <- "Health"
tabla_cargas$Domain[tabla_cargas$Indicator == "animales_ninguno"] <- "Productive"

# 3. Verificar que no queden "Other"
cat("===== VERIFICACIÓN =====\n")
cat("Indicadores con dominio 'Other':\n")
otros <- tabla_cargas %>% filter(Domain == "Other")
if (nrow(otros) == 0) {
  cat("✓ No quedan indicadores con dominio 'Other'\n")
} else {
  print(otros %>% select(Indicator, Domain, Label))
}

# 4. Verificar correcciones específicas
cat("\nVerificación de correcciones:\n")
print(tabla_cargas %>% 
        filter(Indicator %in% c("hay_discapacidad_hogar", 
                                "animales_ninguno",
                                "cultivos_adicionales_sembrados_aguacate")) %>% 
        select(Indicator, Domain, Label))

# 5. Resumen por dominio
cat("\n===== RESUMEN POR DOMINIO =====\n")
resumen_final <- tabla_cargas %>%
  group_by(Domain) %>%
  summarise(
    N = n(),
    Mean_Loading = round(mean(Loading), 3),
    N_Meaningful = sum(Meaningful == "Yes"),
    .groups = "drop"
  )
print(resumen_final)

# 6. Verificar total
cat("\nTotal de indicadores:", nrow(tabla_cargas), "\n")
cat("Indicadores meaningful:", sum(tabla_cargas$Meaningful == "Yes"), "\n")

# 7. Guardar versión final
write.csv(tabla_cargas, "outputs_revisores/TableS6_PCA_loadings_FINAL.csv", 
          row.names = FALSE)
cat("\n✓ TableS6_PCA_loadings_FINAL.csv guardada\n")

# 8. (Opcional) Generar tabla LaTeX automáticamente
cat("\n===== GENERANDO TABLA LATEX =====\n")

# Escapar caracteres especiales para LaTeX
escapar_latex <- function(x) {
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# Generar filas LaTeX
filas_latex <- apply(tabla_cargas, 1, function(row) {
  sprintf("%s & %s & %s & %s & %s \\\\",
          escapar_latex(row["Indicator"]),
          escapar_latex(row["Label"]),
          escapar_latex(row["Domain"]),
          row["Loading"],
          row["Meaningful"])
})

# Escribir tabla completa
tabla_latex <- c(
  "\\begin{longtable}{p{4.5cm}p{4cm}p{1.8cm}cc}",
  "\\caption{Complete PC1 loading table for the 88 indicators retained after prevalence filtering. Loadings with absolute value $\\geq 0.30$ are flagged as meaningful. Complementary binary indicators (e.g., \\texttt{cultivos\\_adicionales\\_sembrados\\_ninguno} vs specific crops) appear with opposite loadings because they encode the same condition with inverted values.}",
  "\\label{tbl:pca_loadings} \\\\",
  "\\toprule",
  "Indicator & Label & Domain & Loading & Meaningful \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Indicator & Label & Domain & Loading & Meaningful \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot",
  filas_latex,
  "\\end{longtable}"
)

writeLines(tabla_latex, "outputs_revisores/TableS6_PCA_loadings.tex")
cat("✓ TableS6_PCA_loadings.tex guardada (lista para pegar en el .tex)\n")

cat("\n========================================\n")
cat("COMPLETADO\n")
cat("========================================\n")
# ============================================================
# CORRECCIÓN FINAL: ELIMINAR TODOS LOS "OTHER"
# ============================================================

# 1. Animales → Productive
tabla_cargas$Domain[tabla_cargas$Indicator %in% c(
  "animales_aves", "animales_bovino", "animales_porcino"
)] <- "Productive"

# 2. hay_ninnos → Health (carga de cuidado)
tabla_cargas$Domain[tabla_cargas$Indicator == "hay_ninnos"] <- "Health"

# 3. es_migrante → Socioeconomic
tabla_cargas$Domain[tabla_cargas$Indicator == "es_migrante"] <- "Socioeconomic"

# 4. Verificar que no queden "Other"
cat("===== VERIFICACIÓN FINAL =====\n")
otros <- tabla_cargas %>% filter(Domain == "Other")
if (nrow(otros) == 0) {
  cat("✓ No quedan indicadores con dominio 'Other'\n")
} else {
  cat("⚠️ Aún quedan", nrow(otros), "indicadores en 'Other':\n")
  print(otros %>% select(Indicator, Label, Loading))
}

# 5. Resumen final
cat("\n===== RESUMEN FINAL POR DOMINIO =====\n")
resumen_final <- tabla_cargas %>%
  group_by(Domain) %>%
  summarise(
    N = n(),
    Mean_Loading = round(mean(Loading), 3),
    N_Meaningful = sum(Meaningful == "Yes"),
    .groups = "drop"
  ) %>%
  arrange(desc(N))
print(resumen_final)

# 6. Verificar total
cat("\nTotal de indicadores:", nrow(tabla_cargas), "\n")
cat("Indicadores meaningful:", sum(tabla_cargas$Meaningful == "Yes"), "\n")

# 7. Guardar versión definitiva
write.csv(tabla_cargas, "outputs_revisores/TableS6_PCA_loadings_FINAL.csv", 
          row.names = FALSE)
cat("\n✓ TableS6_PCA_loadings_FINAL.csv guardada\n")

# 8. Regenerar tabla LaTeX
cat("\n===== REGENERANDO TABLA LATEX =====\n")

escapar_latex <- function(x) {
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("#", "\\\\#", x)
  x
}

filas_latex <- apply(tabla_cargas, 1, function(row) {
  sprintf("%s & %s & %s & %s & %s \\\\",
          escapar_latex(row["Indicator"]),
          escapar_latex(row["Label"]),
          escapar_latex(row["Domain"]),
          row["Loading"],
          row["Meaningful"])
})

tabla_latex <- c(
  "\\begin{longtable}{p{4.5cm}p{4cm}p{1.8cm}cc}",
  "\\caption{Complete PC1 loading table for the 88 indicators retained after prevalence filtering. Loadings with absolute value $\\geq 0.30$ are flagged as meaningful. Complementary binary indicators (e.g., \\texttt{cultivos\\_adicionales\\_sembrados\\_ninguno} vs specific crops) appear with opposite loadings because they encode the same condition with inverted values.}",
  "\\label{tbl:pca_loadings} \\\\",
  "\\toprule",
  "Indicator & Label & Domain & Loading & Meaningful \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Indicator & Label & Domain & Loading & Meaningful \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot",
  filas_latex,
  "\\end{longtable}"
)

writeLines(tabla_latex, "outputs_revisores/TableS6_PCA_loadings.tex")
cat("✓ TableS6_PCA_loadings.tex regenerada\n")

cat("\n========================================\n")
cat("COMPLETADO\n")
cat("========================================\n")
# Cargar el archivo FINAL
final <- read.csv("outputs_revisores/TableS6_PCA_loadings_FINAL.csv", 
                  stringsAsFactors = FALSE)

# Verificar que no hay "Other"
cat("Indicadores en 'Other':", sum(final$Domain == "Other"), "\n")

# Verificar total
cat("Total indicadores:", nrow(final), "\n")
cat("Meaningful:", sum(final$Meaningful == "Yes"), "\n")

# Verificar dominios
print(table(final$Domain))

# Verificar hay_discapacidad_hogar
cat("\nhay_discapacidad_hogar:\n")
print(final %>% filter(Indicator == "hay_discapacidad_hogar") %>% 
        select(Indicator, Domain, Loading, Meaningful))
# ============================================================
# GENERAR TABLA LATEX FINAL
# ============================================================

# 1. Cargar CSV final
final <- read.csv("outputs_revisores/TableS6_PCA_loadings_FINAL.csv", 
                  stringsAsFactors = FALSE)

# 2. Verificar
cat("Indicadores:", nrow(final), "\n")
cat("Other:", sum(final$Domain == "Other"), "\n")
cat("Meaningful:", sum(final$Meaningful == "Yes"), "\n")

# 3. Escapar caracteres LaTeX
escapar_latex <- function(x) {
  x <- gsub("_", "\\\\_", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# 4. Generar filas
filas_latex <- apply(final, 1, function(row) {
  sprintf("%s & %s & %s & %s & %s & %s \\\\",
          escapar_latex(row["Indicator"]),
          escapar_latex(row["Label"]),
          escapar_latex(row["Domain"]),
          row["Loading"],
          row["Direction"],
          row["Meaningful"])
})

# 5. Construir tabla
tabla_latex <- c(
  "\\begin{longtable}{p{4.2cm}p{4.5cm}p{1.6cm}ccc}",
  "\\caption{Complete PC1 loading table for the 88 indicators retained after prevalence filtering. Loadings with absolute value $\\geq 0.30$ are flagged as meaningful. Complementary binary indicators (e.g., \\texttt{cultivos\\_adicionales\\_sembrados\\_ninguno} vs specific crops) appear with opposite loadings because they encode the same condition with inverted values.}",
  "\\label{tbl:pca_loadings} \\\\",
  "\\toprule",
  "Indicator & Label & Domain & Loading & Direction & Meaningful \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Indicator & Label & Domain & Loading & Direction & Meaningful \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot",
  filas_latex,
  "\\end{longtable}"
)

# 6. Guardar
writeLines(tabla_latex, "outputs_revisores/TableS6_PCA_loadings_FINAL.tex")
cat("✓ TableS6_PCA_loadings_FINAL.tex guardada\n")

