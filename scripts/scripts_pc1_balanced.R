# ============================================================
# PC1 BALANCEADO POR DOMINIO (R1-6, R3-4)
# Compara el PC1 original con un PC1 balanceado por dominio
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
library(ggplot2)

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

cat("Variables PCA:", length(vars_pca), "\n")

# ------------------------------------------------------------
# 6. CREAR MATRIZ PCA
# ------------------------------------------------------------
bd_pca <- bd_4747 %>%
  select(all_of(vars_pca)) %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

bd_pca[is.na(bd_pca)] <- 0

cat("Dimensiones bd_pca:", dim(bd_pca), "\n")

# ------------------------------------------------------------
# 7. ASIGNAR DOMINIOS A CADA VARIABLE
# ------------------------------------------------------------
cat("\n===== ASIGNANDO DOMINIOS =====\n")

dominios <- case_when(
  grepl("^(sexo|mujer|estado_civil|edu_baja|max_nivel|situacion_laboral|afiliacion|es_migrante|ha_recibido)", 
        names(bd_pca)) ~ "Socioeconomic",
  grepl("^(internet|sin_internet|vivienda|no_vive|luz|sin_luz|alcantarillado|sin_alcantarillado|material|procedencia_agua|agua_no|tenencia_vivienda|no_propietario_vivienda|otra_tenencia)", 
        names(bd_pca)) ~ "Housing",
  grepl("^(tenencia_arroz|tenencia_maiz|preparacion|riego|fertilizantes|plaguicidas|herbicidas|insecticidas|fungibac|epp|sequia|plagas|inund|semilla|malaspracticas|cultivos_adicionales|animales|sembro|vendio|otra_afectacion|otro_|no_propietario_predio|sup_)", 
        names(bd_pca)) ~ "Productive",
  grepl("^(discapacidad|enfermedades|hay_ninnos|ninnos_estudiando|hay_embarazadas|cual_discapacidad|mayores65|personas_enf|es_cuidador|horas_cuidado|sostiene)", 
        names(bd_pca)) ~ "Health",
  TRUE ~ "Other"
)

tabla_dominios <- data.frame(
  Variable = names(bd_pca),
  Dominio = dominios
)

cat("\nDistribucion por dominio:\n")
print(table(tabla_dominios$Dominio))

# ------------------------------------------------------------
# 8. PCA ORIGINAL (TODAS LAS VARIABLES)
# ------------------------------------------------------------
cat("\n===== PCA ORIGINAL (sin balanceo) =====\n")

# Filtrar variables problemáticas
varianza <- sapply(bd_pca, var, na.rm = TRUE)
prevalencia <- sapply(bd_pca, function(x) {
  p <- mean(x == 1, na.rm = TRUE)
  min(p, 1 - p)
})
vars_ok <- names(varianza[varianza > 0 & prevalencia >= 0.005])
bd_pca_filtrado <- bd_pca %>% select(all_of(vars_ok))

cat("Variables finales:", ncol(bd_pca_filtrado), "\n")

tetra <- tetrachoric(as.matrix(bd_pca_filtrado))
R <- as.matrix(nearPD(tetra$rho, corr = TRUE)$mat)
pca_orig <- principal(R, nfactors = ncol(R), rotate = "none", scores = FALSE)

cargas_orig <- pca_orig$loadings[, 1]
scores_orig <- as.numeric(scale(bd_pca_filtrado) %*% cargas_orig)

var_exp_orig <- pca_orig$values[1] / sum(pca_orig$values) * 100

cat("PC1 varianza:", round(var_exp_orig, 2), "%\n")
cat("Componentes lambda>1:", sum(pca_orig$values > 1), "\n")

# ------------------------------------------------------------
# 9. PCA BALANCEADO POR DOMINIO
# ------------------------------------------------------------
cat("\n===== PCA BALANCEADO POR DOMINIO =====\n")

# Función para calcular PC1 dentro de un dominio
pc1_dominio <- function(data, vars) {
  sub <- data[, vars, drop = FALSE]
  
  if (ncol(sub) < 2) return(NULL)
  
  # Filtrar variables problemáticas
  varianza_sub <- sapply(sub, var, na.rm = TRUE)
  prevalencia_sub <- sapply(sub, function(x) {
    p <- mean(x == 1, na.rm = TRUE)
    min(p, 1 - p)
  })
  vars_ok_sub <- names(varianza_sub[varianza_sub > 0 & prevalencia_sub >= 0.005])
  
  if (length(vars_ok_sub) < 2) return(NULL)
  
  sub <- sub[, vars_ok_sub, drop = FALSE]
  
  # PCA
  tetra_sub <- tetrachoric(as.matrix(sub))
  R_sub <- as.matrix(nearPD(tetra_sub$rho, corr = TRUE)$mat)
  pca_sub <- principal(R_sub, nfactors = ncol(R_sub), rotate = "none", scores = FALSE)
  
  cargas_sub <- pca_sub$loadings[, 1]
  scores_sub <- as.numeric(scale(sub) %*% cargas_sub)
  
  return(list(
    cargas = cargas_sub,
    scores = scores_sub,
    var_exp = pca_sub$values[1]/sum(pca_sub$values)*100,
    n_vars = ncol(sub)
  ))
}

# Calcular PC1 para cada dominio
dominios_unicos <- unique(dominios)
dominios_unicos <- dominios_unicos[dominios_unicos != "Other"]

scores_dominios <- data.frame(row.names = 1:nrow(bd_pca))

for (dom in dominios_unicos) {
  vars_dom <- names(bd_pca)[dominios == dom]
  
  cat("\nDominio:", dom, "\n")
  cat("  Variables:", length(vars_dom), "\n")
  
  resultado <- pc1_dominio(bd_pca, vars_dom)
  
  if (!is.null(resultado)) {
    scores_dominios[[dom]] <- resultado$scores
    cat("  Variables usadas:", resultado$n_vars, "\n")
    cat("  PC1 varianza:", round(resultado$var_exp, 2), "%\n")
  } else {
    cat("  No se pudo calcular PCA\n")
  }
}

# Estandarizar y promediar
scores_std <- scale(scores_dominios)
pc1_balanceado <- rowMeans(scores_std, na.rm = TRUE)

# ------------------------------------------------------------
# 10. COMPARACION
# ------------------------------------------------------------
cat("\n========================================\n")
cat("COMPARACION PC1 ORIGINAL vs BALANCEADO\n")
cat("========================================\n")

cor_balanceado <- cor(scores_orig, pc1_balanceado, use = "complete.obs")

cat("Correlacion:", round(cor_balanceado, 4), "\n")
cat("Tucker congruence > 0.90:", cor_balanceado > 0.90, "\n")

# Estadísticas
cat("\nPC1 original:\n")
cat("  Media:", round(mean(scores_orig), 4), "\n")
cat("  SD:", round(sd(scores_orig), 4), "\n")
cat("  Skewness:", round(psych::skew(scores_orig), 4), "\n")

cat("\nPC1 balanceado:\n")
cat("  Media:", round(mean(pc1_balanceado), 4), "\n")
cat("  SD:", round(sd(pc1_balanceado), 4), "\n")
cat("  Skewness:", round(psych::skew(pc1_balanceado), 4), "\n")

# ------------------------------------------------------------
# 11. GRAFICO DE COMPARACION
# ------------------------------------------------------------
df_comp <- data.frame(
  PC1_original = scores_orig,
  PC1_balanceado = pc1_balanceado
)

fig_comp <- ggplot(df_comp, aes(x = PC1_original, y = PC1_balanceado)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#2471A3") +
  geom_smooth(method = "lm", color = "#C0392B", se = TRUE, linewidth = 1) +
  labs(
    x = "PC1 original",
    y = "PC1 balanceado por dominio",
    title = paste0("Correlation = ", round(cor_balanceado, 3))
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
  )

ggsave("figuras/FigureS3_pc1_comparison.png", fig_comp, 
       width = 140, height = 120, units = "mm", dpi = 300, bg = "white")
cat("\nFigureS3_pc1_comparison.png guardada\n")

# ------------------------------------------------------------
# 12. GUARDAR RESULTADOS
# ------------------------------------------------------------
# Tabla S9
tabla_s9 <- data.frame(
  Especificacion = c("Original PC1", "Domain-balanced PC1"),
  Correlacion_con_original = c(1.000, round(cor_balanceado, 4)),
  Media = c(round(mean(scores_orig), 4), round(mean(pc1_balanceado), 4)),
  SD = c(round(sd(scores_orig), 4), round(sd(pc1_balanceado), 4)),
  Skewness = c(round(psych::skew(scores_orig), 4), round(psych::skew(pc1_balanceado), 4))
)
write.csv(tabla_s9, "figuras/TableS9_pc1_balanced.csv", row.names = FALSE)
cat("TableS9_pc1_balanced.csv guardada\n")

# Guardar scores
write.csv(data.frame(
  HH_ID = sprintf("HH_%04d", 1:nrow(bd_pca)),
  PC1_original = scores_orig,
  PC1_balanceado = pc1_balanceado
), "figuras/TableS9b_pc1_scores.csv", row.names = FALSE)
cat("TableS9b_pc1_scores.csv guardada\n")

# Guardar objetos
saveRDS(list(
  scores_orig = scores_orig,
  scores_balanceado = pc1_balanceado,
  scores_dominios = scores_dominios,
  dominios = dominios,
  correlacion = cor_balanceado
), "datos/pc1_balanced.rds")
cat("pc1_balanced.rds guardado\n")

cat("\n========================================\n")
cat("COMPLETADO\n")
cat("========================================\n")
