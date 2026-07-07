# ==============================================================================
# Evaluación de la distribución nula y robustez del estadístico F bajo 
# diferentes condiciones de ruido (nula, megáfono y outliers).
# ==============================================================================

# 1. Cargar librerías, dependencias y establece tema común para gráficos
if (!require("pacman")) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(ggplot2, dplyr, data.table, tidyr, ggsci, patchwork)
set.seed(220626)
temaFiguras <- theme_classic(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, 
    face = "bold", size = 13, color = "#111111"), plot.subtitle = element_text(hjust = 0.5, 
    size = 10, face = "italic", color = "#444444"), axis.ticks.x = element_blank(), 
    strip.background = element_rect(fill = "#F2F2F2", color = "#999999", 
      linewidth = 0.5), strip.text = element_text(face = "bold", 
      size = 10, color = "#222222"), panel.grid.major.y = element_line(color = "#EAEAEA", 
      linewidth = 0.5), legend.position = "bottom", legend.title = element_text(face = "bold", 
      size = 9), legend.text = element_text(size = 9), 
    panel.spacing.y = unit(1, "lines"))
set_theme(temaFiguras)

configGraficas <- list(
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))),
  coord_cartesian(clip = "on"),
  paletteer::scale_fill_paletteer_d("ggsci::lanonc_lancet"),
  theme_minimal()
)

# 2. Cargar datos y funciones
load("data/allResults.RData") 
source("funcionesGluta.R")
source("funcionesPvalores.R")

# Preparación de datos
matriz <- as.matrix(all_data$BA46_glutamatergic$data)
tods <- all_data$BA46_glutamatergic$time
medianSd <- median(apply(matriz, 1, sd))
n_pacientes <- length(tods)
n_genesNula <- 50000

# 3. Función auxiliar para generar y graficar 
analizar_ruido <- function(matriz_ruido, tods, n_pacientes, titulo_prefijo) {
  
  pvals <- pvalsAnaliticos(matriz_ruido, tods, tau = 2*pi)
  
  # Histograma
  p1 <- ggplot(pvals, aes(x = Fobs)) +
    geom_histogram(aes(y = after_stat(density)), boundary = 0, fill = "#00468BFF", bins = 100, alpha = 0.8) +
    stat_function(color = "#E64B35FF", fun = df, args = list(df1 = 2, df2 = n_pacientes - 3), linewidth = 1) +
    labs(x = "Estadístico F", y = "Densidad", title = titulo_prefijo) +
    configGraficas
  
  # Q-Q Plot
  p2 <- ggplot(pvals, aes(sample = Fobs)) +
    stat_qq(distribution = stats::qf, dparams = list(df1 = 2, df2 = n_pacientes - 3), color = "#00468BFF", alpha = 0.5) +
    stat_qq_line(distribution = stats::qf, dparams = list(df1 = 2, df2 = n_pacientes - 3), color = "#E64B35FF", linewidth = 1) +
    labs(x = "Cuantiles teóricos", y = "Cuantiles observados") +
    configGraficas
  
  # Test estadístico
  ks <- ks.test(pvals$Fobs, function(x) pf(x, df1 = 2, df2 = n_pacientes - 3))
  print(ks)
  
  return(p1 + p2)
}

# 4. Ejecución de escenarios
# Escenario 1: Nula
ruidoNula <- matrix(rnorm(n_genesNula * n_pacientes, 0, medianSd), nrow = n_genesNula)
plot1 <- analizar_ruido(ruidoNula, tods, n_pacientes, "Nula")

# Escenario 2: Megáfono
multiplicadores_sd <- seq(0.5, 3.0, length.out = n_pacientes)
ruidoMegafono <- t(apply(matrix(0, n_genesNula, n_pacientes), 1, function(x) rnorm(n_pacientes, 0, medianSd * multiplicadores_sd)))
plot2 <- analizar_ruido(ruidoMegafono, tods, n_pacientes, "Megáfono")

# Escenario 3: Outliers
ruidoOutliers <- matrix(rnorm(n_genesNula * n_pacientes, 0, medianSd), nrow = n_genesNula)
indicesOutliers <- sample(1:(n_genesNula * n_pacientes), floor(0.02 * n_genesNula * n_pacientes))
ruidoOutliers[indicesOutliers] <- ruidoOutliers[indicesOutliers]+(5 * medianSd)
plot3 <- analizar_ruido(ruidoOutliers, tods, n_pacientes, "Outliers")

# 5. Guardar resultados
if(!dir.exists("outputs")) dir.create("outputs")
ggsave("outputs/nula_distribucion.pdf", plot = plot1, width = 8, height = 4)
ggsave("outputs/megafono_distribucion.pdf", plot = plot2, width = 8, height = 4)
ggsave("outputs/outliers_distribucion.pdf", plot = plot3, width = 8, height = 4)
