# ==============================================================================
# Análisis sobre datos transcriptómicos para diferentes 
# métodos de agregación y cálculo de p-valores (Analítico vs Bootstrap).
# ==============================================================================

# 1. Cargar librerías, dependencias y establece tema común para gráficos
if (!require("pacman")) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(tidyverse, data.table, ggplot2, patchwork, paletteer, parallel, doParallel, doRNG, progress, UpSetR)
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
  paletteer::scale_fill_paletteer_d("ggsci::lanonc_lancet")
)

# 2. Cargar datos y funciones
load("data/allResults.RData") 
source("funcionesGluta.R")
source("funcionesPvalores.R")

# Preparación de datos
matriz <- as.matrix(all_data$BA46_glutamatergic$data)
tods <- all_data$BA46_glutamatergic$time

################################################################################
# BH directo
################################################################################

#Cálculo pvalores analíticos
t0Analitico<-Sys.time()
resultadosAnaliticos<-pvalsAnaliticos(matriz, tods, tau=2*pi)
timeAnalitico<-Sys.time()-t0Analitico
print(timeAnalitico)#Time difference of 0.2663162 secs
resultadosAnaliticos$pvalFAdjust<-p.adjust(resultadosAnaliticos$pvalF, method = "BH")

#Cálculo pvalores bootstrap
nChunks<-50
t0Boot<-Sys.time()
numCores<-parallel::detectCores() - 1
cl<-parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl)
resultadosBootstrap<-pvalsBootstrap(matriz, tods, tau=2*pi, B=10000,n_chunks = nChunks)
parallel::stopCluster(cl)
timeBoot<-Sys.time()-t0Boot
print(timeBoot)#Time difference of 51.72898 secs
resultadosBootstrap$pvalBootstrapAdjust<-p.adjust(resultadosBootstrap$pvalBoot, method = "BH")

df_pval <- bind_rows(
  data.frame(pvalor = resultadosAnaliticos$pvalF, Tipo = "P-valores analíticos"),
  data.frame(pvalor = resultadosBootstrap$pvalBoot, Tipo = "P-valores empíricos")
)

ggplot(df_pval, aes(x = pvalor)) + 
  geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#1A5F96", color = "white", linewidth = 0.3) + 
  facet_grid(~Tipo) + labs(x = "p-valor nominal", y = "Frecuencia") +
  theme(strip.background = element_blank(), panel.spacing = unit(2, "lines")) + configGraficas

#Número de genes rítmicos pre y post ajuste
sum(resultadosAnaliticos$pvalF<0.05)
#1225 (4.446622%)
sum(resultadosAnaliticos$pvalFAdjust<0.05)
#7 (0.02540927%)
sum(resultadosAnaliticos$pvalF<0.1)
#2350 (8.530255%)
sum(resultadosAnaliticos$pvalFAdjust<0.1)
#11 (0.02540927%)
sum(resultadosBootstrap$pvalBoot<0.05)
#1240 (4.52648%)
sum(resultadosBootstrap$pvalBoot<0.1)
#2368
sum(resultadosBootstrap$pvalBootstrapAdjust<0.05)
#0
sum(resultadosBootstrap$pvalBootstrapAdjust<0.1)
#0
rownames(resultadosAnaliticos[resultadosAnaliticos$pvalFAdjust < 0.05, ])
#"AL035252.2" "ARNTL"      "CIART"      "CRY1"       "KCNH4"      "PER3"       "SLC39A14"
genes_BH_directo<-rownames(resultadosAnaliticos[resultadosAnaliticos$pvalFAdjust < 0.05, ])
rownames(resultadosAnaliticos[resultadosAnaliticos$pvalFAdjust < 0.1, ])
#"AFAP1L1"    "AL035252.2" "ARNTL"      "CAMK1G"     "CIART"      "CRY1"       "KCNH4"      "LRRC57"  
#"PER2"       "PER3"       "SLC39A14"
genes_BH_Boot<-rownames(resultadosBootstrap[resultadosBootstrap$pvalBootstrapAdjust < 0.05, ])

################################################################################
# Propuesta metodológica: Mínimo de Bonferroni, Stouffer, Fisher y Truncated t #
################################################################################
numCores<-parallel::detectCores() - 1

t0Split<-Sys.time()
#SPLIT A/B CIRCULAR
sp <- make_split_k(tods, k = 3, screen_mod = 1)
idx_screen <- sp$A
idx_test <- sp$B

#Fold 1
## Suavizado en el screen, sin pesos w = NULL
metrics <- smooth_metrics_all_genes(matriz[,idx_screen],tods[idx_screen], kappa=8)
# Filtro genes: con "buena" SNR
knee <- snr_knee_filter(metrics)
snr_knee <- knee$snr_knee
keep_snr <- knee$keep
sum(keep_snr)
#Nos quedamos con 1697 genes
matriz_fold1<-matriz[keep_snr,idx_test]
tods1<-tods[idx_test]

#Fold 2
# divido en 2 al segundo split
idx_screen <- sp$B[seq(1,length(sp$B),by=2)]
idx_test <- sort(c(sp$B[seq(2,length(sp$B),by=2)],sp$A))

## Suavizado en el screen: mejor sin pesos w = NULL
metrics <- smooth_metrics_all_genes(matriz[,idx_screen],tods[idx_screen], kappa=8)
# Filtro genes: con "buena" SNR
knee <- snr_knee_filter(metrics)
snr_knee <- knee$snr_knee
keep_snr <- knee$keep
sum(keep_snr)
#Nos quedamos con 1973 genes
matriz_fold2<-matriz[keep_snr,idx_test]
tods2<-tods[idx_test]


#Fold 3
# divido en 2 al segundo split
idx_screen <- sp$B[seq(2,length(sp$B),by=2)]
idx_test <- sort(c(sp$B[seq(1,length(sp$B),by=2)],sp$A))
## Suavizado en el screen: sin pesos w = NULL
metrics <- smooth_metrics_all_genes(matriz[,idx_screen],tods[idx_screen], kappa=8)
# Filtro genes: con "buena" SNR
knee <- snr_knee_filter(metrics)
snr_knee <- knee$snr_knee
keep_snr <- knee$keep
sum(keep_snr)
#Nos quedamos con 2523 genes
matriz_fold3<-matriz[keep_snr,idx_test]
tods3<-tods[idx_test]
timeSplit <- as.numeric(difftime(Sys.time(), t0Split, units = "secs"))

#########################
#   Pvals analiticos    #
#########################
t0AnaliticoSplit<-Sys.time()
pvals_fold1_analiticos<-pvalsAnaliticos(matriz_fold1, tods1, tau=2*pi)
matPvals1Analitic<-data.frame(gene=rownames(pvals_fold1_analiticos), pvalF=pvals_fold1_analiticos$pvalF)
pvals_fold2_analiticos<-pvalsAnaliticos(matriz_fold2, tods2, tau=2*pi)
matPvals2Analitic<-data.frame(gene=rownames(pvals_fold2_analiticos), pvalF=pvals_fold2_analiticos$pvalF)
pvals_fold3_analiticos<-pvalsAnaliticos(matriz_fold3, tods3, tau=2*pi)
matPvals3Analitic<-data.frame(gene=rownames(pvals_fold3_analiticos), pvalF=pvals_fold3_analiticos$pvalF)

#JUNTAR P-VALORES
pvalsFoldsAnalitic<-merge(matPvals1Analitic, matPvals2Analitic, by="gene", all=TRUE)
pvalsFoldsAnalitic<-merge(pvalsFoldsAnalitic, matPvals3Analitic, by="gene",all=TRUE)
dim(pvalsFoldsAnalitic)
colnames(pvalsFoldsAnalitic)
timeAnaliticoSplit <- as.numeric(difftime(Sys.time(), t0AnaliticoSplit, units = "secs"))+timeSplit
print(timeAnaliticoSplit)#Time difference of 47.92982 secs

# Bonferroni: min(1, k*min(pj)), donde k es el numero de splits en el que ha sido seleccionado
p_comb_bonfAnalitic <- combine_pvalues(pvalsFoldsAnalitic[,2:4], method="minbonf")
sum(p_comb_bonfAnalitic<0.05)
#92
sum(p_comb_bonfAnalitic<0.1)
#189
q_comb_bonfAnalitic <- p.adjust(p_comb_bonfAnalitic, "BH")
sum(q_comb_bonfAnalitic<0.05)
#4
sum(q_comb_bonfAnalitic<0.1)
#7
pvalsFoldsAnalitic[q_comb_bonfAnalitic<0.1,"gene"]
#"CIART"    "CRY1"     "KCNH4"    "PDZRN3"   "PER2"     "PER3" "SLC39A14"
pvalsFoldsAnalitic[q_comb_bonfAnalitic<0.05,"gene"]
#"CIART" "CRY1"  "PER2"  "PER3" 
genes_Bonferroni<-pvalsFoldsAnalitic[q_comb_bonfAnalitic<0.05,"gene"]
hist(p_comb_bonfAnalitic, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Bonferroni", xlab = "p-valor")

# Fisher: −2*sum(log(pj))∼ chicuadrado 2*k gl
p_comb_fisherAnalitic <- combine_pvalues(pvalsFoldsAnalitic[,2:4], method="fisher")
sum(p_comb_fisherAnalitic<0.05)
#55
sum(p_comb_fisherAnalitic<0.1)
#85
q_comb_fisherAnalitic <- p.adjust(p_comb_fisherAnalitic, "BH")
sum(q_comb_fisherAnalitic<0.05)
#11
sum(q_comb_fisherAnalitic<0.1)
#11
pvalsFoldsAnalitic[q_comb_fisherAnalitic<0.1,"gene"]
#"AFAP1L1"    "AL035252.2" "CIART"      "CRY1"       "KCNH4"      "LRRC57"    
#"MTHFD1L"    "PDZRN3"     "PER2"       "PER3"       "SLC39A14"   
pvalsFoldsAnalitic[q_comb_fisherAnalitic<0.05,"gene"]
#"AFAP1L1"    "AL035252.2" "CIART"      "CRY1"       "KCNH4"      "LRRC57"    
#"MTHFD1L"    "PDZRN3"     "PER2"       "PER3"       "SLC39A14"
genes_Fisher<-pvalsFoldsAnalitic[q_comb_fisherAnalitic<0.05,"gene"]
hist(p_comb_fisherAnalitic, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Fisher", xlab = "p-valor")

# Stouffer / Z-score: convierte p a z y promedio (podria ser ponderado, si los splits fueran de distinto tamaño ponderar por sqrt(nk))
p_comb_zAnalitic <- combine_pvalues(pvalsFoldsAnalitic[,2:4], method="stouffer")
sum(p_comb_zAnalitic<0.05)
#14
sum(p_comb_zAnalitic<0.1)
#14
q_comb_zAnalitic <- p.adjust(p_comb_zAnalitic, "BH")
sum(q_comb_zAnalitic<0.05)
#6
sum(q_comb_zAnalitic<0.1)
#6
pvalsFoldsAnalitic[q_comb_zAnalitic<0.1,"gene"]
#"CIART"      "KCNH4"     "MTHFD1L"    "PDZRN3"     "PER2"       "SLC39A14"
pvalsFoldsAnalitic[q_comb_zAnalitic<0.05,"gene"]
#"CIART"      "KCNH4"     "MTHFD1L"    "PDZRN3"     "PER2"       "SLC39A14"
genes_Stouffer<-pvalsFoldsAnalitic[q_comb_zAnalitic<0.05,"gene"]
hist(p_comb_zAnalitic, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Stouffer", xlab = "p-valor")

#Truncated t
p_comb_truncated_tAnalitic <- combine_pvalues(pvalsFoldsAnalitic[,2:4], method="truncated_t", p0=0.9)
sum(p_comb_truncated_tAnalitic<0.05)
#87
sum(p_comb_truncated_tAnalitic<0.1)
#137
q_comb_truncated_tAnalitic <- p.adjust(p_comb_truncated_tAnalitic, "BH")
sum(q_comb_truncated_tAnalitic<0.05)
#6
sum(q_comb_truncated_tAnalitic<0.1)
#7
pvalsFoldsAnalitic[q_comb_truncated_tAnalitic<0.1,"gene"]
#"CIART"    "CRY1"     "KCNH4"    "PDZRN3"   "PER2"     "PER3"     "SLC39A14"
pvalsFoldsAnalitic[q_comb_truncated_tAnalitic<0.05,"gene"]
#"CIART"    "CRY1"    "KCNH4"  "PER2"     "PER3" "SLC39A14"
genes_Truncated<-pvalsFoldsAnalitic[q_comb_truncated_tAnalitic<0.05,"gene"]
hist(p_comb_truncated_tAnalitic, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Scaled t", xlab = "p-valor")

#######################
#   Pvals empiricos   #
#######################
t0BootSplit<-Sys.time()
cl<-parallel::makeCluster(numCores)
doParallel::registerDoParallel(cl)
pvals_fold1_bootstrap<-pvalsBootstrap(matriz_fold1, tods1, tau=2*pi, B=10000, n_chunks = nChunks)
pvals_fold2_bootstrap<-pvalsBootstrap(matriz_fold2, tods2, tau=2*pi, B=10000, n_chunks = nChunks)
pvals_fold3_bootstrap<-pvalsBootstrap(matriz_fold3, tods3, tau=2*pi, B=10000, n_chunks = nChunks)
parallel::stopCluster(cl)
matPvals1Bootstrap<-data.frame(gene=rownames(pvals_fold1_bootstrap), pvalBoot=pvals_fold1_bootstrap$pvalBoot)
matPvals2Bootstrap<-data.frame(gene=rownames(pvals_fold2_bootstrap), pvalBoot=pvals_fold2_bootstrap$pvalBoot)
matPvals3Bootstrap<-data.frame(gene=rownames(pvals_fold3_bootstrap), pvalBoot=pvals_fold3_bootstrap$pvalBoot)

pvalsFoldsBoot<-merge(matPvals1Bootstrap, matPvals2Bootstrap, by="gene", all=TRUE)
pvalsFoldsBoot<-merge(pvalsFoldsBoot, matPvals3Bootstrap, by="gene",all=TRUE)
dim(pvalsFoldsBoot)
colnames(pvalsFoldsBoot)
timeBootSplit <- as.numeric(difftime(Sys.time(), t0BootSplit, units = "secs"))+timeSplit
print(timeBootSplit)#Time difference of 56.91974 secs

# Bonferroni: min(1, k*min(pj)), donde k es el numero de splits en el que ha sido seleccionado
p_comb_bonfBoot <- combine_pvalues(pvalsFoldsBoot[,2:4], method="minbonf")
sum(p_comb_bonfBoot<0.05)
#95
sum(p_comb_bonfBoot<0.1)
#194
q_comb_bonfBoot <- p.adjust(p_comb_bonfBoot, "BH")
sum(q_comb_bonfBoot<0.05)
#0
sum(q_comb_bonfBoot<0.1)
#0
pvalsFoldsBoot[q_comb_bonfBoot<0.1,"gene"]
#0
genes_BonferroniBoot<-pvalsFoldsBoot[q_comb_bonfBoot<0.05,"gene"]
hist(p_comb_bonfBoot, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Bonferroni", xlab = "p-valor")

# Fisher: −2*sum(log(pj))∼ chicuadrado 2*k gl
p_comb_fisherBoot <- combine_pvalues(pvalsFoldsBoot[,2:4], method="fisher")
sum(p_comb_fisherBoot<0.05)
#57
sum(p_comb_fisherBoot<0.1)
#85
q_comb_fisherBoot <- p.adjust(p_comb_fisherBoot, "BH")
sum(q_comb_fisherBoot<0.05)
#11
sum(q_comb_fisherBoot<0.1)
#12
pvalsFoldsBoot[q_comb_fisherBoot<0.05,"gene"]
#"AFAP1L1"    "AL035252.2"  "CIART"     "CRY1"       "KCNH4"      "LRRC57"     "MTHFD1L"    "PDZRN3"
#"PER2"       "PER3"       "SLC39A14"
genes_FisherBoot<-pvalsFoldsBoot[q_comb_fisherBoot<0.05,"gene"]
pvalsFoldsBoot[q_comb_fisherBoot<0.1,"gene"]
#"AFAP1L1"    "AL035252.2"  "CIART"     "CRY1"       "KCNH4"      "LRRC57"     "MTHFD1L"    "PDZRN3"
#"PER2"       "PER3"       "SLC39A14" "WNT2"
hist(p_comb_fisherBoot, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Fisher", xlab = "p-valor")

# Stouffer / Z-score: convierte p a z y promedio (podria ser ponderado, si los splits fueran de distinto tamaño ponderar por sqrt(nk))
p_comb_zBoot <- combine_pvalues(pvalsFoldsBoot[,2:4], method="stouffer")
sum(p_comb_zBoot<0.05)
#14
sum(p_comb_zBoot<0.1)
#14
q_comb_zBoot <- p.adjust(p_comb_zBoot, "BH")
sum(q_comb_zBoot<0.05)
#6
sum(q_comb_zBoot<0.1)
#6
pvalsFoldsBoot[q_comb_zBoot<0.05,"gene"]
#"CIART"      "KCNH4"     "MTHFD1L"    "PDZRN3"     "PER2"       "SLC39A14"
genes_StoufferBoot<-pvalsFoldsBoot[q_comb_zBoot<0.05,"gene"]
pvalsFoldsBoot[q_comb_zBoot<0.1,"gene"]
#"CIART"      "KCNH4"     "MTHFD1L"    "PDZRN3"     "PER2"       "SLC22A9"
hist(p_comb_zBoot, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Stouffer", xlab = "p-valor")

#Truncated t
p_comb_truncated_tBoot <- combine_pvalues(pvalsFoldsBoot[,2:4], method="truncated_t", p0=0.9)
sum(p_comb_truncated_tBoot<0.05)
#84
sum(p_comb_truncated_tBoot<0.1)
#143
q_comb_truncated_tBoot <- p.adjust(p_comb_truncated_tBoot, "BH")
sum(q_comb_truncated_tBoot<0.05)
#0
sum(q_comb_truncated_tBoot<0.1)
#0
pvalsFoldsAnalitic[q_comb_truncated_tBoot<0.1,"gene"]
pvalsFoldsAnalitic[q_comb_truncated_tBoot<0.05,"gene"]
genes_Scaled_tBoot<-pvalsFoldsBoot[q_comb_truncated_tBoot<0.05,"gene"]
hist(p_comb_truncated_tBoot, breaks = seq(0, 1, by = 0.05), main = "Histograma de p-valores combinados Scaled t", xlab = "p-valor")


##################################################################################
#   Análisis de intersección de los genes detectados por los distintos métodos   #
##################################################################################
lista_combinada <- list(
  "BH (Analítico)" = genes_BH_directo,
  "Bonferroni (Analítico)" = genes_Bonferroni,
  "Fisher (Analítico)" = genes_Fisher,
  "Stouffer (Analítico)" = genes_Stouffer,
  "Truncated t (Analítico)" = genes_Truncated,
  "BH (Empírico)" = genes_BH_Boot,
  "Bonferroni (Empírico)" = genes_BonferroniBoot,
  "Fisher (Empírico)" = genes_FisherBoot,
  "Stouffer (Empírico)" = genes_StoufferBoot,
  "Truncated t (Empírico)" = genes_Scaled_tBoot
)
png(filename = "outputs/upsetAll.png", width = 5400, height = 4200, res = 600, units = "px")
upset(fromList(lista_combinada), 
      order.by = "degree", decreasing = c(TRUE, FALSE),
      sets = c(
        "Truncated t (Empírico)", "Stouffer (Empírico)", "Fisher (Empírico)", "Bonferroni (Empírico)", 
        "Truncated t (Analítico)", "Stouffer (Analítico)", "Fisher (Analítico)", "Bonferroni (Analítico)", 
        "BH (Empírico)", "BH (Analítico)"
      ),
      main.bar.color = "#1A5F96",sets.bar.color = "#333333", 
      matrix.color = "#1A5F96",shade.color = "#F0F0F0",    
      mainbar.y.label = "Número de genes compartidos",
      sets.x.label = "Total de genes por método",
      text.scale = c(1.8, 1.8, 1.5, 1.5, 1.7, 1.6),
      keep.order = TRUE,show.numbers = "yes" 
)
dev.off()

#8 mini graficos de genes destacados
genesInteres<-c("CLOCK", "ARNTL", "PER1", "PER2", "PER3", "CRY1", "CRY2", "CIART")
datosInteres <- as.data.frame(t(matriz[genesInteres, ]))
datosInteres$TOD <- tods
datosInteresLong <- datosInteres %>%
  pivot_longer(cols = -TOD, names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = genesInteres))
resultInteres <- resultadosAnaliticos[genesInteres, c("mesor", "amp", "acro")]
resultInteres$Gene <- rownames(resultInteres)

tiempo <- seq(0, 2*pi, length.out = 200)

datosCurva <- expand.grid(TOD = tiempo, Gene = genesInteres)
datosCurva <- merge(datosCurva, resultInteres, by = "Gene")

datosCurva <- datosCurva %>%
  mutate(Fitted = mesor + amp * cos(TOD + acro)) %>%
  mutate(Gene = factor(Gene, levels = genesInteres))

plotInteres<- ggplot() +
  geom_point(data = datosInteresLong, aes(x = TOD, y = Expression), 
             color = "#333333", alpha = 0.6, size = 1.5) +
  geom_line(data = datosCurva, aes(x = TOD, y = Fitted), 
            color = "#1A5F96", linewidth = 1.2) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4, axes = "all_x") +
  labs(x = "Tiempo Zeitgeber", y = "Nivel de expresión") +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    axis.text = element_text(size = 11.5),  
    axis.title = element_text(size = 14)    
  )
ggsave("outputs/genesInteres.pdf", plot = plotInteres, width = 8, height = 4)


#Tabla con: Snr, pvals analitic y empiric pre y post
tabla_snr<-datosInteresLong%>%
  left_join(resultInteres, by = c("Gene" = "Gene")) %>%
  mutate(Fitted = mesor + amp * cos(TOD + acro)) %>%
  mutate(Residuo = Expression - Fitted) %>%
  group_by(Gene) %>%
  summarize(
    Var_Ruido = var(Residuo),
    Var_Senal = 0.5 * (first(amp)^2),
    SNR = round(Var_Senal / Var_Ruido, 3),
    .groups = "drop"
  )

tablaCompleta <- tabla_snr %>%
  select(Gene, SNR) %>%
  mutate(
    p_Analit = resultadosAnaliticos[as.character(Gene), "pvalF"],
    q_Analit = resultadosAnaliticos[as.character(Gene), "pvalFAdjust"],
    p_Boot   = resultadosBootstrap[as.character(Gene), "pvalBoot"], 
    q_Boot   = resultadosBootstrap[as.character(Gene), "pvalBootstrapAdjust"]
  )

#Histogramas pvals
df_total <- bind_rows(
  # BH
  data.frame(pvalor = resultadosAnaliticos$pvalF, Tipo = "P-valores analíticos", Metodo = "BH_directo"),
  data.frame(pvalor = resultadosBootstrap$pvalBoot,  Tipo = "P-valores empíricos",  Metodo = "BH_directo"),
  
  # Bonferroni
  data.frame(pvalor = p_comb_bonfAnalitic, Tipo = "P-valores analíticos", Metodo = "Bonferroni"),
  data.frame(pvalor = p_comb_bonfBoot,  Tipo = "P-valores empíricos",  Metodo = "Bonferroni"),
  
  # Fisher
  data.frame(pvalor = p_comb_fisherAnalitic, Tipo = "P-valores analíticos", Metodo = "Fisher"),
  data.frame(pvalor = p_comb_fisherBoot,  Tipo = "P-valores empíricos",  Metodo = "Fisher"),
  
  # Stouffer
  data.frame(pvalor = p_comb_zAnalitic, Tipo = "P-valores analíticos", Metodo = "Stouffer"),
  data.frame(pvalor = p_comb_zBoot,  Tipo = "P-valores empíricos",  Metodo = "Stouffer"),
  
  # t trunc
  data.frame(pvalor = p_comb_truncated_tAnalitic, Tipo = "P-valores analíticos", Metodo = "t truncada"),
  data.frame(pvalor = p_comb_truncated_tBoot,  Tipo = "P-valores empíricos",  Metodo = "t truncada")
)


metodos <- c("BH_directo", "Bonferroni", "Fisher", "Stouffer", "t truncada")

lista_p1 <- lapply(metodos, function(m) {
  df_sub <- df_total %>% filter(Tipo == "P-valores analíticos", Metodo == m)
  
  p <- ggplot(df_sub, aes(x = pvalor)) + 
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#1A5F96", color = "white", linewidth = 0.3) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(clip = "on") +
    labs(x = NULL, y = if(m == "BH_directo") "Frecuencia" else NULL) +
    facet_grid(. ~ Metodo)
  
  if (m == "t truncada") {
    df_sub$TipoLabel <- "P-valores analíticos"
    p <- ggplot(df_sub, aes(x = pvalor)) + 
      geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#1A5F96", color = "white", linewidth = 0.3) + 
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(clip = "on") +
      labs(x = NULL, y = NULL) +
      facet_grid(TipoLabel ~ Metodo)
  }
  
  p <- p + theme(
    strip.background = element_rect(fill = "#F2F2F2", color = "#999999", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 11, color = "#222222"),
    panel.border = element_rect(color = "#999999", fill = NA, linewidth = 0.5),
    axis.text = element_text(size = 11.5),  
    axis.title = element_text(size = 12),
    plot.margin = margin(5, 5, 5, 5)
  )
  return(p)
})

lista_p2 <- lapply(metodos, function(m) {
  df_sub <- df_total %>% filter(Tipo == "P-valores empíricos", Metodo == m)
  
  p <- ggplot(df_sub, aes(x = pvalor)) + 
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#1A5F96", color = "white", linewidth = 0.3) + 
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(clip = "on") +
    labs(x = "p-valor nominal", y = if(m == "BH_directo") "Frecuencia" else NULL) +
    theme(strip.background = element_blank(), strip.text = element_blank())
  
  if (m == "t truncada") {
    df_sub$TipoLabel <- "P-valores empíricos"
    p <- ggplot(df_sub, aes(x = pvalor)) + 
      geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#1A5F96", color = "white", linewidth = 0.3) + 
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(clip = "on") +
      labs(x = "p-valor nominal", y = NULL) +
      facet_grid(TipoLabel ~ .) # Caja solo a la derecha
  }
  
  p <- p + theme(
    strip.background = element_rect(fill = "#F2F2F2", color = "#999999", linewidth = 0.5),
    strip.text.y = element_text(face = "bold", size = 11, color = "#222222"),
    panel.border = element_rect(color = "#999999", fill = NA, linewidth = 0.5),
    axis.text = element_text(size = 11.5),  
    axis.title = element_text(size = 12),
    plot.margin = margin(5, 5, 5, 5)
  )
  return(p)
})

fila1 <- wrap_plots(lista_p1, ncol = 5)
fila2 <- wrap_plots(lista_p2, ncol = 5)

composicionHist<-(fila1 / fila2) + 
  plot_layout(heights = c(1, 1)) & 
  theme(
    panel.spacing.x = unit(1.2, "lines"),
    panel.spacing.y = unit(1.5, "lines")
  )

ggsave("outputs/histPvals.pdf", plot = composicionHist, width = 12, height = 6)
