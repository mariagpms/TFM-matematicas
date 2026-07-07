# ==============================================================================
# Simulaciones con evaluación de TPR, FDR y tiempos de cómputo para diferentes 
# métodos de agregación y cálculo de p-valores (Analítico vs Empíricos).
# ==============================================================================

# 1. Cargar librerías, dependencias y establece tema común para gráficos
if (!require("pacman")) install.packages("pacman", repos = "https://cran.rstudio.com/")
pacman::p_load(tidyverse, data.table, ggplot2, patchwork, paletteer, parallel, doParallel, doRNG, progress)
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

# 2. Cargar funciones
source("funcionesCompletoGlutamatergic.R")
source("funcionesPvalores.R")

# 3. Configuración de escenarios y parámetros de simulación
G<-10000
escenarios <- as.data.table(expand.grid(
  n = c(48, 96, 192),
  proporRitmic = c(0.01, 0.05, 0.10),
  A = c(0.408, 0.500, 0.707)
))

numEscenarios<-nrow(escenarios)
numReplicas<-100
numMetodos<-5 #BH directo, Agregación Bonferroni, Agregación Fisher, Agregación Stouffer, Agregación Truncated t
alpha<-0.05

totalFilas<-numEscenarios*numReplicas*numMetodos

#Los resultados se guardan en una lista de data.tables, una para cada forma de calcular los pvalores
#Cada data.table tiene una fila por cada réplica de cada escenario, con columnas para identificar el escenario y TPR y FDR
resultados <- list(
  Analitico = data.table(
    Escenario_ID = character(totalFilas), n = numeric(totalFilas), proporRitmic = numeric(totalFilas),
    A = numeric(totalFilas), Metodo = character(totalFilas), Repeticion = integer(totalFilas),
    TPR = numeric(totalFilas), FDR = numeric(totalFilas), Tiempo = numeric(totalFilas)
  ),
  Boot = data.table(
    Escenario_ID = character(totalFilas), n = numeric(totalFilas), proporRitmic = numeric(totalFilas),
    A = numeric(totalFilas), Metodo = character(totalFilas), Repeticion = integer(totalFilas),
    TPR = numeric(totalFilas), FDR = numeric(totalFilas), Tiempo = numeric(totalFilas)
  )
)

#Función que genera los datos simulados:
#n: número de columans (número de muestras/pacientes)
#ritmicos: porcentaje de genes rítmicos
#A: amplitud de los genes rítmicos
#G: número total de genes (filas) evaluados
generarDatos <- function(n, ritmicos, A, G=10000, tau = 2*pi){
  tods <- runif(n, 0, 2*pi)
  numRitmic <- round(G * ritmicos)
  numNoRitmic <- G - numRitmic
  
  # Rítmicos (Optimizado con outer para evitar sapply)
  faseRitmic <- runif(numRitmic, 0, 2*pi)
  matRitmic <- A * cos(outer(faseRitmic, 2*pi*tods/tau, "-"))
  
  # No rítmicos
  matNoRitmic <- matrix(0, nrow=numNoRitmic, ncol=n)
  
  matriz <- rbind(matRitmic, matNoRitmic)
  ruido <- matrix(rnorm(G*n, mean=0, sd=0.5), nrow=G, ncol=n)
  
  # Asignar nombres fijos a los genes para poder indexar y cruzar los folds
  rownames(matriz) <- paste0("g", seq_len(G))
  verdad<-c(rep(1, numRitmic), rep(0, numNoRitmic))
  names(verdad) <- paste0("g", seq_len(G))
  return(list(
    matriz = matriz + ruido, 
    tods = tods, 
    verdad = verdad
  ))
}

# 4. Generación de escenarios y cálculo de resultados
#Se inicia el cluster para la paralelización del cálculo de los pvalroes empíricos
nChunks<-50
cl<-parallel::makeCluster(parallel::detectCores() - 1)
doParallel::registerDoParallel(cl)
doRNG::registerDoRNG()

contador<-1
for(i in 1:numEscenarios){
  esc_id<-paste0("Escenario_", i)
  escenario <- escenarios[i, ]
  cat("\nIniciando Escenario", i, "de", numEscenarios, ": n =", escenario$n, 
      ", ritmicos =", escenario$proporRitmic, ", A =", escenario$A, "\n")
  
  # Barra de progreso visual por escenario
  pb <- progress_bar$new(
    format = "  Progreso escenario :esc_num [:bar] :percent | ETA: :eta",
    total = numReplicas, clear = FALSE, width = 80
  )
  
  for(rep in 1:numReplicas){
    pb$tick(tokens = list(esc_num = i))
    
    #Generar datos de la réplica
    sim <- generarDatos(n=escenario$n, ritmicos=escenario$proporRitmic, A=escenario$A, tau=2*pi)
    matriz <- sim$matriz
    tods <- sim$tods
    
    #---------------------------------------------------------------------------
    # BH directo
    # --------------------------------------------------------------------------
    #Pvalores analíticos
    t0Analitico<-Sys.time()
    resAnalitico<-pvalsAnaliticos(sim$matriz, sim$tods, tau=2*pi)
    timeAnalitico<-as.numeric(difftime(Sys.time(), t0Analitico, units = "secs"))
    
    pvals_adj_analitico <- p.adjust(resAnalitico$pvalF, method = "BH")
    descubrimientos_analiticos <- sum(pvals_adj_analitico <= alpha)
    
    TPR_analitico <- sum(pvals_adj_analitico <= alpha & sim$verdad == 1) / sum(sim$verdad == 1)
    FDR_analitico <- ifelse(descubrimientos_analiticos > 0,
                            sum(pvals_adj_analitico <= alpha & sim$verdad == 0) / descubrimientos_analiticos,
                            0)
    
    set(resultados$Analitico, contador, 1:ncol(resultados$Analitico),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "BH_directo", rep, TPR_analitico, FDR_analitico, timeAnalitico))
    
    #Pvalores bootstrap
    t0Boot<-Sys.time()
    resBootstrap<-pvalsBootstrap(sim$matriz, sim$tods, tau=2*pi, B=10000, n_chunks = nChunks)
    timeBoot<-as.numeric(difftime(Sys.time(), t0Boot, units = "secs"))
    pvals_adj_boot <- p.adjust(resBootstrap$pvalBoot, method = "BH")
    descubrimientos_boot <- sum(pvals_adj_boot <= alpha)
    
    TPR_boot <- sum(pvals_adj_boot <= alpha & sim$verdad == 1) / sum(sim$verdad == 1)
    FDR_boot <- ifelse(descubrimientos_boot > 0,
                       sum(pvals_adj_boot <= alpha & sim$verdad == 0) / descubrimientos_boot,
                       0)
    set(resultados$Boot, contador, 1:ncol(resultados$Boot),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "BH_directo", rep, TPR_boot, FDR_boot, timeBoot))
    
    # --------------------------------------------------------------------------
    # Propuesta metodológica
    # --------------------------------------------------------------------------
    #Cálculo tiempos parte común
    t0Split <- Sys.time()
    sp <- make_split_k(tods, k = 3, screen_mod = 1)
    
    # --- FOLD 1 ---
    idx_screen_1 <- sp$A
    idx_LRT_1    <- sp$B
    metrics_1  <- smooth_metrics_all_genes(matriz[, idx_screen_1], tods[idx_screen_1], verbose = FALSE)
    keep_snr_1 <- snr_knee_filter(metrics_1)$keep
    matriz_filt_1 <- matriz[keep_snr_1, idx_LRT_1, drop=FALSE]
    
    # --- FOLD 2 ---
    idx_screen_2 <- sp$B[seq(1, length(sp$B), by = 2)]
    idx_LRT_2    <- sort(c(sp$B[seq(2, length(sp$B), by = 2)], sp$A))
    metrics_2  <- smooth_metrics_all_genes(matriz[, idx_screen_2], tods[idx_screen_2], verbose = FALSE)
    keep_snr_2 <- snr_knee_filter(metrics_2)$keep
    matriz_filt_2 <- matriz[keep_snr_2, idx_LRT_2, drop=FALSE]
    
    # --- FOLD 3 ---
    idx_screen_3 <- sp$B[seq(2, length(sp$B), by = 2)]
    idx_LRT_3    <- sort(c(sp$B[seq(1, length(sp$B), by = 2)], sp$A))
    metrics_3  <- smooth_metrics_all_genes(matriz[, idx_screen_3], tods[idx_screen_3], verbose = FALSE)
    keep_snr_3 <- snr_knee_filter(metrics_3)$keep
    matriz_filt_3 <- matriz[keep_snr_3, idx_LRT_3, drop=FALSE]
    
    timeSplit <- as.numeric(difftime(Sys.time(), t0Split, units = "secs"))
    
    #Propuesta analíticos
    t0SplitAnalitic<- Sys.time()
    res_analitic_1 <- pvalsAnaliticos(matriz_filt_1, tods[idx_LRT_1], tau=2*pi)
    mat_fold1_analitic <- data.frame(gene = rownames(matriz_filt_1), p_analitic = res_analitic_1$pvalF, stringsAsFactors = FALSE)
    
    res_analitic_2 <- pvalsAnaliticos(matriz_filt_2, tods[idx_LRT_2], tau=2*pi)
    mat_fold2_analitic <- data.frame(gene = rownames(matriz_filt_2), p_analitic = res_analitic_2$pvalF, stringsAsFactors = FALSE)
    
    res_analitic_3 <- pvalsAnaliticos(matriz_filt_3, tods[idx_LRT_3], tau=2*pi)
    mat_fold3_analitic <- data.frame(gene = rownames(matriz_filt_3), p_analitic = res_analitic_3$pvalF, stringsAsFactors = FALSE)
    
    p_analitic_folds <- merge(mat_fold1_analitic, mat_fold2_analitic, by="gene", all=TRUE)
    p_analitic_folds <- merge(p_analitic_folds, mat_fold3_analitic, by="gene", all=TRUE)
    p_mat_analitic   <- as.matrix(p_analitic_folds[, 2:4])
    rownames(p_mat_analitic) <- p_analitic_folds$gene
    
    timeSplitAnalitic<- as.numeric(difftime(Sys.time(), t0SplitAnalitic, units = "secs"))
    
    #Propuesta bootstrap
    t0SplitBoot<- Sys.time()
    res_boot_1 <- pvalsBootstrap(matriz_filt_1, tods[idx_LRT_1], tau=2*pi, B=10000, n_chunks = nChunks)
    mat_fold1_boot <- data.frame(gene = rownames(matriz_filt_1), p_boot = res_boot_1$pvalBoot, stringsAsFactors = FALSE)
    
    res_boot_2 <- pvalsBootstrap(matriz_filt_2, tods[idx_LRT_2], tau=2*pi, B=10000, n_chunks = nChunks)
    mat_fold2_boot <- data.frame(gene = rownames(matriz_filt_2), p_boot = res_boot_2$pvalBoot, stringsAsFactors = FALSE)
    
    res_boot_3 <- pvalsBootstrap(matriz_filt_3, tods[idx_LRT_3], tau=2*pi, B=10000, n_chunks = nChunks)
    mat_fold3_boot <- data.frame(gene = rownames(matriz_filt_3), p_boot = res_boot_3$pvalBoot, stringsAsFactors = FALSE)
    
    p_boot_folds <- merge(mat_fold1_boot, mat_fold2_boot, by="gene", all=TRUE)
    p_boot_folds <- merge(p_boot_folds, mat_fold3_boot, by="gene", all=TRUE)
    p_mat_boot   <- as.matrix(p_boot_folds[, 2:4])
    rownames(p_mat_boot) <- p_boot_folds$gene
    
    timeSplitBoot<- as.numeric(difftime(Sys.time(), t0SplitBoot, units = "secs"))
    
    #Combinaciones Analítico
    t0Fisher <- Sys.time()
    q_analitic_fisher   <- p.adjust(combine_pvalues(p_mat_analitic, method="fisher"), method="BH")
    timeFisher <- as.numeric(difftime(Sys.time(), t0Fisher, units = "secs"))+timeSplitAnalitic+timeSplit
    t0Stouffer <- Sys.time()
    q_analitic_stouffer <- p.adjust(combine_pvalues(p_mat_analitic, method="stouffer"), method="BH")
    timeStouffer <- as.numeric(difftime(Sys.time(), t0Stouffer, units = "secs"))+timeSplitAnalitic+timeSplit
    t0Bonf <- Sys.time()
    q_analitic_bonf     <- p.adjust(combine_pvalues(p_mat_analitic, method="minbonf"), method="BH")
    timeBonf <- as.numeric(difftime(Sys.time(), t0Bonf, units = "secs"))+timeSplitAnalitic+timeSplit
    t0Trunc <- Sys.time()
    q_analitic_trunc    <- p.adjust(combine_pvalues(p_mat_analitic, method="truncated_t"), method="BH")
    timeTrunc <- as.numeric(difftime(Sys.time(), t0Trunc, units = "secs"))+timeSplitAnalitic+timeSplit
    
    #Combinaciones Bootstrap
    t0FisherBoot <- Sys.time()
    q_boot_fisher   <- p.adjust(combine_pvalues(p_mat_boot, method="fisher"), method="BH")
    timeFisherBoot <- as.numeric(difftime(Sys.time(), t0FisherBoot, units = "secs"))+timeSplitBoot+timeSplit
    t0StoufferBoot <- Sys.time()
    q_boot_stouffer <- p.adjust(combine_pvalues(p_mat_boot, method="stouffer"), method="BH")
    timeStoufferBoot <- as.numeric(difftime(Sys.time(), t0StoufferBoot, units = "secs"))+timeSplitBoot+timeSplit
    t0BonfBoot <- Sys.time()
    q_boot_bonf     <- p.adjust(combine_pvalues(p_mat_boot, method="minbonf"), method="BH")
    timeBonfBoot <- as.numeric(difftime(Sys.time(), t0BonfBoot, units = "secs"))+timeSplitBoot+timeSplit
    t0TruncBoot <- Sys.time()
    q_boot_trunc    <- p.adjust(combine_pvalues(p_mat_boot, method="truncated_t"), method="BH")
    timeTruncBoot <- as.numeric(difftime(Sys.time(), t0TruncBoot, units = "secs"))+timeSplitBoot+timeSplit
    
    # --------------------------------------------------------------------------
    # Cálculo métricas
    # --------------------------------------------------------------------------
    totPositivos <- sum(sim$verdad == 1)
    #Split Analítico + Fisher
    descub_af<-sum(q_analitic_fisher <= alpha, na.rm = TRUE)
    genes_af<-rownames(p_mat_analitic)[which(q_analitic_fisher <= alpha)]
    tp_af<-sum(sim$verdad[genes_af] == 1, na.rm = TRUE)
    fp_af<-sum(sim$verdad[genes_af] == 0, na.rm = TRUE)
    tpr_af<-ifelse(totPositivos > 0, tp_af / totPositivos, 0)
    fdr_af<-ifelse(descub_af > 0, fp_af / descub_af, 0)
    set(resultados$Analitico, contador+1, 1:ncol(resultados$Analitico),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Fisher", rep, tpr_af, fdr_af, timeFisher))
    
    #Split Analítico + Stouffer
    descub_as<- sum(q_analitic_stouffer <= alpha, na.rm = TRUE)    
    genes_as<-rownames(p_mat_analitic)[which(q_analitic_stouffer <= alpha)]
    tp_as<-sum(sim$verdad[genes_as] == 1, na.rm = TRUE)
    fp_as<-sum(sim$verdad[genes_as] == 0, na.rm = TRUE)
    tpr_as<-ifelse(totPositivos > 0, tp_as / totPositivos, 0)
    fdr_as<-ifelse(descub_as > 0, fp_as / descub_as, 0)
    set(resultados$Analitico, contador+2, 1:ncol(resultados$Analitico),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Stouffer", rep, tpr_as, fdr_as, timeStouffer))
    
    #Split Analitico + Bonferroni
    descub_ab <- sum(q_analitic_bonf <= alpha, na.rm = TRUE)
    genes_ab<-rownames(p_mat_analitic)[which(q_analitic_bonf <= alpha)]
    tp_ab<-sum(sim$verdad[genes_ab] == 1, na.rm = TRUE)
    fp_ab<-sum(sim$verdad[genes_ab] == 0, na.rm = TRUE)
    tpr_ab<-ifelse(totPositivos > 0, tp_ab / totPositivos, 0)
    fdr_ab<-ifelse(descub_ab > 0, fp_ab / descub_ab, 0)
    set(resultados$Analitico, contador+3, 1:ncol(resultados$Analitico),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Bonferroni", rep, tpr_ab, fdr_ab, timeBonf))
    
    #Split Analítico + Truncated t
    descub_at <- sum(q_analitic_trunc <= alpha, na.rm = TRUE)
    genes_at<-rownames(p_mat_analitic)[which(q_analitic_trunc <= alpha)]
    tp_at<-sum(sim$verdad[genes_at] == 1, na.rm = TRUE)
    fp_at<-sum(sim$verdad[genes_at] == 0, na.rm = TRUE)
    tpr_at<-ifelse(totPositivos > 0, tp_at / totPositivos, 0)
    fdr_at<-ifelse(descub_at > 0, fp_at / descub_at, 0)
    set(resultados$Analitico, contador+4, 1:ncol(resultados$Analitico),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Truncated_t", rep, tpr_at, fdr_at, timeTrunc))
    
    # Split Bootstrap + Fisher
    descub_bf <- sum(q_boot_fisher <= alpha, na.rm = TRUE)
    genes_bf<-rownames(p_mat_boot)[which(q_boot_fisher <= alpha)]
    tp_bf<-sum(sim$verdad[genes_bf] == 1, na.rm = TRUE)
    fp_bf<-sum(sim$verdad[genes_bf] == 0, na.rm = TRUE)
    tpr_bf<-ifelse(totPositivos > 0, tp_bf / totPositivos, 0)
    fdr_bf<-ifelse(descub_bf > 0, fp_bf / descub_bf, 0)
    set(resultados$Boot, contador+1, 1:ncol(resultados$Boot),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Fisher", rep, tpr_bf, fdr_bf, timeFisherBoot))
    
    #Split Bootstrap + Stouffer
    descub_bs <- sum(q_boot_stouffer <= alpha, na.rm = TRUE)
    genes_bs<-rownames(p_mat_boot)[which(q_boot_stouffer <= alpha)]
    tp_bs<-sum(sim$verdad[genes_bs] == 1, na.rm = TRUE)
    fp_bs<-sum(sim$verdad[genes_bs] == 0, na.rm = TRUE)
    tpr_bs<-ifelse(totPositivos > 0, tp_bs / totPositivos, 0)
    fdr_bs<-ifelse(descub_bs > 0, fp_bs / descub_bs, 0)
    set(resultados$Boot, contador+2, 1:ncol(resultados$Boot),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Stouffer", rep, tpr_bs, fdr_bs, timeStoufferBoot))
    
    #Split Bootstrap + Bonferroni
    descub_bb <- sum(q_boot_bonf <= alpha, na.rm = TRUE)
    genes_bb<-rownames(p_mat_boot)[which(q_boot_bonf <= alpha)]
    tp_bb<-sum(sim$verdad[genes_bb] == 1, na.rm = TRUE)
    fp_bb<-sum(sim$verdad[genes_bb] == 0, na.rm = TRUE)
    tpr_bb<-ifelse(totPositivos > 0, tp_bb / totPositivos, 0)
    fdr_bb<-ifelse(descub_bb > 0, fp_bb / descub_bb, 0)
    set(resultados$Boot, contador+3, 1:ncol(resultados$Boot),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Bonferroni", rep, tpr_bb, fdr_bb, timeBonfBoot))
    
    #Split Bootstrap + Truncated t
    descub_bt <- sum(q_boot_trunc <= alpha, na.rm = TRUE)
    genes_bt<-rownames(p_mat_boot)[which(q_boot_trunc <= alpha)]
    tp_bt<-sum(sim$verdad[genes_bt] == 1, na.rm = TRUE)
    fp_bt<-sum(sim$verdad[genes_bt] == 0, na.rm = TRUE)
    tpr_bt<-ifelse(totPositivos > 0, tp_bt / totPositivos, 0)
    fdr_bt<-ifelse(descub_bt > 0, fp_bt / descub_bt, 0)
    set(resultados$Boot, contador+4, 1:ncol(resultados$Boot),
        list(esc_id, escenario$n, escenario$proporRitmic, escenario$A, "Truncated_t", rep, tpr_bt, fdr_bt, timeTruncBoot))
    
    contador <- contador + numMetodos
  }
  
}
parallel::stopCluster(cl)
save.image("simulaciones.RData")

resultados$Analitico <- resultados$Analitico %>%
  mutate(SNR = case_when(
    A == 0.408 ~ "1/3",
    A == 0.5   ~ "1/2",
    A == 0.707 ~ "1",
    TRUE       ~ as.character(A)
  )) %>%
  mutate(SNR = factor(SNR, levels = c("1/3", "1/2", "1")))
resultados$Analitico[Metodo == "Truncated_t", Metodo := "t truncada"]

resultados$Boot <- resultados$Boot %>%
  mutate(SNR = case_when(
    A == 0.408 ~ "1/3",
    A == 0.5   ~ "1/2",
    A == 0.707 ~ "1",
    TRUE       ~ as.character(A)
  )) %>%
  mutate(SNR = factor(SNR, levels = c("1/3", "1/2", "1")))
resultados$Boot[Metodo == "Truncated_t", Metodo := "t truncada"]

# 5. Visualización 
#Gráfico 1: TPR Analítico (Matriz 3x3)
datos_graficoTPR <- resultados$Analitico %>%
  group_by(n, SNR, proporRitmic, Metodo) %>%
  summarise(Mean_TPR = mean(TPR, na.rm = TRUE), SD_TPR = sd(TPR, na.rm = TRUE), .groups = "drop")

p1<-ggplot(datos_graficoTPR, aes(x = as.factor(proporRitmic), y = Mean_TPR, fill = Metodo)) +
  geom_bar(stat = "identity", alpha=0.8, width=0.85, position = position_dodge(width = 0.85), color = "#333333", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, Mean_TPR - SD_TPR), ymax = Mean_TPR + SD_TPR), 
                position = position_dodge(width = 0.85), width = 0.3, color = "#222222", linewidth = 0.5) +
  facet_grid(n ~ SNR, labeller = label_both) + 
  labs(x = "Porcentaje rítmicos", y = "TPR", fill = "Método") +
  configGraficas
ggsave("outputs/tprAnalitico.pdf", plot = p1, width = 6, height = 4)

#Gráfico 2: FDR Analítico (Matriz 3x3 con Línea en 0.05)
datos_graficoFDR <- resultados$Analitico %>%
  group_by(n, SNR, proporRitmic, Metodo) %>%
  summarise(Mean_FDR = mean(FDR, na.rm = TRUE), SD_FDR = sd(FDR, na.rm = TRUE), .groups = "drop")

p2<-ggplot(datos_graficoFDR, aes(x = as.factor(proporRitmic), y = Mean_FDR, fill = Metodo)) +
  geom_bar(stat = "identity", alpha=0.8, width=0.85, position = position_dodge(width = 0.85), color = "#333333", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, Mean_FDR - SD_FDR), ymax = Mean_FDR + SD_FDR), 
                position = position_dodge(width = 0.85), width = 0.2, color = "#222222", linewidth = 0.4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
  facet_grid(n ~ SNR, labeller = label_both) + 
  labs(x = "Porcentaje rítmicos", y = "FDR", fill = "Método") +
  configGraficas
ggsave("outputs/fdrAnalitico.pdf", plot = p2, width = 6, height = 4)

#Gráfico 3: Tiempo Analítico (Matriz 3x3)
datos_graficoTime <- resultados$Analitico %>%
  group_by(n, SNR, proporRitmic, Metodo) %>%
  summarise(Mean_Time = mean(Tiempo, na.rm = TRUE), SD_Time = sd(Tiempo, na.rm = TRUE), .groups = "drop")

p3<-ggplot(datos_graficoTime, aes(x = as.factor(proporRitmic), y = Mean_Time, fill = Metodo)) +
  geom_bar(stat = "identity", alpha=0.8, width=0.85, position = position_dodge(width = 0.85), color = "#333333", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, Mean_Time - SD_Time), ymax = Mean_Time + SD_Time), 
                position = position_dodge(width = 0.85), width = 0.2, color = "#222222", linewidth = 0.4) +
  facet_grid(n ~ SNR, labeller = label_both) + 
  labs(x = "Porcentaje rítmicos", y = "Tiempo (s)", fill = "Método") +
  configGraficas
ggsave("outputs/timeAnalitico.pdf", plot = p3, width = 6, height = 4)

#Gráfico 4: TPR Bootstrap (Matriz 3x3)
datos_graficoTPRBoot <- resultados$Boot %>%
  group_by(n, SNR, proporRitmic, Metodo) %>%
  summarise(Mean_TPR = mean(TPR, na.rm = TRUE), SD_TPR = sd(TPR, na.rm = TRUE), .groups = "drop")

p4<-ggplot(datos_graficoTPRBoot, aes(x = as.factor(proporRitmic), y = Mean_TPR, fill = Metodo)) +
  geom_bar(stat = "identity", alpha=0.8, width=0.85, position = position_dodge(width = 0.85), color = "#333333", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, Mean_TPR - SD_TPR), ymax = Mean_TPR + SD_TPR), 
                position = position_dodge(width = 0.85), width = 0.2, color = "#222222", linewidth = 0.4) +
  facet_grid(n ~ SNR, labeller = label_both) + 
  labs(x = "Porcentaje rítmicos", y = "TPR", fill = "Método") +
  configGraficas
ggsave("outputs/tprBoot.pdf", plot = p4, width = 6, height = 4)

#Gráfico 5: FDR Bootstrap (Matriz 3x3 con Línea en 0.05)
datos_graficoFDRBoot <- resultados$Boot %>%
  group_by(n, SNR, proporRitmic, Metodo) %>%
  summarise(Mean_FDR = mean(FDR, na.rm = TRUE), SD_FDR = sd(FDR, na.rm = TRUE), .groups = "drop")

p5<-ggplot(datos_graficoFDRBoot, aes(x = as.factor(proporRitmic), y = Mean_FDR, fill = Metodo)) +
  geom_bar(stat = "identity", alpha=0.8, width=0.85, position = position_dodge(width = 0.85), color = "#333333", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, Mean_FDR - SD_FDR), ymax = Mean_FDR + SD_FDR), 
                position = position_dodge(width = 0.85), width = 0.2, color = "#222222", linewidth = 0.4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.6) +
  facet_grid(n ~ SNR, labeller = label_both) + 
  labs(x = "Porcentaje rítmicos", y = "FDR", fill = "Método") +
  configGraficas
ggsave("outputs/fdrBoot.pdf", plot = p5, width = 6, height = 4)

#Gráfico 6: Tiempo Bootstrap (Matriz 3x3)
datos_graficoTimeBoot <- resultados$Boot %>%
  group_by(n, SNR, proporRitmic, Metodo) %>%
  summarise(Mean_Time = mean(Tiempo, na.rm = TRUE), SD_Time = sd(Tiempo, na.rm = TRUE), .groups = "drop")

p6<-ggplot(datos_graficoTimeBoot, aes(x = as.factor(proporRitmic), y = Mean_Time, fill = Metodo)) +
  geom_bar(stat = "identity", alpha=0.8, width=0.85, position = position_dodge(width = 0.85), color = "#333333", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, Mean_Time - SD_Time), ymax = Mean_Time + SD_Time), 
                position = position_dodge(width = 0.85), width = 0.2, color = "#222222", linewidth = 0.4) +
  facet_grid(n ~ SNR, labeller = label_both) + 
  labs(x = "Porcentaje rítmicos", y = "Tiempo (s)", fill = "Método") +
  configGraficas
ggsave("outputs/timeBoot.pdf", plot = p6, width = 6, height = 4)
