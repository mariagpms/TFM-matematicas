pvalsAnaliticos<-function(matriz, tods, tau=2*pi){
  N<-ncol(matriz)
  G<-nrow(matriz)
  
  #Predictores
  xVal<-cos(2*pi*tods/tau)
  zVal<-sin(2*pi*tods/tau)
  Y<-cbind(mesor=1, beta=xVal, gamma=zVal)
  
  #Params lineales
  YtYinvYt<-solve(t(Y)%*%Y)%*%t(Y)
  coeffs<-matriz%*%t(YtYinvYt)
  colnames(coeffs)<-c("mesor", "beta", "gamma")
  
  #Params no lineales
  ampHat<-sqrt(coeffs[, "beta"]^2 + coeffs[, "gamma"]^2)
  phiHat<-atan2(-coeffs[, "gamma"], coeffs[, "beta"])
  
  #F y R2
  YHat<-coeffs%*%t(Y)
  RSS1<-rowSums((matriz-YHat)^2)
  RSS0<-rowSums((matriz-rowMeans(matriz))^2)
  SSR<-RSS0-RSS1
  
  Fobs<-(SSR/2)/(RSS1/(N-3))
  pvalF<-pf(Fobs, 2, N-3, lower.tail = FALSE)
  
  return(data.frame(
      mesor = coeffs[, 1],
      amp = ampHat, 
      acro = phiHat,
      pvalF = pvalF, 
      #r2 = r2Obs, 
      #lrt=lrtObs, 
      Fobs = Fobs
    ))
}


pvalsBootstrap<-function(matriz, tods, tau=2*pi, B=10000, n_chunks = NULL){
  G<-nrow(matriz)
  N<-ncol(matriz)
  
  # Centrar la matriz para simplificar cálculo SSR nulo
  matrizC<-t(scale(t(matriz), scale = FALSE))
  
  # Predictores
  xVal<-cos(2*pi*tods/tau)
  zVal<-sin(2*pi*tods/tau)
  Yc<-scale(cbind(xVal, zVal), scale = FALSE)
  Hc<-solve(t(Yc)%*%Yc)%*%t(Yc)
  YtY<-t(Yc)%*%Yc
  
  # SSR_Observado
  coefsObs<-matrizC%*%t(Hc) 
  SSRobs<-YtY[1,1] * coefsObs[,1]^2 +
    (YtY[1,2] + YtY[2,1]) * coefsObs[,1] * coefsObs[,2] +
    YtY[2,2] * coefsObs[,2]^2
  
  # Params
  mesorHat<-rowMeans(matriz)
  betaHat<-coefsObs[, 1]
  gammaHat<-coefsObs[, 2]
  
  RSS0<-rowSums(matrizC^2)
  RSS1<-RSS0-SSRobs
  Fobs<-(SSRobs/2)/(RSS1/(N-3))
  
  #Si no se especifican bloques, se usa el número de núcleos registrados menos 1
  if(is.null(n_chunks)){
    n_chunks <- foreach::getDoParWorkers()
  }
  
  #Se reparten las B iteraciones de la forma más equitativa posible entre los bloques
  base_b<-B %/% n_chunks
  resto<-B %% n_chunks
  b_chunk <- rep(base_b, n_chunks)
  if(resto > 0) b_chunk[1:resto] <- b_chunk[1:resto] + 1
  
  #Paralelización por bloques
  `%dorng%` <- doRNG::`%dorng%`
  
  #Se itera sobre los bloques creados consiguiendo reducir el overhead
  pCountsTotal <- foreach::foreach(ch = 1:n_chunks, .combine = "+", .inorder = FALSE) %dorng% {
    
    # Inicializamos un vector acumulador local en el nodo para evitar transferencias continuas
    pLocal<-integer(G)
    b_actual<-b_chunk[ch]
    
    for(k in 1:b_actual){
      # Permutación de los tiempos
      idx <- sample.int(N) 
      Yperm <- matrizC[, idx]
      
      # Estimación de coeficientes nulos
      coefsNull <- Yperm %*% t(Hc)
      
      # Cálculo rápido de SSR nulo
      ssrNull <- YtY[1,1] * coefsNull[,1]^2 +
        (YtY[1,2] + YtY[2,1]) * coefsNull[,1] * coefsNull[,2] +
        YtY[2,2] * coefsNull[,2]^2
      
      RSS1Null<-RSS0-ssrNull
      FNull<-(ssrNull/2)/(RSS1Null/(N-3))
      pLocal<-pLocal + (FNull >= Fobs)
    }
    
    pLocal
  }
  
  return(data.frame(
    mesor=mesorHat,
    amp=sqrt(betaHat^2 + gammaHat^2),
    acro=atan2(-gammaHat, betaHat),
    Fobs=Fobs,
    pvalBoot=(pCountsTotal+1)/(B+1)
  ))
}

# pmat es un data.frame con columnas pvalue_split1, pvalue_split2, pvalue_split3 donde NA = “no testeado”.
combine_pvalues <- function(pmat, method = c("minbonf", "tippett", "fisher", "stouffer", "cauchy", "truncated_t"),
                            stouffer_w = NULL, p0 = 0.9) { 
  method <- match.arg(method)
  m <- nrow(pmat)
  out <- rep(NA_real_, m)
  
  for (i in seq_len(m)) {
    pvec <- as.numeric(pmat[i, ])
    pvec[!is.finite(pvec)] <- 1
    if (length(pvec) == 0) { out[i] <- 1; next }
    
    k <- length(pvec)
    
    if (method == "minbonf") {
      out[i] <- min(1, k * min(pvec))
    } else if (method == "tippett") {
      out[i] <- 1 - (1 - min(pvec))^k
    } else if (method == "fisher") {
      pvec2 <- pmax(pvec, 1e-300)
      stat <- -2 * sum(log(pvec2))
      out[i] <- pchisq(stat, df = 2 * k, lower.tail = FALSE)
    } else if (method == "stouffer") {
      pvec2 <- pmin(pmax(pvec, 1e-15), 1 - 1e-15)
      z <- qnorm(1 - pvec2)
      if (is.null(stouffer_w)) {
        w <- rep(1, k)
      } else {
        wfull <- as.numeric(stouffer_w[i, ])
        w <- wfull[is.finite(pmat[i, ])]
      }
      Z <- sum(w * z) / sqrt(sum(w^2))
      out[i] <- 1 - pnorm(Z)
    } else if (method == "cauchy") {
      pvec2 <- pmin(pmax(pvec, 1e-15), 1 - 1e-15)
      stat <- sum(tan((0.5 - pvec2) * pi)) / k
      out[i] <- 0.5 - atan(stat) / pi
    } else if (method == "truncated_t") {
      S<-sum(qt(1-p0*pvec, df=1))
      out[i] <- min(1, (k/p0)*pt(S, df = 1, lower.tail = FALSE)) 
    } 
  }
  out
}
