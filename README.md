# Corrección por comparaciones múltiples en estudios ómicos: evaluación metodológica y aplicaciones a datos de expresión génica.
Este repositorio de GitHub contiene el código empleado para los análisis desarrollados en el Trabajo Fin de Máster.

**Estructura del repositorio**
* [`funcionesPvalores.R`](./funcionesPvalores.R) y [`funcionesGluta.R`](./funcionesGluta.R): Contienen las funciones auxiliares necesarias en el resto de scripts.
* * analisisNula.R: Evaluación de la distribución nula y robustez del estadístico F bajo diferentes condiciones (homocedasticidad, heterocedastiticidad y valores atípicos).
* [`simulaciones.R`](./simulaciones.R): Script encargado de generar escenarios simulados. Evalúa la metodología propuesta mediante métricas de rendimiento (TPR, FDR) y tiempos de cómputo comparando diferentes métodos de agregación y cálculo de _p_-valores (Analítico vs Empíricos).
* [`analisisGluta.R`](./analisisGluta.R): Script para los análisis realizados sobre los datos de expresión de BA46. Incluye el cálculo de _p_-valores analíticos y empíricos, la aplicación de los métodos de agregación (Fisher, Stouffer, Bonferroni, t truncada) y la generación de figuras.
  * _Nota_: Para realizar el análisis de sensibilidad, se debe ajustar el parámetro `kappa` en la función `smooth_metrics_all_genes()`.
