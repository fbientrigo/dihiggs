# File Index

0a0_get_from_lhe
- El profesor Nicolas generó previamente puntos fisicos, en base a una combinación de $m_\phi$ y $m_{12}^2$, est script extrae todos los parametros de input de interés para así obtener datos

0b3_classify_post
- Tras una generacion de los datos y tener los merged en data_batches, puedes usar este script para entrenar modelos de regresion basicos y comprobarlos con una curva ROC, entrenando para positividad, perturbatividad y unitariedad; mediante el uso de RandomForest se obtiene una aproximacion de las dependendicas de los parametros




# TODO

0523: Comprobar si variaciones de $m_{12}^2$ provocan cambios de $\lambda_j$
0523b: usar el calculo de los eigenvalues de las matrices para Tree Unitarity para obtener los eigenvalues, el mayor de los eigenvalues debera ser usado como una función de coste, debe de ser menor que $16 \pi$
