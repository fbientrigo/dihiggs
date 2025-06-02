
# Estudio paramétrico con Machine Learning aplicado al análisis de modelos 2HDM

**Descripción breve:**  
Este repositorio contiene notebooks y scripts en Python para ejecutar simulaciones con el `oracle` en C++, procesar los resultados (parámetros y variables físicas), y entrenar modelos de ML que evalúan la viabilidad de puntos en el espacio paramétrico.  

---

## Índice  
1. [Convenciones](#convenciones)  
2. [Notebooks](#notebooks)  
3. [Estructura de proyecto](#estructura-de-proyecto)  
4. [Referencias](#referencias)

---

## Convenciones  
A continuación se muestra el formato de salida típico al ejecutar el **oracle** con multithreading. Esto ayuda a entender cómo están organizados los datos para entrenar los modelos de ML.

```python
[  
    {  
        'params': array([  
            [1.50000001e+02, 2.99932979e+02, 1.00000000e+00, …, 1.00000000e-01, 0.00000000e+00, 2.24990902e+00],  
            [1.50000001e+02, 3.00088086e+02, 1.00000000e+00, …, … ],  
            [1.50000001e+02, 3.00033035e+02, 1.00000000e+00, …, 1.00000000e-01, 0.00000000e+00, 2.24990903e+00]  
        ]),  
        'results': [  
            {  
                'positivity_ok': 1,  
                'unitarity_ok': 1,  
                'perturbativity_ok': 1,  
                'w_h2_bb': 0.0,  
                'w_h2_tautau': 0.0,  
                'w_h2_uu': 0.0,  
                'w_h2_du': 0.0,  
                'w_h2_ln': 0.0,  
                'w_h2_vv': [5.053510327013813e-11, 0.0, 0.0],  
                'w_h2_gaga': 5.053510327013813e-11,  
                'w_h2_Zga': 1.181312528477372e-11,  
                'w_h2_gg': 0.0,  
                'w_h2_hh': 0.0,  
                'w_total_h2': 6.234822855491185e-11,  
                'w_total_top': 1.338069402459526,  
                'branching_ratio_h2_gaga': 0.8105298970223419,  
                'lambda1': 1.129741358125329,  
                'lambda2': 0.2577337926166231,  
                'lambda3': 2.483237305682109,  
                'lambda4': -1.112766760892846,  
                'lambda5': -1.112766760892846,  
                'lambda6': 0.1,  
                'lambda7': 0.0  
            },  
            {  
                'positivity_ok': 1,  
                'unitarity_ok': 1,  
                'perturbativity_ok': 1,  
                … (otros campos similares)  
            }  
        ]  
    }  
]
````

* **`params`**: arreglo NumPy con los siguientes parámetros físicos (en este orden):

  1. $m_{\phi}$
  2. $m_{A}$
  3. $\sin(\beta - \alpha)$
  4. $\tan\beta$
  5. $\lambda_{6}$
  6. $\lambda_{7}$
  7. $m_{12}^2$

* **`results`**: lista de diccionarios, uno por cada punto calculado. Cada diccionario contiene valores como `positivity_ok`, `unitarity_ok`, `w_h2_bb`, `w_h2_gaga`, etc.

  > **Nota:** Completar la documentación de cada campo en `results` (por ejemplo, `w_h2_bb` = ancho de decaimiento $h_2 \to b\bar b$, `w_total_h2` = ancho total de $h_2$, etc.)

---

## Notebooks

En la carpeta principal (`mlpython/`) se incluyen los siguientes notebooks para distintas etapas del flujo de trabajo:

* **`0a0_get_from_lhe.ipynb`**
  Extrae puntos válidos almacenados en `valid_points_lhe` y prepara los datos iniciales.

* **`0b1_samplings.ipynb`**
  Ejecuta múltiples instancias del `oracle` (en C++) con multithreading para generar puntos en el espacio paramétrico de manera acelerada.

* **`0b2_pipeline_train.ipynb`**
  Orquesta el pipeline de preprocesamiento: carga datos, filtra puntos no físicos y arma los datasets para entrenamiento.

* **`0b3_classify_post.ipynb`**
  Entrena un modelo de clasificación (por ejemplo, Random Forest) para predecir rápidamente la viabilidad de nuevos puntos.

* **`0b4_ML_class.ipynb`**
  Código detallado para entrenar, validar y evaluar un clasificador; incluye gráficos de métricas y matriz de confusión.

* **`0c2_bay_train.ipynb`**
  Entrena un modelo bayesiano (Gaussian Process) para predecir cantidades continuas, como branching ratios, sobre los puntos validados.

* **`02_m12.ipynb`**
  Exploración específica de la dependencia de $m_{12}^2$ y su impacto en restricciones teóricas.

* **`03_wanted_graph.ipynb`**
  Generación de gráficos personalizados (heatmaps, correcciones en 3D, etc.) a partir de datos preprocesados.

* **`04_valid_points.ipynb`**
  Visualización y análisis de puntos válidos extraídos del `oracle`.

* **`05_busquedaPrecisa.ipynb`**
  Implementación de una búsqueda dirigida (search grid fino) para refinar regiones de interés en el espacio paramétrico.

* **`06a_reunion.ipynb`**, **`06wip_genetics.ipynb`**, **`06_reunion.ipynb`**
  Notebooks de trabajo colaborativo/experimentación adicional (WIP).

* **`test_lib.ipynb`**
  Pruebas unitarias y ejemplos de uso de funciones contenidas en `lib/`.

> **Tip:** Si agregas o renombras notebooks, actualiza esta sección para mantenerla sincronizada.

---

## Estructura de proyecto

A modo de referencia, la siguiente salida muestra la jerarquía de archivos en `mlpython/`:

```bash
mlpython/
├── 0a0_get_from_lhe.ipynb
├── 0b1_samplings.ipynb
├── 0b2_pipeline_train.ipynb
├── 0b3_classify_post.ipynb
├── 0b4_ML_class.ipynb
├── 0c2_bay_train.ipynb
├── 02_m12.ipynb
├── 03_wanted_graph.ipynb
├── 04_valid_points.ipynb
├── 05_busquedaPrecisa.ipynb
├── 06a_reunion.ipynb
├── 06wip_genetics.ipynb
├── Analysis/              # Scripts y archivos auxiliares de análisis
├── data_batches/          # Datos intermedios generados por el oracle
├── data_batches_vary/     # Variaciones de data_batches según parámetros
├── data_train/            # Conjuntos de entrenamiento finales
├── docs/                  # Documentación adicional (diagramas, diseños)
├── graphing/              # Funciones de gráficos reutilizables
├── img/                   # Imágenes usadas en notebooks y documentación
├── lib/                   # Módulos Python personalizados (.py)
├── main.py                # Script principal que envuelve el pipeline completo
├── presentaciones/        # PPTs o PDFs para presentaciones internas
├── test_lib.ipynb         # Pruebas de la librería en `lib/`
├── ubenv_requirements.txt # Dependencias de entorno virtual
├── valid_points_lhe/      # Puntos válidos obtenidos desde archivos LHE
├── bitacora.md            # Registro de cambios y notas de desarrollo
├── README.md              # ← Este archivo
└── env/                   # Entorno virtual (no versionar)
```
