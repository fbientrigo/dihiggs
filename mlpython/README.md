# mlpython — Estudio paramétrico con Machine Learning aplicado a modelos 2HDM

**Descripción breve:**
Capa de análisis en Python sobre el binario C++ `dihiggs/app/Oracle`: ejecuta simulaciones, procesa los resultados (parámetros y variables físicas) y explora la viabilidad de puntos en el espacio paramétrico con notebooks y pipelines de datos.

**Nota de arquitectura:** `mlpython/` está completamente desacoplado del pipeline principal del repo (`autoresearch/`, `scripts/`, CI). Nada fuera de esta carpeta lo importa; su única dependencia externa es el ejecutable `Oracle` (compilado con `make` en `dihiggs/`, ver `install.sh` en la raíz del repo).

---

## Índice
1. [Arquitectura](#arquitectura)
2. [Uso de `lib/` desde notebooks](#uso-de-lib-desde-notebooks)
3. [Convenciones del Oracle](#convenciones-del-oracle)
4. [Catálogo de notebooks](#catálogo-de-notebooks)
5. [Datos versionados](#datos-versionados)
6. [Estructura de proyecto](#estructura-de-proyecto)

---

## Arquitectura

`mlpython/` contiene tres islas independientes entre sí:

1. **`lib/`** — paquete Python importable que envuelve al `Oracle`:
   - `oracle.py`: `run_oracle`, `run_oracle_batch`, `safe_run_oracle` y la clase `OracleExecutor` (subprocesos con multithreading). La ruta al ejecutable se deriva automáticamente de la posición del repo — no hay rutas absolutas.
   - `loss.py`: `compute_max_unitarity_eigenvalue` (matrices de scattering 2HDM).
   - `explorer.py`: `explore_points`, `explore_points_parallel` (barridos de grillas).
   - `utils.py`: expansión de espacios de parámetros, estimación de tiempos.
   - `read_merged.py`: lectura/aplanado de batches pickle (`get_data_files`, `collect_all_rows`, …).
   - `branch_ratio.py`: `physpoint`, `plot_branchings`, `plot_branchings_and_lifetime`, …

   Los consumidores de `lib/` son los notebooks de `notebooks/` (ningún script `.py` externo lo importa).

2. **`lake_pipeline/`** (antes `2603/`) — pipeline autónomo de análisis sobre el data lake en Parquet (polars/pandas): consolidación, EDA, plots de ctau/BR/masa y coordenadas paralelas. Tiene su propia documentación (`lake_pipeline/README.md`) y runners (`run_pipeline.sh`, `launch.sh`). No depende de `lib/`.

3. **`Analysis/`** — suite local de notebooks de ML (sklearn: EDA, SMOTE, clasificación, t-SNE/UMAP). **No está versionada** (gitignored); vive solo en copias locales de trabajo.

Los duplicados históricos de código (`lib/pipeline.py`, `lib/advanced/`, `app/`, `main.py`) fueron eliminados; las versiones canónicas viven en `lib/`. Algunos notebooks aún redefinen funciones inline (`compute_max_unitarity_eigenvalue`, `physpoint`, `explore_points`) — al editarlos, preferir las versiones de `lib/`.

## Uso de `lib/` desde notebooks

Convención: cada notebook se ejecuta con el *working directory* en su propia carpeta. Como los notebooks viven en `notebooks/<grupo>/` (profundidad 2), la primera celda agrega `mlpython/` al `sys.path`:

```python
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.getcwd(), "..", "..")))

from lib.oracle import run_oracle, OracleExecutor
```

Requisito previo: compilar el Oracle (`cd dihiggs && make`, o `./install.sh` en la raíz). `lib.oracle.EXECUTABLE_PATH` apunta a `<repo>/dihiggs/app/Oracle` automáticamente.

## Convenciones del Oracle

Formato de salida típico al ejecutar el **oracle** con multithreading:

```python
[
    {
        'params': array([
            [1.50000001e+02, 2.99932979e+02, 1.00000000e+00, …, 1.00000000e-01, 0.00000000e+00, 2.24990902e+00],
            …
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
            …
        ]
    }
]
```

* **`params`**: arreglo NumPy con los parámetros físicos, en este orden:
  1. $m_{\phi}$
  2. $m_{A}$
  3. $\sin(\beta - \alpha)$
  4. $\tan\beta$
  5. $\lambda_{6}$
  6. $\lambda_{7}$
  7. $m_{12}^2$

* **`results`**: lista de diccionarios, uno por punto calculado (`positivity_ok`, `unitarity_ok`, anchos `w_h2_*`, branching ratios, acoplamientos `lambda1..7`).

## Catálogo de notebooks

Agrupados por flujo de trabajo bajo `notebooks/`:

* **`configuration/`**
  * `0_conf_creator.ipynb`, `0b_complementary_conf.ipynb` — generación de archivos de configuración para corridas del scanner (`complementary_config.conf`).

* **`lhe/`**
  * `0a0_get_from_lhe.ipynb` — extrae parámetros de los puntos válidos `.lha` en `valid_points_lhe/` y produce `valid_points_extended.csv`.
  * `07_LHE_study.ipynb` — estudio de los puntos LHE y comparación Oracle vs. valores calculados.

* **`exploration/`**
  * `test_lib.ipynb` — pruebas/ejemplos de uso de `lib/` (`explore_points`, `explore_points_parallel`).
  * `06wip_genetics.ipynb` — búsqueda con algoritmos genéticos (pymoo/NSGA-II) sobre el espacio de parámetros.

* **`valid_points/`** — flujo principal "puntos válidos":
  * `06b_valid_from_data.ipynb` — filtra los batches (`data_batches/merged_*.pkl`) a CSVs de puntos válidos.
  * `06a_get_from_known.ipynb` — recorre puntos conocidos con `OracleExecutor` por lotes.
  * `06a_getm12_search.ipynb`, `06a_reunion.ipynb` — búsqueda de raíces en $m_{12}^2$ vía `compute_max_unitarity_eigenvalue`.
  * `07a_get_m12_from_data.ipynb` — obtiene $m_{12}^2$ desde datos de scans.
  * `06c_construcGraph.ipynb` — construcción de gráficos de branching/lifetime (outputs limpiados; regenerar al ejecutar).
  * `optimized_eigen.ipynb` — optimización del cálculo de autovalores.

* **`graphing/`**
  * `graph.ipynb`, `2_graphing.ipynb` — experimentos de gráficos reutilizables.

> **Nota sobre datos faltantes:** varios notebooks históricos leen `../../data_batches/` y `dihiggs/app/outcsv*` — esos directorios de datos intermedios no están versionados y deben regenerarse con el Oracle antes de re-ejecutar.

## Datos versionados

Entradas que sí viven en git:
* `valid_points_lhe/` — 23 archivos `.lha` (80–300 GeV) + tablas SLHA (`Hbb`, `Haa`, `Haz`, `Hwidth`).
* `notebooks/valid_points/intervals_dict.pkl` — intervalos precomputados para el flujo de puntos válidos.
* `paper_img/` — figuras para el paper.

Todo lo demás (PNGs de diagnóstico en `img/`, pickles de `outs/`, `docs/_build/`, `presentaciones/`, `Analysis/`, outputs de `lake_pipeline/`) está gitignored: se conserva en copias locales pero no se versiona.

## Estructura de proyecto

```bash
mlpython/
├── README.md              # ← este archivo
├── bitacora.md            # registro de notas de desarrollo
├── lib/                   # paquete Python (wrapper del Oracle + análisis)
├── lake_pipeline/         # pipeline Parquet autónomo (ver su README)
├── notebooks/
│   ├── configuration/     # creación de configs
│   ├── lhe/               # extracción y estudio de puntos LHE
│   ├── exploration/       # pruebas de lib + búsqueda genética
│   ├── valid_points/      # flujo principal de puntos válidos
│   └── graphing/          # experimentos de gráficos
├── docs/                  # documentación Sphinx de lib/ (make html; _build/ no versionado)
├── paper_img/             # figuras para el paper (versionadas)
├── valid_points_lhe/      # datos de entrada .lha + tablas SLHA (versionados)
├── Analysis/              # suite local de ML (no versionada)
└── img/, outs/, branching_ratio_plots/, presentaciones/   # artefactos locales (no versionados)
```
