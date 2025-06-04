import numpy as np
from lib.oracle import run_oracle
import time
from itertools import product
from multiprocessing import Pool, cpu_count

# TODO
# en esta branch la logica de run_oracle debe cambiarse
# alpha_range --> sin_ba_range
# beta_range --> tan_beta range

def explore_points(m_phi_range, mA_range, alpha_range, beta_range, lambda6_range, lambda7_range, m12_range):
    """
    Explora un grid de combinaciones de parámetros físicos y devuelve solo los casos físicamente válidos.

    Parámetros:
    -----------
    m_phi_range : iterable
        Valores para el parámetro m_phi.
    mA_range : iterable
        Valores para el parámetro mA.
    sin(beta-alpha)_range : iterable
        Valores para el parámetro alpha.
    tan(beta)_range : iterable
        Valores para el parámetro beta.
    lambda6_range : iterable
        Valores para el parámetro lambda6.
    lambda7_range : iterable
        Valores para el parámetro lambda7.
    m12_range : iterable
        Valores para el parámetro m12.

    Retorna:
    --------
    results : list of dict
        Lista de diccionarios que contienen los parámetros de entrada y los resultados de decays y branching ratios,
        solo si se cumple que:
        - positivity_ok == 1
        - unitarity_ok == 1
        - perturbativity_ok == 1

    Ejemplo de uso:
    ---------------
    >>> from lib.explorer import explore_points
    >>> import numpy as np
    >>> results = explore_points(
    ...     m_phi_range=np.linspace(125, 300, 3),
    ...     mA_range=np.linspace(300, 320, 3),
    ...     alpha_range=np.linspace(0.001, 0.01, 3),
    ...     beta_range=np.linspace(1.55, 1.59, 3),
    ...     lambda6_range=np.linspace(-20, 20, 3),
    ...     lambda7_range=np.linspace(-20, 20, 3),
    ...     m12_range=np.linspace(0.0, 20.0, 3)
    ... )
    >>> print(f"Se encontraron {len(results)} configuraciones válidas.")
    """
    results = []
    total = (len(m_phi_range) * len(mA_range) * len(alpha_range) * len(beta_range) *
             len(lambda6_range) * len(lambda7_range) * len(m12_range))
    print("Total points:", total)
    count = 0
    start = time.time()

    for m_phi in m_phi_range:
        for mA in mA_range:
            for alpha in alpha_range:
                for beta in beta_range:
                    for lambda6 in lambda6_range:
                        for lambda7 in lambda7_range:
                            for m12 in m12_range:
                                count += 1
                                params = [m_phi, mA, alpha, beta, lambda6, lambda7, m12]
                                try:
                                    out = run_oracle(params)
                                    if out.get("positivity_ok") == 1 and out.get("unitarity_ok") == 1 and out.get("perturbativity_ok") == 1:
                                        record = {
                                            "m_phi": m_phi,
                                            "mA": mA,
                                            "alpha": alpha,
                                            "beta": beta,
                                            "lambda6": lambda6,
                                            "lambda7": lambda7,
                                            "m12": m12,
                                            "positivity_ok": out.get("positivity_ok", 0),
                                            "unitarity_ok": out.get("unitarity_ok", 0),
                                            "perturbativity_ok": out.get("perturbativity_ok", 0),
                                            "w_h2_bb": out.get("w_h2_bb"),
                                            "w_h2_tautau": out.get("w_h2_tautau"),
                                            "w_h2_WW": out.get("w_h2_vv", [np.nan, np.nan, np.nan])[0],
                                            "w_h2_ZZ": out.get("w_h2_vv", [np.nan, np.nan, np.nan])[1],
                                            "w_h2_gaga": out.get("w_h2_gaga"),
                                            "w_h2_Zga": out.get("w_h2_Zga"),
                                            "w_h2_gg": out.get("w_h2_gg"),
                                            "w_h2_hh": out.get("w_h2_hh"),
                                            "w_total_h2": out.get("w_total_h2"),
                                            "branching_ratio_h2_gaga": out.get("branching_ratio_h2_gaga")
                                        }
                                        results.append(record)
                                except Exception as e:
                                    print("Error for params:", params, "Error:", e)

                                if count % 1000 == 0:
                                    elapsed = time.time() - start
                                    print(f"Processed {count}/{total} points, elapsed {elapsed:.2f} s")

    return results


import numpy as np
import time
from itertools import product
from multiprocessing import Pool, cpu_count
from .oracle import run_oracle
import numpy as np
from scipy.stats import qmc

def _evaluate_point(params):
    """
    Evalúa un único punto del espacio de parámetros y retorna el resultado si es físicamente válido.
    Esta función se usa internamente con multiprocessing.
    """
    try:
        out = run_oracle(params)
        if out.get("positivity_ok") == 1 and out.get("unitarity_ok") == 1 and out.get("perturbativity_ok") == 1:
            m_phi, mA, alpha, beta, lambda6, lambda7, m12 = params
            return {
                "m_phi": m_phi,
                "mA": mA,
                "alpha": alpha,
                "beta": beta,
                "lambda6": lambda6,
                "lambda7": lambda7,
                "m12": m12,
                "positivity_ok": out.get("positivity_ok", 0),
                "unitarity_ok": out.get("unitarity_ok", 0),
                "perturbativity_ok": out.get("perturbativity_ok", 0),
                "w_h2_bb": out.get("w_h2_bb"),
                "w_h2_tautau": out.get("w_h2_tautau"),
                "w_h2_WW": out.get("w_h2_vv", [np.nan, np.nan, np.nan])[0],
                "w_h2_ZZ": out.get("w_h2_vv", [np.nan, np.nan, np.nan])[1],
                "w_h2_gaga": out.get("w_h2_gaga"),
                "w_h2_Zga": out.get("w_h2_Zga"),
                "w_h2_gg": out.get("w_h2_gg"),
                "w_h2_hh": out.get("w_h2_hh"),
                "w_total_h2": out.get("w_total_h2"),
                "branching_ratio_h2_gaga": out.get("branching_ratio_h2_gaga")
            }
    except Exception as e:
        print("Error for params:", params, "Error:", e)
    return None  # Para puntos inválidos o con error

def explore_points_parallel(m_phi_range, mA_range, alpha_range, beta_range, lambda6_range, lambda7_range, m12_range, n_processes=None):
    """
    Versión paralela de explore_points usando multiprocessing.Pool.

    Parámetros:
    -----------
    Los mismos que explore_points + n_processes opcional.

    Retorna:
    --------
    List of dict con resultados físicamente válidos.
    """
    all_combinations = list(product(m_phi_range, mA_range, alpha_range, beta_range, lambda6_range, lambda7_range, m12_range))
    total = len(all_combinations)
    print(f"Total points: {total}")
    
    start = time.time()
    print(cpu_count())
    with Pool(processes=n_processes or cpu_count()) as pool:
        results_raw = pool.map(_evaluate_point, all_combinations)

    # Filtrar los None (errores o puntos no válidos)
    results = [r for r in results_raw if r is not None]

    elapsed = time.time() - start
    print(f"Exploración completada en {elapsed:.2f} segundos con {len(results)} puntos válidos.")

    return results

vev = 246

def generate_local_variations(
    m_phi_base: float,
    m_A_center: float,
    m12_center: float,
    batch_size: int,
    eps_A: float,
    eps_m12: float,
    seed: int = None
    ) -> np.ndarray:
    """
    Generate `batch_size` points where only m_A and m12_2 vary in a small Latin 
    Hypercube around (m_A_center, m12_center), and all other 5 dims are fixed.

    Returns an array of shape (batch_size, 7) with column order:
      [m_phi, m_A, sin_ba, tan_beta, lambda6, lambda7, m12_2]
    """

    # 1) LatinHypercube in 2D (for m_A and m12)
    sampler = qmc.LatinHypercube(d=2)
    unit = sampler.random(n=batch_size)

    # 2) Scale to [m_A_center±eps_A] × [m12_center±eps_m12]
    bounds = np.array([
        [m_A_center - eps_A, m_A_center + eps_A],
        [m12_center - eps_m12, m12_center + eps_m12]
    ])
    scaled = qmc.scale(unit, bounds[:,0], bounds[:,1])  # shape (batch_size,2)

    # 3) Build the full parameter array, hard-coding the other 5 dims:
    P = np.empty((batch_size, 7), dtype=float)
    P[:, 0] = m_phi_base           # m_phi
    P[:, 1] = scaled[:, 0]         # m_A (varying)
    P[:, 2] = 1.0                  # sin(beta - alpha), fixed
    P[:, 3] = 10000.0              # tan(beta), fixed
    P[:, 4] = 0.1                  # lambda_6, fixed
    P[:, 5] = 0.0                  # lambda_7, fixed
    P[:, 6] = scaled[:, 1]         # m12_2 (varying)

    return P


def get_parameters_from_points(
    csv_path: str,
    batch_size: int,
    eps_A: float,
    eps_m12: float,
    seed: int = None
    ) -> np.ndarray:
    """
    Reads `csv_path` containing base points with columns 'Mh2', 'Mh3', 'm12_2' and
    for each row generates `batch_size` local variations via generate_local_variations.
    Returns a combined array of shape (n_rows * batch_size, 7).
    """
    df = pd.read_csv(csv_path)
    all_batches = []
    for idx, row in df.iterrows():
        m_phi_base   = row["Mh2"]    # heavy CP-even Higgs mass as m_phi
        m_A_center   = row["Mh3"]    # CP-odd Higgs mass as m_A
        m12_center   = row["m12_2"]
        # derive unique seed per batch for reproducibility
        batch_seed = None if seed is None else seed + idx
        P = generate_local_variations(
            m_phi_base, m_A_center, m12_center,
            batch_size, eps_A, eps_m12, seed=batch_seed
        )
        all_batches.append(P)
    # stack all batches into one array
    return np.vstack(all_batches)




def generate_local_variations_phi(
    m_phi_center: float,
    m12_center: float,
    batch_size: int,
    eps_phi: float,
    eps_m12: float,
    m_A_fixed: float = 300.0,
    seed: int = None
) -> np.ndarray:
    """
    Genera `batch_size` puntos variando m_phi y m12_2 en un Latin Hypercube
    alrededor de (m_phi_center, m12_center).
    m_A se fija a `m_A_fixed`. Las otras 4 dimensiones están hardcodeadas.
    Retorna un array de forma (batch_size, 7) con columnas:
    [m_phi, m_A, sin_ba, tan_beta, lambda6, lambda7, m12_2]
    """
    # 1) Muestreo LH en 2D para (m_phi, m12_2)
    sampler = qmc.LatinHypercube(d=2)
    unit = sampler.random(n=batch_size)

    # 2) Escalado a [m_phi_center±eps_phi] × [m12_center±eps_m12]
    bounds = np.array([
        [m_phi_center - eps_phi, m_phi_center + eps_phi],
        [m12_center  - eps_m12,  m12_center  + eps_m12]
    ])
    scaled = qmc.scale(unit, bounds[:,0], bounds[:,1])  # shape (batch_size,2)

    # 3) Construir el array de parámetros
    P = np.empty((batch_size, 7), dtype=float)
    P[:, 0] = scaled[:, 0]       # m_phi (variado)
    P[:, 1] = m_A_fixed           # m_A (fijo)
    P[:, 2] = 1.0                 # sin(beta - alpha)
    P[:, 3] = 10000.0             # tan(beta)
    P[:, 4] = 0.1                 # lambda6
    P[:, 5] = 0.0                 # lambda7
    P[:, 6] = scaled[:, 1]       # m12_2 (variado)

    return P

def get_parameters_from_points_phi(
    csv_path: str,
    batch_size: int,
    eps_phi: float,
    eps_m12: float,
    m_A_fixed: float = 300.0,
    seed: int = None
) -> np.ndarray:
    """
    Lee `csv_path` con columnas 'Mh2' (m_phi_center) y 'm12_2'.
    Para cada fila genera un batch con `generate_local_variations_phi`.
    """
    import pandas as pd
    df = pd.read_csv(csv_path)
    all_batches = []
    for idx, row in df.iterrows():
        m_phi_center = row["Mh2"]
        m12_center   = row["m12_2"]
        batch_seed   = None if seed is None else seed + idx
        P = generate_local_variations_phi(
            m_phi_center,
            m12_center,
            batch_size,
            eps_phi,
            eps_m12,
            m_A_fixed=m_A_fixed,
            seed=batch_seed
        )
        # DEBUG: validación rápida
        # assert np.all(P[:,1] == m_A_fixed), "m_A no está fijo a 300"
        all_batches.append(P)

    return np.vstack(all_batches)


import os
import glob
import time
import pickle
import pandas as pd
from typing import Literal
from .oracle import OracleExecutor

def multiple_runs(
    csv_path: str,
    N_repeat_runs: int,
    batch_size: int,
    eps_A: float,
    eps_m12: float,
    outdir: str,
    executor: OracleExecutor,
    variation_mode: Literal['mA', 'mPhi'] = 'mA',
    eps_phi: float = None,
    base_seed: int = 42
):
    """
    Ejecuta varias corridas locales LH. Dependiendo de `variation_mode`:
      - 'mA': varía (m_A, m12_2), fija m_phi (usa eps_A)
      - 'mPhi': varía (m_phi, m12_2), fija m_A=300, usa eps_phi
    """
    # Validaciones básicas
    if variation_mode == 'mPhi' and eps_phi is None:
        raise ValueError("Para mode='mPhi' debes proporcionar eps_phi")
    
    df = pd.read_csv(csv_path)
    n_runs = len(df)

    # Estimación de tiempo
    time_per_point = 415.7 / 15_000
    total_points = N_repeat_runs * n_runs * batch_size
    pred_mins = total_points * time_per_point / 60
    print(f"Estimado: {pred_mins:.1f} min (~{pred_mins/60:.2f} h) para {total_points} puntos")

    os.makedirs(outdir, exist_ok=True)

    for epoch in range(N_repeat_runs):
        for j, row in df.iterrows():
            m_phi_base = float(row["Mh2"])
            m_A_center = float(row["Mh3"])
            m12_center = float(row["m12_2"])
            # Semilla reproducible única
            seed = base_seed + epoch * n_runs + j

            # Seleccionar la función de variación según el modo
            if variation_mode == 'mA':
                param_list = generate_local_variations(
                    m_phi_base=m_phi_base,
                    m_A_center=m_A_center,
                    m12_center=m12_center,
                    batch_size=batch_size,
                    eps_A=eps_A,
                    eps_m12=eps_m12,
                    seed=seed
                )
            else:  # 'mPhi'
                param_list = generate_local_variations_phi(
                    m_phi_center=m_phi_base,
                    m12_center=m12_center,
                    batch_size=batch_size,
                    eps_phi=eps_phi,
                    eps_m12=eps_m12,
                    m_A_fixed=m_A_center,  # o un valor fijo, p.ej. 300
                    seed=seed
                )

            # Índice de batch basado en archivos existentes
            existing = sorted(glob.glob(os.path.join(outdir, "batch_*.pkl")))
            batch_idx = len(existing) + 1

            # Ejecución
            t0 = time.perf_counter()
            results = executor.map(param_list.tolist(), use_threads=True)
            dt = time.perf_counter() - t0

            # Guardar resultados
            fname = f"batch_{batch_idx}_{int(m_phi_base)}.pkl"
            fout = os.path.join(outdir, fname)
            with open(fout, "wb") as f:
                pickle.dump({"params": param_list, "results": results}, f)

            print(f"[Run {j+1}/{n_runs}] batch {batch_idx} "
                  f"mode={variation_mode} saved in {dt:.1f}s → {fout}")

        print(f"[Epoch {epoch+1}/{N_repeat_runs}] completado")
    print("Todas las corridas finalizadas.")