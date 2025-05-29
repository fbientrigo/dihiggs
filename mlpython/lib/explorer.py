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
from lib.oracle import run_oracle

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
