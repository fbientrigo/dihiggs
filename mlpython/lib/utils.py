import numpy as np

# importing data



# ---
def expand_parameter_space(parameter_space):
    """
    Convierte un diccionario de definiciones de espacio de parámetros con min, max y step
    en un diccionario de arrays de valores a explorar.

    Parámetros
    ----------
    parameter_space : dict
        Diccionario con nombres de parámetros como claves, y valores tipo:
        {
            "param": {"min": float, "max": float, "step": float},
            ...
        }

    Retorna
    -------
    dict
        Diccionario donde cada clave tiene un array de valores generado con np.arange.

    Ejemplo
    -------
    >>> parameter_space = {
    ...     "alpha": {"min": -0.5, "max": 0.5, "step": 0.1},
    ...     "beta": {"min": -1.6, "max": 1.6, "step": 0.5}
    ... }
    >>> ranges = expand_parameter_space(parameter_space)
    >>> print(ranges["alpha"])
    array([-0.5, -0.4, ..., 0.5])
    """
    ranges = {}
    for param, bounds in parameter_space.items():
        pmin, pmax, pstep = bounds["min"], bounds["max"], bounds["step"]
        ranges[f"{param}_range"] = np.round(np.arange(pmin, pmax + pstep, pstep), 6)
    return ranges


def estimate_exploration_time(
    m_phi,
    mA,
    alpha,
    beta,
    lambda6,
    lambda7,
    m12,
    points_per_second=30
):
    """
    Estima la cantidad total de combinaciones de parámetros a explorar y el tiempo aproximado de ejecución.

    Parámetros
    ----------
    m_phi_range : iterable
        Valores del parámetro m_phi.
    mA_range : iterable
        Valores del parámetro mA.
    alpha_range : iterable
        Valores del parámetro alpha.
    beta_range : iterable
        Valores del parámetro beta.
    lambda6_range : iterable
        Valores del parámetro lambda6.
    lambda7_range : iterable
        Valores del parámetro lambda7.
    m12_range : iterable
        Valores del parámetro m12.
    points_per_second : float, opcional
        Velocidad estimada de procesamiento en puntos por segundo (default: 50).

    Retorna
    -------
    dict
        Diccionario con:
        - "n_points": cantidad total de combinaciones
        - "estimated_minutes": tiempo estimado en minutos
        - "estimated_hours": tiempo estimado en horas

    Ejemplo
    -------
    >>> from lib.utils import estimate_exploration_time
    >>> result = estimate_exploration_time(m_phi_range, mA_range, alpha_range, beta_range,
    ...                                    lambda6_range, lambda7_range, m12_range)
    >>> print(result["n_points"], result["estimated_minutes"], "[min]", result["estimated_hours"], "[h]")
    """
    n_points = len(m_phi) * len(mA) * len(alpha) * len(beta) * len(lambda6) * len(lambda7) * len(m12)


    estimated_minutes = n_points / points_per_second / 60.0
    estimated_hours = estimated_minutes / 60.0

    return {
        "n_points": n_points,
        "estimated_minutes": estimated_minutes,
        "estimated_hours": estimated_hours
    }


# Seleccionar una combinación aleatoria
def random_combination(param_ranges, seed=None):
    if seed == None:
        np.random.seed()  # aseguramos aleatoriedad
    else:
        np.random.seed(seed)
    sample = {key: np.random.choice(values) for key, values in param_ranges.items()}
    oracle_input = [
        sample["mphi"],
        sample["mA"],
        sample["alpha"],
        sample["beta"],
        sample["lambda6"],
        sample["lambda7"],
        sample["m12_squared"]
    ]
    return oracle_input


# --- threading-----
# lib/utils.py
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed

def thread_map(func, inputs, max_workers=4):
    """Mapear func(arg) sobre inputs en paralelo con hilos."""
    results = [None] * len(inputs)
    with ThreadPoolExecutor(max_workers=max_workers) as exe:
        fut2i = {exe.submit(func, arg): i for i, arg in enumerate(inputs)}
        for fut in as_completed(fut2i):
            results[fut2i[fut]] = fut.result()
    return results

def process_map(func, inputs, max_workers=4):
    """Mapear func(arg) sobre inputs en paralelo con procesos."""
    results = [None] * len(inputs)
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        fut2i = {exe.submit(func, arg): i for i, arg in enumerate(inputs)}
        for fut in as_completed(fut2i):
            results[fut2i[fut]] = fut.result()
    return results

