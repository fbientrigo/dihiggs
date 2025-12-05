# -*- coding: utf-8 -*-
"""
Funciones para cargar, procesar y filtrar datos de archivos pickle y CSV
generados por la calculadora de modelos 2HDM.

"""

import os
import glob
import pickle
import pandas as pd
from typing import Iterable, Dict, List, Generator, Union


def find_files(pattern: str) -> List[str]:
    """
    Devuelve una lista ordenada de rutas de archivos que coinciden con el patrón dado.
    
    Parámetro:
        pattern (str): Patrón de búsqueda tipo glob, p.ej. "../data_batches/merged_0*.pkl" o "../out/*.csv"
        
    Retorna:
        List[str]: Lista de rutas de archivos ordenada.
    """
    return sorted(glob.glob(pattern))


def flatten_batch(batch: Dict) -> Generator[Dict, None, None]:
    """
    A partir de un diccionario de lote con "params" y "results", genera pares (params, result)
    y empaqueta cada par en un dict intermedio con los campos paramétricos y los resultados.
    
    Parámetros:
        batch (Dict): Diccionario con claves "params" (Nx7 array) y "results" (lista de N dicts).
        
    Genera:
        Dict: Cada diccionario contiene los 7 parámetros (m_phi, m_A, etc.) y las claves de resultado.
    """
    ps = batch["params"]
    rs = batch["results"]
    for p, r in zip(ps, rs):
        yield {
            "m_phi":    p[0],
            "m_A":      p[1],
            "sin_ba":   p[2],
            "tan_beta": p[3],
            "lambda_6": p[4],
            "lambda_7": p[5],
            "m12_2":    p[6],
            **r
        }


def flatten_point(entry: Dict) -> Generator[Dict, None, None]:
    """
    A partir de un diccionario de punto con "params" y "result" o "results", genera un único dict
    que combina los 7 parámetros y el diccionario de resultados (desenvuelto si es lista de tamaño 1).
    
    Parámetros:
        entry (Dict): Diccionario con clave "params" (array de 7) y "result" o "results".
        
    Genera:
        Dict: Un diccionario con los 7 parámetros y las claves de resultado.
    """
    p = entry["params"]
    r = entry.get("result", entry.get("results"))
    if isinstance(r, list):
        r = r[0]
    yield {
        "m_phi":    p[0],
        "m_A":      p[1],
        "sin_ba":   p[2],
        "tan_beta": p[3],
        "lambda_6": p[4],
        "lambda_7": p[5],
        "m12_2":    p[6],
        **r
    }


def _parse_pickle(fp: str) -> Generator[Dict, None, None]:
    with open(fp, "rb") as f:
        data = pickle.load(f)

    if isinstance(data, dict) and "params" in data and "results" in data:
        # Un único batch completo
        yield from flatten_batch(data)

    elif isinstance(data, list):
        # Lista de entradas: cada entry puede ser batch o punto
        for entry in data:
            if not isinstance(entry, dict):
                raise ValueError(f"Entrada inesperada en lista de {fp}: {type(entry)}")
            if "params" in entry and "results" in entry:
                yield from flatten_batch(entry)
            elif "params" in entry:
                yield from flatten_point(entry)
            else:
                raise ValueError(f"Claves inesperadas en entry de {fp}: {getattr(entry, 'keys', lambda: [])()}")
    else:
        raise ValueError(f"Estructura de dato desconocida en {fp}: {type(data)}")


def _parse_csv(fp: str) -> Generator[Dict, None, None]:
    """Lee un CSV y devuelve filas como dicts.
    No asume columnas exactas: emite tal cual vienen en el archivo.
    """
    df = pd.read_csv(fp, low_memory=False)
    for row in df.to_dict(orient="records"):
        yield row


def parse_file(fp: str) -> Generator[Dict, None, None]:
    """
    Procesa un archivo .pkl o .csv y genera dicts aplanados (fila a fila)
    según su estructura interna.
    
    - .pkl → estructuras con "params/results" (batch o punto)
    - .csv → filas tal cual del archivo
    """
    ext = os.path.splitext(fp)[1].lower()
    if ext == ".pkl":
        yield from _parse_pickle(fp)
    elif ext == ".csv":
        yield from _parse_csv(fp)
    else:
        raise ValueError(f"Extensión no soportada: {ext} (archivo: {fp})")


def load_rows(file_paths: List[str]) -> Generator[Dict, None, None]:
    """
    Itera sobre cada archivo y genera filas usando parse_file.
    """
    for fp in file_paths:
        yield from parse_file(fp)


def build_dataframe(rows: Iterable[Dict]) -> pd.DataFrame:
    """Construye un DataFrame a partir de un iterable de dicts."""
    return pd.DataFrame(rows)


def filter_errors(
    df: pd.DataFrame,
    error_col: str = "error",
    prohibited_errors: List[str] = None
) -> pd.DataFrame:
    """
    Filtra el DataFrame para eliminar filas cuyo campo de error esté en la lista de errores prohibidos.
    Luego elimina la columna de error.
    """
    if error_col not in df.columns:
        # Aviso pero devolvemos copia sin tocar
        print(f"[WARN] Columna '{error_col}' no encontrada. Se omite filtrado de errores.")
        return df.copy()

    if prohibited_errors is None:
        prohibited_errors = ["Execution failed", "Invalid parameter set."]
    mask = ~df[error_col].isin(prohibited_errors)
    df_filtered = df.loc[mask].copy()
    return df_filtered.drop(columns=[error_col])


def filter_valid_points(
    df: pd.DataFrame,
    conditions: Dict[str, Union[int, List[int]]]
) -> pd.DataFrame:
    """
    Filtra por condiciones binarias, p.ej. {"positivity_ok": 1, ...}.
    """
    print('\tFilas antes de filtrar:', len(df))
    mask = pd.Series(True, index=df.index)
    for col, valid_vals in conditions.items():
        if col not in df.columns:
            print(f"[WARN] Columna '{col}' no está en el DataFrame; se ignora esa condición.")
            continue
        if isinstance(valid_vals, list):
            mask &= df[col].isin(valid_vals)
        else:
            mask &= (df[col] == valid_vals)
    print('\tfilas despies de filtrar:', len(df.loc[mask]))
    return df.loc[mask].copy()


def save_dataframe(
    df: pd.DataFrame,
    output_path: str,
    format: str = "parquet",
    **kwargs
) -> None:
    """Guarda el DataFrame en parquet o csv."""
    if format == "parquet":
        df.to_parquet(output_path, index=False, **kwargs)
    elif format == "csv":
        df.to_csv(output_path, index=False, **kwargs)
    else:
        raise ValueError(f"Formato de salida no soportado: {format}")


def run_pipeline(
    file_pattern: str,
    prohibited_errors: List[str] = None,
    valid_conditions: Dict[str, Union[int, List[int]]] = None,
    output_path: str = "processed_data.parquet",
    output_format: str = "parquet"
) -> pd.DataFrame:
    """
    Ejecuta la pipeline completa:
      1. Busca archivos por patrón (CSV y/o PKL).
      2. Carga y aplanar/leer todos los datos en un DataFrame.
      3. Filtra filas con errores prohibidos (si existe la columna 'error').
      4. (Opcional) Filtra filas según condiciones binarias.
      5. Guarda el DataFrame resultante en disco.
    """
    files = find_files(file_pattern)
    if not files:
        raise FileNotFoundError(f"No se encontraron archivos con el patrón: {file_pattern}")
    print(files)

    rows_gen = load_rows(files)
    df = build_dataframe(rows_gen)
    print(df.head())

    df_no_errors = filter_errors(df, error_col="error", prohibited_errors=prohibited_errors)

    if valid_conditions:
        df_final = filter_valid_points(df_no_errors, valid_conditions)
    else:
        df_final = df_no_errors

    save_dataframe(df_final, output_path, format=output_format)
    return df_final


def go_by_merged_pattern(search_pattern: str = '', output_file: str = '') -> None:
    # Patrón para seleccionar archivos de entrada (se puede cambiar fácilmente)
    pattern = fr"{search_pattern}"

    # Errores a eliminar
    errores_prohibidos = ["Execution failed", "Invalid parameter set."]

    # Condiciones opcionales: quedarse con filas donde positivity_ok=1 y unitarity_ok=1
    condiciones = {
        "positivity_ok": 1,
        "unitarity_ok": 1,
        "perturbativity_ok": 1,
    }

    # Ejecutar toda la pipeline y guardar en un archivo CSV
    # Nota: si 'search_pattern' incluye '/\\*', el nombre de salida puede ser raro;
    # cámbialo si lo necesitas.
    df_procesado = run_pipeline(
        file_pattern=pattern,
        prohibited_errors=errores_prohibidos,
        valid_conditions=condiciones,
        output_path=output_file,
        output_format="csv"
    )

    print(f"Pipeline completada para {search_pattern}. Filas finales: {len(df_procesado)}")



go_by_merged_pattern(search_pattern='../../dihiggs/app/outcsv/*.csv', output_file='processed_outcsv_all.csv')
go_by_merged_pattern(search_pattern='../../dihiggs/app/outcsv/*/*.csv', output_file='processed_outcsv_subdirs_all.csv')
