





# ======================= extracting from merged ========================
# comes from 06b_valid_from_data
import glob
import pickle
import pandas as pd
from typing import Iterable, Dict, List, Generator, Union


def find_files(pattern: str) -> List[str]:
    """
    Devuelve una lista ordenada de rutas de archivos que coinciden con el patrón dado.
    
    Parámetro:
        pattern (str): Patrón de búsqueda tipo glob, p.ej. "../data_batches/merged_0*.pkl"
        
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
    # "result" simple o "results" que puede ser lista de longitud 1
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


def parse_file(fp: str) -> Generator[Dict, None, None]:
    """
    Procesa un archivo pickle y genera dicts aplanados (fila a fila) según su estructura interna:
      - Caso 1: data es un dict con "params" y "results" → se llama a flatten_batch
      - Caso 2: data es una lista de entries → cada entry puede ser batch-dict o point-dict
    
    Parámetros:
        fp (str): Ruta al archivo pickle.
        
    Genera:
        Dict: Cada diccionario representa una fila con parámetros y resultados.
    """
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
                # Este entry es un batch completo
                yield from flatten_batch(entry)
            elif "params" in entry:
                # Entry de punto simple
                yield from flatten_point(entry)
            else:
                raise ValueError(f"Claves inesperadas en entry de {fp}: {entry.keys()}")

    else:
        raise ValueError(f"Estructura de dato desconocida en {fp}: {type(data)}")


def load_rows(file_paths: List[str]) -> Generator[Dict, None, None]:
    """
    Itera sobre cada archivo y genera filas aplanadas usando parse_file.
    
    Parámetros:
        file_paths (List[str]): Lista de rutas a archivos pickle.
        
    Genera:
        Dict: Cada diccionario representa una fila completa con parámetros y resultados.
    """
    for fp in file_paths:
        yield from parse_file(fp)


def build_dataframe(rows: Iterable[Dict]) -> pd.DataFrame:
    """
    Construye un DataFrame de pandas a partir de un iterable de diccionarios.
    
    Parámetros:
        rows (Iterable[Dict]): Iterable que produce dicts (cada dict = una fila).
        
    Retorna:
        pd.DataFrame: DataFrame con todas las filas cargadas.
    """
    return pd.DataFrame(rows)


def filter_errors(
    df: pd.DataFrame,
    error_col: str = "error",
    prohibited_errors: List[str] = None
) -> pd.DataFrame:
    """
    Filtra el DataFrame para eliminar filas cuyo campo de error esté en la lista de errores prohibidos.
    Luego elimina la columna de error.
    
    Parámetros:
        df (pd.DataFrame): DataFrame original.
        error_col (str): Nombre de la columna que contiene la descripción de error.
        prohibited_errors (List[str]): Lista de mensajes de error a filtrar.
        
    Retorna:
        pd.DataFrame: DataFrame filtrado sin la columna de error.
    """
    # 1) Si no existe la columna, avisar y devolver df intacto
    if error_col not in df.columns:
        print(f"[WARN] Columna '{error_col}' no encontrada. Se omite filtrado de errores.")
        return df.copy()

    if prohibited_errors is None:
        prohibited_errors = ["Execution failed", "Invalid parameter set."]
    # 2) Filtrar filas donde el valor de error NO esté en prohibited_errors
    mask = ~df[error_col].isin(prohibited_errors)
    df_filtered = df.loc[mask].copy()
    # 3) Eliminar la columna de error
    return df_filtered.drop(columns=[error_col])


def filter_valid_points(
    df: pd.DataFrame,
    conditions: Dict[str, Union[int, List[int]]]
) -> pd.DataFrame:
    """
    Filtra el DataFrame para quedarse solo con filas que cumplan ciertas condiciones en columnas binarias.
    Por ejemplo: {"positivity_ok": 1, "unitarity_ok": 1, "perturbativity_ok": [0,1]}.
    
    Parámetros:
        df (pd.DataFrame): DataFrame a filtrar.
        conditions (Dict[str, Union[int, List[int]]]): Diccionario donde cada clave es un nombre de columna,
            y el valor es un entero o lista de enteros que se consideran válidos.
            
    Retorna:
        pd.DataFrame: DataFrame con las filas que satisfacen todas las condiciones.
    """
    print('\tFilas antes de filtrar:', len(df))
    mask = pd.Series(True, index=df.index)
    for col, valid_vals in conditions.items():
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
    """
    Guarda el DataFrame en disco en el formato especificado (parquet o csv).
    
    Parámetros:
        df (pd.DataFrame): DataFrame a guardar.
        output_path (str): Ruta de salida, e.g. "../output/data.parquet" o "../output/data.csv".
        format (str): "parquet" o "csv". Por defecto, "parquet".
        **kwargs: Parámetros adicionales para pandas.DataFrame.to_parquet o to_csv.
    """
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
      1. Busca archivos pickle según patrón.
      2. Carga y aplana todos los datos en un DataFrame.
      3. Filtra filas con errores prohibidos.
      4. (Opcional) Filtra filas según condiciones binarias.
      5. Guarda el DataFrame resultante en disco.
    
    Parámetros:
        file_pattern (str): Patrón glob para localizar archivos de entrada.
        prohibited_errors (List[str]): Lista de mensajes de error a filtrar.
        valid_conditions (Dict[str, Union[int, List[int]]]): Condiciones para filtrar columnas binarias.
        output_path (str): Ruta donde se guardará el DataFrame procesado.
        output_format (str): "parquet" o "csv".
    
    Retorna:
        pd.DataFrame: DataFrame final después de filtrados.
    """
    # 1) Encontrar archivos
    files = find_files(file_pattern)
    if not files:
        raise FileNotFoundError(f"No se encontraron archivos con el patrón: {file_pattern}")
    print(files)

    # 2) Cargar y aplanar datos
    rows_gen = load_rows(files)
    df = build_dataframe(rows_gen)
    print(df.head())

    # 3) Filtrar errores prohibidos
    df_no_errors = filter_errors(df, error_col="error", prohibited_errors=prohibited_errors)

    # 4) Filtrar puntos válidos si se indicaron condiciones
    if valid_conditions:
        df_final = filter_valid_points(df_no_errors, valid_conditions)
    else:
        df_final = df_no_errors

    # 5) Guardar resultado
    save_dataframe(df_final, output_path, format=output_format)
    return df_final


# --->>>>>>> main one
def go_by_merged_pattern(search_pattern = '', output_folder='./output/', output_format = 'csv'):

    # Patrón para seleccionar archivos de entrada (se puede cambiar fácilmente)
    pattern = fr"../data_batches/merged_{search_pattern}*.pkl"

    # Errores a eliminar
    errores_prohibidos = ["Execution failed", "Invalid parameter set."]

    # Condiciones opcionales: quedarse con filas donde positivity_ok=1 y unitarity_ok=1
    condiciones = {
        "positivity_ok": 1,
        "unitarity_ok": 1,
        # "perturbativity_ok": 1,  # Descomentar si se quiere también filtrar por perturbativity_ok
    }

    # Ejecutar toda la pipeline y guardar en un archivo Parquet
    df_procesado = run_pipeline(
        file_pattern=pattern,
        prohibited_errors=errores_prohibidos,
        valid_conditions=condiciones,
        output_path= output_folder + f"filtered_data_{search_pattern}.csv",
        output_format=output_format
    )

    print(f"Pipeline completada para {search_pattern}. Filas finales: {len(df_procesado)}")