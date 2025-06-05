import glob
import pickle
from pathlib import Path
from typing import Any, Dict, Iterator, List, Sequence, Union


def get_data_files(pattern: str = "data_batches/merged_*.pkl") -> List[Path]:
    """
    Retorna la lista de rutas de archivos que coinciden con el patrón dado.
    """
    # Si en el futuro quisieras incluir también archivos batch_*.pkl, podrías hacer:
    # patterns = ["data_batches/merged_*.pkl", "data_batches/batch_*.pkl"]
    # paths = []
    # for p in patterns:
    #     paths.extend(glob.glob(p))
    # return sorted(Path(fp) for fp in paths)
    return sorted(Path(fp) for fp in glob.glob(pattern))


def load_pickle_file(fp: Path) -> Any:
    """
    Carga un objeto Python desde un archivo pickle.
    """
    with fp.open("rb") as f:
        return pickle.load(f)


def flatten_batch(batch: Dict[str, Any]) -> Iterator[tuple]:
    """
    Dado un diccionario de batch con claves "params" (array de N×7) y "results" (lista de N dicts),
    va devolviendo pares (params, result) uno por uno.
    """
    params_array = batch["params"]
    results_list = batch["results"]
    for p_vec, result_dict in zip(params_array, results_list):
        yield p_vec, result_dict


def flatten_point(entry: Dict[str, Any]) -> Iterator[tuple]:
    """
    Dado un diccionario que representa un punto único, extrae (params, result). 
    Puede tener la clave "result" o "results" (lista de 1 elemento).
    """
    p_vec = entry["params"]
    # Prioriza "result"; si no existe, usa "results"
    raw_r = entry.get("result", entry.get("results"))
    # A veces "results" viene como lista de longitud 1
    if isinstance(raw_r, list):
        result_dict = raw_r[0]
    else:
        result_dict = raw_r
    yield p_vec, result_dict


def create_row_from_params(params_vec: Sequence[float], result_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Toma el vector de parámetros (longitud 7) y el diccionario de resultados,
    y construye un nuevo diccionario plano con las claves nombradas.
    """
    # Desempaqueta los 7 valores de parámetros en variables con nombre
    m_phi, m_A, sin_ba, tan_beta, lambda_6, lambda_7, m12_2 = params_vec

    # Crea la fila combinando los campos de parámetros y los resultados
    row: Dict[str, Any] = {
        "m_phi": m_phi,
        "m_A": m_A,
        "sin_ba": sin_ba,
        "tan_beta": tan_beta,
        "lambda_6": lambda_6,
        "lambda_7": lambda_7,
        "m12_2": m12_2,
        **result_dict,
    }
    return row


def extract_rows_from_data(data: Any, source_path: Path) -> Iterator[Dict[str, Any]]:
    """
    Dado el objeto cargado desde pickle (que puede ser un dict de batch,
    o una lista mesclada de batches/puntos), va devolviendo todas las filas
    individuales en forma de diccionarios planos.
    """
    # Caso 1: Un solo batch (dict con "params" y "results")
    if isinstance(data, dict) and "params" in data and "results" in data:
        for p_vec, r_dict in flatten_batch(data):
            yield create_row_from_params(p_vec, r_dict)

    # Caso 2: Lista de entradas (pueden ser batches antiguos o puntos individuales)
    elif isinstance(data, list):
        for entry in data:
            if not isinstance(entry, dict):
                raise ValueError(f"Entrada inesperada en lista de {source_path}: {type(entry)}")

            # Subcaso 2a: este entry es un batch antiguo
            if "params" in entry and "results" in entry:
                for p_vec, r_dict in flatten_batch(entry):
                    yield create_row_from_params(p_vec, r_dict)

            # Subcaso 2b: este entry es un punto único con "params"
            elif "params" in entry:
                for p_vec, r_dict in flatten_point(entry):
                    yield create_row_from_params(p_vec, r_dict)

            else:
                raise ValueError(f"Claves inesperadas en entry de {source_path}: {entry.keys()}")

    else:
        raise ValueError(f"Estructura desconocida en archivo {source_path}: {type(data)}")


def collect_all_rows(file_paths: List[Path]) -> List[Dict[str, Any]]:
    """
    Recorre cada archivo pickle, extrae todas las filas y las acumula en una lista.
    """
    all_rows: List[Dict[str, Any]] = []
    for fp in file_paths:
        data = load_pickle_file(fp)
        for row in extract_rows_from_data(data, fp):
            all_rows.append(row)
    return all_rows

def main() -> None:
    """
    Punto de entrada: obtiene la lista de archivos, procesa cada uno y
    devuelve (o guarda) todas las filas en memoria. Aquí podrías, por ejemplo,
    convertirlas en un DataFrame de pandas si lo deseas:
        import pandas as pd
        df = pd.DataFrame(all_rows)
    """
    # 1) Listar todos los archivos que coinciden con el patrón
    data_files = get_data_files("data_batches/merged_*.pkl")

    # 2) Procesar cada archivo y extraer todas las filas
    rows = collect_all_rows(data_files)

    # 3) Si quisieras convertir a pandas.DataFrame:
    # import pandas as pd
    # df = pd.DataFrame(rows)
    # df.to_csv("tabla_completa.csv", index=False)

    # Por ahora, simplemente mostramos cuántas filas se obtuvieron:
    print(f"Total de filas extraídas: {len(rows)}")


if __name__ == "__main__":
    main()