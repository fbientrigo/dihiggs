import subprocess
import numpy as np
import json
from pathlib import Path

# --------------------------------------------------
# Configuración de ruta al ejecutable
# --------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent
EXECUTABLE_PATH = PROJECT_ROOT / 'dihiggs' / 'app' / 'Oracle'


def run_oracle_batch(param_list, nthreads, executable_path=EXECUTABLE_PATH, debug=False):
    """
    Llama al binario en modo paralelo, entregándole nthreads
    y un array de param_list (lista de listas de 7 floats).
    Devuelve la lista de dicts resultantes.
    """
    flat_args = ["--nthreads", str(nthreads)]
    for params in param_list:
        assert len(params) == 7, "Cada conjunto de parámetros debe tener 7 valores"
        flat_args += list(map(str, params))

    cmd = [str(executable_path)] + flat_args
    if debug:
        print("▶ Ejecutando:", " ".join(cmd))

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True
    )

    if debug:
        print("⮕ returncode:", proc.returncode)
        print("⮕ stderr:\n", proc.stderr.strip())
        print("⮕ stdout (inicio):\n", proc.stdout[:200])

    try:
        outputs = json.loads(proc.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(f"Salida mal formada de Oracle:\n{proc.stdout}")

    return outputs


def run_oracle(params, executable_path=EXECUTABLE_PATH, debug=False):
    """
    Ejecuta el binario Oracle con una lista de parámetros (7 floats) y devuelve el dict JSON.
    """
    cmd = [str(executable_path)] + list(map(str, params))

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        if debug:
            print("▶ Ejecutando:", " ".join(cmd))
            print("⮕ stderr:\n", result.stderr)
            print("⮕ stdout:\n", result.stdout)

        try:
            output = json.loads(result.stdout.strip())
        except json.JSONDecodeError:
            output = {"output": result.stdout.strip()}
        return output

    except subprocess.CalledProcessError as e:
        return {"error": "Execution failed", "stderr": e.stderr.strip(), "params": params}
    except Exception as e:
        return {"error": f"Unexpected error: {str(e)}", "params": params}


def safe_run_oracle(params):
    """
    Llama a run_oracle y devuelve salida uniforme con NaNs en caso de error.
    """
    output = run_oracle(params)
    expected_keys = {
        'positivity_ok': np.nan,
        'unitarity_ok': np.nan,
        'perturbativity_ok': np.nan,
        'w_h2_bb': np.nan,
        'w_h2_tautau': np.nan,
        'w_h2_uu': np.nan,
        'w_h2_du': np.nan,
        'w_h2_ln': np.nan,
        'w_h2_vv': [np.nan, np.nan, np.nan],
        'w_h2_gaga': np.nan,
        'w_h2_Zga': np.nan,
        'w_h2_gg': np.nan,
        'w_h2_hh': np.nan,
        'w_total_h2': np.nan,
        'w_total_top': np.nan,
        'branching_ratio_h2_gaga': np.nan,
        'lambda1': np.nan,
        'lambda2': np.nan,
        'lambda3': np.nan,
        'lambda4': np.nan,
        'lambda5': np.nan,
        'lambda6': np.nan,
        'lambda7': np.nan,
    }
    if "error" in output:
        return {**expected_keys, "error": output["error"]}
    result = expected_keys.copy()
    for key in result:
        result[key] = output.get(key, result[key])
    return result


# --------------------------------------------------
# Main para testing desde línea de comandos
# --------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Test harness para Oracle Python wrapper"
    )
    subparsers = parser.add_subparsers(dest="command")

    # Prueba single
    p1 = subparsers.add_parser("single", help="Test run_oracle con parámetros individuales")
    p1.add_argument("params", nargs=7, type=float, help="7 parámetros para Oracle")
    p1.add_argument("--debug", action="store_true")

    # Prueba batch
    p2 = subparsers.add_parser("batch", help="Test run_oracle_batch con lote de parámetros")
    p2.add_argument("nthreads", type=int, help="Número de hilos OpenMP")
    p2.add_argument("params", nargs='+', type=float, help="Lista plana de parámetros, múltiplo de 7")
    p2.add_argument("--debug", action="store_true")

    args = parser.parse_args()
    if args.command == "single":
        res = run_oracle(args.params, debug=args.debug)
        print(json.dumps(res, indent=2))
    elif args.command == "batch":
        flat = args.params
        if len(flat) % 7 != 0:
            parser.error("El número de parámetros debe ser múltiplo de 7.")
        batches = [flat[i:i+7] for i in range(0, len(flat), 7)]
        res = run_oracle_batch(batches, args.nthreads, debug=True)
        print(json.dumps(res, indent=2))
    else:
        parser.print_help()
