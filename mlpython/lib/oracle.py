import subprocess
import numpy as np
import json

# nota: por ahora es necesario regenerar un ejecutable cada vez
# mas adelante es posible automatizar que el ejecutable quede dentro de lib/
# si se hacen modificaciones con Make, pero prefiero concentrar energia de desarrollo
EXECUTABLE_PATH = "/home/ftrigo/Dihiggs/dihiggs/app/Oracle"

def run_oracle(params, executable_path=EXECUTABLE_PATH):
    """
    Ejecuta el binario Oracle con una lista de parámetros.
    m_phi, mA, alpha, beta, lambda6, lambda7, m12
    
    """
    cmd = [executable_path] + list(map(str, params))
    
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )

        # Intentamos parsear la salida como JSON
        try:
            output = json.loads(result.stdout.strip())
        except json.JSONDecodeError:
            output = {"output": result.stdout.strip()}
        
        return output

    except subprocess.CalledProcessError as e:
        return {
            "error": "Execution failed",
            "stderr": e.stderr.strip(),
            "params": params
        }
    except Exception as e:
        return {
            "error": f"Unexpected error: {str(e)}",
            "params": params
        }




def safe_run_oracle(params):
    """
    Llama a run_oracle y devuelve salida uniforme con NaNs en caso de error.
    
    m_phi, mA, alpha, beta, lambda6, lambda7, m12

    Parameters:
        params (list): Lista de parámetros.

    Returns:
        dict: Diccionario con las claves esperadas y valores o NaNs.
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
    }

    if "error" in output:
        return {**expected_keys, "error": output["error"]}
    
    # Si todo salió bien, rellenar con lo recibido y mantener formato
    result = expected_keys.copy()
    for key in result:
        result[key] = output.get(key, result[key])  # usa valor o NaN por defecto

    return result
