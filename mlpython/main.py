import subprocess

def run_oracle(params):
    executable = "/home/ftrigo/DiHiggsWorking/dihiggs/app/Oracle"  # Ajusta la ruta
    cmd = [executable] + list(map(str, params))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    return result.stdout

# Ejemplo de uso
params = [130.0, 150.0, 0.005, 0.8, 0.1, 0.1, 10000.0]
output = run_oracle(params)
print("Salida de Oracle.cpp:")
print(output)
