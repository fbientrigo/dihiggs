# ML Python
Version de machine learning para el estudio parametrico

## Indice
- Convenciones
- Notebooks
- 
## Convenciones
A lo largo del analisis se utiliza un output tipico de las carpetas del sistema, en el siguiente ejemplo veremos la salida tipica de ejecutar el `oracle` con multithreading, invocando codigo en C++.
Normalmente se obtiene un solo diccionario, pero con outputs de multi hilo tendra el siguiente formato

```python
[{'params': array([[1.50000001e+02, 2.99932979e+02, 1.00000000e+00, ...,
        1.00000000e-01, 0.00000000e+00, 2.24990902e+00],
       [1.50000001e+02, 3.00088086e+02, 1.00000000e+00, ...,
       ...,
       [1.50000001e+02, 3.00033035e+02, 1.00000000e+00, ...,
        1.00000000e-01, 0.00000000e+00, 2.24990903e+00]]), 'results': [{'positivity_ok': 1, 'unitarity_ok': 1, 'perturbativity_ok': 1, 'w_h2_bb': 0.0, 'w_h2_tautau': 0.0, 'w_h2_uu': 0.0, 'w_h2_du': 0.0, 'w_h2_ln': 0.0, 'w_h2_vv': [5.053510327013813e-11, 0.0, 0.0], 'w_h2_gaga': 5.053510327013813e-11, 'w_h2_Zga': 1.181312528477372e-11, 'w_h2_gg': 0.0, 'w_h2_hh': 0.0, 'w_total_h2': 6.234822855491185e-11, 'w_total_top': 1.338069402459526, 'branching_ratio_h2_gaga': 0.8105298970223419, 'lambda1': 1.129741358125329, 'lambda2': 0.2577337926166231, 'lambda3': 2.483237305682109, 'lambda4': -1.112766760892846, 'lambda5': -1.112766760892846, 'lambda6': 0.1, 'lambda7': 0.0}, {'positivity_ok': 1, 'unitarity_ok': 1, 'perturbativity_ok': 1, 'w_h2_bb': 0.0 ...,
        1.00000000e-01, 0.00000000e+00, 2.55990900e+00],
```

los
- `params` son en el formato $[m_\phi, m_A, \sin(\beta-\alpha), \tan \beta, \lambda_6, \lambda_7, m_12^2]$
- `results` requiere escribir toda la documentacion

## Notebooks
- `0b1_samplings.ipynb` permite ejecutar de manera acelerada varias ejecucicones del `oracle`:`oraculo` para generar puntos
- `0a0_get_from_llhe.ipynb` utiliza el archivo `valid_points_lhe`
