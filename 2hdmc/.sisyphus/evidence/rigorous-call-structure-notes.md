# Estructura de llamadas y consistencia de fórmula

## Líneas clave en `src/THDM.cpp`

- Línea 333: `double alpha = -asin(sba)+beta;`
- Línea 341: `lambda[1]=(m_H*m_H*ca2+m_h*m_h*sa2-m12_2*tb)/v2/cb2-1.5*lambda6*tb+0.5*lambda7*tb*tb*tb;`
- Línea 384: `double m12_2 = (m_H*m_H*ca2+m_h*m_h*sa2 ... )/tb;`
- Línea 388: `return set_param_phys(m_h,m_H,m_A,m_Hp,sba,lambda6,lambda7,m12_2,tan_beta);`

## Verificación lógica

1. `set_param_phys_lam1` reutiliza la **misma convención angular** (`alpha = beta - asin(sba)`).
2. Reconstruye `m12_2` invirtiendo la ecuación de `lambda[1]` existente.
3. La inicialización final ocurre por delegación a `set_param_phys`, evitando duplicación del núcleo numérico.

Conclusión: comportamiento append-only con núcleo de cálculo preservado.
