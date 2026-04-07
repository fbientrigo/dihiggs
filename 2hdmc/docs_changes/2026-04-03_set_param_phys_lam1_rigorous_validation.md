# Cambio científico: `set_param_phys_lam1` (entrada con \(\lambda_1\))

## 1) Qué se agregó (append-only)

Se agregó una **nueva API pública** en `THDM`:

```cpp
bool set_param_phys_lam1(double m_h,double m_H, double m_A, double m_Hp,
                         double sba, double lambda1,
                         double lambda6, double lambda7,
                         double tan_beta);
```

Esta API:
- recibe \(\lambda_1\) en lugar de \(m_{12}^2\),
- reconstruye \(m_{12}^2\) con la misma convención interna de 2HDMC,
- y **delega** a `set_param_phys(...)` para no duplicar ni alterar el núcleo.

## 2) Qué NO cambió (integridad del núcleo)

No se tocaron piezas fundamentales del cálculo:

- `set_param_phys(...)` quedó sin refactor funcional.
- No se añadieron nuevas restricciones de estabilidad/unitariedad.
- No se cambiaron convenciones angulares internas (`beta=atan(tan_beta)`, `alpha=-asin(sba)+beta`).

Dif del núcleo (fuente):
- `src/THDM.cpp`: cambio **aditivo** (inserción de método nuevo).
- `src/THDM.h`: declaración + documentación del método nuevo.

## 3) Estructura de llamadas (call structure)

Flujo implementado:

1. `set_param_phys_lam1(...)`
2. validaciones de entrada (mismas guardas base)
3. reconstrucción de \(m_{12}^2\) invirtiendo la fórmula de `lambda[1]` en `set_param_phys`
4. llamada final:

```cpp
return set_param_phys(m_h,m_H,m_A,m_Hp,sba,lambda6,lambda7,m12_2,tan_beta);
```

Esto asegura que la inicialización final del modelo pase por la ruta ya establecida y validada históricamente.

Evidencia adicional de trazabilidad de llamadas:
- `.sisyphus/evidence/rigorous-call-structure-notes.md`

## 4) Fórmula usada (consistente con el source)

\[
m_{12}^2 = \frac{m_H^2\cos^2\alpha + m_h^2\sin^2\alpha - v^2\cos^2\beta\left(\lambda_1 + \frac{3}{2}\lambda_6\tan\beta - \frac{1}{2}\lambda_7\tan^3\beta\right)}{\tan\beta}
\]

con:
- \(\beta = \operatorname{atan}(\tan\beta)\)
- \(\alpha = -\arcsin(s_{\beta-\alpha}) + \beta\)
- \(v^2 = \texttt{THDM::v2}\) (convención interna).

## 5) Verificación rigurosa ejecutada

### 5.1 Dry run de build

- `make -n all` ✅
- `make -n build` ❌ (no existe ese target en este `Makefile`)

Evidencia:
- `.sisyphus/evidence/rigorous-dryrun-all.txt`
- `.sisyphus/evidence/rigorous-dryrun-build.txt`
- `.sisyphus/evidence/rigorous-dryrun-build-exit.txt`
- `.sisyphus/evidence/rigorous-make-build.txt`
- `.sisyphus/evidence/rigorous-make-build-exit.txt`

### 5.2 Build limpio real

- `make clean && make` ✅

Evidencia:
- `.sisyphus/evidence/rigorous-make-clean.txt`
- `.sisyphus/evidence/rigorous-make.txt`
- `.sisyphus/evidence/rigorous-make-clean-rerun-exit.txt`
- `.sisyphus/evidence/rigorous-make-rerun-exit.txt`

### 5.3 Prueba determinística del validador

- `make CalcRoundTrip && ./CalcRoundTrip` ✅
- Resultado:
  - `Case generic: PASS`
  - `Case near-alignment: PASS`
  - `Case invalid-input: PASS`

Evidencia:
- `.sisyphus/evidence/rigorous-make-calcroundtrip.txt`
- `.sisyphus/evidence/rigorous-calcroundtrip.txt`
- `.sisyphus/evidence/rigorous-calcroundtrip-rerun.txt`
- `.sisyphus/evidence/rigorous-calcroundtrip-rerun-exit.txt`
- `.sisyphus/evidence/rigorous-make-calcroundtrip-rerun-exit.txt`
- `.sisyphus/evidence/rigorous-calcroundtrip-postbuild.txt`
- `.sisyphus/evidence/rigorous-calcroundtrip-postbuild-exit.txt`

### 5.4 Prueba ampliada (random scientific points)

Se compiló y ejecutó un verificador adicional con **400 puntos aleatorios válidos** + **3 casos inválidos**:

- `PASS: 400 random round-trip points validated`
- `PASS: invalid-input guards validated (3 cases)`
- `Max abs diff: 2.6193447411060333e-09`
- `Max rel diff: 2.2463382265407876e-12`

Evidencia:
- `.sisyphus/evidence/RigorousRoundTrip.cpp`
- `.sisyphus/evidence/rigorous-random-build.txt`
- `.sisyphus/evidence/rigorous-random-run.txt`
- `.sisyphus/evidence/rigorous-random-rerun.txt`
- `.sisyphus/evidence/rigorous-random-rerun-exit.txt`
- `.sisyphus/evidence/rigorous-random-postbuild.txt`
- `.sisyphus/evidence/rigorous-random-postbuild-exit.txt`
- `.sisyphus/evidence/rigorous-random-rebuild-link-exit.txt`
- `.sisyphus/evidence/rigorous-random-after-relink.txt`
- `.sisyphus/evidence/rigorous-random-after-relink-exit.txt`

## 6) Reproducibilidad (comandos)

```bash
make -n all
make -n build

make clean
make
make CalcRoundTrip
./CalcRoundTrip

g++ .sisyphus/evidence/RigorousRoundTrip.cpp -Isrc -std=c++11 -Wall -O3 \
    -Llib -l2HDMC -lgsl -lgslcblas -lm -o .sisyphus/evidence/RigorousRoundTrip
.sisyphus/evidence/RigorousRoundTrip
```

## 7) Conclusión

El cambio se mantiene en modalidad **append-only**: se añade capacidad nueva sin alterar el flujo de cálculo base.
La evidencia numérica y de compilación confirma consistencia lógica y estabilidad para uso científico.
