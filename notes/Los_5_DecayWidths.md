# 5 Decay Widths

Vamos a explorar las expresiones de la calculadora para así obtener los decaimientos más importantes que aparecen en los calculos,

- h2_b_b
- h2_ta_ta
- h2_ga_ga
- h2_g_g
- h2_Z_ga

## Decaimiento hacia 2 higgs desaparece

En el límite de alineamiento ($\sin(\beta-\alpha) \to 1$), el acoplamiento se vuelve proporcional puramente al parámetro de violación suave de la simetría $Z_2$ en la base de Higgs, dado por $Z_6$:
$$
C_{Hhh}^{\text{align}} = 3 v Z_6
$$


El código implementa $Z_6$ mediante la rotación de los parámetros de entrada:
$$
Z_6 = -\frac{1}{2} \sin(2\beta) \left[ \lambda_1 c_\beta^2 - \lambda_2 s_\beta^2 - \lambda_{345} \cos(2\beta) \right]
+ \cos\beta \cos(3\beta) \lambda_6 + \sin\beta \sin(3\beta) \lambda_7
$$

### Comportamiento asintótico a gran $\tan\beta$

Analizamos la supresión observada numéricamente cuando $\tan\beta \gg 1$. Definiendo $t_\beta \equiv \tan\beta$, tenemos las aproximaciones:
$$
    \sin\beta \approx 1, \quad \cos\beta \approx \frac{1}{t_\beta}, \quad \sin(2\beta) \approx \frac{2}{t_\beta}, \quad \cos(3\beta) \approx \mathcal{O}\left(\frac{1}{t_\beta}\right)
$$

Bajo la condición de entrada $\lambda_7 = 0$ y masa finita para el Higgs pesado, ocurre una cancelación dinámica ("screening") en el término entre corchetes, dejando el comportamiento dominante gobernado por el término explícito de $\lambda_6$. El acoplamiento efectivo escala como:

$$
    Z_6 \xrightarrow{t_\beta \to \infty} \lambda_6 \left( \frac{1}{t_\beta} \right) \left( \frac{3}{t_\beta} \right) \approx \frac{3\lambda_6}{t_\beta^2}
$$

Por lo tanto, la anchura de desintegración experimenta una supresión cuártica respecto a $\tan\beta$:
$$
    \Gamma(H \to hh) \propto |C_{Hhh}|^2 \propto |Z_6|^2 \sim \frac{\lambda_6^2}{\tan^4\beta}
$$

Esto explica el resultado numérico nulo ($\Gamma \sim 10^{-32}$ GeV) obtenido en simulaciones con $\tan\beta = 10^4$ y $\lambda_6 = 0.1$. A pesar de tener un acoplamiento "no nulo" en el lagrangiano, la geometría de la base de Higgs suprime este canal eficazmente, cerrando el espacio de fase dinámico para la producción de pares de Higgs ligeros.

