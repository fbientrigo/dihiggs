import numpy as np
import matplotlib.pyplot as plt
from typing import Iterable, Dict, List, Generator, Union
from .oracle import safe_run_oracle

# Parámetros globales por defecto (pueden modificarse según conveniencia)
lambda6 = 0.1
lambda7 = 0.0
sin_ba_default = 1.0
tan_beta_default = 1e7


def physpoint(
    mphi: np.float64,
    m12_2: np.float64,
    mA: np.float64 = 300.0,
    sin_ba: np.float64 = sin_ba_default,
    tan_beta: np.float64 = tan_beta_default
) -> Union[Dict, None]:
    """
    Envía un punto al “oracle” físico con parámetros fijos y devuelve los
    resultados si la ejecución es exitosa.

    Parámetros:
        mphi (np.float64):
            Masa del estado φ (m_phi) en GeV.
        m12_2 (np.float64):
            Parámtero m12^2 en GeV^2.
        mA (np.float64, opcional):
            Masa del estado A (m_A) en GeV. Por defecto 300.0 GeV.
        sin_ba (np.float64, opcional):
            Valor de sin(β - α). Por defecto, `sin_ba_default`.
        tan_beta (np.float64, opcional):
            Valor de tan(β). Por defecto, `tan_beta_default`.

    Retorna:
        dict o None:
            Un diccionario con resultados del oracle (branching ratios,
            condiciones de chequeo, anchos parciales, etc.), o None si la
            ejecución falla internamente.
    """
    # Construir vector de parámetros en el orden esperado por el oracle:
    # [m_phi, m_A, sin(β - α), tan(β), λ6, λ7, m12^2]
    param_vector = [mphi, mA, sin_ba, tan_beta, lambda6, lambda7, m12_2]
    return safe_run_oracle(param_vector)


def simulate_curve(
    mphi_values: Iterable[np.float64],
    mA: float,
    m12: float,
    lambda6_val: float,
    lambda7_val: float,
    sin_ba: float,
    tan_beta: float
) -> List[Dict]:
    """
    Genera una curva de branching ratios y anchos totales para múltiples valores de m_phi,
    ejecutando `physpoint` en cada punto y recopilando resultados.

    Parámetros:
        mphi_values (Iterable[np.float64]):
            Lista o array de valores de m_phi a evaluar.
        mA (float):
            Masa del estado A (m_A) en GeV (fijo para toda la curva).
        m12 (float):
            Valor de m12^2 en GeV^2 (fijo para toda la curva).
        lambda6_val (float):
            Valor de λ6 a usar (sobrescribe el global).
        lambda7_val (float):
            Valor de λ7 a usar (sobrescribe el global).
        sin_ba (float):
            Valor de sin(β - α) (fijo para toda la curva).
        tan_beta (float):
            Valor de tan(β) (fijo para toda la curva).

    Retorna:
        List[Dict]:
            Lista de diccionarios, uno por cada m_phi, con la forma:
            {
                'mphi': valor de m_phi,
                'br': {
                    'gaga': BR(φ → γγ),
                    'Zga': BR(φ → Zγ),
                    'bb':  BR(φ → bᵦ)
                },
                'total_width': ancho total Γ(φ) en GeV,
                'flags': {
                    'positivity': 0 o 1,
                    'unitarity':  0 o 1,
                    'perturbativity': 0 o 1
                },
                'params': {
                    'mA':      mA,
                    'm12':     m12,
                    'lambda6': lambda6_val,
                    'lambda7': lambda7_val,
                    'sin_ba':  sin_ba,
                    'tan_beta': tan_beta
                }
            }
            Sólo se añade a la lista si `physpoint` devuelve un dict no nulo.
    """
    results: List[Dict] = []

    # Actualizar valores globales de λ6 y λ7 antes de simular
    global lambda6, lambda7
    lambda6 = lambda6_val
    lambda7 = lambda7_val

    for mphi in mphi_values:
        out = physpoint(
            np.float64(mphi),
            np.float64(m12),
            np.float64(mA),
            np.float64(sin_ba),
            np.float64(tan_beta)
        )
        if out is None:
            # Si el oracle devolvió None (ejecución fallida), se omite este punto
            continue

        # Extraer branching ratios; usar NaN o 0 si la clave no existe
        br_gaga = out.get('branching_ratio_h2_gaga', np.nan)
        total_width = out.get('w_total_h2', np.nan)
        br_Zga = (out.get('w_h2_Zga', 0) / total_width) if total_width and total_width != 0 else 0.0
        br_bb = (out.get('w_h2_bb', 0) / total_width) if total_width and total_width != 0 else 0.0

        results.append({
            'mphi': mphi,
            'br': {
                'gaga': br_gaga,
                'Zga': br_Zga,
                'bb':  br_bb
            },
            'total_width': total_width,
            'flags': {
                'positivity':     out.get('positivity_ok', 0),
                'unitarity':      out.get('unitarity_ok', 0),
                'perturbativity': out.get('perturbativity_ok', 0),
            },
            'params': {
                'mA':       mA,
                'm12':      m12,
                'lambda6':  lambda6_val,
                'lambda7':  lambda7_val,
                'sin_ba':   sin_ba,
                'tan_beta': tan_beta
            }
        })

    return results


def make_title(params: Dict[str, float], subplot_tag: str = "(d)") -> str:
    """
    Construye el título del gráfico a partir de los parámetros físicos,
    formateando expresiones LaTeX para mostrar en matplotlib.

    Parámetros:
        params (Dict[str, float]):
            Diccionario con las claves:
                - 'mA'       : Masa m_A en GeV
                - 'm12'      : Valor de m12^2 en GeV^2
                - 'lambda6'  : Valor de λ6
                - 'lambda7'  : Valor de λ7
                - 'sin_ba'   : Valor de sin(β - α)
                - 'tan_beta' : Valor de tan(β)
        subplot_tag (str, opcional):
            Etiqueta corta para indicar subpanel en figura, p.ej. "(a)", "(b)", etc.

    Retorna:
        str:
            Cadena formateada con expresiones LaTeX y valores numéricos
            listos para asignar como título en matplotlib.
    """
    mA_val = params['mA']
    m12_val = params['m12']
    l6 = params['lambda6']
    l7 = params['lambda7']
    sinba = params['sin_ba']
    tanbeta = params['tan_beta']

    title = (
        rf"$\sin(\beta - \alpha) = " + f"{sinba:.2e}" +
        r" \quad \epsilon_A = \frac{1}{" + f"{tanbeta:.2e}" + "}" +
        rf"\quad \Delta\lambda' = {l6:.2f} " +
        rf"\quad m_\phi < m_{{H^\pm}} = m_A = {mA_val:.0f}\,$" +
        f"\n{subplot_tag}"
    )
    return title


def plot_branchings(results: List[Dict], subplot_tag: str = "(d)") -> None:
    """
    Grafica las curvas de branching ratios (φ → γγ, φ → Zγ, φ → b b̄) junto con
    marcadores que resaltan los puntos donde se cumplen las condiciones físicas
    (positivity, unitarity, perturbativity). Además, incluye un subgráfico
    con la evolución del ancho total Γ(φ) en escala log.

    Parámetros:
        results (List[Dict]):
            Lista de diccionarios generados por `simulate_curve`, cada uno con:
                - 'mphi'        : valor de m_phi
                - 'br'          : dict con keys 'gaga', 'Zga', 'bb'
                - 'total_width' : ancho total Γ(φ)
                - 'flags'       : dict con keys 'positivity', 'unitarity', 'perturbativity'
                - 'params'      : diccionario de parámetros (usado para el título)
        subplot_tag (str, opcional):
            Etiqueta corta para indicar subpanel en la figura.

    Retorna:
        None. Muestra la figura con matplotlib.
    """
    # Extraer arrays de m_phi y branching ratios
    mphi_list = [r['mphi'] for r in results]
    brs = {
        'gaga': [r['br']['gaga'] for r in results],
        'Zga':  [r['br']['Zga']  for r in results],
        'bb':   [r['br']['bb']   for r in results]
    }

    fig, ax = plt.subplots(figsize=(5, 5))

    # 1) Marcar puntos válidos donde se cumplen las 3 condiciones físicas
    for channel, color in zip(['gaga', 'Zga', 'bb'], ['r', 'purple', 'b']):
        valid_x = []
        valid_y = []
        for r in results:
            flags = r['flags']
            if flags['positivity'] == 1 and flags['unitarity'] == 1 and flags['perturbativity'] == 1:
                valid_x.append(r['mphi'])
                valid_y.append(r['br'][channel])
        ax.scatter(
            valid_x,
            valid_y,
            color=color,
            marker='o',
            edgecolors='k',
            zorder=5,
            s=30
        )

    # 2) Trazar curvas continuas de BR para cada canal
    ax.plot(mphi_list, brs['gaga'], 'r--',    label=r'$\phi \to \gamma\gamma$')
    ax.plot(mphi_list, brs['Zga'],  'purple', label=r'$\phi \to Z\gamma$')
    ax.plot(mphi_list, brs['bb'],   'b:',     label=r'$\phi \to b\bar{b}$')

    # Configuración de ejes principales
    ax.set_ylim(0, 1)
    ax.set_xlim(min(mphi_list), max(mphi_list))
    ax.set_xlabel(r'$m_\phi$ [GeV]')
    ax.set_ylabel('Branching Ratio')
    ax.legend(loc='upper right')

    # 3) Subgráfico (inset) para ancho total Γ(φ)
    total_widths = [r['total_width'] for r in results]
    inset_ax = ax.inset_axes([0.15, 0.70, 0.35, 0.25])
    inset_ax.plot(mphi_list, total_widths, 'k:', linewidth=1)
    inset_ax.set_yscale('log')

    # Elegir ticks uniformes en m_phi
    x_ticks = np.linspace(min(mphi_list), max(mphi_list), num=4)
    inset_ax.set_xticks(x_ticks)
    inset_ax.set_xticklabels([f"{x:.0f}" for x in x_ticks], fontsize=7)

    # Ticks logarítmicos para ancho total
    min_tw = min([tw for tw in total_widths if tw > 0], default=1e-3)
    max_tw = max(total_widths)
    y_ticks = np.logspace(np.floor(np.log10(min_tw)), np.ceil(np.log10(max_tw)), num=3)
    inset_ax.set_yticks(y_ticks)
    inset_ax.set_yticklabels([f"{y:.1e}" for y in y_ticks], fontsize=7)

    inset_ax.set_title(r'$\Gamma$ [GeV]', fontsize=8)

    # 4) Título principal con parámetros físicos
    title_str = make_title(results[0]['params'], subplot_tag)
    ax.set_title(title_str, fontsize=11)

    plt.tight_layout()
    plt.show()
