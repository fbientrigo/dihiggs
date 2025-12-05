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

    Parameters
    ----------
    mphi : np.float64
        Masa del estado φ (m_phi) en GeV.
    m12_2 : np.float64
        Parámtero m12^2 en GeV^2.
    mA : np.float64, optional
        Masa del estado A (m_A) en GeV. Por defecto 300.0 GeV.
    sin_ba : np.float64, optional
        Valor de sin(β - α). Por defecto, `sin_ba_default`.
    tan_beta : np.float64, optional
        Valor de tan(β). Por defecto, `tan_beta_default`.

    Returns
    -------
    dict or None
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

    Parameters
    ----------
    mphi_values : Iterable[np.float64]
        Lista o array de valores de m_phi a evaluar.
    mA : float
        Masa del estado A (m_A) en GeV (fijo para toda la curva).
    m12 : float
        Valor de m12^2 en GeV^2 (fijo para toda la curva).
    lambda6_val : float
        Valor de λ6 a usar (sobrescribe el global).
    lambda7_val : float
        Valor de λ7 a usar (sobrescribe el global).
    sin_ba : float
        Valor de sin(β - α) (fijo para toda la curva).
    tan_beta : float
        Valor de tan(β) (fijo para toda la curva).

    Returns
    -------
    List[Dict]
        Lista de diccionarios, uno por cada m_phi, con la forma:
        
        .. code-block:: text

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

    Parameters
    ----------
    params : Dict[str, float]
        Diccionario con las claves:
            - 'mA'       : Masa m_A en GeV
            - 'm12'      : Valor de m12^2 en GeV^2
            - 'lambda6'  : Valor de λ6
            - 'lambda7'  : Valor de λ7
            - 'sin_ba'   : Valor de sin(β - α)
            - 'tan_beta' : Valor de tan(β)
    subplot_tag : str, optional
        Etiqueta corta para indicar subpanel en figura, p.ej. "(a)", "(b)", etc.

    Returns
    -------
    str
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

    Parameters
    ----------
    results : List[Dict]
        Lista de diccionarios generados por `simulate_curve`, cada uno con:
            - 'mphi'        : valor de m_phi
            - 'br'          : dict con keys 'gaga', 'Zga', 'bb'
            - 'total_width' : ancho total Γ(φ)
            - 'flags'       : dict con keys 'positivity', 'unitarity', 'perturbativity'
            - 'params'      : diccionario de parámetros (usado para el título)
    subplot_tag : str, optional
        Etiqueta corta para indicar subpanel en la figura.

    Returns
    -------
    None
        Muestra la figura con matplotlib.
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


# get format for plot_branchings_and_lifetime
def get_branching_format_from_df(df_physical,
    total_width_h2_name = "total_width",
    branching_ratio_h2_gaga_name = "br_gaga",
    width_h2_Zga_name = "width_Zga",
    width_h2_bb_name = "width_bb",
    m12_name = "m12"):

    # -- internal ----
    fixed_mA = df_physical["mA"].unique()
    assert len(fixed_mA) == 1, "There should be only one unique value of mA in the DF"
    fixed_mA = fixed_mA[0]  # Get the single value
    fixed_lambda6 = df_physical["lambda6"].unique()
    assert len(fixed_lambda6) == 1, "There should be only one unique"
    fixed_lambda6 = fixed_lambda6[0]  # Get the single value
    fixed_lambda7 = df_physical["lambda7"].unique()
    assert len(fixed_lambda7) == 1, "There should be only one unique"
    fixed_lambda7 = fixed_lambda7[0]  # Get the single value
    fixed_sin_ba = df_physical["sin_ba"].unique()
    assert len(fixed_sin_ba) == 1, "There should be only one unique value of sin_ba in the DF"
    fixed_sin_ba = fixed_sin_ba[0]  # Get the single value
    fixed_tan_beta = df_physical["tan_beta"].unique()
    assert len(fixed_tan_beta) == 1, "There should be only one unique value of tan_beta in the DF"
    fixed_tan_beta = fixed_tan_beta[0]  # Get the single value


    # ---  Ordenar por m_phi y extraer mphi / m12_2 ---
    df_physical = df_physical.sort_values(by='m_phi').reset_index(drop=True)

    # Vector de valores de m_phi (ordenados ascendentemente)
    mphi_values = df_physical['m_phi'].to_numpy()

    # Vector de valores de m12_2 correspondientes a cada m_phi
    m12_values  = df_physical[m12_name].to_numpy()





    # --- Paso  Empaquetar cada fila en el formato esperado por plot_branchings ---
    results = []
    for idx, row in df_physical.iterrows():
        # Extraer anchos y branching ratios
        total_w = row[total_width_h2_name]
        # Si total_width es cero o NaN, evitar división por cero
        if total_w is None or total_w == 0 or np.isnan(total_w):
            br_Zga = 0.0
            br_bb  = 0.0
        else:
            br_Zga = row[width_h2_Zga_name] / total_w
            br_bb  = row[width_h2_bb_name] / total_w

        entry = {
            'mphi': float(row['m_phi']),
            'br': {
                'gaga': float(row[branching_ratio_h2_gaga_name]),
                'Zga':  float(br_Zga),
                'bb':   float(br_bb)
            },
            'total_width': float(total_w),
            'flags': {
                'positivity':     int(row['positivity_ok']),
                'unitarity':      int(row['unitarity_ok']),
                'perturbativity': int(row['perturbativity_ok'])
            },
            'params': {
                'mA':       fixed_mA,
                'm12':      float(row[m12_name]),
                'lambda6':  fixed_lambda6,
                'lambda7':  fixed_lambda7,
                'sin_ba':   fixed_sin_ba,
                'tan_beta': fixed_tan_beta
            }
        }
        results.append(entry)
    return results



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def plot_branchings_and_lifetime(
    results: List[Dict],
    subplot_tag: str = "(d)",
    lifetime_unit: str = "mm",     # "mm", "cm" o "m"
    lifetime_log: bool = True,     # True: escala log en vida media
    width_log: bool = True,         # True: escala log en decay width
    br_log: bool = False        # True: escala log en branching ratios
):
    """
    Plotea:
     - Branching ratios (izq)
     - Decay width Γ(φ) (arriba derecha)
     - Vida media cτ = ħc/Γ (abajo derecha)

    Parámetros adicionales:
      lifetime_unit: unidad para cτ ("mm","cm","m")
      lifetime_log: si True, y-axis de cτ en log
      width_log:    si True, y-axis de Γ en log
      br_log:       si True, y-axis de BR en log
    """
    # --- Extraer datos ----
    mphi       = np.array([r['mphi'] for r in results])
    brs        = {ch: np.array([r['br'][ch] for r in results])
                  for ch in ['gaga','Zga','bb']}
    total_width= np.array([r['total_width'] for r in results])

    # --- Vida media base en mm ---
    hbar = 6.582119569e-25    # GeV·s
    c    = 2.99792458e11      # mm/s
    c_tau_mm = hbar * c / total_width  # mm

    # --- Conversión de unidades ---
    unit_factors = {"mm": 1.0, "cm": 0.1, "m": 1e-3}
    if lifetime_unit not in unit_factors:
        raise ValueError(f"Unknown lifetime_unit '{lifetime_unit}'")
    c_tau = c_tau_mm * unit_factors[lifetime_unit]
    y_label_tau = rf"$c\tau\,[{lifetime_unit}]$"

    # --- Preparar figura ---
    global_title = make_title(results[0]['params'])
    fig = plt.figure(figsize=(12, 6))
    fig.suptitle(global_title, fontsize=14, y=0.98)
    plt.subplots_adjust(top=0.85)

    gs = GridSpec(2, 2, figure=fig,
                  width_ratios=[1,1], height_ratios=[1,1],
                  wspace=0.3, hspace=0.4)

    # 1) Branching ratios (col 0)
    ax_br = fig.add_subplot(gs[:,0])
    styles = {'gaga':'--','Zga':'-','bb':':'}
    colors = {'gaga':'r','Zga':'purple','bb':'b'}
    labels = {
        'gaga': r'$\phi\to \gamma\gamma$',
        'Zga':  r'$\phi\to Z\gamma$',
        'bb':   r'$\phi\to b\bar{b}$'
    }
    for ch in ['gaga','Zga','bb']:
        ax_br.plot(
            mphi, brs[ch],
            linestyle=styles[ch],
            color=colors[ch],
            label=labels[ch]
        )
        valid = [(r['mphi'], r['br'][ch]) for r in results
                 if all(r['flags'][f]==1 for f in ['positivity','unitarity','perturbativity'])]
        if valid:
            vx, vy = zip(*valid)
            ax_br.scatter(vx, vy,
                          color=colors[ch],
                          edgecolors='k',
                          s=30, zorder=5)
    ax_br.set_xlabel(r'$m_\phi\,[\mathrm{GeV}]$')
    ax_br.set_ylabel('Branching Ratio')
    ax_br.legend(loc='upper left', fontsize=8)
    ax_br.set_title(subplot_tag, loc='left', fontsize=12)
    if br_log:
        ax_br.set_yscale('log')
    else:
        ax_br.set_ylim(0,1)


    # 2) Decay width (fila 0,col 1)
    ax_w = fig.add_subplot(gs[0,1])
    ax_w.plot(mphi, total_width, linestyle='-', color='k', linewidth=1.2)
    if width_log:
        ax_w.set_yscale('log')
    ax_w.set_xlabel(r'$m_\phi\,[\mathrm{GeV}]$')
    ax_w.set_ylabel(r'$\Gamma(\phi)\,[\mathrm{GeV}]$')
    ax_w.set_title('Decay Width', fontsize=10)

    # 3) Vida media (fila 1,col 1)
    ax_tau = fig.add_subplot(gs[1,1])
    ax_tau.plot(mphi, c_tau, linestyle='-', color='g', linewidth=1.2)
    if lifetime_log:
        ax_tau.set_yscale('log')
    ax_tau.set_xlabel(r'$m_\phi\,[\mathrm{GeV}]$')
    ax_tau.set_ylabel(y_label_tau)
    ax_tau.set_title(r'Lifetime $c\tau = \hbar c / \Gamma$', fontsize=10)

    plt.show()