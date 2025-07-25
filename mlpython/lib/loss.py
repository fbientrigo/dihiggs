import numpy as np

def compute_max_unitarity_eigenvalue(lambdas):
    """
    Given the tuple or list `lambdas` = (λ1, λ2, λ3, λ4, λ5, λ6, λ7),
    construct the 2HDM scattering matrices S_(2,1), S_(2,0), S_(0,1), S_(0,0)
    (each multiplied by 16π) and return the largest-magnitude eigenvalue
    among them (i.e. max |eigenvalue|).
    
    According to unitarity, each eigenvalue must satisfy |a_i| < 8π.
    We compute eigenvalues of each block and return max absolute value.
    """
    λ1, λ2, λ3, λ4, λ5, λ6, λ7 = lambdas
    sixteen_pi = 16.0 * np.pi

    # 1) S_(2,1): 3×3
    # 16π·S_(2,1) = 
    sqrt2 = np.sqrt(2.0)
    S21 = np.array([
        [λ1,           λ5,        sqrt2 * λ6],
        [λ5,           λ2,        sqrt2 * λ7],
        [sqrt2 * λ6,   sqrt2 * λ7, λ3 + λ4  ]
    ])
    # Rescale: actual matrix = S21 (here we will compare eigenvalues of 16π·S21
    # S21 *= sixteen_pi

    # 2) S_(2,0): scalar
    # 16π·S_(2,0) = λ3 - λ4  (i.e. 1×1 matrix)
    S20 = np.array([[ (λ3 - λ4) * sixteen_pi ]])

    # 3) S_(0,1): 4×4
    # 16π·S_(0,1) = 
    S01 = np.array([
        [λ1,  λ4,  λ6,  λ6],
        [λ4,  λ2,  λ7,  λ7],
        [λ6,  λ7,  λ3,  λ5],
        [λ6,  λ7,  λ5,  λ3]
    ])
    # S01 *= sixteen_pi

    # S_(0,0): 4×4
    # 16π·S_(0,0) 
    S00 = np.array([
        [3*λ1,         2*λ3 + λ4,  3*λ6,         3*λ6      ],
        [2*λ3 + λ4,    3*λ2,       3*λ7,         3*λ7      ],
        [3*λ6,         3*λ7,       λ3 + 2*λ4,    3*λ5      ],
        [3*λ6,         3*λ7,       3*λ5,         λ3 + 2*λ4 ]
    ])
    # S00 *= sixteen_pi

    # Compute eigenvalues
    # For small symmetric matrices, `np.linalg.eigvalsh` is more efficient & accurate.
    ev_S21 = np.linalg.eigvalsh(S21)
    ev_S20 = np.linalg.eigvalsh(S20)      # just one value
    ev_S01 = np.linalg.eigvalsh(S01)
    ev_S00 = np.linalg.eigvalsh(S00)

    # Collect all eigenvalues, take absolute values, and return the maximum
    all_evs = np.hstack([ev_S21, ev_S20, ev_S01, ev_S00])
    max_abs_ev = np.max(np.abs(all_evs))
    return max_abs_ev
