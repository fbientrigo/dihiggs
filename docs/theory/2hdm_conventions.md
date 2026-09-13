# Canonical 2HDM conventions

This page is a compact code-facing summary of the project's canonical
theoretical appendix. When a derivation or sign is in doubt, the canonical
appendix is the primary source of truth; the repository physics authority is
the machine-facing contract for production metadata.

## Generic basis

The project uses two hypercharge-\(1/2\) scalar doublets with

\[
v^2=v_1^2+v_2^2,\qquad
\tan\beta=\frac{v_2}{v_1},\qquad
v_1=v c_\beta,\qquad v_2=v s_\beta .
\]

For the CP-conserving branch,

\[
\begin{aligned}
V={}&m_{11}^2\Phi_1^\dagger\Phi_1
+m_{22}^2\Phi_2^\dagger\Phi_2
-\left[m_{12}^2\Phi_1^\dagger\Phi_2+\mathrm{h.c.}\right] \\
&+\frac12\lambda_1(\Phi_1^\dagger\Phi_1)^2
+\frac12\lambda_2(\Phi_2^\dagger\Phi_2)^2
+\lambda_3(\Phi_1^\dagger\Phi_1)(\Phi_2^\dagger\Phi_2)\\
&+\lambda_4(\Phi_1^\dagger\Phi_2)(\Phi_2^\dagger\Phi_1)
+\left\{\frac12\lambda_5(\Phi_1^\dagger\Phi_2)^2
+\left[\lambda_6(\Phi_1^\dagger\Phi_1)+
\lambda_7(\Phi_2^\dagger\Phi_2)\right]\Phi_1^\dagger\Phi_2
+\mathrm{h.c.}\right\}.
\end{aligned}
\]

Only after fixing this quadratic convention define

\[
\boxed{M^2=\frac{m_{12}^2}{s_\beta c_\beta}}.
\]

Therefore `M2`, `m12_sq`, and `m22_sq` are distinct quantities.

## CP-even state convention and exact alignment

The Davidson--Haber sign convention is

\[
h=-s_\alpha\rho_1+c_\alpha\rho_2,\qquad
\phi\equiv H=c_\alpha\rho_1+s_\alpha\rho_2.
\]

With

\[
\rho_v=c_\beta\rho_1+s_\beta\rho_2,\qquad
\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2,
\]

exact alignment \(s_{\beta-\alpha}=1\) gives

\[
\boxed{h=\rho_v,\qquad \phi=-\rho_\perp}.
\]

The minus sign on \(\phi\) is part of the convention and must be retained in
odd-\(\phi\) interactions.

## Higgs basis

The fixed Higgs-basis rotation is

\[
\boxed{H_1=c_\beta\Phi_1+s_\beta\Phi_2,\qquad
H_2=-s_\beta\Phi_1+c_\beta\Phi_2}.
\]

Stationarity gives

\[
Y_1=-\frac12 Z_1v^2,\qquad
Y_3=-\frac12 Z_6v^2,
\]

and the CP-even Higgs-basis mass matrix is

\[
\mathcal M^2_{\rm even,HB}=
\begin{pmatrix}
Z_1v^2 & Z_6v^2\\
Z_6v^2 & m_A^2+Z_5v^2
\end{pmatrix}.
\]

Away from an exactly degenerate CP-even eigenspace,

\[
\boxed{\text{exact alignment}\Longleftrightarrow Z_6=0}.
\]

The common quartic shift is

\[
\Delta_Z=
\frac14\sin^22\beta(\lambda_1+\lambda_2-2\lambda_{345})
-\sin2\beta\cos2\beta(\lambda_6-\lambda_7),
\]

so

\[
\boxed{Z_3=\lambda_3+\Delta_Z,\quad
Z_4=\lambda_4+\Delta_Z,\quad
Z_5=\lambda_5+\Delta_Z}.
\]

## \(Y_2\) and the physical \(h\phi\phi\) coupling

On the exact-alignment branch,

\[
\boxed{Y_2=M^2-\frac12m_h^2}.
\]

The physical trilinear is defined by

\[
\mathcal L_{\rm int}\supset
-\frac12 g_{h\phi\phi}^{\rm phys}h\phi^2,
\]

with

\[
\boxed{g_{h\phi\phi}^{\rm phys}
=v(Z_3+Z_4+Z_5)
=\frac{m_h^2+2m_\phi^2-2M^2}{v}},
\]

and Feynman rule \(-i g_{h\phi\phi}^{\rm phys}\).

The almost-inert \(m_2^2\) notation used by Gao--Neill must not be identified by
name with the generic-basis coefficient \(m_{22}^2\).
