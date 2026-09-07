# 06 — Alignment and gauge couplings from the scalar kinetic terms

Scope: Phase 6 of issue #81.  
Input convention: the scalar fields and VEVs fixed in Phases 1–4 and the CP-even mass-eigenstate signs frozen in Phase 4.  
Primary source used **after** the independent derivation: Davidson–Haber (DH05), especially the gauge-coupling table `littletable`; BFLRS11 gauge-coupling discussion; GHOO18 exact-alignment gauge projection.  
Implementation cross-check used only after the analytic result: active 2HDMC `THDM::get_qki` and `THDM::get_coupling_vvh`.

This phase does not use the scalar potential. It derives the tree-level neutral-scalar couplings to `W^+W^-` and `ZZ` entirely from gauge invariance, canonical kinetic terms and the already-fixed scalar mixing convention.

The fields are

\[
\Phi_i=
\begin{pmatrix}
\phi_i^+\\[1mm]
\dfrac{v_i+\rho_i+i\eta_i}{\sqrt2}
\end{pmatrix},
\qquad
v_1=v c_\beta,\quad v_2=v s_\beta,
\qquad
v^2=v_1^2+v_2^2.
\]

The CP-even mass eigenstates are

\[
\boxed{
\begin{pmatrix}h\\\phi\end{pmatrix}
=
\begin{pmatrix}
-s_\alpha & c_\alpha\\
c_\alpha & s_\alpha
\end{pmatrix}
\begin{pmatrix}\rho_1\\\rho_2\end{pmatrix}}
\]

and therefore

\[
\boxed{\rho_1=-s_\alpha h+c_\alpha\phi},
\qquad
\boxed{\rho_2=c_\alpha h+s_\alpha\phi}.
\]

---

## What is established

[DERIVED] The gauge-boson masses depend only on the VEV norm

\[
\boxed{v^2=v_1^2+v_2^2},
\]

with

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},
\qquad
\boxed{m_Z^2=\frac{(g^2+g'^2)v^2}{4}}.
\]

[DERIVED] The interactions linear in a CP-even neutral fluctuation and quadratic in massive gauge fields depend only on

\[
\boxed{v_1\rho_1+v_2\rho_2=v\rho_v},
\qquad
\boxed{\rho_v\equiv c_\beta\rho_1+s_\beta\rho_2}.
\]

Thus the CP-even direction parallel to the VEV is the unique neutral scalar direction with a tree-level linear `VV` coupling.

The orthogonal direction

\[
\boxed{\rho_\perp=-s_\beta\rho_1+c_\beta\rho_2}
\]

has no tree-level interaction of the form `rho_perp W^+W^-` or `rho_perp ZZ`.

[DERIVED] In the mass basis,

\[
\boxed{\rho_v=s_{\beta-\alpha}h+c_{\beta-\alpha}\phi}.
\]

Defining the Lagrangian modifiers by

\[
\boxed{
\mathcal L_{SVV}
=\kappa_V^S\left[
\frac{2m_W^2}{v}S W_\mu^+W^{-\mu}
+\frac{m_Z^2}{v}S Z_\mu Z^\mu
\right],
}
\]

we obtain

\[
\boxed{\kappa_V^h=\sin(\beta-\alpha)},
\qquad
\boxed{\kappa_V^\phi=\cos(\beta-\alpha)}.
\]

On the project exact-alignment branch,

\[
s_{\beta-\alpha}=+1,\qquad c_{\beta-\alpha}=0,
\]

so

\[
\boxed{\kappa_V^h=1},
\qquad
\boxed{\kappa_V^\phi=0}.
\]

Therefore claim C4 is **VERIFIED** in the frozen project convention.

[DERIVED] The CP-odd fluctuations enter the gauge-boson mass operator only quadratically. Hence the physical CP-odd scalar `A` has no tree-level vertex of the form `A W^+W^-` or `A Z Z`.

[DERIVED] The CP-even quadratic gauge interaction satisfies

\[
\rho_1^2+\rho_2^2=h^2+\phi^2,
\]

so the orthogonal CP-even rotation produces no `h phi VV` cross term from this norm. This is an independent structural check of the rotation.

---

## Derivation

### 1. Covariant derivative and hypercharge convention

Use the common modern normalization

\[
Q=T^3+Y,
\qquad
Y(\Phi_i)=\frac12,
\]

so

\[
\boxed{
D_\mu=\partial_\mu+i g\frac{\sigma^a}{2}W_\mu^a+i g'\frac12 B_\mu.
}
\]

This is equivalent to the DH05 convention `Q=T3+Y/2` with scalar-doublet hypercharge `Y=1`; only the hypercharge normalization differs.

Define

\[
W_\mu^\pm\equiv\frac{W_\mu^1\mp iW_\mu^2}{\sqrt2},
\qquad
g_Z\equiv\sqrt{g^2+g'^2},
\]

and

\[
Z_\mu\equiv\frac{gW_\mu^3-g'B_\mu}{g_Z}.
\]

For the neutral lower component, write

\[
\Phi_i^0=\frac{x_i}{\sqrt2},
\qquad
x_i\equiv v_i+\rho_i+i\eta_i.
\]

The gauge part of the covariant derivative acting on the neutral configuration is

\[
D_\mu\Phi_i\supset
\begin{pmatrix}
 i\,\dfrac{g}{2}W_\mu^+ x_i\\[2mm]
 -i\,\dfrac{gW_\mu^3-g'B_\mu}{2\sqrt2}x_i
\end{pmatrix}.
\]

The relative signs are irrelevant after taking the norm, but the factors of `1/2` and `1/sqrt(2)` are not.

### 2. Gauge-boson terms from one doublet

Squaring the charged upper component gives

\[
\boxed{
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset
\frac{g^2}{4}|x_i|^2 W_\mu^+W^{-\mu}.
}
\]

Since

\[
|x_i|^2=(v_i+\rho_i)^2+\eta_i^2,
\]

we have

\[
\frac{g^2}{4}\left[v_i^2+2v_i\rho_i+\rho_i^2+\eta_i^2\right]W^+W^-.
\]

For the neutral gauge field,

\[
\boxed{
(D_\mu\Phi_i)^\dagger(D^\mu\Phi_i)
\supset
\frac{g_Z^2}{8}|x_i|^2 Z_\mu Z^\mu.
}
\]

Thus

\[
\frac{g_Z^2}{8}\left[v_i^2+2v_i\rho_i+\rho_i^2+\eta_i^2\right]Z^2.
\]

The photon does not appear because the neutral VEV component has electric charge zero.

### 3. Sum over the two doublets: masses and linear interactions

Summing `i=1,2`, the mass terms are

\[
\mathcal L_{\rm mass}
\supset
\frac{g^2}{4}(v_1^2+v_2^2)W^+W^-
+\frac{g_Z^2}{8}(v_1^2+v_2^2)Z^2.
\]

Comparing with

\[
\mathcal L_{\rm mass}
\supset m_W^2W^+W^-+\frac12m_Z^2Z^2
\]

gives

\[
\boxed{m_W^2=\frac{g^2v^2}{4}},
\qquad
\boxed{m_Z^2=\frac{g_Z^2v^2}{4}}.
\]

The terms linear in `rho_i` are

\[
\boxed{
\mathcal L_{VV,\rm linear}
=
(v_1\rho_1+v_2\rho_2)
\left[
\frac{g^2}{2}W^+W^-+\frac{g_Z^2}{4}ZZ
\right].
}
\]

Using

\[
v_1\rho_1+v_2\rho_2=v\rho_v
\]

and the mass relations,

\[
\boxed{
\mathcal L_{VV,\rm linear}
=
\rho_v\left[
\frac{2m_W^2}{v}W^+W^-+
\frac{m_Z^2}{v}ZZ
\right].
}
\]

### 4. Why the `ZZ` coefficient in the Lagrangian differs by a factor two from the vertex

For `W^+W^-`, the two vector fields are distinct, so

\[
\mathcal L\supset \frac{2m_W^2}{v}S W^+W^-
\]

corresponds directly to the Feynman rule

\[
\boxed{i\frac{2m_W^2}{v}g_{\mu\nu}}.
\]

For two identical `Z` fields,

\[
\mathcal L\supset \frac{m_Z^2}{v}S Z_\mu Z^\mu
\]

has a factor of two upon functional differentiation with respect to the two identical `Z` fields, giving

\[
\boxed{i\frac{2m_Z^2}{v}g_{\mu\nu}}.
\]

Thus the physical SM-normalized modifier is the same for `WW` and `ZZ`. Confusing the Lagrangian coefficient with the vertex coefficient would create a spurious factor-of-two discrepancy in the neutral channel.

### 5. Projection onto the mass eigenstates

From Phase 4,

\[
\rho_1=-s_\alpha h+c_\alpha\phi,
\qquad
\rho_2=c_\alpha h+s_\alpha\phi.
\]

Therefore

\[
\begin{aligned}
\rho_v
&=c_\beta\rho_1+s_\beta\rho_2\\
&=(-c_\beta s_\alpha+s_\beta c_\alpha)h
 +(c_\beta c_\alpha+s_\beta s_\alpha)\phi\\
&=\sin(\beta-\alpha)h+\cos(\beta-\alpha)\phi.
\end{aligned}
\]

Substitution gives

\[
\boxed{
\mathcal L_{VV,\rm linear}
=
\left[s_{\beta-\alpha}h+c_{\beta-\alpha}\phi\right]
\left[
\frac{2m_W^2}{v}W^+W^-+
\frac{m_Z^2}{v}ZZ
\right].
}
\]

Hence

\[
\boxed{\kappa_V^h=s_{\beta-\alpha}},
\qquad
\boxed{\kappa_V^\phi=c_{\beta-\alpha}}.
\]

### 6. Exact alignment acquires its gauge-coupling meaning

Phase 4 had established only the mixing geometry. On the project branch,

\[
s_{\beta-\alpha}=1
\quad\Longrightarrow\quad
h=\rho_v,
\qquad
\phi=-\rho_\perp.
\]

The kinetic-term derivation now establishes the physical content:

\[
\boxed{h=\rho_v\Longrightarrow hVV\text{ is exactly SM-like at tree level}},
\]

\[
\boxed{\phi\perp\rho_v\Longrightarrow \phi VV=0\text{ at tree level}}.
\]

This closes C4 without using the scalar potential or a literature coupling table.

It does **not** yet establish `Z6=0`. That statement belongs to the Higgs-basis mass matrix and remains for Phase 7.

### 7. CP-odd sector

Because

\[
|x_i|^2=(v_i+\rho_i)^2+\eta_i^2,
\]

there is no term linear in `eta_i` multiplying `W^+W^-` or `ZZ`. Consequently, any orthogonal linear combination of the `eta_i`, including the physical `A`, has no tree-level `AVV` vertex of this form.

### 8. Quartic-gauge cross-check

The CP-even quadratic part is proportional to

\[
\rho_1^2+\rho_2^2.
\]

Since the CP-even mass rotation is orthogonal,

\[
\boxed{\rho_1^2+\rho_2^2=h^2+\phi^2}.
\]

There is therefore no `h phi W^+W^-` or `h phi ZZ` cross term from this norm. This check is independent of the trigonometric identities used for the linear coupling.

---

## Convention map

| Object | Project / DH-reference convention | BFLRS11 | GHOO18 | Active 2HDMC | Status |
|---|---|---|---|---|---|
| hypercharge | `Q=T3+Y`, `Y(Phi)=1/2` in this derivation | equivalent | equivalent | SM convention internally | `[TRANSLATED]` |
| vacuum direction | `rho_v=c_beta rho1+s_beta rho2` | implicit in `hVV/HVV` discussion | aligned state is VEV direction | first `q_ki` component | `[DERIVED][SOURCE]` |
| `kappa_V^h` | `sin(beta-alpha)` | same physical modifier after sign caution | aligned state SM-like | `q_11=sba` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| `kappa_V^phi` | `cos(beta-alpha)` | heavy-state modifier `cos(alpha-beta)` | orthogonal state vanishes in exact AL | `q_21=cba` | `[DERIVED][SOURCE][IMPLEMENTATION-CHECKED]` |
| exact alignment | `h=rho_v`, `phi=-rho_perp` | compare only after field-sign translation | `e_aligned=v`, others zero | `sba=1,cba=0` | `[DERIVED][TRANSLATED]` |
| `AVV` | zero at tree level | zero | zero for CP-odd neutral state | `q_31=0` | `[DERIVED][IMPLEMENTATION-CHECKED]` |

### BFLRS11 sign caution retained

The early BFLRS11 CP-even field display uses the global negatives of the DH05 fields, while nearby gauge-coupling prose quotes the standard positive `sin(beta-alpha)` and `cos(beta-alpha)` modifiers. We do not resolve that internal review inconsistency by choosing one quoted sign. The project signs come from the kinetic-term derivation above in the already-frozen Phase-4 DH convention.

---

## What was checked against the source

The source comparison was performed only after the analytic result above was frozen.

1. **DH05:** the audited source inventory locates the CP-even gauge-coupling angle table at `littletable` (around TeX line 1617), with `h_DH VV` proportional to `sin(beta-alpha)` and `H_DH VV` proportional to `cos(beta-alpha)`.
2. **BFLRS11:** its gauge-coupling prose gives the same physical `sin(beta-alpha)`/`cos(alpha-beta)` pattern, but its early displayed scalar-field signs are not used as normalization evidence.
3. **GHOO18:** exact alignment is described geometrically by the mass eigenstate aligned with the VEV; the aligned state carries the SM `VV` coupling and the orthogonal neutral states have zero `VV` coupling at exact alignment.
4. **2HDMC implementation check:** `THDM::get_qki` returns first components `(sba,cba,0,i)` for `(h,H,A,H+)`, and `THDM::get_coupling_vvh` multiplies `Re(q_{k1})` by `i g M_W` or `i gM_Z/c_W`. Using `M_W=gv/2` and `M_Z=g_Zv/2`, these are exactly the Feynman-rule coefficients `i(2m_V^2/v) kappa_V` derived above.

### Algebra audit

The following identities are sufficient for a symbolic audit independent of the gauge-group matrix multiplication:

```python
rho1 = -sin(alpha)*h + cos(alpha)*phi
rho2 =  cos(alpha)*h + sin(alpha)*phi
rho_v = cos(beta)*rho1 + sin(beta)*rho2
expand_trig(rho_v)
# -> sin(beta-alpha)*h + cos(beta-alpha)*phi

assert simplify((rho1**2 + rho2**2) - (h**2 + phi**2)) == 0
```

The group-theory factors are derived explicitly above rather than delegated to CAS.

---

## What remains uncertain

1. C4 (`kappa_V^phi=0` in exact alignment) is now **VERIFIED**.
2. The finite-mass relation between exact alignment and `Z6=0` is **not yet promoted**. Phase 7 must rotate the scalar potential explicitly into the Higgs basis and identify the CP-even off-diagonal mass entry as `Z6 v^2` in the selected sign convention.
3. The residual Higgs-basis sign `H2 -> -H2` still changes `Y3,Z6,Z7` and all odd-`H2` scalar interactions. Phase 7 must freeze
   \[
   H_1=c_\beta\Phi_1+s_\beta\Phi_2,
   \qquad
   H_2=-s_\beta\Phi_1+c_\beta\Phi_2.
   \]
4. The statement `phi VV=0` is a tree-level statement about the linear neutral-scalar gauge vertex. It does not by itself establish any loop-induced `phi -> gamma gamma` or `phi -> Z gamma` amplitude; loop amplitudes remain deferred.
5. A potential implementation sign layer involving 2HDMC Higgs-basis aliases (`Lambda_i`, local `l_i`, and `Z_i`) remains to be audited before any Phase-9 trilinear sign is promoted.

---

## Next smallest validation

Proceed to Phase 7 and define the Higgs basis explicitly,

\[
\boxed{H_1=c_\beta\Phi_1+s_\beta\Phi_2},
\qquad
\boxed{H_2=-s_\beta\Phi_1+c_\beta\Phi_2}.
\]

Then:

1. verify directly that only `H1` has a VEV;
2. rewrite the complete scalar potential in `Y_i,Z_i` form;
3. rederive the Higgs-basis stationarity conditions, especially the sign of `Y3`;
4. derive the CP-even Higgs-basis mass matrix and identify its off-diagonal element;
5. only after that decide the finite-mass exact-alignment statement involving `Z6`.
