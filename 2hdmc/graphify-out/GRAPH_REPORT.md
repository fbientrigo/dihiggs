# Graph Report - .  (2026-04-11)

## Corpus Check
- Corpus is ~41,146 words - fits in a single context window. You may not need a graph.

## Summary
- 269 nodes · 569 edges · 12 communities detected
- Extraction: 100% EXTRACTED · 0% INFERRED · 0% AMBIGUOUS
- Token cost: 0 input · 0 output

## God Nodes (most connected - your core abstractions)
1. `write_LesHouches()` - 24 edges
2. `print_decays()` - 16 edges
3. `clear_lookup()` - 16 edges
4. `get_gammatot_h()` - 13 edges
5. `thdmc_set_param()` - 13 edges
6. `get_param_phys()` - 10 edges
7. `oblique_param()` - 8 edges
8. `set_yukawas_type()` - 8 edges
9. `get_qki()` - 8 edges
10. `get_cba()` - 8 edges

## Surprising Connections (you probably didn't know these)
- `Constraints()` --calls--> `set_THDM()`  [EXTRACTED]
  src/Constraints.h → src/Constraints.cpp
- `DecayTable()` --calls--> `set_model()`  [EXTRACTED]
  src/DecayTable.h → src/DecayTable.cpp

## Communities

### Community 0 - "Decay width engine"
Cohesion: 0.09
Nodes (53): BHp(), br(), DecayTable(), DHm(), DHp(), F_0(), F_1(), F_pf() (+45 more)

### Community 1 - "SM parameter utilities"
Cohesion: 0.07
Nodes (41): clear_lookup(), get_CKM_element(), get_costw(), get_dmass_MSbar(), get_dmass_pole(), get_e(), get_g(), get_gamma_top() (+33 more)

### Community 2 - "Calculator frontends"
Cohesion: 0.08
Nodes (25): main(), pc(), check_lep(), check_masses(), check_perturbativity(), check_positivity(), check_stability(), check_unitarity() (+17 more)

### Community 3 - "THDM couplings"
Cohesion: 0.16
Nodes (22): get_alpha(), get_cba(), get_coupling_hdd(), get_coupling_hdu(), get_coupling_hll(), get_coupling_hln(), get_coupling_huu(), get_coupling_vhh() (+14 more)

### Community 4 - "FCNC helpers"
Cohesion: 0.1
Nodes (2): htb_fcn1(), sqrtlambda()

### Community 5 - "Yukawa model setup"
Cohesion: 0.15
Nodes (13): get_rho_lepton(), get_yukawas_lepton(), init(), set_inert(), set_kappa(), set_kappa_D(), set_kappa_L(), set_kappa_U() (+5 more)

### Community 6 - "HBHS integration"
Cohesion: 0.18
Nodes (6): charged_input(), check(), effective_couplings(), neutral_input(), nonSM_branching_ratios(), zeroIfNaN()

### Community 7 - "Parameter getters"
Cohesion: 0.22
Nodes (10): get_coupling_hhh(), get_m12_2(), get_param_gen(), get_param_HHG(), get_param_higgs(), get_param_hybrid(), print_param_gen(), print_param_HHG() (+2 more)

### Community 8 - "Parameter setters"
Cohesion: 0.27
Nodes (10): read_LesHouches(), set_hMSSM(), set_MSSM(), set_param_gen(), set_param_HHG(), set_param_higgs(), set_param_phys(), set_param_phys_lam1() (+2 more)

### Community 9 - "Stability and hybrid"
Cohesion: 0.25
Nodes (9): calc_unitarity(), check_stability(), check_unitarity(), get_SM(), recalc_tan_beta(), set_param_hybrid(), set_param_hybrid_sba(), set_SM() (+1 more)

### Community 10 - "Round-trip tests"
Cohesion: 0.83
Nodes (3): main(), nearly_equal(), run_case()

### Community 11 - "Perturbativity checks"
Cohesion: 0.67
Nodes (3): calc_perturbativity(), check_perturbativity(), get_coupling_hhhh()

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `write_LesHouches()` connect `THDM couplings` to `Stability and hybrid`, `Perturbativity checks`, `Yukawa model setup`, `Parameter getters`?**
  _High betweenness centrality (0.003) - this node is a cross-community bridge._
- **Should `Decay width engine` be split into smaller, more focused modules?**
  _Cohesion score 0.09 - nodes in this community are weakly interconnected._
- **Should `SM parameter utilities` be split into smaller, more focused modules?**
  _Cohesion score 0.07 - nodes in this community are weakly interconnected._
- **Should `Calculator frontends` be split into smaller, more focused modules?**
  _Cohesion score 0.08 - nodes in this community are weakly interconnected._
- **Should `FCNC helpers` be split into smaller, more focused modules?**
  _Cohesion score 0.1 - nodes in this community are weakly interconnected._