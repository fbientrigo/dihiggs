# Decisions — set_param_phys_lam1

## Source-consistent m12² inversion formula
```
m12² = (m_H² cos²α + m_h² sin²α 
        - v² cos²β (λ1 + 3/2 λ6 tanβ - 1/2 λ7 tan³β)) / tanβ
```
Where:
- v² = THDM::v2 (from SM::get_v2(), NOT v²/2)
- α = β − asin(sba)
- β = atan(tan_beta)

## Implementation decisions
- Delegate to existing set_param_phys after computing m12_2
- Replicate guard checks from set_param_phys (tan_beta>0, |sba|<=1, masses>=0)
- NO refactoring of existing code
- NO new physics constraints beyond existing guards
