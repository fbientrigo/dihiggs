#include "THDM.h"
#include <cmath>
#include <iostream>

int main() {
  THDM model;
  if (!model.set_param_gen_exact(0.,0.,0.,0.,1.,0.,0.,1.,2.)) return 1;
  double s20[1],s21[3],s01[4],s00[4];
  const double maximum=model.calc_unitarity_spectrum(s20,s21,s01,s00);
  if (maximum != model.calc_unitarity()) return 2;
  if (!(s21[0] <= s21[1] && s21[1] <= s21[2])) return 3;
  for (int i=1;i<4;i++) if (s01[i-1] > s01[i] || s00[i-1] > s00[i]) return 4;
  // The lambda5-only fixture has exact +/-1 ties across all four blocks.
  // The legacy signed maximum is selected in GSL's unsorted iteration order.
  if (std::fabs(maximum-3.0) > 1e-12) {
    std::cerr << "signed tie changed: " << maximum << "\n";
    return 5;
  }
  std::cout << "complete spectrum ordering and legacy signed tie PASS\n";
  return 0;
}
