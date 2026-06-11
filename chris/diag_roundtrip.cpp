// Diagnostic (throwaway): trace what happens to lambda6, lambda7, m12_sq,
// and the Higgs-basis Lambda6 across the hybrid round trip for point 1.
#include "THDM.h"
#include "SM.h"
#include <cmath>
#include <iomanip>
#include <iostream>

int main() {
  std::cout.setf(std::ios::scientific);
  std::cout << std::setprecision(12);

  // Point 1
  double mh=125, mH=197.3684, mA=300, mHp=300, sba=1;
  double l6=9.993333333333333e-05, l7=0, lambda1_target=1, tanb=1.5e7;

  SM sm; double v2 = sm.get_v2();
  double beta = std::atan(tanb);
  double cb = std::cos(beta), tb = std::tan(beta);
  double alpha = -std::asin(sba)+beta;
  double sa = std::sin(alpha), ca = std::cos(alpha);
  double M11 = mH*mH*ca*ca + mh*mh*sa*sa;
  double m12sq = (M11 - v2*cb*cb*(lambda1_target + 1.5*l6*tb - 0.5*l7*tb*tb*tb))/tb;

  THDM model; model.set_SM(sm);
  model.set_param_phys(mh,mH,mA,mHp,sba,l6,l7,m12sq,tanb);

  double L1g,L2g,L3g,L4g,L5g,L6g,L7g,m12g,tbg;
  model.get_param_gen(L1g,L2g,L3g,L4g,L5g,L6g,L7g,m12g,tbg);

  double H1,H2,H3,H4,H5,H6,H7,Hmhp;
  model.get_param_higgs(H1,H2,H3,H4,H5,H6,H7,Hmhp);

  std::cout << "ORIGINAL model (generic basis):\n";
  std::cout << "  lambda1_gen = " << L1g << "\n";
  std::cout << "  lambda6_gen = " << L6g << "  (input was " << l6 << ")\n";
  std::cout << "  lambda7_gen = " << L7g << "  (input was " << l7 << ")\n";
  std::cout << "  m12sq_gen   = " << m12g << "  (computed " << m12sq << ")\n";
  std::cout << "  Higgs-basis Lambda6 = " << H6 << "   <-- this is DROPPED by get_param_hybrid\n";
  std::cout << "  Higgs-basis Lambda7 = " << H7 << "  (-> becomes hybrid Z7)\n";

  double M11_check = v2*cb*cb*(L1g + 1.5*L6g*tbg - 0.5*L7g*tbg*tbg*tbg) + m12g*tbg;
  std::cout << "  M11 invariant = " << M11_check << "  (vs analytic M11=" << M11 << ")\n\n";

  double hyb_mh,hyb_mH,hyb_cba,Z4,Z5,Z7,hyb_tb;
  model.get_param_hybrid(hyb_mh,hyb_mH,hyb_cba,Z4,Z5,Z7,hyb_tb);
  std::cout << "HYBRID extraction:\n";
  std::cout << "  mh="<<hyb_mh<<" mH="<<hyb_mH<<" cba="<<hyb_cba
            << " Z4="<<Z4<<" Z5="<<Z5<<" Z7="<<Z7<<" tanb="<<hyb_tb<<"\n\n";

  THDM rt; rt.set_SM(sm);
  rt.set_param_hybrid(hyb_mh,hyb_mH,hyb_cba,Z4,Z5,Z7,hyb_tb);

  double R1,R2,R3,R4,R5,R6,R7,Rm12,Rtb;
  rt.get_param_gen(R1,R2,R3,R4,R5,R6,R7,Rm12,Rtb);
  std::cout << "ROUND-TRIP model (generic basis):\n";
  std::cout << "  lambda1_rt = " << R1 << "\n";
  std::cout << "  lambda6_rt = " << R6 << "  (vs original lambda6_gen=" << L6g << ", input l6=" << l6 << ")\n";
  std::cout << "  lambda7_rt = " << R7 << "  (vs original lambda7_gen=" << L7g << ")\n";
  std::cout << "  m12sq_rt   = " << Rm12 << "  (vs original m12sq_gen=" << m12g << ")\n";

  double M11_rt = v2*cb*cb*(R1 + 1.5*R6*Rtb - 0.5*R7*Rtb*Rtb*Rtb) + Rm12*Rtb;
  std::cout << "  M11 invariant (rt) = " << M11_rt << "  (vs original " << M11_check << ")\n\n";

  std::cout << "Decomposition check: lambda6_rt * tanb = " << R6*Rtb
            << "   lambda1_rt - lambda1_gen = " << (R1-L1g) << "\n";
  std::cout << "1.5*(lambda6_rt - lambda6_gen)*tanb = " << 1.5*(R6-L6g)*Rtb << "\n";

  return 0;
}
