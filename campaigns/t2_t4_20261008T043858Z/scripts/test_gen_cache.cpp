#include "THDM.h"
#include <cmath>
#include <iostream>

int main() {
  THDM model;
  SM sm; model.set_SM(sm);
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 1;
  double a,b,c,d,e,f,g,m12,tb;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12!=.007 || tb!=1.e6) return 2;
  // A file-open failure leaves the model untouched, so preserve the cache.
  if (model.read_LesHouches("/tmp/t2_t4_missing_lha_input.slha")) return 13;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12!=.007 || tb!=1.e6) return 14;
  // A failed setter must invalidate the exact cache too.
  if (model.set_param_HHG(1,1,1,1,1,1,-1)) return 3;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 4;
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 5;
  if (!model.set_param_phys(125.,200.,250.,250.,1.,.01,.02,.1,4.)) return 6;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 7;
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 8;
  if (!model.set_param_higgs(.3,.4,.2,.1,-.05,.01,.02,250.)) return 9;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 10;
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 11;
  SM changed=sm; changed.set_GF(SM::GF*1.001); model.set_SM(changed);
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 12;
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 15;
  if (!model.set_inert(125.,200.,210.,220.,.5,.2)) return 16;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 17;
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 18;
  if (model.set_MSSM(-1.,1.e6)) return 19;
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 20;
  if (!model.set_param_gen_exact(.3,.4,.2,.1,-.05,.01,.02,.007,1.e6)) return 21;
  model.recalc_tan_beta(2.);
  model.get_param_gen(a,b,c,d,e,f,g,m12,tb);
  if (m12==.007 || tb==1.e6) return 22;
  std::cout << "exact generic cache preservation and setter invalidation PASS\n";
  return 0;
}
