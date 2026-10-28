// Isolated T4 calculator bridge. Reads one parameter set from command line.
#include "DecayTable.h"
#include "THDM.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

static void number(const char *name, double x, bool comma=true) {
  std::cout << "\"" << name << "\":";
  if (std::isfinite(x)) std::cout << std::setprecision(17) << x;
  else std::cout << "null";
  if (comma) std::cout << ",";
}

static void print_array(const char *name, const double *v, size_t n, bool comma=true) {
  std::cout << "\"" << name << "\":[";
  for(size_t i=0;i<n;++i) { if(i) std::cout << ","; std::cout << std::setprecision(17) << v[i]; }
  std::cout << "]"; if(comma) std::cout << ",";
}

int main(int argc, char **argv) {
  if (argc != 16) {
    std::cerr << "mode mh mH mA mHp tanb l1 l2 l3 l4 l5 l6 l7 m12 requested_v\n";
    return 2;
  }
  const std::string mode=argv[1];
  double p[14]; for(int i=2;i<16;++i) p[i-2]=std::stod(argv[i]);
  const double mh=p[0],mH=p[1],mA=p[2],mHp=p[3],tb=p[4];
  const double l1=p[5],l2=p[6],l3=p[7],l4=p[8],l5=p[9],l6=p[10],l7=p[11];
  const double m12=p[12],requested_v=p[13];
  THDM model;
  SM sm;
  if (mode=="explicit" || mode=="analytic")
    sm.set_GF(1.0/(std::sqrt(2.0)*requested_v*requested_v));
  model.set_SM(sm);
  bool ok=false;
  if (mode=="legacy") ok=model.set_param_phys(mh,mH,mA,mHp,1.0,l6,l7,m12,tb);
  else if (mode=="explicit") ok=model.set_param_phys_lam1(mh,mH,mA,mHp,1.0,l1,l6,l7,tb);
  else if (mode=="analytic") ok=model.set_param_gen_exact(l1,l2,l3,l4,l5,l6,l7,m12,tb);
  if (!ok) {
    std::cout << "{\"mode\":\"" << mode << "\",\"constructor_ok\":false}\n";
    return 0;
  }
  model.set_yukawas_type(1);
  double a,b,c,d,e,f,g,soft,bet;
  model.get_param_gen(a,b,c,d,e,f,g,soft,bet);
  double z1,z2,z3,z4,z5,z6,z7,hpc;
  model.get_param_higgs(z1,z2,z3,z4,z5,z6,z7,hpc);
  double s20[1],s21[3],s01[4],s00[4];
  const double max_unitarity=model.calc_unitarity_spectrum(s20,s21,s01,s00);
  DecayTable decay(model);
  const double width=decay.get_gammatot_h(2);
  const double diphoton=decay.get_gamma_hgaga(2);
  const double width_gg=decay.get_gamma_hgg(2);
  const double width_Zga=decay.get_gamma_hZga(2);
  const double width_bb=decay.get_gamma_hdd(2,3,3);
  const double width_cc=decay.get_gamma_huu(2,2,2);
  const double width_tautau=decay.get_gamma_hll(2,3,3);
  const double width_mumu=decay.get_gamma_hll(2,2,2);
  const double width_WW=decay.get_gamma_hvv(2,3);
  const double width_ZZ=decay.get_gamma_hvv(2,2);
  const double width_hh=decay.get_gamma_hhh(2,1,1);
  std::complex<double> scalar,pseudo;
  model.get_coupling_huu(2,3,3,scalar,pseudo);
  SM sm_copy=model.get_SM();
  const double mt_run=sm_copy.run_qmass_MSbar(sm_copy.get_umass_MSbar(3),
      sm_copy.get_umass_MSbar(3),model.get_hmass(2),sm_copy.get_qmass_pole(6),sm_copy.get_qmass_pole(5));
  std::complex<double> hphi,phiHp;
  model.get_coupling_hhh(1,2,2,hphi);
  model.get_coupling_hhh(2,4,4,phiHp);
  std::cout << "{\"mode\":\"" << mode << "\",\"constructor_ok\":true,";
  number("requested_v_GeV",requested_v);
  number("backend_v_GeV",model.get_SM().get_v());
  number("lambda1",a);number("lambda2",b);number("lambda3",c);number("lambda4",d);
  number("lambda5",e);number("lambda6",f);number("lambda7",g);
  number("m12_sq_GeV2",soft);number("tan_beta_backend",bet);
  number("Z1",z1);number("Z2",z2);number("Z3",z3);number("Z4",z4);
  number("Z5",z5);number("Z6",z6);number("Z7",z7);
  number("m_h_backend_GeV",model.get_hmass(1));number("m_H_backend_GeV",model.get_hmass(2));
  number("m_A_backend_GeV",model.get_hmass(3));number("m_Hp_backend_GeV",model.get_hmass(4));
  number("calc_unitarity_legacy_max_abs_Lambda",max_unitarity);
  print_array("backend_eigenvalues_Y2_sigma1",s21,3);print_array("backend_eigenvalues_Y2_sigma0",s20,1);
  print_array("backend_eigenvalues_Y0_sigma1",s01,4);print_array("backend_eigenvalues_Y0_sigma0",s00,4);
  number("width_total_phi_GeV",width);number("width_diphoton_GeV",diphoton);
  number("width_gg_GeV",width_gg);number("width_Zgamma_GeV",width_Zga);
  number("width_bb_GeV",width_bb);number("width_cc_GeV",width_cc);
  number("width_tautau_GeV",width_tautau);number("width_mumu_GeV",width_mumu);
  number("width_WW_GeV",width_WW);number("width_ZZ_GeV",width_ZZ);
  number("width_hh_GeV",width_hh);
  number("BR_diphoton",width>0?diphoton/width:std::numeric_limits<double>::quiet_NaN());
  number("ctau_mm",width>0?1.973269804e-13/width:std::numeric_limits<double>::quiet_NaN());
  number("phi_top_scalar_coupling_real",scalar.real());number("phi_top_scalar_coupling_imag",scalar.imag());
  number("phi_top_pseudoscalar_coupling_real",pseudo.real());number("phi_top_pseudoscalar_coupling_imag",pseudo.imag());
  number("phi_top_running_mass_GeV",mt_run);
  // The 2HDMC h2 fermion vertex is +i*rho/sqrt(2) in this convention;
  // the signed scalar modifier is -Im(cs)*v/m_run.
  number("kappa_phi_top_backend_signed",mt_run>0?-scalar.imag()*model.get_SM().get_v()/mt_run:std::numeric_limits<double>::quiet_NaN());
  number("sin_beta_minus_alpha_backend",model.get_sba());number("cos_beta_minus_alpha_backend",model.get_cba());
  number("g_hphi_phi_backend_GeV",-hphi.imag());
  number("g_phi_HpHm_backend_GeV",-phiHp.imag(),false);
  std::cout << "}\n";
  return 0;
}
