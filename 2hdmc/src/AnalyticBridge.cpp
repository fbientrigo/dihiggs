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

static void print_complex(const std::complex<double> &z) {
  std::cout << "{\"real\":" << std::setprecision(17) << z.real()
            << ",\"imag\":" << std::setprecision(17) << z.imag() << "}";
}

static int sign_nonzero(double x) { return x < 0. ? -1 : 1; }

static void print_hgamma_diagnostics(const DecayTable::HgammaDiagnostics &d, int phase) {
  std::cout << "\"diphoton_native_diagnostics\":{";
  number("mass_GeV",d.mass_GeV); number("vev_GeV",d.vev_GeV);
  number("alpha0",d.alpha0); number("mW_GeV",d.mW_GeV); number("mHp_GeV",d.mHp_GeV);
  number("alpha_s_at_mass",d.alpha_s_at_mass);
  number("loop_width_GeV",d.loop_width_GeV);
  std::cout << "\"qcd_enabled\":" << (d.qcd_enabled?"true":"false") << ",";
  std::cout << "\"heavy_damping_enabled\":" << (d.heavy_damping_enabled?"true":"false") << ",";
  std::cout << "\"state_phase\":" << phase << ",\"W_amplitude_raw\":"; print_complex(d.W_amplitude);
  std::cout << ",\"charged_higgs_amplitude_raw\":"; print_complex(d.charged_higgs_amplitude);
  std::cout << ",\"W_coupling_raw\":"; print_complex(d.W_coupling_raw);
  std::cout << ",\"charged_higgs_coupling_raw\":"; print_complex(d.charged_higgs_coupling_raw);
  std::cout << ",\"scalar_amplitude_raw\":"; print_complex(d.scalar_amplitude);
  std::cout << ",\"pseudoscalar_amplitude_raw\":"; print_complex(d.pseudoscalar_amplitude);
  std::cout << ",\"scalar_amplitude_mapped\":"; print_complex(static_cast<double>(phase)*d.scalar_amplitude);
  std::cout << ",\"pseudoscalar_amplitude_mapped\":"; print_complex(static_cast<double>(phase)*d.pseudoscalar_amplitude);
  std::cout << ",\"fermions\":[";
  for (size_t i=0;i<d.fermions.size();++i) {
    const DecayTable::HgammaFermionDiagnostic &f=d.fermions[i];
    if(i) std::cout << ",";
    std::cout << "{\"family\":\"" << f.family << "\",\"generation\":" << f.generation << ",";
    number("msbar_mass_GeV",f.msbar_mass_GeV); number("pole_mass_GeV",f.pole_mass_GeV);
    number("heavy_damping_mass_GeV",f.heavy_damping_mass_GeV);
    number("vertex_running_mass_GeV",f.vertex_running_mass_GeV); number("tau",f.tau);
    number("q_init_GeV",f.q_init_GeV);
    number("running_mass_at_half_scale_GeV",f.running_mass_at_half_scale_GeV);
    number("running_mass_at_pole_scale_GeV",f.running_mass_at_pole_scale_GeV);
    number("alpha_s",f.alpha_s);
    std::cout << "\"included\":" << (f.included?"true":"false") << ",\"scalar_coupling\":";
    print_complex(f.scalar_coupling); std::cout << ",\"pseudoscalar_coupling\":"; print_complex(f.pseudoscalar_coupling);
    std::cout << ",\"scalar_amplitude_before_qcd\":"; print_complex(f.scalar_amplitude_before_qcd);
    std::cout << ",\"pseudoscalar_amplitude_before_qcd\":"; print_complex(f.pseudoscalar_amplitude_before_qcd);
    std::cout << ",\"scalar_amplitude_after_qcd\":"; print_complex(f.scalar_amplitude_after_qcd);
    std::cout << ",\"pseudoscalar_amplitude_after_qcd\":"; print_complex(f.pseudoscalar_amplitude_after_qcd);
    std::cout << "}";
  }
  std::cout << "]},";
}

int main(int argc, char **argv) {
  if (argc != 16 && argc != 17) {
    std::cerr << "mode mh mH mA mHp tanb l1 l2 l3 l4 l5 l6 l7 m12 requested_v [qcd_on]\n";
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
  const double sorted_m1=model.get_hmass(1), sorted_m2=model.get_hmass(2);
  const double cost_as_is=std::fabs(sorted_m1-mh)+std::fabs(sorted_m2-mH);
  const double cost_swapped=std::fabs(sorted_m2-mh)+std::fabs(sorted_m1-mH);
  const int h_index=(cost_as_is<=cost_swapped)?1:2;
  const int phi_index=(h_index==1)?2:1;
  const bool ordering_ok=(mh!=mH)&&std::fabs(cost_as_is-cost_swapped)>1e-8;
  const double sba=model.get_sba(), cba=model.get_cba();
  const double h_vev_projection=(h_index==1)?sba:cba;
  const double phi_vev_projection=(phi_index==1)?sba:cba;
  const int phase_h=sign_nonzero(h_vev_projection);
  const int phase_phi=(h_index==1)?phase_h:-phase_h;
  const double h_v_mapped=phase_h*h_vev_projection;
  const double phi_v_mapped=phase_phi*phi_vev_projection;
  double a,b,c,d,e,f,g,soft,bet;
  model.get_param_gen(a,b,c,d,e,f,g,soft,bet);
  double z1,z2,z3,z4,z5,z6,z7,hpc;
  model.get_param_higgs(z1,z2,z3,z4,z5,z6,z7,hpc);
  double s20[1],s21[3],s01[4],s00[4];
  const double max_unitarity=model.calc_unitarity_spectrum(s20,s21,s01,s00);
  double width=std::numeric_limits<double>::quiet_NaN();
  double diphoton=width,width_gg=width,width_Zga=width,width_bb=width,width_cc=width;
  double width_tautau=width,width_mumu=width,width_WW=width,width_ZZ=width,width_hh=width;
  DecayTable decay(model);
  if (argc==17) decay.set_qcd(std::string(argv[16])!="0");
  DecayTable::HgammaDiagnostics hgamma;
  if (ordering_ok) {
    width=decay.get_gammatot_h(phi_index);
    diphoton=decay.get_gamma_hgaga(phi_index);
    decay.get_hgaga_diagnostics(phi_index,hgamma);
    width_gg=decay.get_gamma_hgg(phi_index);
    width_Zga=decay.get_gamma_hZga(phi_index);
    width_bb=decay.get_gamma_hdd(phi_index,3,3);
    width_cc=decay.get_gamma_huu(phi_index,2,2);
    width_tautau=decay.get_gamma_hll(phi_index,3,3);
    width_mumu=decay.get_gamma_hll(phi_index,2,2);
    width_WW=decay.get_gamma_hvv(phi_index,3);
    width_ZZ=decay.get_gamma_hvv(phi_index,2);
    width_hh=decay.get_gamma_hhh(phi_index,h_index,h_index);
  }
  std::complex<double> scalar,pseudo;
  model.get_coupling_huu(phi_index,3,3,scalar,pseudo);
  SM sm_copy=model.get_SM();
  const double mt_run=sm_copy.run_qmass_MSbar(sm_copy.get_umass_MSbar(3),
      sm_copy.get_umass_MSbar(3),model.get_hmass(phi_index),sm_copy.get_qmass_pole(6),sm_copy.get_qmass_pole(5));
  std::complex<double> hphi,phiHp;
  model.get_coupling_hhh(h_index,phi_index,phi_index,hphi);
  model.get_coupling_hhh(phi_index,4,4,phiHp);
  std::cout << "{\"mode\":\"" << mode << "\",\"constructor_ok\":true,";
  std::cout << "\"physical_state_ordering\":\"" << (ordering_ok?"MAPPED":"INCONCLUSIVE") << "\",";
  std::cout << "\"exact_failure_reason\":\"" << (ordering_ok?"":"CP-even mass ordering is degenerate or ambiguous at backend precision") << "\",";
  number("requested_h_index_backend",h_index); number("requested_phi_index_backend",phi_index);
  number("requested_h_phase_factor",phase_h); number("requested_phi_phase_factor",phase_phi);
  number("raw_h_vev_projection",h_vev_projection); number("raw_phi_vev_projection",phi_vev_projection);
  number("mapped_h_vev_projection",h_v_mapped); number("mapped_phi_vev_projection",phi_v_mapped);
  number("m_h_sorted_backend_GeV",sorted_m1); number("m_H_sorted_backend_GeV",sorted_m2);
  number("requested_v_GeV",requested_v);
  number("backend_v_GeV",model.get_SM().get_v());
  number("lambda1",a);number("lambda2",b);number("lambda3",c);number("lambda4",d);
  number("lambda5",e);number("lambda6",f);number("lambda7",g);
  number("m12_sq_GeV2",soft);number("tan_beta_backend",bet);
  number("Z1",z1);number("Z2",z2);number("Z3",z3);number("Z4",z4);
  number("Z5",z5);number("Z6",z6);number("Z7",z7);
  number("m_h_backend_GeV",model.get_hmass(h_index));number("m_H_backend_GeV",model.get_hmass(phi_index));
  number("m_A_backend_GeV",model.get_hmass(3));number("m_Hp_backend_GeV",model.get_hmass(4));
  number("calc_unitarity_legacy_max_abs_Lambda",max_unitarity);
  print_array("backend_eigenvalues_Y2_sigma1",s21,3);print_array("backend_eigenvalues_Y2_sigma0",s20,1);
  print_array("backend_eigenvalues_Y0_sigma1",s01,4);print_array("backend_eigenvalues_Y0_sigma0",s00,4);
  std::cout << "\"widths_computed\":" << (ordering_ok?"true":"false") << ",";
  if (ordering_ok) print_hgamma_diagnostics(hgamma,phase_phi);
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
  number("kappa_phi_top_backend_signed",mt_run>0?phase_phi*(-scalar.imag()*model.get_SM().get_v()/mt_run):std::numeric_limits<double>::quiet_NaN());
  number("kappa_phi_top_backend_raw_signed",mt_run>0?-scalar.imag()*model.get_SM().get_v()/mt_run:std::numeric_limits<double>::quiet_NaN());
  number("sin_beta_minus_alpha_backend",h_v_mapped);number("cos_beta_minus_alpha_backend",phi_v_mapped);
  number("sin_beta_minus_alpha_raw",model.get_sba());number("cos_beta_minus_alpha_raw",model.get_cba());
  number("g_hphi_phi_backend_GeV",-phase_h*hphi.imag());
  number("g_hphi_phi_backend_raw_GeV",-hphi.imag());
  number("g_phi_HpHm_backend_GeV",-phase_phi*phiHp.imag());
  number("g_phi_HpHm_backend_raw_GeV",-phiHp.imag(),false);
  std::cout << "}\n";
  return 0;
}
