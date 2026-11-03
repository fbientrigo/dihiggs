// Campaign-only accounting audit. This enumerates the same native DecayTable
// terms as get_gammatot_h without changing the calculator library or bridge.
#include "DecayTable.h"
#include "THDM.h"
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct Term { std::string family; int i; int j; double width; };

static void number(const char *name, double value, bool comma=true) {
  std::cout << "\"" << name << "\":";
  if (std::isfinite(value)) std::cout << std::setprecision(17) << value;
  else std::cout << "null";
  if (comma) std::cout << ",";
}

static void add_term(std::vector<Term> &terms, const std::string &family,
                     int i, int j, double width, double &ordered, double &family_total) {
  terms.push_back(Term{family,i,j,width});
  ordered += width;
  family_total += width;
}

int main(int argc, char **argv) {
  if (argc != 16) {
    std::cerr << "usage: lambda1..lambda7 m12_sq tan_beta requested_v h_index phi_index expected_mh expected_mphi expected_backend_v\n";
    return 2;
  }
  try {
    double lambda[7];
    for (int i=0;i<7;++i) lambda[i]=std::stod(argv[i+1]);
    const double m12=std::stod(argv[8]);
    const double tan_beta=std::stod(argv[9]);
    const double requested_v=std::stod(argv[10]);
    const int h_index=std::stoi(argv[11]);
    const int phi_index=std::stoi(argv[12]);
    const double expected_mh=std::stod(argv[13]);
    const double expected_mphi=std::stod(argv[14]);
    const double expected_backend_v=std::stod(argv[15]);
    for (double x : lambda) if (!std::isfinite(x)) throw std::invalid_argument("nonfinite quartic");
    if (!std::isfinite(m12) || !std::isfinite(tan_beta) || !std::isfinite(requested_v) ||
        !std::isfinite(expected_mh) || !std::isfinite(expected_mphi) || !std::isfinite(expected_backend_v) ||
        tan_beta<=0. || requested_v<=0. || (h_index!=1 && h_index!=2) ||
        (phi_index!=1 && phi_index!=2) || h_index==phi_index)
      throw std::invalid_argument("invalid audit input domain");

    THDM model;
    SM sm;
    sm.set_GF(1.0/(std::sqrt(2.0)*requested_v*requested_v));
    model.set_SM(sm);
    const bool constructed=model.set_param_gen_exact(lambda[0],lambda[1],lambda[2],lambda[3],
        lambda[4],lambda[5],lambda[6],m12,tan_beta);
    if (!constructed) {
      std::cout << "{\"constructor_ok\":false,\"reason\":\"set_param_gen_exact returned false\"}\n";
      return 0;
    }
    model.set_yukawas_type(1);
    const double mass_h=model.get_hmass(h_index);
    const double mass_phi=model.get_hmass(phi_index);
    const double backend_v=model.get_SM().get_v();
    if (std::fabs(mass_h-expected_mh)>1e-8 || std::fabs(mass_phi-expected_mphi)>1e-8 ||
        std::fabs(backend_v-expected_backend_v)>1e-10) {
      std::cout << "{\"constructor_ok\":true,\"state_identity_ok\":false,\"reason\":\"reconstructed state masses or electroweak scale differ from frozen bridge record\",\"m_h_index_GeV\":"
                << std::setprecision(17) << mass_h << ",\"m_phi_index_GeV\":" << mass_phi
                << ",\"backend_v_GeV\":" << backend_v << "}\n";
      return 0;
    }

    DecayTable decay(model);
    const double total=decay.get_gammatot_h(phi_index);
    const double diphoton=decay.get_gamma_hgaga(phi_index);
    const double ctau=total>0. ? 1.973269804e-13/total : std::numeric_limits<double>::quiet_NaN();
    std::vector<Term> terms;
    terms.reserve(78);
    double ordered_sum=0.;
    double hdd=0., huu=0., hdu=0., hll=0., hln=0.;
    for (int i=1;i<4;++i) for (int j=1;j<4;++j) {
      add_term(terms,"hdd",i,j,decay.get_gamma_hdd(phi_index,i,j),ordered_sum,hdd);
      add_term(terms,"huu",i,j,decay.get_gamma_huu(phi_index,i,j),ordered_sum,huu);
      add_term(terms,"hdu",i,j,decay.get_gamma_hdu(phi_index,i,j),ordered_sum,hdu);
      add_term(terms,"hll",i,j,decay.get_gamma_hll(phi_index,i,j),ordered_sum,hll);
      add_term(terms,"hln",i,j,decay.get_gamma_hln(phi_index,i,j),ordered_sum,hln);
    }
    double hvv=0.;
    for (int i=1;i<4;++i) add_term(terms,"hvv",i,0,decay.get_gamma_hvv(phi_index,i),ordered_sum,hvv);
    const double zgamma=decay.get_gamma_hZga(phi_index);
    double zgamma_total=0.;
    add_term(terms,"hZga",0,0,zgamma,ordered_sum,zgamma_total);
    double gluons=0.;
    add_term(terms,"hgg",0,0,decay.get_gamma_hgg(phi_index),ordered_sum,gluons);
    double hvh=0.;
    for (int i=1;i<4;++i) for (int j=1;j<=4;++j)
      add_term(terms,"hvh",i,j,decay.get_gamma_hvh(phi_index,i,j),ordered_sum,hvh);
    double hhh=0.;
    for (int i=1;i<5;++i) for (int j=1;j<5;++j)
      add_term(terms,"hhh",i,j,decay.get_gamma_hhh(phi_index,i,j),ordered_sum,hhh);

    std::cout << "{\"constructor_ok\":true,\"state_identity_ok\":true,";
    number("h_index",h_index); number("phi_index",phi_index);
    number("m_h_index_GeV",mass_h); number("m_phi_index_GeV",mass_phi);
    number("requested_v_GeV",requested_v); number("backend_v_GeV",backend_v);
    number("total_width_GeV",total); number("ordered_term_sum_GeV",ordered_sum);
    number("diphoton_width_GeV",diphoton); number("ctau_mm",ctau);
    std::cout << "\"family_totals_GeV\":{";
    number("hdd",hdd); number("huu",huu); number("hdu",hdu); number("hll",hll);
    number("hln",hln); number("hvv",hvv); number("hZga",zgamma_total); number("hgg",gluons);
    number("hvh",hvh); number("hhh",hhh,false);
    std::cout << "},\"terms\":[";
    for (size_t i=0;i<terms.size();++i) {
      if (i) std::cout << ",";
      const Term &term=terms[i];
      std::cout << "{\"family\":\"" << term.family << "\",\"i\":" << term.i << ",\"j\":" << term.j << ",\"width_GeV\":";
      if (std::isfinite(term.width)) std::cout << std::setprecision(17) << term.width;
      else std::cout << "null";
      std::cout << "}";
    }
    std::cout << "]}\n";
  } catch (const std::exception &error) {
    std::cerr << "decay-accounting audit input error: " << error.what() << "\n";
    return 2;
  }
  return 0;
}
