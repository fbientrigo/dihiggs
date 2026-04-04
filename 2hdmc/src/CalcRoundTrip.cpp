// src/CalcRoundTrip.cpp

#include "THDM.h"

#include <algorithm>
#include <cmath>
#include <iostream>

static bool nearly_equal(double a, double b, double rtol, double atol) {
  const double diff = std::fabs(a - b);
  const double scale = std::max(std::fabs(a), std::fabs(b));
  return diff <= (atol + rtol * scale);
}

static int run_case(const char* name,
                    double m_h, double m_H, double m_A, double m_Hp,
                    double sba, double lambda6, double lambda7,
                    double m12_2, double tan_beta,
                    double rtol, double atol) {

  THDM m1;
  if (!m1.set_param_phys(m_h, m_H, m_A, m_Hp, sba, lambda6, lambda7, m12_2, tan_beta)) {
    std::cout << "Case " << name << ": FAIL (set_param_phys returned false)\n";
    return 1;
  }

  double l1,l2,l3,l4,l5,l6,l7,m12_ref,tb_ref;
  m1.get_param_gen(l1,l2,l3,l4,l5,l6,l7,m12_ref,tb_ref);

  THDM m2;
  if (!m2.set_param_phys_lam1(m_h, m_H, m_A, m_Hp, sba, l1, lambda6, lambda7, tan_beta)) {
    std::cout << "Case " << name << ": FAIL (set_param_phys_lam1 returned false)\n";
    return 1;
  }

  double L1,L2,L3,L4,L5,L6,L7,m12_new,tb_new;
  m2.get_param_gen(L1,L2,L3,L4,L5,L6,L7,m12_new,tb_new);

  bool ok = true;
  ok &= nearly_equal(m12_new, m12_ref, rtol, atol);
  ok &= nearly_equal(tb_new, tb_ref, rtol, atol);
  ok &= nearly_equal(L1, l1, rtol, atol);
  ok &= nearly_equal(L2, l2, rtol, atol);
  ok &= nearly_equal(L3, l3, rtol, atol);
  ok &= nearly_equal(L4, l4, rtol, atol);
  ok &= nearly_equal(L5, l5, rtol, atol);
  ok &= nearly_equal(L6, l6, rtol, atol);
  ok &= nearly_equal(L7, l7, rtol, atol);

  for (int h = 1; h <= 4; ++h) {
    ok &= nearly_equal(m2.get_hmass(h), m1.get_hmass(h), rtol, atol);
  }

  if (ok) {
    std::cout << "Case " << name << ": PASS\n";
    return 0;
  }

  std::cout << "Case " << name << ": FAIL\n";
  std::cout << "  m12_2: ref=" << m12_ref << " new=" << m12_new << "\n";
  std::cout << "  tanb : ref=" << tb_ref  << " new=" << tb_new  << "\n";
  return 1;
}

int main() {
  const double rtol = 1e-10;
  const double atol = 1e-12;

  int rc = 0;

  // Generic point
  rc |= run_case("generic",
                 125.0, 300.0, 350.0, 360.0,
                 0.80,
                 0.10, -0.20,
                 20000.0,
                 2.0,
                 rtol, atol);

  // Near-alignment point
  rc |= run_case("near-alignment",
                 125.0, 300.0, 350.0, 360.0,
                 0.9999,
                 0.10, -0.20,
                 20000.0,
                 2.0,
                 rtol, atol);

  // Invalid-input guard (expected failure)
  {
    THDM m_bad;
    const bool ok = !m_bad.set_param_phys_lam1(125.0, 300.0, 350.0, 360.0,
                                               0.80,
                                               1.0,
                                               0.10, -0.20,
                                               -1.0);
    if (ok) {
      std::cout << "Case invalid-input: PASS\n";
    } else {
      std::cout << "Case invalid-input: FAIL\n";
      rc |= 1;
    }
  }

  if (rc == 0) {
    std::cout << "All tests passed\n";
  }
  return rc;
}
