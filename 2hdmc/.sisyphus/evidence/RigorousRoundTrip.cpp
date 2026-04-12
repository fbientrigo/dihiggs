#include "THDM.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

static bool nearly_equal(double a, double b, double rtol, double atol) {
	const double diff = std::fabs(a - b);
	const double scale = std::max(std::fabs(a), std::fabs(b));
	return diff <= (atol + rtol * scale);
}

struct Point {
	double m_h;
	double m_H;
	double m_A;
	double m_Hp;
	double sba;
	double lambda6;
	double lambda7;
	double m12_2;
	double tan_beta;
};

static bool check_finite(double x) {
	return std::isfinite(x);
}

static bool run_point(const Point &p, int idx, double rtol, double atol,
			      double &max_abs, double &max_rel, std::string &failure) {
	THDM m_ref;
	if (!m_ref.set_param_phys(p.m_h, p.m_H, p.m_A, p.m_Hp,
				  p.sba, p.lambda6, p.lambda7, p.m12_2, p.tan_beta)) {
		std::ostringstream oss;
		oss << "Point " << idx << ": set_param_phys returned false";
		failure = oss.str();
		return false;
	}

	double l1, l2, l3, l4, l5, l6, l7, m12_ref, tb_ref;
	m_ref.get_param_gen(l1, l2, l3, l4, l5, l6, l7, m12_ref, tb_ref);

	if (!check_finite(l1) || !check_finite(l2) || !check_finite(l3) ||
	    !check_finite(l4) || !check_finite(l5) || !check_finite(l6) ||
	    !check_finite(l7) || !check_finite(m12_ref) || !check_finite(tb_ref)) {
		std::ostringstream oss;
		oss << "Point " << idx << ": non-finite values in reference model";
		failure = oss.str();
		return false;
	}

	THDM m_lam1;
	if (!m_lam1.set_param_phys_lam1(p.m_h, p.m_H, p.m_A, p.m_Hp,
					 p.sba, l1, p.lambda6, p.lambda7, p.tan_beta)) {
		std::ostringstream oss;
		oss << "Point " << idx << ": set_param_phys_lam1 returned false";
		failure = oss.str();
		return false;
	}

	double L1, L2, L3, L4, L5, L6, L7, m12_new, tb_new;
	m_lam1.get_param_gen(L1, L2, L3, L4, L5, L6, L7, m12_new, tb_new);

	auto check_pair = [&](const char *name, double a, double b) {
		const double abs_diff = std::fabs(a - b);
		const double scale = std::max(std::fabs(a), std::fabs(b));
		const double rel_diff = (scale > 0.0) ? (abs_diff / scale) : abs_diff;
		max_abs = std::max(max_abs, abs_diff);
		max_rel = std::max(max_rel, rel_diff);

		if (!check_finite(a) || !check_finite(b)) {
			std::ostringstream oss;
			oss << "Point " << idx << " " << name << ": non-finite values";
			failure = oss.str();
			return false;
		}

		if (!nearly_equal(a, b, rtol, atol)) {
			std::ostringstream oss;
			oss << "Point " << idx << " " << name
			    << " mismatch: ref=" << std::setprecision(17) << a
			    << " new=" << std::setprecision(17) << b
			    << " abs_diff=" << abs_diff
			    << " rel_diff=" << rel_diff;
			failure = oss.str();
			return false;
		}

		return true;
	};

	if (!check_pair("m12_2", m12_ref, m12_new) ||
	    !check_pair("tan_beta", tb_ref, tb_new) ||
	    !check_pair("lambda1", l1, L1) ||
	    !check_pair("lambda2", l2, L2) ||
	    !check_pair("lambda3", l3, L3) ||
	    !check_pair("lambda4", l4, L4) ||
	    !check_pair("lambda5", l5, L5) ||
	    !check_pair("lambda6", l6, L6) ||
	    !check_pair("lambda7", l7, L7)) {
		return false;
	}

	for (int h = 1; h <= 4; ++h) {
		double m_ref_h = m_ref.get_hmass(h);
		double m_new_h = m_lam1.get_hmass(h);
		std::ostringstream name;
		name << "hmass(" << h << ")";
		if (!check_pair(name.str().c_str(), m_ref_h, m_new_h)) {
			return false;
		}
	}

	return true;
}

int main() {
	const double rtol = 1e-10;
	const double atol = 1e-12;
	const int n_points = 400;

	std::mt19937_64 rng(20260403ULL);
	std::uniform_real_distribution<double> dist_mH(150.0, 1200.0);
	std::uniform_real_distribution<double> dist_mA(130.0, 1200.0);
	std::uniform_real_distribution<double> dist_mHp(130.0, 1200.0);
	std::uniform_real_distribution<double> dist_sba(-0.999999, 0.999999);
	std::uniform_real_distribution<double> dist_l6l7(-2.0, 2.0);
	std::uniform_real_distribution<double> dist_m12(-8.0e5, 8.0e5);
	std::uniform_real_distribution<double> dist_tb(0.15, 60.0);

	double max_abs = 0.0;
	double max_rel = 0.0;

	for (int i = 0; i < n_points; ++i) {
		Point p;
		p.m_h = 125.0;
		p.m_H = dist_mH(rng);
		p.m_A = dist_mA(rng);
		p.m_Hp = dist_mHp(rng);
		p.sba = dist_sba(rng);
		p.lambda6 = dist_l6l7(rng);
		p.lambda7 = dist_l6l7(rng);
		p.m12_2 = dist_m12(rng);
		p.tan_beta = dist_tb(rng);

		std::string failure;
		if (!run_point(p, i, rtol, atol, max_abs, max_rel, failure)) {
			std::cout << "FAIL: " << failure << "\n";
			return 2;
		}
	}

	THDM bad;
	bool guard_ok = true;
	guard_ok &= !bad.set_param_phys_lam1(125.0, 300.0, 350.0, 360.0,
					 0.8, 1.0, 0.1, -0.2, -1.0);
	guard_ok &= !bad.set_param_phys_lam1(125.0, 300.0, 350.0, 360.0,
					 1.2, 1.0, 0.1, -0.2, 2.0);
	guard_ok &= !bad.set_param_phys_lam1(300.0, 125.0, 350.0, 360.0,
					 0.8, 1.0, 0.1, -0.2, 2.0);

	if (!guard_ok) {
		std::cout << "FAIL: invalid-input guards were not enforced\n";
		return 3;
	}

	std::cout << "PASS: " << n_points << " random round-trip points validated\n";
	std::cout << "PASS: invalid-input guards validated (3 cases)\n";
	std::cout << "Max abs diff: " << std::setprecision(17) << max_abs << "\n";
	std::cout << "Max rel diff: " << std::setprecision(17) << max_rel << "\n";

	return 0;
}
