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

struct Metric {
	double max_abs;
	double max_rel;
	int idx_abs;
	int idx_rel;
	double ref_abs;
	double new_abs;
	double ref_rel;
	double new_rel;

	Metric() : max_abs(0.0), max_rel(0.0), idx_abs(-1), idx_rel(-1),
		   ref_abs(0.0), new_abs(0.0), ref_rel(0.0), new_rel(0.0) {}
};

static bool check_finite(double x) {
	return std::isfinite(x);
}

static bool compare_pair(const char *name, int metric_idx,
			 double ref_v, double new_v,
			 const Point &p, int idx,
			 double rtol, double atol,
			 Metric metrics[],
			 std::string &failure) {
	const double abs_diff = std::fabs(ref_v - new_v);
	const double scale = std::max(std::fabs(ref_v), std::fabs(new_v));
	const double rel_diff = (scale > 0.0) ? (abs_diff / scale) : abs_diff;

	if (abs_diff > metrics[metric_idx].max_abs) {
		metrics[metric_idx].max_abs = abs_diff;
		metrics[metric_idx].idx_abs = idx;
		metrics[metric_idx].ref_abs = ref_v;
		metrics[metric_idx].new_abs = new_v;
	}

	if (rel_diff > metrics[metric_idx].max_rel) {
		metrics[metric_idx].max_rel = rel_diff;
		metrics[metric_idx].idx_rel = idx;
		metrics[metric_idx].ref_rel = ref_v;
		metrics[metric_idx].new_rel = new_v;
	}

	if (!check_finite(ref_v) || !check_finite(new_v)) {
		std::ostringstream oss;
		oss << "Point " << idx << " " << name << ": non-finite values"
		    << " [mH=" << p.m_H << ", mA=" << p.m_A
		    << ", mHp=" << p.m_Hp << ", sba=" << p.sba
		    << ", l6=" << p.lambda6 << ", l7=" << p.lambda7
		    << ", m12=" << p.m12_2 << ", tb=" << p.tan_beta << "]";
		failure = oss.str();
		return false;
	}

	if (!nearly_equal(ref_v, new_v, rtol, atol)) {
		std::ostringstream oss;
		oss << "Point " << idx << " " << name
		    << " mismatch: ref=" << std::setprecision(17) << ref_v
		    << " new=" << std::setprecision(17) << new_v
		    << " abs_diff=" << abs_diff
		    << " rel_diff=" << rel_diff
		    << " [mH=" << p.m_H << ", mA=" << p.m_A
		    << ", mHp=" << p.m_Hp << ", sba=" << p.sba
		    << ", l6=" << p.lambda6 << ", l7=" << p.lambda7
		    << ", m12=" << p.m12_2 << ", tb=" << p.tan_beta << "]";
		failure = oss.str();
		return false;
	}

	return true;
}

static bool run_point(const Point &p, int idx,
		      double rtol, double atol,
		      Metric metrics[], std::string &failure) {
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

	if (!compare_pair("m12_2", 0, m12_ref, m12_new, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("tan_beta", 1, tb_ref, tb_new, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda1", 2, l1, L1, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda2", 3, l2, L2, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda3", 4, l3, L3, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda4", 5, l4, L4, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda5", 6, l5, L5, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda6", 7, l6, L6, p, idx, rtol, atol, metrics, failure) ||
	    !compare_pair("lambda7", 8, l7, L7, p, idx, rtol, atol, metrics, failure)) {
		return false;
	}

	for (int h = 1; h <= 4; ++h) {
		double m_ref_h = m_ref.get_hmass(h);
		double m_new_h = m_lam1.get_hmass(h);
		std::ostringstream name;
		name << "hmass(" << h << ")";
		const int metric_idx = 8 + h;
		if (!compare_pair(name.str().c_str(), metric_idx, m_ref_h, m_new_h,
				  p, idx, rtol, atol, metrics, failure)) {
			return false;
		}
	}

	return true;
}

int main() {
	const double rtol = 1e-10;
	const double atol = 1e-12;
	const int n_points = 400;

	const char *names[] = {
		"m12_2", "tan_beta",
		"lambda1", "lambda2", "lambda3", "lambda4", "lambda5", "lambda6", "lambda7",
		"hmass(1)", "hmass(2)", "hmass(3)", "hmass(4)"
	};
	Metric metrics[13];

	std::mt19937_64 rng(20260403ULL);
	std::uniform_real_distribution<double> dist_mH(150.0, 1200.0);
	std::uniform_real_distribution<double> dist_mA(130.0, 1200.0);
	std::uniform_real_distribution<double> dist_mHp(130.0, 1200.0);
	std::uniform_real_distribution<double> dist_sba(-0.999999, 0.999999);
	std::uniform_real_distribution<double> dist_l6l7(-2.0, 2.0);
	std::uniform_real_distribution<double> dist_m12(-8.0e5, 8.0e5);
	std::uniform_real_distribution<double> dist_tb(0.15, 60.0);

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
		if (!run_point(p, i, rtol, atol, metrics, failure)) {
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

	double global_max_abs = 0.0;
	double global_max_rel = 0.0;
	int global_abs_idx = -1;
	int global_rel_idx = -1;
	for (int i = 0; i < 13; ++i) {
		if (metrics[i].max_abs > global_max_abs) {
			global_max_abs = metrics[i].max_abs;
			global_abs_idx = i;
		}
		if (metrics[i].max_rel > global_max_rel) {
			global_max_rel = metrics[i].max_rel;
			global_rel_idx = i;
		}
	}

	std::cout << "PASS: " << n_points << " random round-trip points validated\n";
	std::cout << "PASS: invalid-input guards validated (3 cases)\n";
	std::cout << "Global max abs diff: " << std::setprecision(17) << global_max_abs
		  << " in " << names[global_abs_idx]
		  << " at point " << metrics[global_abs_idx].idx_abs << "\n";
	std::cout << "Global max rel diff: " << std::setprecision(17) << global_max_rel
		  << " in " << names[global_rel_idx]
		  << " at point " << metrics[global_rel_idx].idx_rel << "\n";

	std::cout << "\nPer-observable maxima:\n";
	for (int i = 0; i < 13; ++i) {
		std::cout << "- " << names[i]
			  << ": max_abs=" << std::setprecision(17) << metrics[i].max_abs
			  << " (point " << metrics[i].idx_abs << ")"
			  << ", max_rel=" << std::setprecision(17) << metrics[i].max_rel
			  << " (point " << metrics[i].idx_rel << ")"
			  << "\n";
	}

	return 0;
}
