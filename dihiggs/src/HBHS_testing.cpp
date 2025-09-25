#include <Higgs/Predictions.hpp>
#include <Higgs/Signals.hpp>
#include <Higgs/Bounds.hpp>

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>
#include <complex>

// Conveniencia de nombres
namespace HP = Higgs::predictions;
namespace HB = Higgs::bounds;
namespace HS = Higgs::signals;

// -----------------------------------------------------------------------------
// Higgs Predictions, toy model
// Permite la creación de un conjunto mínimo de predicciones
// Es un ejemplo sencillo de lo que puede hacer 2HDMC en memoria
// Disponible
// -----------------------------------------------------------------------------
// Construye un conjunto mínimo de predicciones:
//  - h(125) SM-like (para HS)
//  - H(300) CP-par con acoplamientos efectivos reescalados (para HB/HS)
// -----------------------------------------------------------------------------
static Higgs::Predictions make_toy_predictions() {
    Higgs::Predictions pred;

    // h ~125 GeV SM-like
    auto &h = pred.addParticle(HP::BsmParticle("h", HP::ECharge::neutral, HP::CP::even));
    h.setMass(125.0);
    // SM-like mediante effectiveCouplingInput (usa tablas de referencia internas)
    HP::effectiveCouplingInput(
        h,
        HP::smLikeEffCouplings,
        HP::ReferenceModel::SMHiggsEW, // precisión alta alrededor de 125 GeV
        /*calcggH=*/true,
        /*calcHgamgam=*/true
    );

    // H más pesado (ej. 300 GeV), CP-par neutro con reescalado sencillo
    auto &H = pred.addParticle(HP::BsmParticle("H", HP::ECharge::neutral, HP::CP::even));
    H.setMass(300.0);

    // los acoples se almacenan en una estructura dedicada:
    HP::NeutralEffectiveCouplings c{};

    // Aqui escribí acoplamientos de fermiones y bosones vectoriales a mano.
    // solo a modo de ejemplo pedagógico.
    // estos mismos acoplamientos los calcularía 2HDMC en un caso real.
    
    c.tt = std::complex<double>(0.2, 0.0);
    c.bb = std::complex<double>(0.5, 0.0);
    c.tautau = std::complex<double>(0.5, 0.0);
    c.uu = std::complex<double>(0.2, 0.0);
    c.dd = std::complex<double>(0.5, 0.0);
    c.WW = 0.3;
    c.ZZ = 0.3;


    // Dejar gg y gamgam en 0 para que se calculen con los bucles efectivos de top/W si se desea:
    HP::effectiveCouplingInput(
        H, c, HP::ReferenceModel::SMHiggsInterp, /*calcggH=*/true, /*calcHgamgam=*/true
    );

    // Ejemplo de acceso a BR/anchos de H:
    // double br_H_gaga = H.br(HP::Decay::gamgam);
    // double wtot_H    = H.totalWidth();

    return pred;
}

// -----------------------------------------------------------------------------
// Pretty-print de resultados de HiggsBounds
// -----------------------------------------------------------------------------
static void print_hb_result(const HB::HBResult &res) {
    std::cout << "\n=== HiggsBounds ===\n";
    std::cout << "allowed (95% CL): " << (res.allowed ? "YES" : "NO") << "\n";

    // Seleccionados por sensibilidad esperada (uno por partícula)
    std::cout << "selectedLimits (per particle):\n";
    for (const auto &kv : res.selectedLimits) {
        const std::string &pid = kv.first;
        const HB::AppliedLimit &app = kv.second;
        auto lim = app.limit();

        std::cout << "  - particle: " << pid
                  << " | obsRatio=" << std::fixed << std::setprecision(3) << app.obsRatio()
                  << " | expRatio=" << app.expRatio() << "\n"
                  << "    ref: " << lim->reference()
                  << " | citeKey: " << lim->citeKey() << "\n"
                  << "    process: " << lim->processDesc() << "\n"
                  << "    extent : " << lim->extentDesc() << "\n";

        // Quién contribuye (útil si hay clustering de masas)
        const auto &parts = app.contributingParticles();
        if (!parts.empty()) {
            std::cout << "    contributors: ";
            for (size_t i = 0; i < parts.size(); ++i) {
                std::cout << parts[i] << (i + 1 < parts.size() ? " " : "");
            }
            std::cout << "\n";
        }
    }

    // Los aplicados (ordenados por sensibilidad esperada). Muestra los 3 primeros.
    std::cout << "top applied limits by expected sensitivity (expRatio):\n";
    for (size_t i = 0; i < std::min<size_t>(3, res.appliedLimits.size()); ++i) {
        const auto &app = res.appliedLimits[i];
        auto lim = app.limit();
        std::cout << "  [" << i << "] expRatio=" << app.expRatio()
                  << " | obsRatio=" << app.obsRatio()
                  << " | " << lim->reference() << "\n";
    }
}

// -----------------------------------------------------------------------------
// Pretty-print de resultados de HiggsSignals
// -----------------------------------------------------------------------------
static void run_and_print_hs(HS::Signals &signals, const Higgs::Predictions &pred) {
    double chi2 = signals(pred);
    std::cout << "\n=== HiggsSignals ===\n";
    std::cout << "chi^2 total = " << std::fixed << std::setprecision(3) << chi2
              << "  (observables cargados: " << signals.observableCount() << ")\n";
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main() {
    try {
        const char* hb_dir = std::getenv("HB_DATASET");
        const char* hs_dir = std::getenv("HS_DATASET");
        if (!hb_dir || !hs_dir) {
            std::cerr << "Debe exportar HB_DATASET y HS_DATASET (rutas al dataset JSON).\n";
            return 1;
        }

        // 1) Predicciones mínimas con HiggsTools usando la herramienta de Predictions
        // es equivalente a 2HDMC en memoria
        auto pred = make_toy_predictions();
        std::cout << "Particles in Predictions: ";
        for (const auto &id : pred.particleIds()) std::cout << id << " ";
        std::cout << "\n";

        // 2) HiggsBounds
        // Comprobación de las predicciones contra los límites
        HB::Bounds bounds(hb_dir);
        std::cout << "HB: limits loaded = " << bounds.limits().size() << "\n";
        auto hbRes = bounds(pred);
        print_hb_result(hbRes);

        // 3) HiggsSignals
        // Cálculo del chi2 de las predicciones contra las medidas
        HS::Signals signals(hs_dir);
        run_and_print_hs(signals, pred);

        std::cout << "\nOK: Test básico de HiggsTools completado.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}
