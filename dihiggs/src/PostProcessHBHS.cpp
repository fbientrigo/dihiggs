#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>

// HiggsTools
#include <Higgs/Predictions.hpp>
#include <Higgs/Signals.hpp>
#include <Higgs/Bounds.hpp>

// Conveniencia de nombres (idéntico a HBHS_testing.cpp)
namespace HP = Higgs::predictions;
namespace HB = Higgs::bounds;
namespace HS = Higgs::signals;

using namespace std;

struct Row {
    // columnas mínimas para reconstruir el modelo
    double m_phi, mA, alpha, beta, lambda6, lambda7, m12;
    double sin_ba, tan_beta;
    // fila original completa
    vector<string> original;
};

Row parse_row(const string &line, const vector<string> &headers) {
    Row r;
    r.original.clear();

    string token;
    stringstream ss(line);
    int col = 0;
    while (getline(ss, token, ',')) {
        r.original.push_back(token);
        switch (col) {
            case 0: r.m_phi    = stod(token); break;
            case 1: r.mA       = stod(token); break;
            case 2: r.alpha    = stod(token); break;
            case 3: r.beta     = stod(token); break;
            case 4: r.lambda6  = stod(token); break;
            case 5: r.lambda7  = stod(token); break;
            case 6: r.m12      = stod(token); break;
            case 7: r.sin_ba   = stod(token); break;
            case 8: r.tan_beta = stod(token); break;
        }
        col++;
    }
    return r;
}

// Construcción mínima de Predictions a partir de la fila del CSV
static Higgs::Predictions build_predictions_from_row(const Row &r) {
    Higgs::Predictions pred;

    // Partículas
    auto& h  = pred.addParticle(HP::BsmParticle("h",  HP::ECharge::neutral, HP::CP::even));
    auto& H  = pred.addParticle(HP::BsmParticle("H",  HP::ECharge::neutral, HP::CP::even));
    auto& A  = pred.addParticle(HP::BsmParticle("A",  HP::ECharge::neutral, HP::CP::odd));
    auto& Hp = pred.addParticle(HP::BsmParticle("H+", HP::ECharge::single, HP::CP::even));

    double mh = 125.0;
    h.setMass(mh);
    H.setMass(r.m_phi);
    A.setMass(r.mA);
    Hp.setMass(r.mA); // aprox: mHp = mA

    // ⚠️ Importante: aquí deberías poblar acoplamientos o anchos parciales.
    // Por ahora solo fijamos masas. En un pipeline real, deberías reconstruir
    // los HP::NeutralEffectiveCouplings desde 2HDMC para HB/HS consistentes.

    // Ejemplo: darle a h los acoplamientos SM-like para que HS pueda evaluarlo:
    HP::effectiveCouplingInput(
        h,
        HP::smLikeEffCouplings,
        HP::ReferenceModel::SMHiggsEW,
        /*calcggH=*/true,
        /*calcHgamgam=*/true
    );

    return pred;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Uso: " << argv[0] << " input.csv output.csv\n";
        return 1;
    }

    string input_csv  = argv[1];
    string output_csv = argv[2];

    const char* hbdir = getenv("HB_DATASET");
    const char* hsdir = getenv("HS_DATASET");
    if (!hbdir || !hsdir) {
        cerr << "ERROR: Debes exportar HB_DATASET y HS_DATASET\n";
        return 1;
    }

    HB::Bounds bounds(hbdir);
    HS::Signals signals(hsdir);

    // número de observables cargados
    size_t nObs = signals.observableCount();
    cerr << "HiggsSignals cargó " << nObs << " observables.\n";

    ifstream in(input_csv);
    ofstream out(output_csv);

    if (!in.is_open() || !out.is_open()) {
        cerr << "ERROR: no puedo abrir input o output\n";
        return 1;
    }

    string header_line;
    getline(in, header_line);
    vector<string> headers;
    {
        stringstream ss(header_line);
        string token;
        while (getline(ss, token, ',')) headers.push_back(token);
    }

    // cabecera extendida
    out << header_line << ",hb_allowed,hs_chi2,hs_ndof, hs_chi2_ndof\n";

    string line;
    int nrows = 0;
    while (getline(in, line)) {
        Row r = parse_row(line, headers);
        auto pred = build_predictions_from_row(r);

        auto hbRes = bounds(pred);
        double chi2 = signals(pred);

        // escribir fila original
        for (size_t i=0; i<r.original.size(); ++i) {
            out << r.original[i];
            if (i+1 < r.original.size()) out << ",";
        }
        // resultados HB/HS + NDOF global
        out << "," << (hbRes.allowed ? 1 : 0)
            << "," << chi2
            << "," << nObs
            << "," << (nObs > 0 ? chi2/nObs : 0.0)
            << "\n";

        if (++nrows % 100 == 0) {
            cerr << "Procesadas " << nrows << " filas...\r";
        }
    }
    cerr << "\nTerminado. Filas="<<nrows<<"\n";

    return 0;
}