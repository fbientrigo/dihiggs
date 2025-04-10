void plot_h_br() {
    TFile *file = TFile::Open("Demo_out.root");
    TTree *tree = (TTree*)file->Get("DecayTable");

    // Crear un histograma
    TH1F *h_br = new TH1F("h_br", "Branching Ratios for h", 20, 0, 1);
    tree->Draw("BR >> h_br", "Particle == 'h'");
    h_br->GetXaxis()->SetTitle("Branching Ratio (BR)");
    h_br->GetYaxis()->SetTitle("Counts");
    h_br->Draw();

    // Guardar como imagen
    TCanvas *c = new TCanvas();
    h_br->Draw();
    c->SaveAs("h_br_histogram.png");
}
