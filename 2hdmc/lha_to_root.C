#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void lha_to_root() {
    // Abrir el archivo Les Houches
    std::ifstream infile("Demo_out.lha");
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open LHA file." << std::endl;
        return;
    }

    // Crear el archivo ROOT y el TTree
    TFile *outfile = new TFile("Demo_out.root", "RECREATE");
    TTree *tree = new TTree("DecayTable", "Decay Table from LHA");

    // Variables para almacenar datos
    std::string particle = "";
    std::string channel = "";
    float width = 0.0;
    float br = 0.0;

    // Conectar las variables al TTree
    tree->Branch("Particle", &particle);
    tree->Branch("Width", &width, "Width/F");
    tree->Branch("BR", &br, "BR/F");
    tree->Branch("DecayChannel", &channel);

    // Leer el archivo línea por línea
    std::string line;
    while (std::getline(infile, line)) {
        std::cout << "Processing line: " << line << std::endl;

        // Buscar bloques de decaimiento
        if (line.find("DECAY") != std::string::npos) {
            std::istringstream iss(line);
            std::string dummy;
            iss >> dummy >> particle >> width; // Leer el PDG y la anchura
            std::cout << "Found particle: " << particle << ", Width: " << width << std::endl;
        } else if (!line.empty() && isdigit(line[0])) {
            std::istringstream iss(line);
            iss >> br; // Leer BR
            int nda;
            iss >> nda; // Leer NDA (número de partículas finales), ignorar
            std::getline(iss, channel); // Leer canal de decaimiento
            channel.erase(0, channel.find_first_not_of(" \t")); // Eliminar espacios al inicio
            std::cout << "BR: " << br << ", Channel: " << channel << std::endl;
            tree->Fill(); // Llenar el árbol
        }
    }

    // Guardar el TTree en el archivo ROOT
    tree->Write();
    outfile->Close();

    std::cout << "Decay table saved to Demo_out.root" << std::endl;
}
