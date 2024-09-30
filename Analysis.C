#include "zdcTree.h"

void Analysis(){
    TFile *file = TFile::Open("../NTUP.root");

    if (!file || !file->IsOpen()) {
        std::cerr << "Failed to open file" << std::endl;
        return;
    }
    
    TTree *tree;
    file->GetObject("zdcTree", tree);
    if (!tree) {
        std::cerr << "Failed to get tree from file" << std::endl;
        return;
    }
    // Create an instance of the zdcTree class
    zdcTree myTree(tree);
    const int side = 1;
    std::vector<zdcTree::Event> V;
    myTree.Init(tree);
    std::cout << "Creating V" << "\n";
    myTree.Loop(V);
    std::cout << "Finished creating V" << "\n";
    file->Close();
    file->Clear();
    delete file;
    gStyle->SetOptStat(0);
    gSystem->Exec("mkdir -p Results");
    // myTree.TruthRecoDifByRows(V,side);
    myTree.NumberOfEventsAsActiveChannels(V,side);
    return;
}