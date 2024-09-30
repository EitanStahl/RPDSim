void Macro() {
    gROOT->ProcessLine(".L zdcTree.h");         // Load header file
    gROOT->ProcessLine(".L zdcTree.C");       // Load function implementations
    gROOT->ProcessLine(".L Analysis.C");// Load script containing analysis function
    gROOT->ProcessLine("Analysis()");
    // gROOT->ProcessLine(".q");
}
