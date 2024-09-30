#define zdcTree_cxx
#include "zdcTree.h"
#include <TH2.h>
#include <TStyle.h>
// #include <TCanvas.h>

zdcTree::zdcTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../NTUP.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../NTUP.root");
      }
      f->GetObject("zdcTree",tree);

   }
   Init(tree);
}

zdcTree::~zdcTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t zdcTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t zdcTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void zdcTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   zdc_ZdcTruthParticlePosx = 0;
   zdc_ZdcTruthParticlePosy = 0;
   zdc_ZdcTruthParticlePosz = 0;
   zdc_ZdcTruthParticleTime = 0;
   zdc_ZdcTruthParticlePx = 0;
   zdc_ZdcTruthParticlePy = 0;
   zdc_ZdcTruthParticlePz = 0;
   zdc_ZdcTruthParticleEnergy = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("bunchGroup", &bunchGroup, &b_bunchGroup);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("avgIntPerCrossing", &avgIntPerCrossing, &b_avgIntPerCrossing);
   fChain->SetBranchAddress("actIntPerCrossing", &actIntPerCrossing, &b_actIntPerCrossing);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("trigger_TBP", &trigger_TBP, &b_trigger_TBP);
   fChain->SetBranchAddress("tbp", tbp, &b_tbp);
   fChain->SetBranchAddress("tav", tav, &b_tav);
   fChain->SetBranchAddress("passBits", &passBits, &b_passBits);
   fChain->SetBranchAddress("extendedLevel1ID", &extendedLevel1ID, &b_extendedLevel1ID);
   fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   fChain->SetBranchAddress("timeStampNSOffset", &timeStampNSOffset, &b_timeStampNSOffset);
   fChain->SetBranchAddress("zdcEventInfoError", &zdcEventInfoError, &b_zdcEventInfoError);
   fChain->SetBranchAddress("zdcEventInfoErrorWord", &zdcEventInfoErrorWord, &b_zdcEventInfoErrorWord);
   fChain->SetBranchAddress("zdc_raw", zdc_raw, &b_zdc_raw);
   fChain->SetBranchAddress("rpd_raw", rpd_raw, &b_rpd_raw);
   fChain->SetBranchAddress("zdc_ZdcAmp", zdc_ZdcAmp, &b_zdc_ZdcAmp);
   fChain->SetBranchAddress("zdc_ZdcAmpErr", zdc_ZdcAmpErr, &b_zdc_ZdcAmpErr);
   fChain->SetBranchAddress("zdc_ZdcEnergy", zdc_ZdcEnergy, &b_zdc_ZdcEnergy);
   fChain->SetBranchAddress("zdc_ZdcEnergyErr", zdc_ZdcEnergyErr, &b_zdc_ZdcEnergyErr);
   fChain->SetBranchAddress("zdc_ZdcTime", zdc_ZdcTime, &b_zdc_ZdcTime);
   fChain->SetBranchAddress("zdc_ZdcStatus", zdc_ZdcStatus, &b_zdc_ZdcStatus);
   fChain->SetBranchAddress("zdc_ZdcTrigEff", zdc_ZdcTrigEff, &b_zdc_ZdcTrigEff);
   fChain->SetBranchAddress("zdc_ZdcModuleMask", &zdc_ZdcModuleMask, &b_zdc_ZdcModuleMask);
   fChain->SetBranchAddress("zdc_ZdcLucrodTriggerSideAmp", zdc_ZdcLucrodTriggerSideAmp, &b_zdc_ZdcLucrodTriggerSideAmp);
   fChain->SetBranchAddress("zdc_ZdcModuleAmp", zdc_ZdcModuleAmp, &b_zdc_ZdcModuleAmp);
   fChain->SetBranchAddress("zdc_ZdcModuleTime", zdc_ZdcModuleTime, &b_zdc_ZdcModuleTime);
   fChain->SetBranchAddress("zdc_ZdcModuleFitAmp", zdc_ZdcModuleFitAmp, &b_zdc_ZdcModuleFitAmp);
   fChain->SetBranchAddress("zdc_ZdcModuleFitT0", zdc_ZdcModuleFitT0, &b_zdc_ZdcModuleFitT0);
   fChain->SetBranchAddress("zdc_ZdcModuleStatus", zdc_ZdcModuleStatus, &b_zdc_ZdcModuleStatus);
   fChain->SetBranchAddress("zdc_ZdcModuleChisq", zdc_ZdcModuleChisq, &b_zdc_ZdcModuleChisq);
   fChain->SetBranchAddress("zdc_ZdcModuleCalibAmp", zdc_ZdcModuleCalibAmp, &b_zdc_ZdcModuleCalibAmp);
   fChain->SetBranchAddress("zdc_ZdcModuleCalibTime", zdc_ZdcModuleCalibTime, &b_zdc_ZdcModuleCalibTime);
   fChain->SetBranchAddress("zdc_ZdcModuleBkgdMaxFraction", zdc_ZdcModuleBkgdMaxFraction, &b_zdc_ZdcModuleBkgdMaxFraction);
   fChain->SetBranchAddress("zdc_ZdcModuleAmpError", zdc_ZdcModuleAmpError, &b_zdc_ZdcModuleAmpError);
   fChain->SetBranchAddress("zdc_ZdcModuleMinDeriv2nd", zdc_ZdcModuleMinDeriv2nd, &b_zdc_ZdcModuleMinDeriv2nd);
   fChain->SetBranchAddress("zdc_ZdcModulePresample", zdc_ZdcModulePresample, &b_zdc_ZdcModulePresample);
   fChain->SetBranchAddress("zdc_ZdcModulePreSampleAmp", zdc_ZdcModulePreSampleAmp, &b_zdc_ZdcModulePreSampleAmp);
   fChain->SetBranchAddress("zdc_ZdcLucrodTriggerAmp", zdc_ZdcLucrodTriggerAmp, &b_zdc_ZdcLucrodTriggerAmp);
   fChain->SetBranchAddress("zdc_ZdcModuleMaxADC", zdc_ZdcModuleMaxADC, &b_zdc_ZdcModuleMaxADC);
   fChain->SetBranchAddress("zdc_ZdcModuleTruthTotal", zdc_ZdcModuleTruthTotal, &b_zdc_ZdcModuleTruthTotal);
   fChain->SetBranchAddress("zdc_ZdcModuleTruthInvisible", zdc_ZdcModuleTruthInvisible, &b_zdc_ZdcModuleTruthInvisible);
   fChain->SetBranchAddress("zdc_ZdcModuleTruthEM", zdc_ZdcModuleTruthEM, &b_zdc_ZdcModuleTruthEM);
   fChain->SetBranchAddress("zdc_ZdcModuleTruthNonEM", zdc_ZdcModuleTruthNonEM, &b_zdc_ZdcModuleTruthNonEM);
   fChain->SetBranchAddress("zdc_ZdcModuleTruthEscaped", zdc_ZdcModuleTruthEscaped, &b_zdc_ZdcModuleTruthEscaped);
   fChain->SetBranchAddress("zdc_ZdcModuleTruthNphotons", zdc_ZdcModuleTruthNphotons, &b_zdc_ZdcModuleTruthNphotons);
   fChain->SetBranchAddress("zdc_RpdModuleTruthNphotons", zdc_RpdModuleTruthNphotons, &b_zdc_RpdModuleTruthNphotons);
   fChain->SetBranchAddress("zdc_ZdcTruthTotal", zdc_ZdcTruthTotal, &b_zdc_ZdcTruthTotal);
   fChain->SetBranchAddress("zdc_ZdcTruthInvisible", zdc_ZdcTruthInvisible, &b_zdc_ZdcTruthInvisible);
   fChain->SetBranchAddress("zdc_ZdcTruthEM", zdc_ZdcTruthEM, &b_zdc_ZdcTruthEM);
   fChain->SetBranchAddress("zdc_ZdcTruthNonEM", zdc_ZdcTruthNonEM, &b_zdc_ZdcTruthNonEM);
   fChain->SetBranchAddress("zdc_ZdcTruthEscaped", zdc_ZdcTruthEscaped, &b_zdc_ZdcTruthEscaped);
   fChain->SetBranchAddress("zdc_ZdcTruthParticlePosx", &zdc_ZdcTruthParticlePosx, &b_zdc_ZdcTruthParticlePosx);
   fChain->SetBranchAddress("zdc_ZdcTruthParticlePosy", &zdc_ZdcTruthParticlePosy, &b_zdc_ZdcTruthParticlePosy);
   fChain->SetBranchAddress("zdc_ZdcTruthParticlePosz", &zdc_ZdcTruthParticlePosz, &b_zdc_ZdcTruthParticlePosz);
   fChain->SetBranchAddress("zdc_ZdcTruthParticleTime", &zdc_ZdcTruthParticleTime, &b_zdc_ZdcTruthParticleTime);
   fChain->SetBranchAddress("zdc_ZdcTruthParticlePx", &zdc_ZdcTruthParticlePx, &b_zdc_ZdcTruthParticlePx);
   fChain->SetBranchAddress("zdc_ZdcTruthParticlePy", &zdc_ZdcTruthParticlePy, &b_zdc_ZdcTruthParticlePy);
   fChain->SetBranchAddress("zdc_ZdcTruthParticlePz", &zdc_ZdcTruthParticlePz, &b_zdc_ZdcTruthParticlePz);
   fChain->SetBranchAddress("zdc_ZdcTruthParticleEnergy", &zdc_ZdcTruthParticleEnergy, &b_zdc_ZdcTruthParticleEnergy);
   fChain->SetBranchAddress("zdc_RpdChannelBaseline", zdc_RpdChannelBaseline, &b_zdc_RpdChannelBaseline);
   fChain->SetBranchAddress("zdc_RpdChannelPileupExpFitParams", zdc_RpdChannelPileupExpFitParams, &b_zdc_RpdChannelPileupExpFitParams);
   fChain->SetBranchAddress("zdc_RpdChannelPileupStretchedExpFitParams", zdc_RpdChannelPileupStretchedExpFitParams, &b_zdc_RpdChannelPileupStretchedExpFitParams);
   fChain->SetBranchAddress("zdc_RpdChannelPileupExpFitParamErrs", zdc_RpdChannelPileupExpFitParamErrs, &b_zdc_RpdChannelPileupExpFitParamErrs);
   fChain->SetBranchAddress("zdc_RpdChannelPileupStretchedExpFitParamErrs", zdc_RpdChannelPileupStretchedExpFitParamErrs, &b_zdc_RpdChannelPileupStretchedExpFitParamErrs);
   fChain->SetBranchAddress("zdc_RpdChannelPileupExpFitMSE", zdc_RpdChannelPileupExpFitMSE, &b_zdc_RpdChannelPileupExpFitMSE);
   fChain->SetBranchAddress("zdc_RpdChannelPileupStretchedExpFitMSE", zdc_RpdChannelPileupStretchedExpFitMSE, &b_zdc_RpdChannelPileupStretchedExpFitMSE);
   fChain->SetBranchAddress("zdc_RpdChannelAmplitude", zdc_RpdChannelAmplitude, &b_zdc_RpdChannelAmplitude);
   fChain->SetBranchAddress("zdc_RpdChannelAmplitudeCalib", zdc_RpdChannelAmplitudeCalib, &b_zdc_RpdChannelAmplitudeCalib);
   fChain->SetBranchAddress("zdc_RpdChannelMaxADC", zdc_RpdChannelMaxADC, &b_zdc_RpdChannelMaxADC);
   fChain->SetBranchAddress("zdc_RpdChannelMaxADCCalib", zdc_RpdChannelMaxADCCalib, &b_zdc_RpdChannelMaxADCCalib);
   fChain->SetBranchAddress("zdc_RpdChannelMaxSample", zdc_RpdChannelMaxSample, &b_zdc_RpdChannelMaxSample);
   fChain->SetBranchAddress("zdc_RpdChannelStatus", zdc_RpdChannelStatus, &b_zdc_RpdChannelStatus);
   fChain->SetBranchAddress("zdc_RpdChannelPileupFrac", zdc_RpdChannelPileupFrac, &b_zdc_RpdChannelPileupFrac);
   fChain->SetBranchAddress("zdc_RpdSideStatus", zdc_RpdSideStatus, &b_zdc_RpdSideStatus);
   fChain->SetBranchAddress("zdc_centroidEventValid", &zdc_centroidEventValid, &b_zdc_centroidEventValid);
   fChain->SetBranchAddress("zdc_centroidStatus", zdc_centroidStatus, &b_zdc_centroidStatus);
   fChain->SetBranchAddress("zdc_RPDChannelSubtrAmp", zdc_RPDChannelSubtrAmp, &b_zdc_RPDChannelSubtrAmp);
   fChain->SetBranchAddress("zdc_RPDSubtrAmpSum", zdc_RPDSubtrAmpSum, &b_zdc_RPDSubtrAmpSum);
   fChain->SetBranchAddress("zdc_xCentroidPreGeomCorPreAvgSubtr", zdc_xCentroidPreGeomCorPreAvgSubtr, &b_zdc_xCentroidPreGeomCorPreAvgSubtr);
   fChain->SetBranchAddress("zdc_yCentroidPreGeomCorPreAvgSubtr", zdc_yCentroidPreGeomCorPreAvgSubtr, &b_zdc_yCentroidPreGeomCorPreAvgSubtr);
   fChain->SetBranchAddress("zdc_xCentroidPreAvgSubtr", zdc_xCentroidPreAvgSubtr, &b_zdc_xCentroidPreAvgSubtr);
   fChain->SetBranchAddress("zdc_yCentroidPreAvgSubtr", zdc_yCentroidPreAvgSubtr, &b_zdc_yCentroidPreAvgSubtr);
   fChain->SetBranchAddress("zdc_xCentroid", zdc_xCentroid, &b_zdc_xCentroid);
   fChain->SetBranchAddress("zdc_yCentroid", zdc_yCentroid, &b_zdc_yCentroid);
   fChain->SetBranchAddress("zdc_xRowCentroid", zdc_xRowCentroid, &b_zdc_xRowCentroid);
   fChain->SetBranchAddress("zdc_yColCentroid", zdc_yColCentroid, &b_zdc_yColCentroid);
   fChain->SetBranchAddress("zdc_reactionPlaneAngle", zdc_reactionPlaneAngle, &b_zdc_reactionPlaneAngle);
   fChain->SetBranchAddress("zdc_cosDeltaReactionPlaneAngle", &zdc_cosDeltaReactionPlaneAngle, &b_zdc_cosDeltaReactionPlaneAngle);
   Notify();
}

Bool_t zdcTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void zdcTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t zdcTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

zdcTree::Event::Event(Long64_t jentry, zdcTree* tree){
   m_jentry = jentry;
   m_ActiveChannels = tree->GetNumOfActiveChannels(jentry);
   m_TruthCentroid = tree->GetTruthCentroid(jentry);
   m_UnsubCentroid = tree->GetUnSubtrCentroid(jentry);
   m_SubCentroid.resize(2,std::vector<Float_t>(2,0.0f));
   for(int side=0; side<2; side++){
      m_SubCentroid[side][0] = tree->zdc_xCentroid[side];
      m_SubCentroid[side][1] = tree->zdc_yCentroid[side];
      m_TruthUnsubCentroidDiff[side][0] = (m_TruthCentroid[side][0] - x_RPD[side]) - m_UnsubCentroid[side][0];
      m_TruthUnsubCentroidDiff[side][1] = (m_TruthCentroid[side][1] - y_RPD[side]) - m_UnsubCentroid[side][1];
      m_TruthSubCentroidDiff[side][0] = (m_TruthCentroid[side][0] - x_RPD[side]) - m_SubCentroid[side][0];
      m_TruthSubCentroidDiff[side][1] = (m_TruthCentroid[side][1] - y_RPD[side]) - m_SubCentroid[side][1];
      for(int ch=0; ch<number_channels; ch++){
         m_UnsbtrAmplitude[side][ch] = tree->zdc_RpdChannelAmplitude[side][ch];
      }
   }
}

zdcTree::Event::~Event() {
}

std::array<Long64_t,2> zdcTree::GetNumOfActiveChannels(Long_t jentry){
   // zdcTree::b_zdc_RpdChannelAmplitude->GetEntry(entry);
   std::array<Long64_t,2> countCh = {0,0};
   for(int side=0; side<2; side++){
      for(int ch=0; ch<number_channels; ch++){
         if(zdc_RpdChannelAmplitude[side][ch] != 0) countCh[side]++;
      }
   }
   // std::cout << "countCh " << countCh << "\n";
   return countCh;
}

std::vector<std::vector<Float_t>> zdcTree::GetTruthCentroid(Long64_t jentry){
   // std::cout << "Calling GetTruthCentroid\n";
   // BranchGetEntry(jentry);
   std::vector<std::vector<Float_t>> TruthCentroid(2, std::vector<Float_t>(2,0.0f));

   Float_t px_C = zdc_ZdcTruthParticlePx->at(2);
   Float_t py_C = zdc_ZdcTruthParticlePy->at(2);
   Float_t pz_C = zdc_ZdcTruthParticlePz->at(2);

   Float_t px_A = zdc_ZdcTruthParticlePx->at(3);
   Float_t py_A = zdc_ZdcTruthParticlePy->at(3);
   Float_t pz_A = zdc_ZdcTruthParticlePz->at(3);

   Float_t x_C = (px_C/pz_C)*z_RPD[0];
   Float_t y_C = (py_C/pz_C)*z_RPD[0];
   Float_t x_A = (px_A/pz_A)*z_RPD[1];
   Float_t y_A = (py_A/pz_A)*z_RPD[1];

   TruthCentroid[0] = std::vector<Float_t>{x_C,y_C};
   TruthCentroid[1] = std::vector<Float_t>{x_A,y_A};

   // std::cout << "GetTruthCentroid is done\n";
   
   return TruthCentroid;
}

std::vector<std::vector<Float_t>> zdcTree::GetUnSubtrCentroid(Long64_t jentry){
   // std::cout << "Calling GetUnSubtrCentroid\n";
   std::vector<std::vector<Float_t>> UnSubtrCentroid(2, std::vector<Float_t>(2,0.0f));

   Float_t x_Centroid = 0;
   Float_t y_Centroid = 0;
   Float_t sum_Amplitude = 0;

   // BranchGetEntry(jentry);
   for(int side=0; side<2; side++){
      for(int ch=0; ch<number_channels; ch++){
         x_Centroid += m_Channel_center[side][ch][0]*zdc_RpdChannelAmplitude[side][ch];
         y_Centroid += m_Channel_center[side][ch][1]*zdc_RpdChannelAmplitude[side][ch];
         sum_Amplitude += zdc_RpdChannelAmplitude[side][ch];
      }
      x_Centroid /= sum_Amplitude;
      y_Centroid /= sum_Amplitude;
      UnSubtrCentroid[side] = std::vector<Float_t>{x_Centroid,y_Centroid};
   }
   // std::cout << "GetUnSubtrCentroid is done\n";
   return UnSubtrCentroid;
}

bool zdcTree::IsRowActive(const Event& event,const int& side, const int& row){
   for(int col=0; col<4; col++){
      if(event.m_UnsbtrAmplitude[side][m_ChannelsByRow[row][col]] != 0) return true;
   }
   return false;
}

void zdcTree::TruthRecoDifByRows(const std::vector<Event>& V, const int& side){
   std::cout << "Entering TruthRecoDifByRows\n";
   TLegend* legend = new TLegend(0.1,0.9,0.4,0.95);
   std::cout << "Creating legend\n";
   TString SavePath, hist_name;
   TH1F* TruthUnsubtXDiff [16][4];
   TH1F* TruthUnsubtYDiff [16][4];
   TH1F* TruthSubtXDiff [16][4];
   TH1F* TruthSubtYDiff [16][4];
   std::cout << "Created histogram array\n";

   for(int active_channels=1; active_channels<17; active_channels++){
      std::cout << "active_channels: " << active_channels << "\n";
      SavePath = Form("Results/%d_Side/%d-Active_Channels",side,active_channels);
      gSystem->Exec(Form("mkdir -p %s", SavePath.Data()));
      for(int row=0; row<4; row++){
         hist_name = Form("TruthUnsubtXDiff_%d_%d_row",active_channels,row);
         TruthUnsubtXDiff[active_channels - 1][row] = new TH1F(hist_name.Data(),Form("%s;Truth X_{Centroid}-Unsubt X_{Centroid}[mm];counts",hist_name.Data()),30,-60,60);
         TruthUnsubtXDiff[active_channels - 1][row]->SetFillColor(kRed);
         TruthUnsubtXDiff[active_channels - 1][row]->SetFillStyle(3001);
         hist_name = Form("TruthUnsubtYDiff_%d_%d_row",active_channels,row);
         TruthUnsubtYDiff[active_channels - 1][row] = new TH1F(hist_name.Data(),Form("%s;Truth Y_{Centroid}-Unsubt Y_{Centroid}[mm];counts",hist_name.Data()),30,-60,60);
         TruthUnsubtYDiff[active_channels - 1][row]->SetFillColor(kRed);
         TruthUnsubtYDiff[active_channels - 1][row]->SetFillStyle(3001);
         hist_name = Form("TruthSubtXDiff_%d_%d_row",active_channels,row);
         TruthSubtXDiff[active_channels - 1][row] = new TH1F(hist_name.Data(),Form("%s;Truth X_{Centroid}-Subt X_{Centroid}[mm];counts",hist_name.Data()),30,-60,60);
         TruthSubtXDiff[active_channels - 1][row]->SetLineColor(kBlue);
         hist_name = Form("TruthSubtYDiff_%d_%d_row",active_channels,row);
         TruthSubtYDiff[active_channels - 1][row] = new TH1F(hist_name.Data(),Form("%s;Truth Y_{Centroid}-Subt Y_{Centroid}[mm];counts",hist_name.Data()),30,-60,60);
         TruthSubtYDiff[active_channels - 1][row]->SetLineColor(kBlue);
      }
   }
   std::cout << "Finished creatnig histograms\n";
   for(const auto& event : V){
      for(int active_channels=1; active_channels<17; active_channels++){
         if(event.m_ActiveChannels[side] != active_channels) continue;
            for(int row=0; row<4; row++){
               if(!IsRowActive(event, side, row)) continue;
               TruthUnsubtXDiff[active_channels - 1][row]->Fill(event.m_TruthUnsubCentroidDiff[side][0]);
               TruthUnsubtYDiff[active_channels - 1][row]->Fill(event.m_TruthUnsubCentroidDiff[side][1]);
               TruthSubtXDiff[active_channels - 1][row]->Fill(event.m_TruthSubCentroidDiff[side][0]);
               TruthSubtYDiff[active_channels - 1][row]->Fill(event.m_TruthSubCentroidDiff[side][1]);  
            }
      }
   }
   std::cout << "Finished filling histograms\n";

   TFile* file;
   TString filePath;
   std::cout << "Creating TCanvas\n";
   TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
   std::cout << "TCanvas Created\n";
   for(int active_channels=0; active_channels<number_channels; active_channels++){
      SavePath = Form("Results/%d_Side/%d-Active_Channels",side,active_channels + 1);
      for(int row=0; row<4; row++){
         filePath = Form("%s/%d_row.root",SavePath.Data(),row);
         file = new TFile(filePath.Data(),"RECREATE");
         if(TruthUnsubtXDiff[active_channels][row]) TruthUnsubtXDiff[active_channels][row]->Write();
         if(TruthUnsubtYDiff[active_channels][row]) TruthUnsubtYDiff[active_channels][row]->Write();
         if(TruthSubtXDiff[active_channels][row]) TruthSubtXDiff[active_channels][row]->Write();
         if(TruthSubtYDiff[active_channels][row]) TruthSubtYDiff[active_channels][row]->Write();

         if(TruthUnsubtXDiff[active_channels][row] && TruthSubtXDiff[active_channels][row]){
            c1->cd();
            c1->SetName("X");
            legend->AddEntry(TruthUnsubtXDiff[active_channels][row],"TruthUnsubtXDiff");
            legend->AddEntry(TruthSubtXDiff[active_channels][row],"TruthSubtXDiff");
            if(TruthUnsubtXDiff[active_channels][row]->GetMaximum() > TruthSubtXDiff[active_channels][row]->GetMaximum()){
               TruthUnsubtXDiff[active_channels][row]->SetTitle("");
               TruthUnsubtXDiff[active_channels][row]->GetXaxis()->SetTitle("Truth X_{Centroid} - Reco X_{Centroid}[mm]");
               TruthUnsubtXDiff[active_channels][row]->Draw(); 
               TruthSubtXDiff[active_channels][row]->Draw("SAME");
            }
            else{
               TruthSubtXDiff[active_channels][row]->SetTitle("");
               TruthSubtXDiff[active_channels][row]->GetXaxis()->SetTitle("Truth X_{Centroid} - Reco X_{Centroid}[mm]");
               TruthSubtXDiff[active_channels][row]->Draw();
               TruthUnsubtXDiff[active_channels][row]->Draw("SAME"); 
            }
            legend->Draw();
            c1->Write();
            legend->Clear();
            c1->Clear();
         }
         
         if(TruthUnsubtYDiff[active_channels][row] && TruthSubtYDiff[active_channels][row]){
            c1->cd();
            c1->SetName("Y");      
            legend->AddEntry(TruthUnsubtYDiff[active_channels][row],"TruthUnsubtYDiff");
            legend->AddEntry(TruthSubtYDiff[active_channels][row],"TruthSubtYDiff");
            if(TruthUnsubtYDiff[active_channels][row]->GetMaximum() > TruthSubtYDiff[active_channels][row]->GetMaximum()){
               TruthUnsubtYDiff[active_channels][row]->SetTitle("");
               TruthUnsubtYDiff[active_channels][row]->GetXaxis()->SetTitle("Truth Y_{Centroid} - Reco Y_{Centroid}[mm]");
               TruthUnsubtYDiff[active_channels][row]->Draw();
               TruthSubtYDiff[active_channels][row]->Draw("SAME");
            }
            else{
               TruthSubtYDiff[active_channels][row]->SetTitle("");
               TruthSubtYDiff[active_channels][row]->GetXaxis()->SetTitle("Truth Y_{Centroid} - Reco Y_{Centroid}[mm]");
               TruthSubtYDiff[active_channels][row]->Draw();
               TruthUnsubtYDiff[active_channels][row]->Draw("SAME");   
            }
            legend->Draw();
            c1->Write();
            legend->Clear();
            c1->Clear();
         }

         file->Close();
         // file->Clear();
         std::cout << "Finised active_channels " << active_channels << "\n";
      }
   }
   std::cout << "Finished TruthRecoDifByRows\n";
}

void zdcTree::EventDisplay(const Event& event, const int& side){
   std::cout << "Entering EventDisplay\n";
   TFile* file;
   TString SavePath = Form("Results/%d_Side/EventDisplay/%lld_Active_Channels",side,event.m_ActiveChannels[side]);
   TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
   TH2F* EventDisplayUnSubtr = new TH2F("EventDisplayUnSubtr","Unsubtracted Amplitude Heat Map;x[mm];y[mm]",4,-22.8,22.8,4,-22.8,22.8);
   TMarker* TruthCentroidMark = new TMarker(event.m_TruthCentroid[side][0] - x_RPD[side],event.m_TruthCentroid[side][1] - y_RPD[side],20);
   TMarker* UnsubtCentroidMark= new TMarker(event.m_UnsubCentroid[side][0],event.m_UnsubCentroid[side][1],20);
   TLegend* legend = new TLegend(0.1, 0.9, 0.4, 0.95);
   gSystem->Exec(Form("mkdir -p %s", SavePath.Data()));
   // TH2F* EventDisplayPostSubtr = new TH2F("EventDisplayPostSubtr","Post subtraction Amplitude Heat Map;x[mm];y[mm]",4,-22.8,22.8,4,-22.8,22.8);
   for(int ch=0; ch<number_channels; ch++){
      EventDisplayUnSubtr->Fill(zdcTree::m_Channel_center[side][ch][0],zdcTree::m_Channel_center[side][ch][1],event.m_UnsbtrAmplitude[side][ch]);
   }
   
   file = new TFile(Form("%s/%lld_Unsubtr.root",SavePath.Data(),event.m_jentry),"RECREATE");
   if(!file || file->IsZombie()){
      std::cerr << "Could Not Open File\n";
      return;
   }   
   c1->cd();

   TruthCentroidMark->SetMarkerColor(kRed); // Set truth color
   TruthCentroidMark->SetMarkerSize(1.5); // Set truth size
   UnsubtCentroidMark->SetMarkerColor(kBlack);
   UnsubtCentroidMark->SetMarkerSize(1.5);
   legend->AddEntry(EventDisplayUnSubtr,Form("Entry number: %lld", event.m_jentry));
   legend->AddEntry(TruthCentroidMark,Form("Truth centroid (%.2f, %.2f)", event.m_TruthCentroid[side][0] - x_RPD[side],event.m_TruthCentroid[side][1] - y_RPD[side]),"p");
   legend->AddEntry(UnsubtCentroidMark,Form("Unsubtracted centroid (%.2f, %.2f)", event.m_UnsubCentroid[side][0],event.m_UnsubCentroid[side][1]), "p");
   EventDisplayUnSubtr->Draw("TEXT COL");
   TruthCentroidMark->Draw("SAME");
   UnsubtCentroidMark->Draw("SAME");
   DrawLines();
   legend->Draw();

   c1->Write();
   // legend->Clear();
   c1->Clear();
   file->Close();
   file->Clear();
   delete EventDisplayUnSubtr;
   delete legend;
   delete c1;
   delete file;
   std::cout << "Finished EventDisplay\n";
}

void zdcTree::DrawLines(){
   TLine *lines[6];
   lines[0] = new TLine(-11.4,22.8,-11.4,-22.8);
   lines[1] = new TLine(0,22.8,0,-22.8);
   lines[2] = new TLine(11.4,22.8,11.4,-22.8);
   lines[3] = new TLine(-22.8,-11.4,22.8,-11.4);
   lines[4] = new TLine(-22.8,0,22.8,0);
   lines[5] = new TLine(-22.8,11.4,22.8,11.4);

   // Set line attributes (optional)
   for (int i = 0; i < 6; ++i) {
      lines[i]->SetLineColor(kBlack); // Different color for each line
      lines[i]->SetLineWidth(2);        // Line width
   }

   // Draw the lines
   for (int i = 0; i < 6; ++i) {
      lines[i]->Draw("same");
   }
}

void zdcTree::NumberOfEventsAsActiveChannels(const std::vector<Event>& V, const int& side){
   std::cout << "Entering NumberOfEventsAsActiveChannels\n";
   TString SavePath = Form("Results/%d_Side", side);
   TCanvas* c1 = new TCanvas("c1","c1",800,800);
   gSystem->Exec(Form("mkdir -p %s",SavePath.Data()));
   Long64_t number_events_active_channels[number_channels+1] = {0};
   TH1F* NumberOfEventsAsActiveChannels = new TH1F("NumberOfEventsAsActiveChannels","Number of events as a function of active channels; number of active channels; N/N_{tot}",number_channels+1,-0.5,16.5);
   TFile* file = new TFile(Form("%s/NumberOfEventsAsActiveChannels.root",SavePath.Data()),"RECREATE");
   for(const auto& event : V){
      number_events_active_channels[event.m_ActiveChannels[side]]++;
      // NumberOfEventsAsActiveChannels->Fill(event.m_ActiveChannels[side]);
   }
   for(const auto a : number_events_active_channels){
      std::cout << a << "\n";
   }
   for(int bin=0; bin<NumberOfEventsAsActiveChannels->GetNbinsX(); bin++){
      NumberOfEventsAsActiveChannels->SetBinContent(bin+1,number_events_active_channels[bin]/1e5);
   }
   NumberOfEventsAsActiveChannels->SetMarkerStyle(3);
   NumberOfEventsAsActiveChannels->Draw("P");
   // NumberOfEventsAsActiveChannels->Write();
   c1->Write("NumberOfEventsAsActiveChannels");
   file->Close();
   c1->Clear();
   // delete NumberOfEventsAsActiveChannels;
   // delete c1;
   // delete file;
   std::cout << "Finished NumberOfEventsAsActiveChannels\n";
}

void zdcTree::Loop(std::vector<zdcTree::Event>& V)
{
//   In a ROOT session, you can do:
//      root> .L zdcTree.C
//      root> zdcTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Event event(jentry, this);
      V.push_back(event);
   }
}

