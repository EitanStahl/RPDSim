//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 17 15:23:35 2024 by ROOT version 6.30/08
// from TTree zdcTree/ZDC Tree
// found on file: ../NTUP.root
//////////////////////////////////////////////////////////

#ifndef zdcTree_h
#define zdcTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
// Header file for the classes stored in the TTree if any.
#include "c++/11/vector"
using namespace std;

const Float_t z_RPD[2] = {-141489.5,141459.0}; //mm
const Float_t x_RPD[2] = {-2.0,1.8}; //mm
const Float_t y_RPD[2] = {21.4,21.3}; //mm
const Float_t x_BP_ATL[2] = {6,-5}; //mm
const Float_t y_BP_ATL[2] = {-1,1}; //mm
const int number_rows = 4;
const int number_columns = 4;
const int number_channels = number_rows*number_columns;


class zdcTree {
public :


   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          eventNumber;
   UInt_t          lumiBlock;
   UChar_t         bunchGroup;
   UInt_t          bcid;
   Float_t         avgIntPerCrossing;
   Float_t         actIntPerCrossing;
   ULong64_t       trigger;
   UInt_t          trigger_TBP;
   UInt_t          tbp[16];
   UInt_t          tav[16];
   UInt_t          passBits;
   UInt_t          extendedLevel1ID;
   UInt_t          timeStamp;
   UInt_t          timeStampNSOffset;
   UChar_t         zdcEventInfoError;
   UInt_t          zdcEventInfoErrorWord;
   UShort_t        zdc_raw[2][4][2][2][24];
   UShort_t        rpd_raw[2][16][24];
   Float_t         zdc_ZdcAmp[2];
   Float_t         zdc_ZdcAmpErr[2];
   Float_t         zdc_ZdcEnergy[2];
   Float_t         zdc_ZdcEnergyErr[2];
   Float_t         zdc_ZdcTime[2];
   Short_t         zdc_ZdcStatus[2];
   Float_t         zdc_ZdcTrigEff[2];
   UInt_t          zdc_ZdcModuleMask;
   Short_t         zdc_ZdcLucrodTriggerSideAmp[2];
   Float_t         zdc_ZdcModuleAmp[2][4];
   Float_t         zdc_ZdcModuleTime[2][4];
   Float_t         zdc_ZdcModuleFitAmp[2][4];
   Float_t         zdc_ZdcModuleFitT0[2][4];
   UInt_t          zdc_ZdcModuleStatus[2][4];
   Float_t         zdc_ZdcModuleChisq[2][4];
   Float_t         zdc_ZdcModuleCalibAmp[2][4];
   Float_t         zdc_ZdcModuleCalibTime[2][4];
   Float_t         zdc_ZdcModuleBkgdMaxFraction[2][4];
   Float_t         zdc_ZdcModuleAmpError[2][4];
   Float_t         zdc_ZdcModuleMinDeriv2nd[2][4];
   Float_t         zdc_ZdcModulePresample[2][4];
   Float_t         zdc_ZdcModulePreSampleAmp[2][4];
   Short_t         zdc_ZdcLucrodTriggerAmp[2][4];
   Float_t         zdc_ZdcModuleMaxADC[2][4];
   Float_t         zdc_ZdcModuleTruthTotal[2][7];
   Float_t         zdc_ZdcModuleTruthInvisible[2][7];
   Float_t         zdc_ZdcModuleTruthEM[2][7];
   Float_t         zdc_ZdcModuleTruthNonEM[2][7];
   Float_t         zdc_ZdcModuleTruthEscaped[2][7];
   UInt_t          zdc_ZdcModuleTruthNphotons[2][7];
   UInt_t          zdc_RpdModuleTruthNphotons[2][16];
   Float_t         zdc_ZdcTruthTotal[2];
   Float_t         zdc_ZdcTruthInvisible[2];
   Float_t         zdc_ZdcTruthEM[2];
   Float_t         zdc_ZdcTruthNonEM[2];
   Float_t         zdc_ZdcTruthEscaped[2];
   vector<float>   *zdc_ZdcTruthParticlePosx;
   vector<float>   *zdc_ZdcTruthParticlePosy;
   vector<float>   *zdc_ZdcTruthParticlePosz;
   vector<float>   *zdc_ZdcTruthParticleTime;
   vector<float>   *zdc_ZdcTruthParticlePx;
   vector<float>   *zdc_ZdcTruthParticlePy;
   vector<float>   *zdc_ZdcTruthParticlePz;
   vector<float>   *zdc_ZdcTruthParticleEnergy;
   Float_t         zdc_RpdChannelBaseline[2][16];
   Float_t         zdc_RpdChannelPileupExpFitParams[2][16][2];
   Float_t         zdc_RpdChannelPileupStretchedExpFitParams[2][16][3];
   Float_t         zdc_RpdChannelPileupExpFitParamErrs[2][16][2];
   Float_t         zdc_RpdChannelPileupStretchedExpFitParamErrs[2][16][3];
   Float_t         zdc_RpdChannelPileupExpFitMSE[2][16];
   Float_t         zdc_RpdChannelPileupStretchedExpFitMSE[2][16];
   Float_t         zdc_RpdChannelAmplitude[2][16];
   Float_t         zdc_RpdChannelAmplitudeCalib[2][16];
   Float_t         zdc_RpdChannelMaxADC[2][16];
   Float_t         zdc_RpdChannelMaxADCCalib[2][16];
   UInt_t          zdc_RpdChannelMaxSample[2][16];
   UInt_t          zdc_RpdChannelStatus[2][16];
   Float_t         zdc_RpdChannelPileupFrac[2][16];
   UInt_t          zdc_RpdSideStatus[2];
   Char_t          zdc_centroidEventValid;
   UInt_t          zdc_centroidStatus[2];
   Float_t         zdc_RPDChannelSubtrAmp[2][16];
   Float_t         zdc_RPDSubtrAmpSum[2];
   Float_t         zdc_xCentroidPreGeomCorPreAvgSubtr[2];
   Float_t         zdc_yCentroidPreGeomCorPreAvgSubtr[2];
   Float_t         zdc_xCentroidPreAvgSubtr[2];
   Float_t         zdc_yCentroidPreAvgSubtr[2];
   Float_t         zdc_xCentroid[2];
   Float_t         zdc_yCentroid[2];
   Float_t         zdc_xRowCentroid[2][4];
   Float_t         zdc_yColCentroid[2][4];
   Float_t         zdc_reactionPlaneAngle[2];
   Float_t         zdc_cosDeltaReactionPlaneAngle;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_bunchGroup;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_avgIntPerCrossing;   //!
   TBranch        *b_actIntPerCrossing;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_trigger_TBP;   //!
   TBranch        *b_tbp;   //!
   TBranch        *b_tav;   //!
   TBranch        *b_passBits;   //!
   TBranch        *b_extendedLevel1ID;   //!
   TBranch        *b_timeStamp;   //!
   TBranch        *b_timeStampNSOffset;   //!
   TBranch        *b_zdcEventInfoError;   //!
   TBranch        *b_zdcEventInfoErrorWord;   //!
   TBranch        *b_zdc_raw;   //!
   TBranch        *b_rpd_raw;   //!
   TBranch        *b_zdc_ZdcAmp;   //!
   TBranch        *b_zdc_ZdcAmpErr;   //!
   TBranch        *b_zdc_ZdcEnergy;   //!
   TBranch        *b_zdc_ZdcEnergyErr;   //!
   TBranch        *b_zdc_ZdcTime;   //!
   TBranch        *b_zdc_ZdcStatus;   //!
   TBranch        *b_zdc_ZdcTrigEff;   //!
   TBranch        *b_zdc_ZdcModuleMask;   //!
   TBranch        *b_zdc_ZdcLucrodTriggerSideAmp;   //!
   TBranch        *b_zdc_ZdcModuleAmp;   //!
   TBranch        *b_zdc_ZdcModuleTime;   //!
   TBranch        *b_zdc_ZdcModuleFitAmp;   //!
   TBranch        *b_zdc_ZdcModuleFitT0;   //!
   TBranch        *b_zdc_ZdcModuleStatus;   //!
   TBranch        *b_zdc_ZdcModuleChisq;   //!
   TBranch        *b_zdc_ZdcModuleCalibAmp;   //!
   TBranch        *b_zdc_ZdcModuleCalibTime;   //!
   TBranch        *b_zdc_ZdcModuleBkgdMaxFraction;   //!
   TBranch        *b_zdc_ZdcModuleAmpError;   //!
   TBranch        *b_zdc_ZdcModuleMinDeriv2nd;   //!
   TBranch        *b_zdc_ZdcModulePresample;   //!
   TBranch        *b_zdc_ZdcModulePreSampleAmp;   //!
   TBranch        *b_zdc_ZdcLucrodTriggerAmp;   //!
   TBranch        *b_zdc_ZdcModuleMaxADC;   //!
   TBranch        *b_zdc_ZdcModuleTruthTotal;   //!
   TBranch        *b_zdc_ZdcModuleTruthInvisible;   //!
   TBranch        *b_zdc_ZdcModuleTruthEM;   //!
   TBranch        *b_zdc_ZdcModuleTruthNonEM;   //!
   TBranch        *b_zdc_ZdcModuleTruthEscaped;   //!
   TBranch        *b_zdc_ZdcModuleTruthNphotons;   //!
   TBranch        *b_zdc_RpdModuleTruthNphotons;   //!
   TBranch        *b_zdc_ZdcTruthTotal;   //!
   TBranch        *b_zdc_ZdcTruthInvisible;   //!
   TBranch        *b_zdc_ZdcTruthEM;   //!
   TBranch        *b_zdc_ZdcTruthNonEM;   //!
   TBranch        *b_zdc_ZdcTruthEscaped;   //!
   TBranch        *b_zdc_ZdcTruthParticlePosx;   //!
   TBranch        *b_zdc_ZdcTruthParticlePosy;   //!
   TBranch        *b_zdc_ZdcTruthParticlePosz;   //!
   TBranch        *b_zdc_ZdcTruthParticleTime;   //!
   TBranch        *b_zdc_ZdcTruthParticlePx;   //!
   TBranch        *b_zdc_ZdcTruthParticlePy;   //!
   TBranch        *b_zdc_ZdcTruthParticlePz;   //!
   TBranch        *b_zdc_ZdcTruthParticleEnergy;   //!
   TBranch        *b_zdc_RpdChannelBaseline;   //!
   TBranch        *b_zdc_RpdChannelPileupExpFitParams;   //!
   TBranch        *b_zdc_RpdChannelPileupStretchedExpFitParams;   //!
   TBranch        *b_zdc_RpdChannelPileupExpFitParamErrs;   //!
   TBranch        *b_zdc_RpdChannelPileupStretchedExpFitParamErrs;   //!
   TBranch        *b_zdc_RpdChannelPileupExpFitMSE;   //!
   TBranch        *b_zdc_RpdChannelPileupStretchedExpFitMSE;   //!
   TBranch        *b_zdc_RpdChannelAmplitude;   //!
   TBranch        *b_zdc_RpdChannelAmplitudeCalib;   //!
   TBranch        *b_zdc_RpdChannelMaxADC;   //!
   TBranch        *b_zdc_RpdChannelMaxADCCalib;   //!
   TBranch        *b_zdc_RpdChannelMaxSample;   //!
   TBranch        *b_zdc_RpdChannelStatus;   //!
   TBranch        *b_zdc_RpdChannelPileupFrac;   //!
   TBranch        *b_zdc_RpdSideStatus;   //!
   TBranch        *b_zdc_centroidEventValid;   //!
   TBranch        *b_zdc_centroidStatus;   //!
   TBranch        *b_zdc_RPDChannelSubtrAmp;   //!
   TBranch        *b_zdc_RPDSubtrAmpSum;   //!
   TBranch        *b_zdc_xCentroidPreGeomCorPreAvgSubtr;   //!
   TBranch        *b_zdc_yCentroidPreGeomCorPreAvgSubtr;   //!
   TBranch        *b_zdc_xCentroidPreAvgSubtr;   //!
   TBranch        *b_zdc_yCentroidPreAvgSubtr;   //!
   TBranch        *b_zdc_xCentroid;   //!
   TBranch        *b_zdc_yCentroid;   //!
   TBranch        *b_zdc_xRowCentroid;   //!
   TBranch        *b_zdc_yColCentroid;   //!
   TBranch        *b_zdc_reactionPlaneAngle;   //!
   TBranch        *b_zdc_cosDeltaReactionPlaneAngle;   //!

   class Event{
      public:
         Long64_t m_jentry;
         Float_t m_UnsbtrAmplitude [2][16];
         Float_t m_TruthUnsubCentroidDiff [2][2];
         Float_t m_TruthSubCentroidDiff [2][2];
         bool  m_ActiveRows [4] = {false};
         bool  m_ActiveColumns [4] = {false};
         bool m_GoodCentroid[2] = {false};

         std::array<Long64_t,2> m_ActiveChannels;
         std::vector<std::vector<Float_t>> m_TruthCentroid;
         std::vector<std::vector<Float_t>> m_UnsubCentroid;
         std::vector<std::vector<Float_t>> m_SubCentroid;

         Event(Long64_t jentry, zdcTree* tree);
         ~Event();
   };

   zdcTree(TTree *tree=0);
   virtual ~zdcTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(std::vector<zdcTree::Event>& V);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   std::array<Long64_t,2> GetNumOfActiveChannels(Long_t jentry);
   std::vector<std::vector<Float_t>> GetTruthCentroid (Long64_t jentry);
   std::vector<std::vector<Float_t>> GetUnSubtrCentroid(Long64_t jentry);
   bool IsRowActive(const Event& event,const int& side, const int& row);
   void TruthRecoDifByRows(const std::vector<Event>& V, const int& side);
   void EventDisplay(const Event& event, const int& side);
   void DrawLines();
   void NumberOfEventsAsActiveChannels(const std::vector<Event>& V, const int& side);

   const int m_ChannelsByRow [4][4] = {{12,13,14,15},{8,9,10,11},{4,5,6,7},{0,1,2,3}};
   const Float_t m_Channel_center[2][16][2] =
                     {{{-17.1,17.1}    ,{-5.7,17.1}   ,{5.7,17.1}    ,{17.1,17.1},
                       {-17.1,5.7}     ,{-5.7,5.7}    ,{5.7,5.7}     ,{17.1,5.7},
                       {-17.1,-5.7}    ,{-5.7,-5.7}   ,{5.7,-5.7}    ,{17.1,-5.7},
                       {-17.1,-17.1}   ,{-5.7,-17.1}  ,{5.7,-17.1}   ,{17.1,-17.1}},
                       
                       {{17.1,17.1}    ,{5.7,17.1}    ,{-5.7,17.1}   ,{-17.1,17.1},
                       {17.1,5.7}      ,{5.7,5.7}     ,{-5.7,5.7}    ,{-17.1,5.7},
                       {17.1,-5.7}     ,{5.7,-5.7}    ,{-5.7,-5.7}   ,{-17.1,-5.7},
                       {17.1,-17.1}    ,{5.7,-17.1}   ,{-5.7,-17.1}  ,{-17.1,-17.1}}}; 
};



#endif

#ifdef zdcTree_cxx
// zdcTree::zdcTree(TTree *tree) : fChain(0) 
// {
// // if parameter tree is not specified (or zero), connect the file
// // used to generate this class and read the Tree.
//    if (tree == 0) {
//       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../NTUP.root");
//       if (!f || !f->IsOpen()) {
//          f = new TFile("../NTUP.root");
//       }
//       f->GetObject("zdcTree",tree);

//    }
//    Init(tree);
// }

// zdcTree::~zdcTree()
// {
//    if (!fChain) return;
//    delete fChain->GetCurrentFile();
// }

// Int_t zdcTree::GetEntry(Long64_t entry)
// {
// // Read contents of entry.
//    if (!fChain) return 0;
//    return fChain->GetEntry(entry);
// }

// Long64_t zdcTree::LoadTree(Long64_t entry)
// {
// // Set the environment to read one entry
//    if (!fChain) return -5;
//    Long64_t centry = fChain->LoadTree(entry);
//    if (centry < 0) return centry;
//    if (fChain->GetTreeNumber() != fCurrent) {
//       fCurrent = fChain->GetTreeNumber();
//       Notify();
//    }
//    return centry;
// }

// void zdcTree::Init(TTree *tree)
// {
//    // The Init() function is called when the selector needs to initialize
//    // a new tree or chain. Typically here the branch addresses and branch
//    // pointers of the tree will be set.
//    // It is normally not necessary to make changes to the generated
//    // code, but the routine can be extended by the user if needed.
//    // Init() will be called many times when running on PROOF
//    // (once per file to be processed).

//    // Set object pointer
//    zdc_ZdcTruthParticlePosx = 0;
//    zdc_ZdcTruthParticlePosy = 0;
//    zdc_ZdcTruthParticlePosz = 0;
//    zdc_ZdcTruthParticleTime = 0;
//    zdc_ZdcTruthParticlePx = 0;
//    zdc_ZdcTruthParticlePy = 0;
//    zdc_ZdcTruthParticlePz = 0;
//    zdc_ZdcTruthParticleEnergy = 0;
//    // Set branch addresses and branch pointers
//    if (!tree) return;
//    fChain = tree;
//    fCurrent = -1;
//    fChain->SetMakeClass(1);

//    fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
//    fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
//    fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
//    fChain->SetBranchAddress("bunchGroup", &bunchGroup, &b_bunchGroup);
//    fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
//    fChain->SetBranchAddress("avgIntPerCrossing", &avgIntPerCrossing, &b_avgIntPerCrossing);
//    fChain->SetBranchAddress("actIntPerCrossing", &actIntPerCrossing, &b_actIntPerCrossing);
//    fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
//    fChain->SetBranchAddress("trigger_TBP", &trigger_TBP, &b_trigger_TBP);
//    fChain->SetBranchAddress("tbp", tbp, &b_tbp);
//    fChain->SetBranchAddress("tav", tav, &b_tav);
//    fChain->SetBranchAddress("passBits", &passBits, &b_passBits);
//    fChain->SetBranchAddress("extendedLevel1ID", &extendedLevel1ID, &b_extendedLevel1ID);
//    fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
//    fChain->SetBranchAddress("timeStampNSOffset", &timeStampNSOffset, &b_timeStampNSOffset);
//    fChain->SetBranchAddress("zdcEventInfoError", &zdcEventInfoError, &b_zdcEventInfoError);
//    fChain->SetBranchAddress("zdcEventInfoErrorWord", &zdcEventInfoErrorWord, &b_zdcEventInfoErrorWord);
//    fChain->SetBranchAddress("zdc_raw", zdc_raw, &b_zdc_raw);
//    fChain->SetBranchAddress("rpd_raw", rpd_raw, &b_rpd_raw);
//    fChain->SetBranchAddress("zdc_ZdcAmp", zdc_ZdcAmp, &b_zdc_ZdcAmp);
//    fChain->SetBranchAddress("zdc_ZdcAmpErr", zdc_ZdcAmpErr, &b_zdc_ZdcAmpErr);
//    fChain->SetBranchAddress("zdc_ZdcEnergy", zdc_ZdcEnergy, &b_zdc_ZdcEnergy);
//    fChain->SetBranchAddress("zdc_ZdcEnergyErr", zdc_ZdcEnergyErr, &b_zdc_ZdcEnergyErr);
//    fChain->SetBranchAddress("zdc_ZdcTime", zdc_ZdcTime, &b_zdc_ZdcTime);
//    fChain->SetBranchAddress("zdc_ZdcStatus", zdc_ZdcStatus, &b_zdc_ZdcStatus);
//    fChain->SetBranchAddress("zdc_ZdcTrigEff", zdc_ZdcTrigEff, &b_zdc_ZdcTrigEff);
//    fChain->SetBranchAddress("zdc_ZdcModuleMask", &zdc_ZdcModuleMask, &b_zdc_ZdcModuleMask);
//    fChain->SetBranchAddress("zdc_ZdcLucrodTriggerSideAmp", zdc_ZdcLucrodTriggerSideAmp, &b_zdc_ZdcLucrodTriggerSideAmp);
//    fChain->SetBranchAddress("zdc_ZdcModuleAmp", zdc_ZdcModuleAmp, &b_zdc_ZdcModuleAmp);
//    fChain->SetBranchAddress("zdc_ZdcModuleTime", zdc_ZdcModuleTime, &b_zdc_ZdcModuleTime);
//    fChain->SetBranchAddress("zdc_ZdcModuleFitAmp", zdc_ZdcModuleFitAmp, &b_zdc_ZdcModuleFitAmp);
//    fChain->SetBranchAddress("zdc_ZdcModuleFitT0", zdc_ZdcModuleFitT0, &b_zdc_ZdcModuleFitT0);
//    fChain->SetBranchAddress("zdc_ZdcModuleStatus", zdc_ZdcModuleStatus, &b_zdc_ZdcModuleStatus);
//    fChain->SetBranchAddress("zdc_ZdcModuleChisq", zdc_ZdcModuleChisq, &b_zdc_ZdcModuleChisq);
//    fChain->SetBranchAddress("zdc_ZdcModuleCalibAmp", zdc_ZdcModuleCalibAmp, &b_zdc_ZdcModuleCalibAmp);
//    fChain->SetBranchAddress("zdc_ZdcModuleCalibTime", zdc_ZdcModuleCalibTime, &b_zdc_ZdcModuleCalibTime);
//    fChain->SetBranchAddress("zdc_ZdcModuleBkgdMaxFraction", zdc_ZdcModuleBkgdMaxFraction, &b_zdc_ZdcModuleBkgdMaxFraction);
//    fChain->SetBranchAddress("zdc_ZdcModuleAmpError", zdc_ZdcModuleAmpError, &b_zdc_ZdcModuleAmpError);
//    fChain->SetBranchAddress("zdc_ZdcModuleMinDeriv2nd", zdc_ZdcModuleMinDeriv2nd, &b_zdc_ZdcModuleMinDeriv2nd);
//    fChain->SetBranchAddress("zdc_ZdcModulePresample", zdc_ZdcModulePresample, &b_zdc_ZdcModulePresample);
//    fChain->SetBranchAddress("zdc_ZdcModulePreSampleAmp", zdc_ZdcModulePreSampleAmp, &b_zdc_ZdcModulePreSampleAmp);
//    fChain->SetBranchAddress("zdc_ZdcLucrodTriggerAmp", zdc_ZdcLucrodTriggerAmp, &b_zdc_ZdcLucrodTriggerAmp);
//    fChain->SetBranchAddress("zdc_ZdcModuleMaxADC", zdc_ZdcModuleMaxADC, &b_zdc_ZdcModuleMaxADC);
//    fChain->SetBranchAddress("zdc_ZdcModuleTruthTotal", zdc_ZdcModuleTruthTotal, &b_zdc_ZdcModuleTruthTotal);
//    fChain->SetBranchAddress("zdc_ZdcModuleTruthInvisible", zdc_ZdcModuleTruthInvisible, &b_zdc_ZdcModuleTruthInvisible);
//    fChain->SetBranchAddress("zdc_ZdcModuleTruthEM", zdc_ZdcModuleTruthEM, &b_zdc_ZdcModuleTruthEM);
//    fChain->SetBranchAddress("zdc_ZdcModuleTruthNonEM", zdc_ZdcModuleTruthNonEM, &b_zdc_ZdcModuleTruthNonEM);
//    fChain->SetBranchAddress("zdc_ZdcModuleTruthEscaped", zdc_ZdcModuleTruthEscaped, &b_zdc_ZdcModuleTruthEscaped);
//    fChain->SetBranchAddress("zdc_ZdcModuleTruthNphotons", zdc_ZdcModuleTruthNphotons, &b_zdc_ZdcModuleTruthNphotons);
//    fChain->SetBranchAddress("zdc_RpdModuleTruthNphotons", zdc_RpdModuleTruthNphotons, &b_zdc_RpdModuleTruthNphotons);
//    fChain->SetBranchAddress("zdc_ZdcTruthTotal", zdc_ZdcTruthTotal, &b_zdc_ZdcTruthTotal);
//    fChain->SetBranchAddress("zdc_ZdcTruthInvisible", zdc_ZdcTruthInvisible, &b_zdc_ZdcTruthInvisible);
//    fChain->SetBranchAddress("zdc_ZdcTruthEM", zdc_ZdcTruthEM, &b_zdc_ZdcTruthEM);
//    fChain->SetBranchAddress("zdc_ZdcTruthNonEM", zdc_ZdcTruthNonEM, &b_zdc_ZdcTruthNonEM);
//    fChain->SetBranchAddress("zdc_ZdcTruthEscaped", zdc_ZdcTruthEscaped, &b_zdc_ZdcTruthEscaped);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticlePosx", &zdc_ZdcTruthParticlePosx, &b_zdc_ZdcTruthParticlePosx);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticlePosy", &zdc_ZdcTruthParticlePosy, &b_zdc_ZdcTruthParticlePosy);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticlePosz", &zdc_ZdcTruthParticlePosz, &b_zdc_ZdcTruthParticlePosz);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticleTime", &zdc_ZdcTruthParticleTime, &b_zdc_ZdcTruthParticleTime);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticlePx", &zdc_ZdcTruthParticlePx, &b_zdc_ZdcTruthParticlePx);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticlePy", &zdc_ZdcTruthParticlePy, &b_zdc_ZdcTruthParticlePy);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticlePz", &zdc_ZdcTruthParticlePz, &b_zdc_ZdcTruthParticlePz);
//    fChain->SetBranchAddress("zdc_ZdcTruthParticleEnergy", &zdc_ZdcTruthParticleEnergy, &b_zdc_ZdcTruthParticleEnergy);
//    fChain->SetBranchAddress("zdc_RpdChannelBaseline", zdc_RpdChannelBaseline, &b_zdc_RpdChannelBaseline);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupExpFitParams", zdc_RpdChannelPileupExpFitParams, &b_zdc_RpdChannelPileupExpFitParams);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupStretchedExpFitParams", zdc_RpdChannelPileupStretchedExpFitParams, &b_zdc_RpdChannelPileupStretchedExpFitParams);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupExpFitParamErrs", zdc_RpdChannelPileupExpFitParamErrs, &b_zdc_RpdChannelPileupExpFitParamErrs);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupStretchedExpFitParamErrs", zdc_RpdChannelPileupStretchedExpFitParamErrs, &b_zdc_RpdChannelPileupStretchedExpFitParamErrs);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupExpFitMSE", zdc_RpdChannelPileupExpFitMSE, &b_zdc_RpdChannelPileupExpFitMSE);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupStretchedExpFitMSE", zdc_RpdChannelPileupStretchedExpFitMSE, &b_zdc_RpdChannelPileupStretchedExpFitMSE);
//    fChain->SetBranchAddress("zdc_RpdChannelAmplitude", zdc_RpdChannelAmplitude, &b_zdc_RpdChannelAmplitude);
//    fChain->SetBranchAddress("zdc_RpdChannelAmplitudeCalib", zdc_RpdChannelAmplitudeCalib, &b_zdc_RpdChannelAmplitudeCalib);
//    fChain->SetBranchAddress("zdc_RpdChannelMaxADC", zdc_RpdChannelMaxADC, &b_zdc_RpdChannelMaxADC);
//    fChain->SetBranchAddress("zdc_RpdChannelMaxADCCalib", zdc_RpdChannelMaxADCCalib, &b_zdc_RpdChannelMaxADCCalib);
//    fChain->SetBranchAddress("zdc_RpdChannelMaxSample", zdc_RpdChannelMaxSample, &b_zdc_RpdChannelMaxSample);
//    fChain->SetBranchAddress("zdc_RpdChannelStatus", zdc_RpdChannelStatus, &b_zdc_RpdChannelStatus);
//    fChain->SetBranchAddress("zdc_RpdChannelPileupFrac", zdc_RpdChannelPileupFrac, &b_zdc_RpdChannelPileupFrac);
//    fChain->SetBranchAddress("zdc_RpdSideStatus", zdc_RpdSideStatus, &b_zdc_RpdSideStatus);
//    fChain->SetBranchAddress("zdc_centroidEventValid", &zdc_centroidEventValid, &b_zdc_centroidEventValid);
//    fChain->SetBranchAddress("zdc_centroidStatus", zdc_centroidStatus, &b_zdc_centroidStatus);
//    fChain->SetBranchAddress("zdc_RPDChannelSubtrAmp", zdc_RPDChannelSubtrAmp, &b_zdc_RPDChannelSubtrAmp);
//    fChain->SetBranchAddress("zdc_RPDSubtrAmpSum", zdc_RPDSubtrAmpSum, &b_zdc_RPDSubtrAmpSum);
//    fChain->SetBranchAddress("zdc_xCentroidPreGeomCorPreAvgSubtr", zdc_xCentroidPreGeomCorPreAvgSubtr, &b_zdc_xCentroidPreGeomCorPreAvgSubtr);
//    fChain->SetBranchAddress("zdc_yCentroidPreGeomCorPreAvgSubtr", zdc_yCentroidPreGeomCorPreAvgSubtr, &b_zdc_yCentroidPreGeomCorPreAvgSubtr);
//    fChain->SetBranchAddress("zdc_xCentroidPreAvgSubtr", zdc_xCentroidPreAvgSubtr, &b_zdc_xCentroidPreAvgSubtr);
//    fChain->SetBranchAddress("zdc_yCentroidPreAvgSubtr", zdc_yCentroidPreAvgSubtr, &b_zdc_yCentroidPreAvgSubtr);
//    fChain->SetBranchAddress("zdc_xCentroid", zdc_xCentroid, &b_zdc_xCentroid);
//    fChain->SetBranchAddress("zdc_yCentroid", zdc_yCentroid, &b_zdc_yCentroid);
//    fChain->SetBranchAddress("zdc_xRowCentroid", zdc_xRowCentroid, &b_zdc_xRowCentroid);
//    fChain->SetBranchAddress("zdc_yColCentroid", zdc_yColCentroid, &b_zdc_yColCentroid);
//    fChain->SetBranchAddress("zdc_reactionPlaneAngle", zdc_reactionPlaneAngle, &b_zdc_reactionPlaneAngle);
//    fChain->SetBranchAddress("zdc_cosDeltaReactionPlaneAngle", &zdc_cosDeltaReactionPlaneAngle, &b_zdc_cosDeltaReactionPlaneAngle);
//    Notify();
// }

// Bool_t zdcTree::Notify()
// {
//    // The Notify() function is called when a new file is opened. This
//    // can be either for a new TTree in a TChain or when when a new TTree
//    // is started when using PROOF. It is normally not necessary to make changes
//    // to the generated code, but the routine can be extended by the
//    // user if needed. The return value is currently not used.

//    return kTRUE;
// }

// void zdcTree::Show(Long64_t entry)
// {
// // Print contents of entry.
// // If entry is not specified, print current entry
//    if (!fChain) return;
//    fChain->Show(entry);
// }

// Int_t zdcTree::Cut(Long64_t entry)
// {
// // This function may be called from Loop.
// // returns  1 if entry is accepted.
// // returns -1 otherwise.
//    return 1;
// }
#endif // #ifdef zdcTree_cxx
