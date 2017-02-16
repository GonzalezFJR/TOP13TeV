#ifndef AnalMiniTree_h
#define AnalMiniTree_h

#include <TROOT.h>
#include <TChain.h>
#include <THStack.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <TLorentzVector.h>
#include "Histo.h"



//############################################################################
//# Definition of constants
//############################################################################
//const TString path = "/mnt_pool/fanae105/user/juanr/stop/TOP13TeV/outputFiles/_feb11/";
const TString path = "/mnt_pool/fanae105/user/juanr/stop/TOP13TeV/outputFiles/feb13/dilepton/";
//const TString path = "/mnt_pool/fanae105/user/juanr/stop/TOP13TeV/outputFiles/feb13/";
const TString prefix = "Tree_"; const TString sufix = ".root";
const TString treeName = "sTopTree";

const Int_t nBkgs = 13;
const TString Bkgs[nBkgs] = {
"DYJetsToLL_M10to50", "DYJetsToLL_M50_aMCatNLO",  // DY
"TTbar_Powheg", "TbarW", "TW",                    // ttbar, tW
"TTWToLNu", "TTWToQQ", "TTZToLLNuNu", "TTZToQQ",  // TTX
"WW", "WZ", "ZZ",                                 // Dibosons
"WJetsToLNu_aMCatNLO"                             // WJets
};
const Int_t nDataSamples = 5;
const TString DataSamples[nDataSamples] = {"DoubleMuon", "DoubleEG", "MuonEF", "SigleElectron", "SingleMuon"};
const Int_t nSignals = 10;
const TString Signals[nSignals] = {"T2tt_mStop175_mLsp1", "T2tt_mStop183_mLsp1", "T2tt_mStop192_mLsp25", "T2tt_mStop200_mLsp25", "T2tt_mStop208_mLsp25", "T2tt_mStop217_mLsp50", "T2tt_mStop225_mLsp50", "T2tt_mStop242_mLsp75", "T2tt_mStop250_mLsp75", "T2tt_mStop258_mLsp75"};

//############################################################################
//# Global variables
//############################################################################

class AnalMiniTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *tree;
   Histo* histo;
   TString theSample;
   TString chan;
   TString sys = "0";

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
	 Float_t         TWeight;
	 Float_t TWeight_LepEffUp;
	 Float_t  TWeight_LepEffDown;
	 Float_t  TWeight_TrigUp;
	 Float_t  TWeight_TrigDown;
	 Float_t  TWeight_FSUp;
	 Float_t  TWeight_FSDown;
	 Float_t  TWeight_PUUp;
	 Float_t  TWeight_PUDown;
   Int_t           TChannel;
   Int_t           TNJets;
   Int_t           TNVert;
   Int_t           TNJetsJESUp;
   Int_t           TNJetsJESDown;
   Int_t           TNJetsJER;
   Int_t           TNJetsBtag;
   Int_t           TNJetsBtagUp;
   Int_t           TNJetsBtagDown;
   Int_t           TNJetsBtagMisTagUp;
   Int_t           TNJetsBtagMisTagDown;
   Float_t         TMET;
   Float_t         TGenMET;
   Float_t         TMETJESUp;
   Float_t         TMETJESDown;
   Float_t         TMT2ll;
   Float_t         TMT2llJESUp;
   Float_t         TMT2llJESDown;
   Float_t         TMT2jj;
   Float_t         TMT2lblb;
   Float_t         TMll;
   Float_t         TMeff;
   Float_t         TPtllb;
   Float_t         THT;
   Float_t         THTJESUp;
   Float_t         THTJESDown;
   Float_t         TMET_Phi;
   Float_t         TdPhiPtllbMET;
   Float_t         TMinDPhiMetJets;
   Float_t         TdPhiJetMet;
   Float_t         TdPhiLepMet;
   Float_t         TdPhiLepJet;
   Float_t         TdPhill;
   Float_t         TLep1_Px;
   Float_t         TLep1_Py;
   Float_t         TLep1_Pz;
   Float_t         TLep1_E;
   Float_t         TLep1_Carge;
   Float_t         TLep2_Px;
   Float_t         TLep2_Py;
   Float_t         TLep2_Pz;
   Float_t         TLep2_E;
   Float_t         TLep2_Carge;
   Int_t           TJet_isBJet[13];   //[TNJets]
   Float_t         TJet_Px[13];   //[TNJets]
   Float_t         TJet_Py[13];   //[TNJets]
   Float_t         TJet_Pz[13];   //[TNJets]
   Float_t         TJet_E[13];   //[TNJets]
   Float_t         TJetJESUp_Pt[13];   //[TNJets]
   Float_t         TJetJESDown_Pt[13];   //[TNJets]

	 // List of branches
	 TBranch        *b_TWeight;   //!
	 TBranch        *b_TChannel;   //!
	 TBranch		*b_TWeight_LepEffUp;
	 TBranch	*b_TWeight_LepEffDown;
	 TBranch		*b_TWeight_TrigUp;
	 TBranch		*b_TWeight_TrigDown;
	 TBranch		*b_TWeight_FSUp;
	 TBranch		*b_TWeight_FSDown;
	 TBranch		*b_TWeight_PUUp;
	 TBranch		*b_TWeight_PUDown;

   TBranch        *b_TIsDoubleMuon;   //!
   TBranch        *b_TIsDoubleElec;   //!
   TBranch        *b_TIsElMu;   //!
   TBranch        *b_TNJets;   //!
   TBranch        *b_TNVert;   //!
	 TBranch *b_TNJetsJESUp;
	 TBranch *b_TNJetsJESDown;
	 TBranch *b_TNJetsJER;
	 TBranch *b_TNJetsBtag;
	 TBranch *b_TNJetsBtagUp;
	 TBranch *b_TNJetsBtagDown;
	 TBranch *b_TNJetsBtagMisTagUp;
	 TBranch *b_TNJetsBtagMisTagDown;
   TBranch        *b_TMET;   //!
   TBranch        *b_TGenMET;   //!
   TBranch        *b_TMETJESUp;   //!
   TBranch        *b_TMETJESDown;   //!
   TBranch        *b_TMT2ll;   //!
   TBranch        *b_TMT2llJESUp;   //!
   TBranch        *b_TMT2llJESDown;   //!
   TBranch        *b_TMT2jj;   //!
   TBranch        *b_TMT2lblb;   //!
   TBranch        *b_TMll;   //!
   TBranch        *b_TMeff;   //!
   TBranch        *b_TPtllb;   //!
   TBranch        *b_THT;   //!
   TBranch        *b_THTJESUp;   //!
   TBranch        *b_THTJESDown;   //!
   TBranch        *b_TMET_Phi;   //!
   TBranch        *b_TdPhiPtllbMET;   //!
   TBranch        *b_TMinDPhiMetJets;   //!
   TBranch        *b_TdPhiJetMet;   //!
   TBranch        *b_TdPhiLepMet;   //!
   TBranch        *b_TdPhiLepJet;   //!
   TBranch        *b_TdPhill;   //!
   TBranch        *b_TLep1_Px;   //!
   TBranch        *b_TLep1_Py;   //!
   TBranch        *b_TLep1_Pz;   //!
   TBranch        *b_TLep1_E;   //!
   TBranch        *b_TLep1_Charge;   //!
   TBranch        *b_TLep2_Px;   //!
   TBranch        *b_TLep2_Py;   //!
   TBranch        *b_TLep2_Pz;   //!
   TBranch        *b_TLep2_E;   //!
   TBranch        *b_TLep2_Charge;   //!
   TBranch        *b_TJet_isBJet;   //!
   TBranch        *b_TJet_Px;   //!
   TBranch        *b_TJet_Py;   //!
   TBranch        *b_TJet_Pz;   //!
   TBranch        *b_TJet_E;   //!
   TBranch        *b_TJetJESUp_Pt;   //!
   TBranch        *b_TJetJESDown_Pt;   //!

	 AnalMiniTree(TString sample = "TTbar_Powheg", TString ch = "ElMu", TString thesys = "0"){
		 chan = ch;
		 theSample = sample;
     sys = thesys;
		 TString Path = path;
		 if(sample.BeginsWith("T2tt")) Path = path + "Susy/";
		 TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Path + prefix + sample + sufix);
		 if (!f || !f->IsOpen()) {
			 f = new TFile(Path + prefix + sample + sufix);
		 }
		 f->GetObject(treeName,tree);

		 Init(tree);
	 }
	 virtual ~AnalMiniTree();
	 virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   void     Loop(TString plot);
};


AnalMiniTree::~AnalMiniTree(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnalMiniTree::GetEntry(Long64_t entry){
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnalMiniTree::LoadTree(Long64_t entry){ 
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void AnalMiniTree::Init(TTree *tree){
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("TWeight", &TWeight, &b_TWeight);
   fChain->SetBranchAddress("TChannel", &TChannel, &b_TChannel);
	 fChain->SetBranchAddress("TWeight_LepEffUp", &TWeight_LepEffUp, &b_TWeight_LepEffUp);
	 fChain->SetBranchAddress("TWeight_LepEffDown", &TWeight_LepEffDown, &b_TWeight_LepEffDown);
	 fChain->SetBranchAddress("TWeight_FSUp", &TWeight_FSUp, &b_TWeight_FSUp);
	 fChain->SetBranchAddress("TWeight_FSDown", &TWeight_FSDown, &b_TWeight_FSDown);
	 fChain->SetBranchAddress("TWeight_TrigUp", &TWeight_TrigUp, &b_TWeight_TrigUp);
	 fChain->SetBranchAddress("TWeight_TrigDown", &TWeight_TrigDown, &b_TWeight_TrigDown);
	 fChain->SetBranchAddress("TWeight_PUUp", &TWeight_PUUp, &b_TWeight_PUUp);
	 fChain->SetBranchAddress("TWeight_PUDown", &TWeight_PUDown, &b_TWeight_PUDown);
   fChain->SetBranchAddress("TNVert", &TNVert, &b_TNVert);
   fChain->SetBranchAddress("TNJets", &TNJets, &b_TNJets);
   fChain->SetBranchAddress("TNJetsJESUp", &TNJetsJESUp, &b_TNJetsJESUp);
   fChain->SetBranchAddress("TNJetsJESDown", &TNJetsJESDown, &b_TNJetsJESDown);
   fChain->SetBranchAddress("TNJetsJER", &TNJetsJER, &b_TNJetsJER);
   fChain->SetBranchAddress("TNJetsBtag", &TNJetsBtag, &b_TNJetsBtag);
   fChain->SetBranchAddress("TNJetsBtagUp", &TNJetsBtagUp, &b_TNJetsBtagUp);
   fChain->SetBranchAddress("TNJetsBtagDown", &TNJetsBtagDown, &b_TNJetsBtagDown);
   fChain->SetBranchAddress("TNJetsBtagMisTagUp", &TNJetsBtagMisTagUp, &b_TNJetsBtagMisTagUp);
   fChain->SetBranchAddress("TNJetsBtagMisTagDown", &TNJetsBtagMisTagDown, &b_TNJetsBtagMisTagDown);
   fChain->SetBranchAddress("TMET", &TMET, &b_TMET);
   fChain->SetBranchAddress("TGenMET", &TGenMET, &b_TGenMET);
   fChain->SetBranchAddress("TMETJESUp", &TMETJESUp, &b_TMETJESUp);
   fChain->SetBranchAddress("TMETJESDown", &TMETJESDown, &b_TMETJESDown);
   fChain->SetBranchAddress("TMT2ll", &TMT2ll, &b_TMT2ll);
   fChain->SetBranchAddress("TMT2llJESUp", &TMT2llJESUp, &b_TMT2llJESUp);
   fChain->SetBranchAddress("TMT2llJESDown", &TMT2llJESDown, &b_TMT2llJESDown);
   fChain->SetBranchAddress("TMT2jj", &TMT2jj, &b_TMT2jj);
   fChain->SetBranchAddress("TMT2lblb", &TMT2lblb, &b_TMT2lblb);
   fChain->SetBranchAddress("TMll", &TMll, &b_TMll);
   fChain->SetBranchAddress("TMeff", &TMeff, &b_TMeff);
   fChain->SetBranchAddress("TPtllb", &TPtllb, &b_TPtllb);
   fChain->SetBranchAddress("THT", &THT, &b_THT);
   fChain->SetBranchAddress("THTJESUp", &THTJESUp, &b_THTJESUp);
   fChain->SetBranchAddress("THTJESDown", &THTJESDown, &b_THTJESDown);
   fChain->SetBranchAddress("TMET_Phi", &TMET_Phi, &b_TMET_Phi);
   fChain->SetBranchAddress("TdPhiPtllbMET", &TdPhiPtllbMET, &b_TdPhiPtllbMET);
   fChain->SetBranchAddress("TMinDPhiMetJets", &TMinDPhiMetJets, &b_TMinDPhiMetJets);
   fChain->SetBranchAddress("TdPhiJetMet", &TdPhiJetMet, &b_TdPhiJetMet);
   fChain->SetBranchAddress("TdPhiLepMet", &TdPhiLepMet, &b_TdPhiLepMet);
   fChain->SetBranchAddress("TdPhiLepJet", &TdPhiLepJet, &b_TdPhiLepJet);
   fChain->SetBranchAddress("TdPhill", &TdPhill, &b_TdPhill);
   fChain->SetBranchAddress("TLep1_Px", &TLep1_Px, &b_TLep1_Px);
   fChain->SetBranchAddress("TLep1_Py", &TLep1_Py, &b_TLep1_Py);
   fChain->SetBranchAddress("TLep1_Pz", &TLep1_Pz, &b_TLep1_Pz);
   fChain->SetBranchAddress("TLep1_E", &TLep1_E, &b_TLep1_E);
   fChain->SetBranchAddress("TLep1_Carge", &TLep1_Carge, &b_TLep1_Charge);
   fChain->SetBranchAddress("TLep2_Px", &TLep2_Px, &b_TLep2_Px);
   fChain->SetBranchAddress("TLep2_Py", &TLep2_Py, &b_TLep2_Py);
   fChain->SetBranchAddress("TLep2_Pz", &TLep2_Pz, &b_TLep2_Pz);
   fChain->SetBranchAddress("TLep2_E", &TLep2_E, &b_TLep2_E);
   fChain->SetBranchAddress("TLep2_Carge", &TLep2_Carge, &b_TLep2_Charge);
   fChain->SetBranchAddress("TJet_isBJet", TJet_isBJet, &b_TJet_isBJet);
   fChain->SetBranchAddress("TJet_Px", TJet_Px, &b_TJet_Px);
   fChain->SetBranchAddress("TJet_Py", TJet_Py, &b_TJet_Py);
   fChain->SetBranchAddress("TJet_Pz", TJet_Pz, &b_TJet_Pz);
   fChain->SetBranchAddress("TJet_E", TJet_E, &b_TJet_E);
   fChain->SetBranchAddress("TJetJESUp_Pt",  TJetJESUp_Pt, &b_TJetJESUp_Pt);
   fChain->SetBranchAddress("TJetJESDown_Pt",  TJetJESDown_Pt, &b_TJetJESDown_Pt);
}
#endif
