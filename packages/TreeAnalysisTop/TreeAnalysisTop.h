#pragma once
// Includes
// + PAF
#include "PAFChainItemSelector.h"

// + Packages
//#include "GlobalVariables.h"
#include "mt2.h"
#include "PUWeight.h"
#include "BTagSFUtil.h"
#include "LeptonSF.h"

// + ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TRandom3.h"

// + C++
#include <vector>

const int nWeights = 248;
const int nGenb = 0;
const double pi = 3.1415926535897932384;

enum gChannel{
  channels_begin,
  Muon = channels_begin,
  Elec,
  ElMu,
  gNCHANNELS,
};
const TString gChanLabel[gNCHANNELS] = {"Muon","Elec","ElMu"};
enum gFPSwitch{
  SigSup,
  ZDecay,
  Sig
};
enum iCut{
  iDilepton, 
  iZVeto, 
  iMET, 
  i2jets, 
  i1btag, 
  iExact1btag,
  iExact2btag,
  iNCUTS
};
const TString sCut[iNCUTS] = {"dilepton", "ZVeto", "MET", "2jets", "1btag","Exact1btag","Exact2btag"};
enum gSystFlag{
  Norm,
  BtagUp,
  BtagDown,
  MisTagUp,
  MisTagDown,
  JESUp,
  JESDown,
  JER,
  LESUp,
  LESDown,
  /*  LepUp,
      LepDown,
      TrigUp,
      TrigDown,
  */
  PUUp,
  PUDown,
  TopPt,
  gNSYST
};
const TString SystName[gNSYST] = {
  "Normal",
  "BtagUp",
  "BtagDown",
  "MisTagUp",
  "MisTagDown",
  "JESUp",
  "JESDown",
  "JER",
  "LESUp",
  "LESDown",
  /*  "LepUp",
  "LepDown",
  "TrigUp",
  "TrigDown",*/
//  "METUp",
//  "METDown",
  "PUUp",
  "PUDown",
  "TopPt",
};
enum FakeSource{
  HF_mu,
  Other_mu,
  HF_el,
  Conv_el,
  Other_el,
  RightSign,
  WrongSign,
  gNFAKESOURCE
};

enum gNLOWeight{
  muR1muF1,
  muR1muF2,
  muR1muF05,
  muR2muF1,
  muR2muF2,
  muR2muF05,
  muR05muF1,
  muR05muF2,
  muR05muF05,
  gNWEIGHT
};

const TString WeiName[gNWEIGHT] = {
  "muR1muF1",
  "muR1muF2",
  "muR1muF05",
  "muR2muF1",
  "muR2muF2",
  "muR2muF05",
  "muR05muF1",
  "muR05muF2",
  "muR05muF05"
};

class lepton{
 public:
  //lepton(){}
  //lepton(const lepton &l): p(l.p), charge(l.charge), type(l.type), index(l.index){ };
  lepton(TLorentzVector vec = TLorentzVector(0,0,0,0), int ch = 0, int ty = -1, int ind = -1){
    p = vec;
    charge = ch;
    type = ty;
    index = ind;
  }
  TLorentzVector p;
  int charge;
  int type; // -1(unknown), 0(mu), 1(ele)
  int index;
};

class jet{
 public:
  jet(){};
  jet(TLorentzVector vec, bool btag, int ind){
    p = vec;
    isbtag = btag;
    index = ind;
  };
  TLorentzVector p;
  bool isbtag;
  int index;
};
 
 
const int gNMuFPtBins = 6;
const int gNMuPPtbins = 10;
const int gNMuEtabins = 5;
const int gNElFPtBins = 8;
const int gNElPPtbins = 10;
const int gNElEtabins = 5;
const int gNElCMIdbins = 2;
 
// Muon Binning
const double gMuFPtBins[gNMuFPtBins+1] = {20., 25., 30., 35., 40., 50., 60.}; 
const double gMuPPtbins[gNMuPPtbins+1] = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.};
const double gMuEtabins[gNMuEtabins+1] = {0., 0.5, 1.0, 1.479, 2.0, 2.5};
 
// Electron Binning //////////////////////////////////////////////////////////////
const double gElFPtBins[gNElFPtBins+1]   = {20., 25., 30., 40., 50., 60., 70., 80., 100.};
const double gElPPtbins[gNElPPtbins+1]   = {20., 25., 30., 35., 40., 50., 60., 70., 80., 90., 100.}; 
const double gElEtabins[gNElEtabins+1]   = {0., 0.5, 1.0, 1.479, 2.0, 2.5};

 
int     getNFPtBins(gChannel chan);
const double *getFPtBins (gChannel chan);
int     getNPPtBins(gChannel chan);
const double *getPPtBins (gChannel chan);
int     getNEtaBins(gChannel chan);
const double *getEtaBins (gChannel chan);

// Your selector class
class TreeAnalysisTop : public PAFChainItemSelector 
{
 public:
  // Constructor and destructor
  TreeAnalysisTop();
  virtual ~TreeAnalysisTop() {}
  
  // Mandatory methods
  virtual void Initialise();
  virtual void InsideLoop();
    virtual void Summary();

  // Frequantly used variables in the trees
  Int_t   nLepGood;
  Int_t   ngenLep;
  Int_t   ngenLepFromTau;
  Int_t   nJet;
  Float_t genWeight;
  Float_t LepGood_px[30];
  Float_t LepGood_py[30];
  Float_t LepGood_pz[30];
  Float_t LepGood_energy[30];
  Float_t LepGood_pt[30];
  Float_t LepGood_etaSc[30];
  Float_t LepGood_eta[30];
  Float_t LepGood_dxy[30];
  Float_t LepGood_dz[30];
  Float_t LepGood_relIso03[30];
  Float_t LepGood_relIso04[30];
  Int_t   LepGood_charge[30];
  Int_t   LepGood_pdgId[30];
  Float_t Jet_px[50];
  Float_t Jet_py[50];
  Float_t Jet_pz[50];
  Float_t Jet_energy[50];
  Float_t Jet_eta[50];
  Float_t Jet_btagCSV[50];
  Float_t genLep_pt[50];
  Float_t genLep_eta[50];
  Float_t genLep_phi[50];
  Float_t genLep_mass[50];
  Int_t   genLep_pdgId[50];
  Int_t   genLepFromTau_pdgId[50];

  // You would tipically add here your protected methods
  // and data members
  virtual void InitialiseYieldsHistos();
  virtual void InitialiseKinematicHistos();
  virtual void InitialiseDYHistos();
  virtual void InitialiseGenHistos();
  virtual void InitialiseSystematicHistos();

  // Saving Histograms.
  void WriteHistos();
  void WriteValidationsHistos(){};

  void  GetParameters();
  void  GetTreeVariables();
//   Int_t SelectedVertexIndex();
  
  bool PassTriggerMuMu();
  bool PassTriggerEE();
  bool PassTriggerEMu();
  bool PassesZVeto();
  bool PassesNJetsCut();
  bool PassesMETCut();
  bool PassesNBtagCut();
  bool PassesMllVeto();
  bool Passes3rdLeptonVeto();
  bool PassesMuonEta2p1(gChannel);
  bool PassesTopDCut();

  int   getNJets();
  int   getNBTags();
  int   getLeadingJetbTag();
  float getDRClosestJet(TLorentzVector);
  float getDPhiClosestJet(TLorentzVector);
  //  int   getNBTagsMed();
  void  setMET(float);
  float getMET();
  float getMETPhi();
  //float getMT(int, gChannel);
  float getHT();
  float getJetPtIndex(unsigned int);
  float getJetEtaIndex(unsigned int);
  float getBtagJetPtIndex(unsigned int);
  float getErrPt(float,float);
  float getJERScaleUp(int);
  float getJERScale(int);
  float getJERScaleDown(int);
  float getSF(gChannel);
  float getLeptonError(gChannel);
  float getTriggerError(gChannel);
  float getTopPtSF();
  float getTopD();
  float getDeltaPhillJet();
  float weightNvtx(int);
  float getMT(gChannel);
  float getMT2(gChannel);
  float getMeff();
  float getDPhiLepJet();
  float getDelPhill();
  float getDPhiJetMet();
  float getDPhiLepMet();
  float getDPhibMet();

  TLorentzVector getPtllb();

  bool METFilter();
  // Lepton selection methods
  int  getSelectedLeptons();
//   bool IsVetoMuon(unsigned int, float ptcut=20.);
  bool IsTightMuon(unsigned int, float ptcut=20.);
  float getMuonIso(int);
//   bool IsVetoElectron(unsigned int,float ptcut=20.);
//   bool IsMVAIDElectron(unsigned int);
  bool IsTightElectron(unsigned int,float ptcut=20.);
  float getElecIso(int);
  float getEACorrection(float);
  std::vector<lepton> SortLeptonsByPt(std::vector<lepton>&);
  
  int getSelectedJets();
  bool IsGoodJet(unsigned int, float ptcut=30.);
  std::vector<int> CleanedJetIndices(float);
  void SmearJetPts(int);
  void propagateMET(TLorentzVector,TLorentzVector);
  void ScaleMET(int);
  void ScaleLeptons(int);
  
  //void GetGenMuon();
  //void GetGenElec();
  void SelectedGenLepton();
  
  ///////////////////////////////////////////////////////////////////////////// 
  //    Selecting Methods
  ///////////////////////////////////////////////////////////////////////////// 
  int  IsDileptonEvent();
  bool IsMuMuEvent();
  bool IsElMuEvent();
  bool IsElElEvent();
  ///////////////////////////////////////////////////////////////////////////// 
  //    Filling Methods
  ///////////////////////////////////////////////////////////////////////////// 
  void FillYieldsHistograms(gChannel, iCut, gSystFlag);
  void FillYields(gSystFlag sys=Norm);
  void FillDYHistograms();
  void FillKinematicHistos(gChannel,iCut);
  
  ///////////////////////////////////////////////////////////////////////////// 
  //    Set/Reset methods
  ///////////////////////////////////////////////////////////////////////////// 
  void SetOriginalObjects();
  void ResetOriginalObjects();
  void SetEventObjects();
  void ResetHypLeptons();

 protected:

  // You would tipically add here your protected methods
  // and data members

  // Input parameters
  //----------------------------------------------------------------------------
  TString gSampleName;
  TString gfileSuffix;
  Float_t gWeight;
  Float_t gLumiForPU;
  Float_t gTotalLumi;
  Int_t   gSysSource;
  Int_t   gSysDirection;
  Bool_t  gDoSystStudies;
  Bool_t  gIsData;
  Bool_t  gUseCSVM;
  Bool_t  gDoSF;
  Bool_t  gDoDF;
  Int_t   gSelection;
  Bool_t  gIsMCatNLO;

  PUWeight *fPUWeight;      //The PU weight utility
  PUWeight *fPUWeightUp;    //The PU weight utility
  PUWeight *fPUWeightDown;  //The PU weight utility
  //BTagSFUtil *fBTagSF[5]; //The new BTag SF utility 
  BTagSFUtil *fBTagSFnom ;
  BTagSFUtil *fBTagSFbUp ;
  BTagSFUtil *fBTagSFbDo ;
  BTagSFUtil *fBTagSFlUp ;
  BTagSFUtil *fBTagSFlDo ;
  LeptonSF *fLeptonSF;
  TRandom3 *fRand3;

  // EventWeight
  //----------------------------------------------------------------------------
  float EventWeight;
  float PUSF;
  bool  fChargeSwitch;


  //////////////////////////////////////////////////////////////////////////////
  //               Data members
  //////////////////////////////////////////////////////////////////////////////
  // HISTOGRAMS
  
  //++ Yields
  TH1F* fHDummy;
  TH1F* hWeight;
  TH1F* fHyields     [gNCHANNELS][gNSYST];
  TH1F* fHWeightyield[gNCHANNELS][gNWEIGHT];
  TH1F* fHSSyields   [gNCHANNELS][gNSYST];
  TH1F* fHTopPtWeight;
  TH1F* fHLepSys[gNCHANNELS][iNCUTS];
  TH1F* fHTrigSys[gNCHANNELS][iNCUTS];
  TH1F* fHnGenEle;
  TH1F* fHnGenMuo;
  TH1F* fHnGenLep0;
  TH1F* fHGenElePt;
  TH1F* fHGenMuoPt;

  TH2F* fHDY_InvMassVsNPV   [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsMET   [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsNjets [gNCHANNELS][iNCUTS];
  TH2F* fHDY_InvMassVsNbtags[gNCHANNELS][iNCUTS];
  TH1F* fHDY_InvMass        [gNCHANNELS][iNCUTS];
  
  //++ Origin Histos
//  TH2F* fHSSOrigins[gNCHANNELS][iNCUTS];
//  TH2F* fHOrigins[gNCHANNELS][iNCUTS];
  
  //++ Kinematic  
  TH1F* fHLHEweights[gNCHANNELS][iNCUTS];
  TH1F* fHMET[gNCHANNELS][iNCUTS];       
  TH1F* fHLep0Eta[gNCHANNELS][iNCUTS];    
  TH1F* fHLep1Eta[gNCHANNELS][iNCUTS];    
  TH1F* fHDelLepPhi[gNCHANNELS][iNCUTS]; 
  TH1F* fHHT[gNCHANNELS][iNCUTS];        
  TH1F* fHHT2[gNCHANNELS][iNCUTS];        
  TH1F* fHHT3[gNCHANNELS][iNCUTS];        
  TH1F* fHHT4[gNCHANNELS][iNCUTS];        
  TH1F* fHHT5[gNCHANNELS][iNCUTS];        
  TH1F* fHJet0Eta[gNCHANNELS][iNCUTS];    
  TH1F* fHJet1Eta[gNCHANNELS][iNCUTS];    
  TH1F* fHBtagJet0Pt[gNCHANNELS][iNCUTS];

  TH1F* fHMT[gNCHANNELS][iNCUTS];       
  TH1F* fHMT2[gNCHANNELS][iNCUTS];       
  TH1F* fHPtllb[gNCHANNELS][iNCUTS];       
  TH1F* fHMeff[gNCHANNELS][iNCUTS];       
  TH1F* fHDelPhiLepMet[gNCHANNELS][iNCUTS]; 
  TH1F* fHDelPhiJetMet[gNCHANNELS][iNCUTS]; 
  TH1F* fHDelPhiPllbMet[gNCHANNELS][iNCUTS]; 
  TH1F* fHDelPhiLepJet[gNCHANNELS][iNCUTS]; 

  TH1F* fHDiLepPt[gNCHANNELS][iNCUTS][gNSYST];   
  TH1F* fHLep0Pt[gNCHANNELS][iNCUTS][gNSYST];    
  TH1F* fHLep1Pt[gNCHANNELS][iNCUTS][gNSYST];    
  TH1F* fHJet0Pt[gNCHANNELS][iNCUTS][gNSYST];    
  TH1F* fHJet1Pt[gNCHANNELS][iNCUTS][gNSYST];    
  TH1F* fHNJets[gNCHANNELS][iNCUTS][gNSYST];     
  TH1F* fHNBtagJets[gNCHANNELS][iNCUTS][gNSYST]; 

  TH1F* fHInvMass[gNCHANNELS][iNCUTS][gNSYST];   
  TH1F* fHInvMass2[gNCHANNELS][iNCUTS][gNSYST];   
  TH1F* fHSSInvMass[gNCHANNELS][iNCUTS][gNSYST];   
  TH1F* fHNBtagsNJets[gNCHANNELS][iNCUTS][gNSYST]; 
  TH1F* fHSSNBtagsNJets[gNCHANNELS][iNCUTS][gNSYST]; 
  TH1F* fHCSVTag[gNCHANNELS][iNCUTS]; 
  TH1F* fHJetCSV[gNCHANNELS][iNCUTS];
  TH1F* fHTopD[gNCHANNELS][iNCUTS];
  TH1F* fHDelPhillJet[gNCHANNELS][iNCUTS];

  TH1F* fHDRLep[gNCHANNELS][iNCUTS];
  TH1F* fHDRLep0Jet[gNCHANNELS][iNCUTS];
  TH1F* fHDPhiLep0Jet[gNCHANNELS][iNCUTS];
  TH1F* fHLep0Iso[gNCHANNELS][iNCUTS];
  TH1F* fHDRLep1Jet[gNCHANNELS][iNCUTS];
  TH1F* fHDPhiLep1Jet[gNCHANNELS][iNCUTS];
  TH1F* fHLep1Iso[gNCHANNELS][iNCUTS];
  
  /// STOP
  //TH1F* fHAbsDelPhiLep[gNCHANNELS][iNCUTS];
  TH1F* fHminDelRJetsLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHSSminDelRJetsLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHdelPhi2LeadJets[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHSSdelPhi2LeadJets[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHAbsDelPhiLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHSSAbsDelPhiLeps[gNCHANNELS][iNCUTS][gNSYST];
  TH1F* fHStopMass[gNCHANNELS][iNCUTS];
  TH1F* fHChi0Mass[gNCHANNELS][iNCUTS];
  TH2F* fHChi0StopMass[gNCHANNELS][iNCUTS];
  TH1F* fHvertices[gNCHANNELS][iNCUTS];
  TH1F* fHgoodvertices[gNCHANNELS][iNCUTS];
	TH1F* fnGenLep[gNCHANNELS][iNCUTS];
  
  //++ Gen Info
  TH1F* fHDeltaRLepJet[gNCHANNELS-1];

  lepton fHypLepton1;
  lepton fHypLepton2;
  
  std::vector<Double_t>       Gen_Muon_Charge;
  std::vector<Double_t>       Gen_Elec_Charge;
  std::vector<TLorentzVector> Gen_Muon;
  std::vector<TLorentzVector> Gen_Elec;
 
  std::vector<Int_t>          NGen_Jet;
  std::vector<Int_t>          NGen_b;
  
  std::vector<Double_t>       PtGen_Jet;
  std::vector<Double_t>       PtGen_b;

  Int_t nGenElec;
  Int_t nGenMuon;
  Int_t nGenTau;
  Int_t nGenLepton;
  Int_t nTauElec;
  Int_t nTauMuon;
  Int_t nSGenMuon;
  Int_t nSGenElec;
  
  Int_t nGoodVertex;
  Int_t nVertex;
  Int_t nBtags;
  Int_t nJets;
  Int_t nMuon;
  Int_t nElec;
  Int_t nLeptons;

  ///// OBJECTS
  std::vector<lepton> Lepton;
  std::vector<jet>    Jet;
  //  std::vector<jet>    Jet15;
  
  std::vector<float> JetEt;
  std::vector<float> JetPt;
  std::vector<float> JetPhi;
  std::vector<float> MuPx;
  std::vector<float> MuPy;
  std::vector<float> MuPz;
  std::vector<float> MuEnergy;
  std::vector<float> ElPx;
  std::vector<float> ElPy;
  std::vector<float> ElPz;
  std::vector<float> ElEnergy;
  float MET;
  float MET_Phi;
  
  ClassDef(TreeAnalysisTop,0);
};
