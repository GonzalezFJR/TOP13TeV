#pragma once
// Includes
// + PAF
#include "PAFChainItemSelector.h"

// + Packages
//#include "GlobalVariables.h"
//#include "PUWeight.h"
//#include "BTagSFUtil.h"
//#include "LeptonSF.h"

// + ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TRandom3.h"

// + C++
#include <vector>

using namespace std;

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
  ElecUp,
  ElecDown,
  TrigUp,
  TrigDown,
  MuonUp,
  MuonDown,
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
  "ElecUp",
  "ElecDown",
  "TrigUp",
  "TrigDown",
  "MuonUp",
  "MuonDown",
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
 
 
const Int_t maxdim = 60;
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
class TOP5TeVAnalyzer : public PAFChainItemSelector 
{
 public:
  // Constructor and destructor
  TOP5TeVAnalyzer();
  virtual ~TOP5TeVAnalyzer() {}
  
  // Mandatory methods
  virtual void Initialise();
  virtual void InsideLoop();
    virtual void Summary();

  // Frequantly used variables in the trees
  
  Int_t nElec;

  Float_t ElecPt[maxdim];
  Float_t ElecPx[maxdim];
  Float_t ElecPy[maxdim];
  Float_t ElecPz[maxdim];
  Float_t ElecEnergy[maxdim];
  Float_t ElecEta[maxdim];
  Int_t ElecCharge[maxdim];
  Int_t ElecIDVeto[maxdim];
  Int_t ElecIDLoose[maxdim];
  Int_t ElecIDMedium[maxdim];
  Int_t ElecIDTight[maxdim];

  Int_t nMuon;
  Float_t MuonPt[maxdim];
  Float_t MuonPx[maxdim];
  Float_t MuonPy[maxdim];
  Float_t MuonPz[maxdim];
  Float_t MuonEnergy[maxdim];
  Float_t MuonEta[maxdim];
  Int_t MuonCharge[maxdim];
  Float_t MuonDxy[maxdim];
  Float_t MuonDz[maxdim];
  Float_t MuonChi2NDF[maxdim];
  Int_t MuonPixelHits[maxdim];
  Int_t MuonStations[maxdim];
  Int_t MuonTrkLayers[maxdim];
  Int_t MuonTrkQuality[maxdim];

  Int_t ngenLep;
  Int_t genLep_pdgId[maxdim];
  Float_t genLep_pt[maxdim];
  Float_t genLep_eta[maxdim];
  Float_t genLep_phi[maxdim];
  Float_t genLep_energy[maxdim];
  Int_t genLep_status[maxdim];
  Int_t genLep_MomPID[maxdim];
  Int_t genLep_GMomPID[maxdim];

  Int_t nJet;
  Float_t Jet_px[maxdim];
  Float_t Jet_py[maxdim];
  Float_t Jet_pz[maxdim];
  Float_t Jet_energy[maxdim];
  Float_t Jet_eta[maxdim];
  Float_t Jet_btagCSV[maxdim];
  Int_t Jet_mcFlavour[maxdim];

  Float_t  JetPfCHF[maxdim];
  Float_t  JetPfNHF[maxdim];
  Float_t  JetPfCEF[maxdim];
  Float_t  JetPfNEF[maxdim];
  Float_t  JetPfMUF[maxdim];
  Int_t    chargedMult[maxdim];
  Int_t    totalMult[maxdim];

  Int_t   ngenJet;
  Float_t genJet_pt[maxdim];
  Float_t genJet_eta[maxdim];
  Float_t genJet_phi[maxdim];
  Float_t genJet_m[maxdim];
  Int_t   genJet_matchId[maxdim];

  Int_t evt;
  Int_t run;
  Int_t lum;
  Float_t met_pt;
  Float_t met_phi;
  Float_t genWeight;
  Int_t HLT_HIL2Mu15_v1;
  Int_t HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1;

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
  float getSF(gChannel, gSystFlag sys = Norm);
  float getLeptonError(gChannel);
  float getTriggerError(gChannel);
  float getTopD();
  float getDeltaPhillJet();
  float getMT(gChannel);
  float getDelPhill();
  Int_t getJet_id(int);

  TLorentzVector getPtllb();

  // Lepton selection methods
  int  getSelectedLeptons();
//   bool IsVetoMuon(unsigned int, float ptcut=20.);
  bool IsTightMuon(unsigned int, float ptcut=18.);
	float getMuonIso(unsigned int iMuon);
//   bool IsVetoElectron(unsigned int,float ptcut=20.);
//   bool IsMVAIDElectron(unsigned int);
  void CoutEvent(long unsigned int en = 0, TString t = " ");
  bool IsTightElectron(unsigned int,float ptcut=20.);
  std::vector<lepton> SortLeptonsByPt(std::vector<lepton>&);
  
  int getSelectedJets();
  bool IsGoodJet(unsigned int, float ptcut=25.);
  std::vector<int> CleanedJetIndices(float);
  void SmearJetPts(int);
	float getJetSF(int index, bool up = false);
	float getJerSF(float eta);
	float getJerSFerror(float eta);
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

  //PUWeight *fPUWeight;      //The PU weight utility
  //PUWeight *fPUWeightUp;    //The PU weight utility
  //PUWeight *fPUWeightDown;  //The PU weight utility
  //BTagSFUtil *fBTagSF[5]; //The new BTag SF utility 
  //BTagSFUtil *fBTagSFnom ;
  //BTagSFUtil *fBTagSFbUp ;
  //BTagSFUtil *fBTagSFbDo ;
  //BTagSFUtil *fBTagSFlUp ;
  //BTagSFUtil *fBTagSFlDo ;
  //LeptonSF *fLeptonSF;
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
	TH1F* fHDelPhillJet[gNCHANNELS][iNCUTS];

	TH1F* fMuonIsoCharged[gNCHANNELS][iNCUTS];
	TH1F* fMuonIsoNeutral[gNCHANNELS][iNCUTS];
	TH1F* fMuonIsoPhotons[gNCHANNELS][iNCUTS];
	TH1F* fMuonIsoPU[gNCHANNELS][iNCUTS];       
	TH1F* fMuonIso[gNCHANNELS][iNCUTS];          

	TH1F* fHDRLep[gNCHANNELS][iNCUTS];
	TH1F* fHDRLep0Jet[gNCHANNELS][iNCUTS];
	TH1F* fHDPhiLep0Jet[gNCHANNELS][iNCUTS];
	TH1F* fHDRLep1Jet[gNCHANNELS][iNCUTS];
  TH1F* fHDPhiLep1Jet[gNCHANNELS][iNCUTS];
  
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
  Int_t nMuons;
  Int_t nElecs;
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
  
  ClassDef(TOP5TeVAnalyzer,0);
};
