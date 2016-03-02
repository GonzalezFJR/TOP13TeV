#ifndef TOPPLOTTER_H
#define TOPPLOTTER_H

#include<iostream>
#include <fstream>
#include<iomanip>      
#include<sstream>

#include "TH1F.h"
#include "THStack.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"

#include "TSystem.h"
#include "TROOT.h"
#include "tdrstyle.h"

#include "TH2F.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

enum gChannel{
  channels_begin,
  Muon = channels_begin,
  Elec,
  ElMu,
  gNCHANNELS,
};  
TString gChanLabel[gNCHANNELS+1] = {"Muon","Elec","ElMu","EEMM"};
enum iCut{
  iDilepton, 
  iZVeto, 
  iMET, 
  i2jets, 
  i1btag, 
  iNCUTS
};
TString sCut[iNCUTS] = {"dilepton", "ZVeto", "MET", "2jets", "1btag"}; 
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
  PUUp,
  PUDown,
  TopPt,
  gNSYST,
  Q2ScaleUp=gNSYST,
  Q2ScaleDown,
  MatchingUp,
  MatchingDown,
  gNALLSYST
};
TString sysname[gNALLSYST]={
  "",
  "__btag__plus",
  "__btag__minus",
  "__mistag__plus",
  "__mistag__minus",
  "__jes__plus",
  "__jes__minus",
  "__jer",
  "__les__plus",
  "__les__minus",
  "__pu__plus",
  "__pu__minus",
  "__toppt",
  "__q2__minus",
  "__q2__plus",
  "__mpts__plus",
  "__mpts__minus"
};
TString SystName[gNALLSYST] = {
  "Nominal",
  "BtagUp",
  "BtagDown",
  "MisTagUp",
  "MisTagDown",
  "JESUp",
  "JESDown",
  "JER",
  "LESUp",
  "LESDown",
  "PUUp",
  "PUDown",
  "TopPt",
  "Q2ScaleUp",
  "Q2ScaleDown",
  "MatchingUp",
  "MatchingDown"
};
enum gSystErrType{
  SFIDISO,
  SFTrig,
  les,
  jes,
  jer,
  btag,
  mistag,
  pu,
  toppt,
  had,
  pdf,
  cr,
  Q2,
  Matching,
  gNSYSTERRTypes,
  stop=gNSYSTERRTypes,
  vv,
  dy,
  fake,
  rare,
  gNSYSTERRTypesALL
};
TString SystErrLabel[gNSYSTERRTypesALL] = {"SFIDISO", 
					   "SFTrig",
					   "LES",
					   "JES",
					   "JER",
					   "btag",
					   "mistag",
					   "pu"
					   "toppt",
					   "had",
					   "pdf"
					   "cr",
					   "Q2",
					   "Matching",
					   "STop",
					   "VV",
					   "DrellYan",
					   "Fake"
					   "Rare"};
enum Samples{
  DoubleEG        ,
  DoubleMuon		,
  MuEG  	        ,
  TTJets_MadSpin        ,
  ////TTJetsFullLeptMGtauola,
  TTJetsSemiLeptMGtauola,
  // TTbar systematics
  ////TTJets_matchingup     ,
  ////TTJets_matchingdown   ,
  TTJets_scaleup        ,
  TTJets_scaledown      ,
  ////TTJets_MadSpinPDF     ,
  ////TTbar_Powheg          ,
  ////TTbar_Powheg_Herwig   ,
  ////TTJetsFullLeptMGTuneP11,
  ////TTJetsFullLeptMGTuneP11noCR,
  TbarWDilep		,
  TWDilep		,
  DYJets_Madgraph	,
  ZJets_Madgraph	,
  Wbb_Madgraph		,
  ////WgammaToLNuG		,
  WWTo2L2Nu_Madgraph	,
  WZ			,
  ZZ			,
  //  TTGJets		,
  ////TTWJets		,
  ////TTWWJets		,
  ////TTZJets		,
  ////WWWJets		,
  ////WWZJets		,
  ////WZZJets		,
  ////ZZZJets               ,
  ////T2tt_150to250LSP1to100_LeptonFilter,
  gNSAMPLES
};
TString SampleName[gNSAMPLES] = { 
  "DoubleEGsum"  ,
  "DoubleMuonSum",
  "MuonEGsum"    ,
  "TTbar_Powheg", //"TTbar_Powheg"    , 
  "TTbarSemi_Powheg", 
  // TTbar systematics
  ////"TTJets_matchingup"     ,
  ////"TTJets_matchingdown"   ,
  "TTbar_Powheg_scaleUp"      ,    
  "TTbar_Powheg_scaleDown_ext"    ,   
  ////"TTJets_MadSpinPDF"     ,
  ////"TTbar_Powheg"          ,
  ////"TTbar_Powheg_Herwig"   ,
  ////"TTJetsFullLeptMGTuneP11",
  ////"TTJetsFullLeptMGTuneP11noCR",
  "TbarW"                      , 
  "TW"                         , 
  "DYJetsToLL_M10to50_aMCatNLO_ext", 
  "DYJetsToLL_M50_aMCatNLO"    , 
  "WJetsToLNu_aMCatNLO"        , 
  "WW"                         ,
  "WZ"	                       ,
  "ZZ"			  
};
//++ categories
enum iVar{
  MET,     
  InvMass,
  DiLepPt,
  Lep0Pt,
  Lep1Pt,
  DelLepPhi,
  NJets,
  NBtagJets,
  Jet0Pt,
  Jet1Pt,
  NBTagsNJets,
  DelPhillJet,
  AbsDelPhiLeps,
  delPhi2LeadJets,
  minDelRJetsLeps,
  Vtx,
  gNVARS
};
TString KinVarName[gNVARS] = {"MET","InvMass","DiLepPt", "Lep0Pt","Lep1Pt","DelLepPhi","NJets",
			     
"NBtagJets","Jet0Pt","Jet1Pt","NBtagsNJets","DelPhillJet","AbsDelPhiLeps", "delPhi2LeadJets", "minDelRJetsLeps", "Vtx"}; 
TString KinAxisLabel[gNVARS] = {"E_{T}^{miss} (GeV)",
				"M_{ll} (GeV)",
				"Di-Lepton p_{T} (GeV)",
				"Leading Lepton p_{T} (GeV)",
				"Subleading Lepton p_{T} (GeV)",
				"#Delta #phi (ll)",
				"Number of Jets (E_{T} > 30 GeV)",
				"Number of b-jets",
				"Leading Jet p_{T} (GeV)",
				"Subleading Jet p_{T} (GeV)",
				"Number of b-jets (N_{jets})",
				//				"CSV Tagger of the Leading Jet",
				//				"Top Discrimator",
				"#Delta #phi (ll,jet)",
				"|#Delta #phi (l,l)|/#Pi",
				"|#Delta #phi (j1,j2)|",
				"min(#Delta R(j1,l), #Delta R(j2,l))",
				"Number of vertices",
				};
class TopPlotter {
 public:
  TopPlotter();
  ~TopPlotter(){};
  
  void Init(TString);
  void Loop();
  void LoadSamples(TString);
  //  void LoadCategory(Categories &, TString);
  void LoadCategories();
  void ResetDataMembers();
  void ResetSystematicErrors();
  struct Categories{
    TString name;
    float Yields       [gNCHANNELS]  [iNCUTS];
    float Yields_syst  [gNCHANNELS]  [iNCUTS][gNSYST]; 
    float Yields_stat  [gNCHANNELS]  [iNCUTS];
    float SSYields     [gNCHANNELS]  [iNCUTS];
    float SSYields_stat[gNCHANNELS]  [iNCUTS];
    TH1F* KinHistos    [gNCHANNELS+1][iNCUTS][gNVARS];
    TH1F* SysHistos    [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* SSSysHistos  [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* NBtagsNJets  [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* SSNBtagsNJets[gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* InvMass      [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* SSInvMass    [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* AbsDelPhiLeps      [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* SSAbsDelPhiLeps    [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* delPhi2LeadJets    [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* SSdelPhi2LeadJets  [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* minDelRJetsLeps    [gNCHANNELS][iNCUTS][gNALLSYST]; 
    TH1F* SSminDelRJetsLeps  [gNCHANNELS][iNCUTS][gNALLSYST]; 
    float SystError    [gNCHANNELS][gNSYSTERRTypesALL];
    TH1F* MllHistos    [gNCHANNELS][iNCUTS];
    TH1F* pdfWeights;
    TH1F* pdfWeightsSum;
  };
  
  struct XSection {
    float xsec     [gNCHANNELS];
    float xsec_syst[gNCHANNELS];
    float xsec_stat[gNCHANNELS];
    float xsec_lumi[gNCHANNELS];
    float acc      [gNCHANNELS];
    float acc_stat [gNCHANNELS];
    float acc_syst [gNCHANNELS];
    
    float err_VV   [gNCHANNELS];
    float err_DY   [gNCHANNELS];
    float err_STop [gNCHANNELS];
    float err_Fake [gNCHANNELS];
    float err_Rare [gNCHANNELS];
    float err_IDIso[gNCHANNELS];
    float err_Trig [gNCHANNELS];
    float err_LES  [gNCHANNELS];
    float err_JES  [gNCHANNELS];
    float err_JER  [gNCHANNELS];
    float err_Btag [gNCHANNELS];
    float err_mtag [gNCHANNELS];
    float err_TopPt[gNCHANNELS];
    float err_cr   [gNCHANNELS];
    float err_pdf  [gNCHANNELS];
    float err_had  [gNCHANNELS];
    float err_PU   [gNCHANNELS];
    float err_Q2   [gNCHANNELS];
    float err_Match[gNCHANNELS];
  };
  void PrintYieldsWithMC();
  void PrintYieldsWithDD();

  void PrintSystematicErrors();
  void CalculateSystematicErrors(Categories&, Int_t);
  void CalculateNonWZLeptonsBkg();
  void CalculateDYBkg();
  void CalculateCrossSection(Bool_t DD = false);
  void CalculateSystematicErrorsWithXSec(Categories&, Int_t);
  float GetPDFUncertainty();

  void DrawKinematicPlots(Bool_t DD,Int_t onechan = -1,Int_t onevar = -1, Int_t onecut =-1);
  //void DrawKinematicPlotsWithDD(Int_t onechan = -1,Int_t onevar = -1, Int_t onecut =-1);
  void DrawNbjetsNjets(bool);
  void SaveHistosForLH(bool);

  void SetOutputDir   (TString dir){ fOutputDir    = dir; };
  void SetOutputSubDir(TString dir){ fOutputSubDir = dir; };
  void SetVerbose(Int_t v)         { fVerbose      = v;   };
  
  void SetupDraw(TH1F*,int, int);
  TH1F* GetHisto1D(TFile*, TString);
  TH2F* GetHisto2D(TFile*, TString);
  Int_t GetRebin(TH1F*, Int_t);
  void DrawTopLine(Int_t, Float_t=0.93);

 private:
  Categories S[gNSAMPLES];

  Categories Data ;
  Categories TTbar;
  Categories STop ;
  Categories DY   ;
  Categories VV   ;
  Categories Rare ;
  Categories Fake ;
  Categories Total;
  Categories SUSYstop;
  
  // Data Driven backgrounds.
  Categories DD_DY  ;
  Categories DD_NonW;
  Float_t    DY_SF[gNCHANNELS][iNCUTS];
  
  XSection ttbar;
  
  float ttbar_TLWG;
  
  // Other helper 
  TString fOutputDir;
  TString fOutputSubDir;
  Int_t   fVerbose;
  
  Bool_t DoDF;
  Bool_t DoSF;
  
  float toppt_weight;

  std::ofstream fOUTSTREAM, fOUTSTREAM2, fOUTSTREAM3;

};
#endif
