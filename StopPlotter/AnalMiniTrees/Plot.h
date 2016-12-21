#ifndef Plot_h
#define Plot_h 1


#include "AnalMiniTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TString.h"
#include <iostream>

//##################################################
//# Configuration
//##################################################
const bool doSys       = false;
const bool doData      = false;
const bool doYields    = true;
const bool doSingleLep = false;
const bool doStackOverflow = true;
const bool doSignal = true;
const bool save = false;

const TString Chan = "ElMu";
const TString ThePlot = "MT2";
const TString plotfolder = "/nfs/fanae/user/juanr/stop/TOP13TeV/StopPlotter/AnalMiniTrees/";
const TString LimitFolder = "/nfs/fanae/user/juanr/CMSSW_7_4_7_patch1/src/StopDatacards/";
//const float Lumi = 12.85; //fb-1
const float Lumi = 35; //fb-1

class Plot {
public:
	Histo* httbar; Histo* hDY; Histo* hWJets; Histo* hVV; Histo* htW; Histo* httV;
	Histo* hData; Histo* hElMu; Histo* hMuon; Histo* hElec; Histo* hSignal; Histo* AllBkg;
  Bool_t  doSetLogy   = true;
  THStack* hStack;
  Int_t nBkgs = 0;
	TPad* plot; TPad* pratio;
	TLegend* leg;
	TLatex* texlumi;
	TLatex* texcms;
	TLatex* texchan;
	TH1F* hratio;
  std::vector<Histo*> VBkgs;
  std::vector<Histo*> VSignals;
  std::vector<Histo*> VSyst;
  std::vector<TString> VSystLabel;
  TString sys = "0";
  Float_t* TotalSysUp;
  Float_t* TotalSysDown;

  Float_t sys_lumi = 0.10;

  virtual ~Plot(){};            // Destructor
	//Plot(TString v = "MET", TString ch = "ElMu", Bool_t logy = 1, TString signal = "T2tt_mStop183_mLsp1");
	Plot(TString v = "MET", TString ch = "ElMu", Bool_t logy = 1, TString sig = "T2tt_mStop183_mLsp1"){ // Constructor
		var = v; 
		chan = ch;
		signal = sig;
		SetAllProcesses();
		SetAllBkg();
		SetData();
		SetSignal();
    //SetOtherSignal();
	  //AddAllSystematics();
	}


  // ######### Methods ########
  Histo* GetH(TString sample);
  TCanvas *SetCanvas();
  void SetLegend(bool doyi);
  void SetTexChan(TString cuts);
  void SetHRatio();

	void SetTTbar();
	void SetDY();
	void SetWJets();
	void SetVV();
	void SettW();
	void SetttV();
	void SetData();
	void SetAllProcesses();
	void SetAllBkg();
	void SetStack();
  void SetSignal();
  void SetPlotStyle();

	void DrawStack(TString tag, bool sav);
  void SaveHistograms(TString tag);
  TString GetStatUncDatacard();
	TString GetShapeUncDatacard();
  void MakeDatacard(TString tag);
  void MakeDatacardAllBins(TString tag);
  void MakeDatacardBin(Int_t bin, TString tag);

  TString getVar();
  TString getChan();
  TString getSignal();
  void setVar(TString v);
  void setChan(TString c);
  void setSignal(TString s);
	void SetOtherSignal(TString signal2 = "T2tt_mStop250_mLsp75", Int_t ccol = kGreen+4);

  void AddToHistos(Histo* p);
	void AddSystematic(TString s);
	void AddAllSystematics();

protected: 
	TString var;
  TString chan;
  TString signal;
};

#endif
