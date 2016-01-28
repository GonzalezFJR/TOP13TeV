#include <iostream>
#include "TH1F.h"
#include "THStack.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"

float fLumiNorm = 19664.225;

enum gChannel{
  channels_begin,
  Muon = channels_begin,
  Elec,
  ElMu,
  gNCHANNELS,
};
TString gChanLabel[gNCHANNELS] = {"Muon","Elec","ElMu"};
enum iCut{
  iDilepton, 
  iZVeto, 
  i2jets, 
  iMET, 
  i1btag, 
  i2btag,
  iNCUTS
};
TString sCut[iNCUTS] = {"dilepton", "ZVeto", "2jets", "MET", "1btag", "2btag"};
enum iVar{
  MET,     
  InvMass,
  Lep0Pt,
  Lep1Pt,
  DelLepPhi,
  NJets,
  NBtagJets,
  Jet0Pt,
  gNVARS
};
TString sVar[gNVARS]     = {"MET","InvMass","Lep0Pt","Lep1Pt","DelLepPhi","NJets","NBtagJets","Jet0Pt"};

enum Samples{
  DoubleElectron        ,
  DoubleMu   	        ,
  MuEG   		,
  TTJets_MadSpin,
  TTJetsSemiLeptMGtauola,
  TbarWDilep		,
  TWDilep		,
  DYJets_Madgraph	,
  ZJets_Madgraph	,
  Wbb_Madgraph		,
  WgammaToLNuG		,
  WWTo2L2Nu_Madgraph	,
  WZ			,
  ZZ			,
  TTGJets		,
  TTWJets		,
  TTWWJets		,
  TTZJets		,
  WWWJets		,
  WWZJets		,
  WZZJets		,
  ZZZJets               ,
  gNSAMPLES
};

TString SampleName[gNSAMPLES] = { 
  "DoubleElectron"        ,
  "DoubleMu"  	    	  ,
  "MuEG"    	      	  ,
  //  "TTJets_MadSpin",
  "TTJetsFullLeptMGtauola",
  "TTJetsSemiLeptMGtauola",
  "TbarWDilep"		  ,
  "TWDilep"		  ,
  "DYJets_Madgraph"	  ,
  "ZJets_Madgraph"	  ,
  "Wbb_Madgraph"	  	,
  "WgammaToLNuG"	  	,
  "WWTo2L2Nu_Madgraph"	  ,
  "WZ"			  ,
  "ZZ"			  ,
  "TTGJets"		  ,
  "TTWJets"		  ,
  "TTWWJets"		  ,
  "TTZJets"		  ,
  "WWWJets"  		  ,
  "WWZJets"		  ,
  "WZZJets"		  ,
  "ZZZJets"             
};
TH1F *Histo[gNSAMPLES][gNCHANNELS];

float   Yields[gNSAMPLES][gNCHANNELS][iNCUTS];
float SSYields[gNSAMPLES][gNCHANNELS][iNCUTS];

//++ categories
float Data [gNCHANNELS][iNCUTS];
float TTbar[gNCHANNELS][iNCUTS];
float STop [gNCHANNELS][iNCUTS];
float DY   [gNCHANNELS][iNCUTS];
float VV   [gNCHANNELS][iNCUTS];
float Rare [gNCHANNELS][iNCUTS];
float Fake [gNCHANNELS][iNCUTS];
float Total[gNCHANNELS][iNCUTS];

float SS_Data [gNCHANNELS][iNCUTS];
float SS_TTbar[gNCHANNELS][iNCUTS];
float SS_STop [gNCHANNELS][iNCUTS];
float SS_DY   [gNCHANNELS][iNCUTS];
float SS_VV   [gNCHANNELS][iNCUTS];
float SS_Rare [gNCHANNELS][iNCUTS];
float SS_Fake [gNCHANNELS][iNCUTS];
float SS_Total[gNCHANNELS][iNCUTS];

void MakeYields(TString pathtofiles="TopTrees/ReReco_CSVM_Feb11/"){
  
  TFile *_file ;
  TH1F *hyields;
  TH1F *hSSyields;
  
  for (size_t sample=0; sample<gNSAMPLES; sample++){
    _file = new TFile(pathtofiles + "Tree_Legacy_"+SampleName[sample]+".root");
    
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      hyields   = (TH1F*) _file->Get("H_Yields_"  +gChanLabel[chan]);
      hSSyields = (TH1F*) _file->Get("H_SSYields_"+gChanLabel[chan]);
      
      for (size_t cut=0; cut<iNCUTS; cut++){
	if (sample == DoubleElectron || sample == DoubleMu || sample == MuEG){ 
	  Yields  [sample][chan][cut] =   hyields->GetBinContent(cut+1);
	  SSYields[sample][chan][cut] = hSSyields->GetBinContent(cut+1);
	}
	else {
	  Yields  [sample][chan][cut] =   hyields->GetBinContent(cut+1) * fLumiNorm;
	  SSYields[sample][chan][cut] = hSSyields->GetBinContent(cut+1) * fLumiNorm;
	}
	if (sample == 0){
	  Data [chan][cut] = 0.;  SS_Data [chan][cut] = 0.;  
	  TTbar[chan][cut] = 0.;  SS_TTbar[chan][cut] = 0.;
	  STop [chan][cut] = 0.;  SS_STop [chan][cut] = 0.;
	  DY   [chan][cut] = 0.;  SS_DY   [chan][cut] = 0.;
	  VV   [chan][cut] = 0.;  SS_VV   [chan][cut] = 0.;
	  Rare [chan][cut] = 0.;  SS_Rare [chan][cut] = 0.;
	  Fake [chan][cut] = 0.;  SS_Fake [chan][cut] = 0.;
	  Total[chan][cut] = 0.;  SS_Total[chan][cut] = 0.;
	}
      }
    }
  }
  
  for (size_t cut=0; cut<iNCUTS-1; cut++){
    Data [Muon][cut] = Yields[DoubleMu]      [Muon][cut];
    Data [Elec][cut] = Yields[DoubleElectron][Elec][cut];
    Data [ElMu][cut] = Yields[MuEG]          [ElMu][cut];

    SS_Data [Muon][cut] = SSYields[DoubleMu]      [Muon][cut];
    SS_Data [Elec][cut] = SSYields[DoubleElectron][Elec][cut];
    SS_Data [ElMu][cut] = SSYields[MuEG]          [ElMu][cut];
    
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      TTbar[chan][cut]  = Yields[TTJets_MadSpin][chan][cut];
      
      STop [chan][cut]  = Yields[TbarWDilep][chan][cut];
      STop [chan][cut] += Yields[TWDilep]   [chan][cut];
            
      DY   [chan][cut]  = Yields[ZJets_Madgraph][chan][cut];
      DY   [chan][cut] += Yields[DYJets_Madgraph][chan][cut];
      
      VV   [chan][cut]  = Yields[WZ]                [chan][cut];
      VV   [chan][cut] += Yields[ZZ]                [chan][cut];
      VV   [chan][cut] += Yields[WWTo2L2Nu_Madgraph][chan][cut];
      
      Rare [chan][cut]  = Yields[TTWJets] [chan][cut];
      Rare [chan][cut] += Yields[TTZJets] [chan][cut];
      Rare [chan][cut] += Yields[TTGJets] [chan][cut];
      Rare [chan][cut] += Yields[TTWWJets][chan][cut];
      Rare [chan][cut] += Yields[WWWJets] [chan][cut];
      Rare [chan][cut] += Yields[WWZJets] [chan][cut];
      Rare [chan][cut] += Yields[WZZJets] [chan][cut];
      Rare [chan][cut] += Yields[ZZZJets] [chan][cut];

      Fake [chan][cut]  = Yields[TTJetsSemiLeptMGtauola][chan][cut];
      Fake [chan][cut] += Yields[Wbb_Madgraph]          [chan][cut];
      Fake [chan][cut] += Yields[WgammaToLNuG]          [chan][cut];

      Total[chan][cut]  = STop[chan][cut];
      Total[chan][cut] += DY  [chan][cut];
      Total[chan][cut] += VV  [chan][cut];
      Total[chan][cut] += Rare[chan][cut];
      Total[chan][cut] += Fake[chan][cut];
      
      // SAME SIGN YIELDS
      SS_TTbar[chan][cut]  = SSYields[TTJets_MadSpin][chan][cut];

      SS_STop [chan][cut]  = SSYields[TbarWDilep][chan][cut];
      SS_STop [chan][cut] += SSYields[TWDilep]   [chan][cut];

      SS_DY   [chan][cut]  = SSYields[ZJets_Madgraph][chan][cut];
      SS_DY   [chan][cut] += SSYields[DYJets_Madgraph][chan][cut];

      SS_VV   [chan][cut]  = SSYields[WZ]                [chan][cut];
      SS_VV   [chan][cut] += SSYields[ZZ]                [chan][cut];
      SS_VV   [chan][cut] += SSYields[WWTo2L2Nu_Madgraph][chan][cut];

      SS_Rare [chan][cut]  = SSYields[TTWJets] [chan][cut];
      SS_Rare [chan][cut] += SSYields[TTZJets] [chan][cut];
      SS_Rare [chan][cut] += SSYields[TTGJets] [chan][cut];
      SS_Rare [chan][cut] += SSYields[TTWWJets][chan][cut];
      SS_Rare [chan][cut] += SSYields[WWWJets] [chan][cut];
      SS_Rare [chan][cut] += SSYields[WWZJets] [chan][cut];
      SS_Rare [chan][cut] += SSYields[WZZJets] [chan][cut];
      SS_Rare [chan][cut] += SSYields[ZZZJets] [chan][cut];

      SS_Fake [chan][cut]  = SSYields[TTJetsSemiLeptMGtauola][chan][cut];
      SS_Fake [chan][cut] += SSYields[Wbb_Madgraph]          [chan][cut];
      SS_Fake [chan][cut] += SSYields[WgammaToLNuG]          [chan][cut];

      SS_Total[chan][cut]  = SS_STop[chan][cut];
      SS_Total[chan][cut] += SS_DY  [chan][cut];
      SS_Total[chan][cut] += SS_VV  [chan][cut];
      SS_Total[chan][cut] += SS_Rare[chan][cut];
      SS_Total[chan][cut] += SS_Fake[chan][cut];
    }
  }
  for (size_t cut=0; cut<iNCUTS-1; cut++){
    cout << "Cut level: " << sCut[cut] << endl;
    cout << "=========================================" << endl;
    cout << " Source           |  ElEl  |  MuMu  |  ElMu  |" << endl;
    cout << "-----------------------------------------" << endl;
    cout << Form(" Drell-Yan        | %5.1f | %5.1f | %5.1f |", 
		 DY[Elec][cut], DY[Muon][cut], DY[ElMu][cut] )<< endl;
    cout << Form(" Non W/Z leptons  | %5.1f | %5.1f | %5.1f |",  
		 Fake[Elec][cut], Fake[Muon][cut], Fake[ElMu][cut] )<< endl;
    cout << Form(" Single top quark | %5.1f | %5.1f | %5.1f |",  
		 STop[Elec][cut], STop[Muon][cut], STop[ElMu][cut] )<< endl;
    cout << Form(" Dibosons         | %5.1f | %5.1f | %5.1f |",  
		 VV[Elec][cut], VV[Muon][cut], VV[ElMu][cut] )<< endl;
    cout << Form(" Rare             | %5.1f | %5.1f | %5.1f |",  
		 Rare[Elec][cut], Rare[Muon][cut], Rare[ElMu][cut] )<< endl;
    cout << " ------------------------------------------" << endl;
    cout << Form(" Total background | %5.1f | %5.1f | %5.1f |",  
		 Total[Elec][cut], Total[Muon][cut], Total[ElMu][cut] )<< endl;
    cout << Form(" TTbar dilepton   | %5.1f | %5.1f | %5.1f |",  
		 TTbar[Elec][cut], TTbar[Muon][cut], TTbar[ElMu][cut] )<< endl;
    cout << " ------------------------------------------" << endl;
    cout << Form(" Data             | %5.0f | %5.0f | %5.0f |", 
		 Data[Elec][cut], Data[Muon][cut], Data[ElMu][cut] )<< endl;
    cout << endl;
    cout << " SAME SIGN: "                              << endl;
    cout << "=========================================" << endl;
    cout << " Source           |  ElEl  |  MuMu  |  ElMu  |" << endl;
    cout << "-----------------------------------------" << endl;
    cout << Form(" Drell-Yan        | %5.1f | %5.1f | %5.1f |", 
		 SS_DY[Elec][cut], SS_DY[Muon][cut], SS_DY[ElMu][cut] )<< endl;
    cout << Form(" Non W/Z leptons  | %5.1f | %5.1f | %5.1f |",  
		 SS_Fake[Elec][cut], SS_Fake[Muon][cut], SS_Fake[ElMu][cut] )<< endl;
    cout << Form(" Single top quark | %5.1f | %5.1f | %5.1f |",  
		 SS_STop[Elec][cut], SS_STop[Muon][cut], SS_STop[ElMu][cut] )<< endl;
    cout << Form(" Dibosons         | %5.1f | %5.1f | %5.1f |",  
		 SS_VV[Elec][cut], SS_VV[Muon][cut], SS_VV[ElMu][cut] )<< endl;
    cout << Form(" Rare             | %5.1f | %5.1f | %5.1f |",  
		 SS_Rare[Elec][cut], SS_Rare[Muon][cut], SS_Rare[ElMu][cut] )<< endl;
    cout << " ------------------------------------------" << endl;
    cout << Form(" Total background | %5.1f | %5.1f | %5.1f |",  
		 SS_Total[Elec][cut], SS_Total[Muon][cut], SS_Total[ElMu][cut] )<< endl;
    cout << Form(" TTbar dilepton   | %5.1f | %5.1f | %5.1f |",  
		 SS_TTbar[Elec][cut], SS_TTbar[Muon][cut], SS_TTbar[ElMu][cut] )<< endl;
    cout << " ------------------------------------------" << endl;
    cout << Form(" Data             | %5.0f | %5.0f | %5.0f |", 
		 SS_Data[Elec][cut], SS_Data[Muon][cut], SS_Data[ElMu][cut] )<< endl;
    cout << endl;
  }
}
void MakePlots(TString pathtofiles, iVar histovar, iCut cut){

  // for each histogram print mm (0), ee (1), DF (2) and  SF (3) channel
  TH1F* Data [4];
  TH1F* TTbar[4];
  TH1F* STop [4];
  TH1F* DY   [4];
  TH1F* VV   [4];
  TH1F* Rare [4];
  TH1F* Fake [4];
  TH1F* Total[4];

  THStack* MC[4];
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  gROOT->LoadMacro("~folgueras/TOP/TopCode/tdrstyle.C"); 
  setTDRStyle();
  
  cout << "Get Histograms from files... ";
  for (size_t sample=0; sample<gNSAMPLES; sample++){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      TString histoname = "H_" + sVar[histovar] + "_" + gChanLabel[chan] + "_" + sCut[cut];
      Histo[sample][chan] = GetHisto(pathtofiles+"Tree_Legacy_"+SampleName[sample]+".root", histoname);
    }
  }
  cout << "DONE!" << endl;
  
  // First pile-up DF histograms.
  for (size_t chan=0; chan<gNCHANNELS+1; chan++){
    if (chan == gNCHANNELS){
      Data [chan] = (TH1F*)Histo[DoubleMu]              [Muon]->Clone();
      Data [chan] ->Add(   Histo[DoubleElectron]        [Elec]);
      TTbar[chan] = (TH1F*)Histo[TTJets_MadSpin]        [Muon]->Clone();
      TTbar[chan] ->Add(   Histo[TTJets_MadSpin]        [Elec]);
      STop [chan] = (TH1F*)Histo[TbarWDilep]            [Muon]->Clone();
      STop [chan] ->Add(   Histo[TbarWDilep]            [Elec]);
      STop [chan] ->Add(   Histo[TWDilep]               [Muon]);
      STop [chan] ->Add(   Histo[TWDilep]               [Elec]);
      DY   [chan] = (TH1F*)Histo[DYJets_Madgraph]       [Muon]->Clone();
      DY   [chan] ->Add(   Histo[ZJets_Madgraph]        [Muon]);
      DY   [chan] ->Add(   Histo[DYJets_Madgraph]       [Elec]);
      DY   [chan] ->Add(   Histo[ZJets_Madgraph]        [Elec]);
      VV   [chan] = (TH1F*)Histo[WWTo2L2Nu_Madgraph]    [Muon]->Clone();
      VV   [chan] ->Add(   Histo[WZ]                    [Muon]);
      VV   [chan] ->Add(   Histo[ZZ]                    [Muon]);
      VV   [chan] ->Add(   Histo[WWTo2L2Nu_Madgraph]    [Elec]);
      VV   [chan] ->Add(   Histo[WZ]                    [Elec]);
      VV   [chan] ->Add(   Histo[ZZ]                    [Elec]);
      Fake [chan] = (TH1F*)Histo[TTJetsSemiLeptMGtauola][Muon]->Clone();  
      Fake [chan] ->Add(   Histo[WgammaToLNuG]          [Muon]);
      Fake [chan] ->Add(   Histo[Wbb_Madgraph]          [Muon]);
      Fake [chan] ->Add(   Histo[TTJetsSemiLeptMGtauola][Elec]);
      Fake [chan] ->Add(   Histo[WgammaToLNuG]          [Elec]);
      Fake [chan] ->Add(   Histo[Wbb_Madgraph]          [Elec]);
      Rare [chan] = (TH1F*)Histo[TTGJets]               [Muon]->Clone();  
      Rare [chan] ->Add(   Histo[TTWJets]               [Muon]);
      Rare [chan] ->Add(   Histo[TTWWJets]              [Muon]);
      Rare [chan] ->Add(   Histo[TTZJets]               [Muon]);
      Rare [chan] ->Add(   Histo[WWWJets]               [Muon]);
      Rare [chan] ->Add(   Histo[WWZJets]               [Muon]);
      Rare [chan] ->Add(   Histo[WZZJets]               [Muon]);
      Rare [chan] ->Add(   Histo[ZZZJets]               [Muon]);
      Rare [chan] ->Add(   Histo[TTGJets]               [Elec]);
      Rare [chan] ->Add(   Histo[TTWJets]               [Elec]);
      Rare [chan] ->Add(   Histo[TTWWJets]              [Elec]);
      Rare [chan] ->Add(   Histo[TTZJets]               [Elec]);
      Rare [chan] ->Add(   Histo[WWWJets]               [Elec]);
      Rare [chan] ->Add(   Histo[WWZJets]               [Elec]);
      Rare [chan] ->Add(   Histo[WZZJets]               [Elec]);
      Rare [chan] ->Add(   Histo[ZZZJets]               [Elec]);
    }
    else {
      if (chan == Muon){
	Data [chan] = (TH1F*)Histo[DoubleMu]              [chan]->Clone();
      }
      else if (chan == Elec){
	Data [chan] = (TH1F*)Histo[DoubleElectron]        [chan]->Clone();
      }
      else if (chan == ElMu){
	Data [chan] = (TH1F*)Histo[MuEG]                  [chan]->Clone();
      }
      TTbar[chan] = (TH1F*)Histo[TTJets_MadSpin]        [chan]->Clone();
      STop [chan] = (TH1F*)Histo[TbarWDilep]            [chan]->Clone();
      STop [chan] ->Add(   Histo[TWDilep]               [chan]);
      DY   [chan] = (TH1F*)Histo[DYJets_Madgraph]       [chan]->Clone();
      DY   [chan] ->Add(   Histo[ZJets_Madgraph]        [chan]);
      VV   [chan] = (TH1F*)Histo[WWTo2L2Nu_Madgraph]    [chan]->Clone();
      VV   [chan] ->Add(   Histo[WZ]                    [chan]);
      VV   [chan] ->Add(   Histo[ZZ]                    [chan]);
      Fake [chan] = (TH1F*)Histo[TTJetsSemiLeptMGtauola][chan]->Clone();  
      Fake [chan] ->Add(   Histo[WgammaToLNuG]          [chan]);
      Fake [chan] ->Add(   Histo[Wbb_Madgraph]          [chan]);
      Rare [chan] = (TH1F*)Histo[TTGJets]               [chan]->Clone();  
      Rare [chan] ->Add(   Histo[TTWJets]               [chan]);
      Rare [chan] ->Add(   Histo[TTWWJets]              [chan]);
      Rare [chan] ->Add(   Histo[TTZJets]               [chan]);
      Rare [chan] ->Add(   Histo[WWWJets]               [chan]);
      Rare [chan] ->Add(   Histo[WWZJets]               [chan]);
      Rare [chan] ->Add(   Histo[WZZJets]               [chan]);
      Rare [chan] ->Add(   Histo[ZZZJets]               [chan]);
    }
  }
 
  // Now change the style of each of the histos...
  // TTbar, kRed+1
  // STop , kPink-3
  // Fake , kGreen-3
  // VV   , kOrange-3
  // DY   , kAzure-2
  // Rare , kYellow
  // Data , kBlack
  for (size_t i=0; i<4; i++){
    SetupDraw(Data [i], kBlack   , histovar ,  cut); 
    SetupDraw(TTbar[i], kRed+1   , histovar ,  cut);  TTbar[i]->Scale(fLumiNorm);
    SetupDraw(STop [i], kPink-3  , histovar ,  cut);  STop [i]->Scale(fLumiNorm);
    SetupDraw(DY   [i], kAzure-2 , histovar ,  cut);  DY   [i]->Scale(fLumiNorm);
    SetupDraw(VV   [i], kOrange-3, histovar ,  cut);  VV   [i]->Scale(fLumiNorm);
    SetupDraw(Rare [i], kYellow  , histovar ,  cut);  Rare [i]->Scale(fLumiNorm);
    SetupDraw(Fake [i], kGreen-3 , histovar ,  cut);  Fake [i]->Scale(fLumiNorm);
    
    MC[i]->Add(TTbar[i]);  
    MC[i]->Add(STop [i]);  
    MC[i]->Add(DY   [i]);  
    MC[i]->Add(VV   [i]);  
    MC[i]->Add(Rare [i]);  
    MC[i]->Add(Fake [i]);  

    Total[i] = (TH1F*)TTbar[i]->Clone(); 
    Total[i]->Add(STop [i]);
    Total[i]->Add(DY   [i]);
    Total[i]->Add(VV   [i]);
    Total[i]->Add(Rare [i]);
    Total[i]->Add(Fake [i]);
  }

  TLegend *leg = new TLegend(0.73,0.58,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(Data [0], "Data", "PL");
  leg->AddEntry(Rare [0], "t#bar{t}V, VVV", "F");
  leg->AddEntry(VV   [0], "VV", "F");
  leg->AddEntry(Fake [0], "Non W/Z", "F");
  leg->AddEntry(STop [0], "Single Top", "F");
  leg->AddEntry(DY   [0], "Z/\\gamma* l^{+}l^{-}", "F");
  leg->AddEntry(TTbar[0],"t#bar{t}", "F");
   
  cout << "Drawing everything.. " << endl;
  TCanvas *c1 = new TCanvas("c1","Plot");
  c1->Divide(1,2);
   
  //Plot Pad
  TPad *plot = (TPad*)c1->GetPad(1); 
  plot->SetPad(0.01, 0.23, 0.99, 0.99);
  plot->SetTopMargin(0.1);
  plot->SetRightMargin(0.04);
  
  TPad *ratio = (TPad*)c1->GetPad(2);
  ratio->SetPad(0.01, 0.02, 0.99, 0.3);
  ratio->SetGridx();
  ratio->SetGridy();
  ratio->SetTopMargin(0.05);
  ratio->SetBottomMargin(0.4);
  ratio->SetRightMargin(0.04);
  
  for (size_t i=0; i<4; i++){
    plot->cd();
    MC[i]->Draw("hist");
    
    MC[i]->GetYaxis()->SetTitle("Events");
    MC[i]->GetYaxis()->SetTitleOffset(1.2);
    MC[i]->GetYaxis()->SetTitleSize(0.07);
    MC[i]->GetYaxis()->SetLabelSize(0.055);
    MC[i]->GetYaxis()->SetNdivisions(607);
    
    MC[i]->GetXaxis()->SetLabelSize(0.0);
    MC[i]->GetXaxis()->SetTitle("");
     
    Data[i]->Sumw2();
    Data[i]->SetMarkerStyle(20);
    Data[i]->Draw("SAME");
    
    leg->Draw("SAME");
    
    float maxh = Data[i]->GetMaximum();
    if (maxh < MC[i]->GetMaximum()) maxh = MC[i]->GetMaximum();
    MC[i]->SetMaximum(1.7*maxh);
    MC[i]->SetMinimum(50);

    // set logY
    if (cut == iDilepton) plot->SetLogy();
    if (cut == iZVeto   ) plot->SetLogy();
    
    ratio->cd();
     
    TH1F *H_Ratio;
    H_Ratio = (TH1F*)Data[i]->Clone();
    H_Ratio->Divide(Total[i]);
    
    H_Ratio->SetFillStyle(1001);
    H_Ratio->SetLineWidth(1);
    H_Ratio->SetFillColor(  kGray+1);
    H_Ratio->SetLineColor(  kGray+1);
    H_Ratio->SetMarkerColor(kGray+1);
    H_Ratio->SetMarkerStyle(20);
    //	  H_Ratio->SetMarkerColor(1);
    //    H_Ratio->SetLineColor(0);
    H_Ratio->SetMaximum(2);
    H_Ratio->SetMinimum(0);
    H_Ratio->SetTitle("");

    H_Ratio->GetYaxis()->SetTitle("Obs/Exp");
    H_Ratio->GetYaxis()->CenterTitle();
    H_Ratio->GetYaxis()->SetTitleOffset(0.45);
    H_Ratio->GetYaxis()->SetTitleSize(0.16);
    H_Ratio->GetYaxis()->SetLabelSize(0.15);
    H_Ratio->GetYaxis()->SetNdivisions(402);
    H_Ratio->GetXaxis()->SetTitleOffset(1.1);
    H_Ratio->GetXaxis()->SetLabelSize(0.15);
    H_Ratio->GetXaxis()->SetTitleSize(0.16);

    H_Ratio->SetMinimum(0.6);
    H_Ratio->SetMaximum(1.4);

    H_Ratio->DrawCopy("E2");
    TGraphErrors *thegraphRatio = new TGraphErrors(H_Ratio);
    thegraphRatio->SetFillStyle(3001);
    thegraphRatio->SetFillColor(kGray+1);
//   
    thegraphRatio->Draw("E2 SAME");
    
    /***********************
           CMS Legend
    **********************/
    c1->cd();
    plot->cd();
    
    TString htitleCMSChannel;
    if(i==0) htitleCMSChannel="#mu#mu";
    if(i==1) htitleCMSChannel="ee";
    if(i==2) htitleCMSChannel="e#mu";
    if(i==3) htitleCMSChannel="ee#mu#mu";
    
    title  = new TLatex(-20.,50.,"CMS Preliminary, #sqrt{s}=8TeV, 19.5 fb^{-1}");
    title->SetNDC();
    title->SetTextAlign(12);
    title->SetX(0.16);
    title->SetY(0.93);
    title->SetTextFont(42);
    title->SetTextSize(0.04);
    title->SetTextSizePixels(22);
    title->Draw("SAME");

    chtitle  = new TLatex(-20.,50.,htitleCMSChannel);
    chtitle->SetNDC();
    chtitle->SetTextAlign(12);
    chtitle->SetX(0.85);
    chtitle->SetY(0.93);
    chtitle->SetTextFont(42);
    chtitle->SetTextSize(0.04);
    chtitle->SetTextSizePixels(22);
    chtitle->Draw("SAME");
   
    TString dirfigname_pdf= pathtofiles + "figures/";
    gSystem->mkdir(dirfigname_pdf,       kTRUE);
    
    TString channel = "";
    if (i==0) channel = "_MM";
    if (i==1) channel = "_EE";
    if (i==2) channel = "_DF";
    if (i==3) channel = "_SF";
    
    dirfigname_pdf += sVar[histovar] + channel + "_" + sCut[cut];
    c1->SaveAs(dirfigname_pdf + ".pdf" );
    c1->SaveAs(dirfigname_pdf + ".png" );
    c1->SaveAs(dirfigname_pdf + ".root");
  }
  cout << "DONE!"<< endl;
  
};
TH1F* GetHisto(TString filename, TString histoname) {
  
  TFile* file  = TFile::Open(filename);
  if (!file) {
    std::cerr << "ERROR: Could not load file" << std::endl
	      << "                        " << filename << std::endl;
    return 0;
  }
  TH1F* h = (TH1F*) file->Get(histoname)->Clone(histoname);
  if (!h) {
    std::cerr << "ERROR[TopSelector]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
  }
  h->SetDirectory(0);
  file->Close();
  
  return h;
};
void SetupDraw(TH1F* h,int color, iVar var, iCut cut){
  h->Rebin(GetRebin(h,var,cut));
  
  TString xaxis = "";
  if (var==MET      ) xaxis = "E_{T}^{miss} (GeV)";
  if (var==InvMass  ) xaxis = "M_{ll} (GeV)"                   ; //2
  if (var==Lep0Pt   ) xaxis = "Leading Lepton p_{T} (GeV)"     ; //3
  if (var==Lep1Pt   ) xaxis = "Leading Lepton p_{T} (GeV)"     ; //4
  if (var==DelLepPhi) xaxis = "#Delta #phi(ll)"                ; //5
  if (var==NJets    ) xaxis = "Number of Jets (E_{T} > 30 GeV)"; //6
  if (var==NBtagJets) xaxis = "Number of b-jets (CSVM)"        ; //7
  if (var==Jet0Pt   ) xaxis = "Leading Jet p_{T} (GeV)"        ; //8

  h->SetTitle(h->GetTitle());
  h->GetXaxis()->SetTitle(xaxis);
  h->SetLineColor(1);
  h->SetFillColor(color);
  h->SetFillStyle(1001);
  
  h->SetBinContent(h->GetNbinsX(),(h->GetBinContent(h->GetNbinsX()+1)+h->GetBinContent(h->GetNbinsX())));
  h->SetBinContent(h->GetNbinsX()+1,0);
  
  return;
}
int GetRebin(TH1F* h, iVar var, iCut cut){
  TString histo = h->GetTitle();
  
  int nbins   = h->GetNbinsX();
  float xmax  = h->GetXaxis()->GetBinUpEdge(nbins);
  float xmin  = h->GetXaxis()->GetBinLowEdge(1);
  

  int rebin = nbins/(xmax-xmin); // normalize to 1GeV?
  if (var==MET      ) return 5 * rebin;
  if (var==InvMass  ) return 5 * rebin;
  if (var==Lep0Pt   ) return 5 * rebin;
  if (var==Lep1Pt   ) return 5 * rebin;
  if (var==DelLepPhi) return 1.;
  if (var==NJets    ) return rebin;
  if (var==NBtagJets) return rebin;
  if (var==Jet0Pt   ) return 5 * rebin;

  
  return rebin;
}
