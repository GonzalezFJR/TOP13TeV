#include <iomanip>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLine.h"
#include "TChain.h"
#include "TLatex.h"
#include "stdio.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <fstream>
#include <string>

using namespace std;


const TString path = "/mnt_pool/fanae105/user/juanr/TOP13TeV/temp/";
//const TString path = "/mnt_pool/fanae105/user/juanr/TOP13TeV/Trees5TeV/may11/";

const Int_t nVars=23;
const Int_t nLevels=6;
const Int_t nChannels = 5;
const Int_t nSamples = 7;
const Int_t nProcesses = 5+1;
const Int_t nSyst = 10; //12;
const Int_t nTSyst = 6;
const Int_t nBkgs = 12;
//const float Lumi = 1264; // pb^-1
const float Lumi = 25.8; // pb^-1
const TString syst[nSyst] = {"BtagUp", "BtagDown", "MisTagUp", "MisTagDown", "LESUp", "LESDown", "PUUp", "PUDown", "JER", "TopPt"}; // "JESUp", "JESDown", 
const TString sample[nSamples] = {
  "TTJets",
  "TW","TbarW",
  "DYJetsToLL_M50_aMCatNLO",
  "WW","WZ",
  //"WJetsToLNu_aMCatNLO",
  "Data_SingleMu"
};
//const TString bkg[nBkgs] = { "TTJets", "TW","TbarW", "DYJetsToLL_M10to50_aMCatNLO", "DYJetsToLL_M50_aMCatNLO", "WW","WZ","ZZ", "TTWToLNu", "TTWToQQ", "TTZToQQ", "WJetsToLNu_aMCatNLO"};
const TString bkg[nBkgs] = { "TTJets", "TW","TbarW", "DYJetsToLL_M5to50_MLM", "DYJetsToLL_M50_MLM", "WW","WZ","ZZ", "TTWToLNu", "TTWToQQ", "TTZToQQ", "WJetsToLNu_aMCatNLO"};

const TString level[nLevels] = {"dilepton","ZVeto","MET","2jets","1btag","DYVeto"};
const TString process[nProcesses] = {"ttbar", "tW", "DY", "VV", "NonW", "data"}; 
const TString chan[nChannels] = { "ElMu", "Elec", "Muon", "All", "sameF"};
const TString var[nVars] = {
  "DiLepPt", "Lep0Pt", "Jet0Pt", "Jet1Pt",  "Lep1Pt", // 0, 1, 2, 3, 4
  "NJets", "NBtagJets", "NBtagsNJets",                // 5, 6 , 7
  "InvMass", "MT2", "MT2bb", "MT2lblb",               // 8, 9, 10, 11
  "MET", "Ptllb", "Meff", "HT",                            // 12, 13, 14, 15
  "DelLepPhi", "DelPhiLepMet", "DelPhiJetMet", "DelPhiLepJet", "DelPhiPllbMet", // 16, 17, 18, 19, 20
  "METHT", "MinDPhiMetJets"}; // 21, 22  
const int nUsefulvars = 9;
const TString usefulvars[nUsefulvars] = {"DelLepPhi", "DiLepPt", "Lep0Pt", "Jet0Pt", "Jet1Pt",  "Lep1Pt", "NJets", "InvMass", "HT"}; // MET NBtagJets NBtagsNJets

TPad* plot; TPad* pratio;
TLegend* leg;
TLatex* texlumi;
TLatex* texcms;
TLatex* texchan;
TH1F* hratio;

TCanvas *SetCanvas();
void SetLegend();
void SetTexChan(TString chan, TString level);
void SetHRatio(TString var);

TCanvas *SetCanvas(){
  TCanvas* c= new TCanvas("c","c",10,10,800,600);
  c->Divide(1,2);

  plot = (TPad*)c->GetPad(1);
  plot->SetPad(0.0, 0.23, 1.0, 1.0);
  plot->SetTopMargin(0.06); plot->SetRightMargin(0.02);

  pratio = (TPad*)c->GetPad(2);
  pratio->SetPad(0.0, 0.0, 1.0, 0.29);
  pratio->SetGridy();// pratio->SetGridx();
  pratio->SetTopMargin(0.03); pratio->SetBottomMargin(0.4); pratio->SetRightMargin(0.02);

  texlumi = new TLatex(-20.,50., Form("%4.0f pb^{-1}, #sqrt{s} = 5 TeV", Lumi));  
  texlumi->SetNDC();
  texlumi->SetTextAlign(12);
  texlumi->SetX(0.72);
  texlumi->SetY(0.97);
  texlumi->SetTextFont(42);
  texlumi->SetTextSize(0.045);
  texlumi->SetTextSizePixels(22);

  texcms = new TLatex(0.,0., "CMS Preliminary");
  texcms->SetNDC();
  texcms->SetTextAlign(12);
  texcms->SetX(0.15);
  texcms->SetY(0.9);
//  texcms->SetTextFont(61);
  texcms->SetTextSize(0.052);
  texcms->SetTextSizePixels(23);
  return c;
}

void SetLegend(){
  leg = new TLegend(0.80,0.70,0.93,0.93);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
}

void SetTexChan(TString chan, TString level){
  TString t = "";
  if (chan == "ElMu") t += "e^{#pm}#mu^{#mp}";
  else if (chan == "Elec") t += "e^{+}e^{-}";
  else if (chan == "Muon") t += "#mu^{+}#mu^{-}";
  else if (chan == "All") t += "#mu^{+}#mu^{-} + e^{+}e^{-} + e^{#pm}#mu^{#mp}";
  else if (chan == "sameF") t += "#mu^{+}#mu^{-} + e^{+}e^{-}";
  if (level == "MET")   t += ", MET > 80";
	if (level == "2jets") t += ", #geq 2 jets";
	if (level == "1btag") t += ", #geq 2 jets, #geq 1 b-jet";
	if (level == "DYVeto"){
    TString tt = TString(t);
    tt += ", MET > 80, #geq 2 jets, #geq 1 b-jet";
		t = "#splitline{" + tt + "}{+ DYVeto}";
	}
  texchan = new TLatex(-20, 50, t);
  texchan->SetNDC();
  texchan->SetTextAlign(12);
  texchan->SetX(0.15);
  texchan->SetY(0.84);
  texchan->SetTextFont(42);
  texchan->SetTextSize(0.05);
  texchan->SetTextSizePixels(22);
}

void SetHRatio(TString var){
  hratio->SetTitle("");
  hratio->GetYaxis()->SetTitle("Data/MC");
  hratio->GetXaxis()->SetTitleSize(0.05);
  hratio->GetYaxis()->CenterTitle();
  hratio->GetYaxis()->SetTitleOffset(0.25);
  hratio->GetYaxis()->SetTitleSize(0.1);
  hratio->GetYaxis()->SetLabelSize(0.1);
  hratio->GetYaxis()->SetNdivisions(402);
  hratio->GetXaxis()->SetTitleOffset(0.9);
  hratio->GetXaxis()->SetLabelSize(0.13);
  hratio->GetXaxis()->SetTitleSize(0.16);
  if (var == "NBtagsNJets") {  //change bin labels
    hratio->GetXaxis()->SetBinLabel( 1, "(0, 0)");
    hratio->GetXaxis()->SetBinLabel( 2, "(1, 0)");
    hratio->GetXaxis()->SetBinLabel( 3, "(1, 1)");
    hratio->GetXaxis()->SetBinLabel( 4, "(2, 0)");
    hratio->GetXaxis()->SetBinLabel( 5, "(2, 1)");
    hratio->GetXaxis()->SetBinLabel( 6, "(2, 2)");
    hratio->GetXaxis()->SetBinLabel( 7, "(3, 0)");
    hratio->GetXaxis()->SetBinLabel( 8, "(3, 1)");
    hratio->GetXaxis()->SetBinLabel( 9, "(3, 2)");
    hratio->GetXaxis()->SetBinLabel(10, "(3, 3)");
    hratio->GetXaxis()->SetBinLabel(11, "(4, 0)");
    hratio->GetXaxis()->SetBinLabel(12, "(4, 1)");
    hratio->GetXaxis()->SetBinLabel(13, "(4, 2)");
    hratio->GetXaxis()->SetBinLabel(14, "(4, 3)");
    hratio->GetXaxis()->SetBinLabel(15, "(4, 4)");
    hratio->GetXaxis()->SetLabelSize(0.14);
    hratio->GetXaxis()->SetLabelOffset(0.02);
  } else  if (var == "Yields") {  //change bin labels
    hratio->GetXaxis()->SetBinLabel( 1, "Dilepton");
    hratio->GetXaxis()->SetBinLabel( 2, "Z veto");
    hratio->GetXaxis()->SetBinLabel( 3, "MET");
    hratio->GetXaxis()->SetBinLabel( 4, "#geq2 jets");
    hratio->GetXaxis()->SetBinLabel( 5, "#geq1 b tags");
    hratio->GetXaxis()->SetBinLabel( 6, "=1 b tag");
    hratio->GetXaxis()->SetBinLabel( 7, "=2 b tags");
    hratio->GetXaxis()->SetLabelSize(0.19);
    hratio->GetXaxis()->SetLabelOffset(0.04);
  } 
  TString xvar = "";
    if(var == "DiLepPt") xvar = "P_{T}^{ll}";
    else if(var == "Jet0Pt") xvar = "P_{T}^{jet0}"; 
    else if(var == "Jet1Pt") xvar = "P_{T}^{jet1}"; 
    else if(var == "Lep0Pt") xvar = "P_{T}^{lep0}"; 
    else if(var == "Lep1Pt") xvar = "P_{T}^{lep1}"; 
    else if(var == "NJets") xvar = "Jet Multiplicity"; 
    else if(var == "NBtagJets") xvar = "b-jet multiplicity"; 
    else if(var == "InvMass") xvar = "M_{ll}"; 
    else if(var == "HT") xvar = "HT"; 
    else if(var == "METHT") xvar = "MET/#sqrt{HT}"; 
    else if(var == "MinDPhiMetJets") xvar = "Min[#Delta#phi_{jet-MET}]"; 
    else if(var == "MT") xvar = "M_{T}"; 
    else if(var == "MT2") xvar = "M_{T2}"; 
    else if(var == "MET") xvar = "E_{T}^{miss}"; 
    else if(var == "Ptllb") xvar = "P_{T}^{llb}"; 
    else if(var == "Meff") xvar = "M_{eff}"; 
    else if(var == "DelLepPhi") xvar = "#Delta#phi_{ll}"; 
    else if(var == "DelLepPhi") xvar = "#Delta#phi_{ll}"; 
    else if(var == "DelPhiLepMet") xvar = "#Delta#phi_{l-MET}"; 
    else if(var == "DelPhiJetMet") xvar = "#Delta#phi_{jet-MET}"; 
    else if(var == "DelPhiLepJet") xvar = "#Delta#phi_{lep-jet}"; 
    else if(var == "DelPhiPllbMet") xvar = "#Delta#phi_{P_{T}^{llb}-MET}"; 
    else {xvar = var;}
  hratio->GetXaxis()->SetTitle(xvar);
  hratio->SetMinimum(0.2);
  hratio->SetMaximum(1.8);
}
