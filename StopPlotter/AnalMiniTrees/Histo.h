#ifndef Histo_h
#define Histo_h

#include <TROOT.h>
#include <TChain.h>
#include <THStack.h>
#include <TH1F.h>
#include <TFile.h>
#include <TMath.h>
#include <TLegend.h>
#include <iostream>

class Histo : public TH1F{
public:
  Int_t type; // Type = 0: background, Type =1: signal, Type = 2: data. 
  Int_t color;
  TString tag = "0"; TString process="0"; TString cuts="0"; TString xlabel="0";
  Float_t syst;
  Double_t yield;
  Double_t max;
  Bool_t IsStackOverflow = true;

  Histo(const char *name, const char *title, Int_t nbins, Double_t xlow, Double_t xup);
  Histo(const char *name, const char *title, Int_t nbins, const Float_t* xbins);
  Histo(const TH1F &h, Int_t tipo = 0, Int_t c = 1){
    ((Histo&)h).Copy(*this);
    cuts = "cut"; xlabel = "[GeV]"; syst = 0.25;
    SetType(tipo);
    SetColor(c);
    SetStyle();
  }
  ~Histo(){};

  void SetType(Int_t tipo = 0);
  void SetColor(Int_t c);
  void SetStyle();

  void StackOverflow(Bool_t doStackOverflow = 1);
  void SetTag(TString p, TString t="", TString x = "", TString c = "");
  void SetTitles(TString x, TString c = "");
  void SetSyst(Float_t s);
  void SetParams(Int_t tipo = 0, Int_t color = 1, Float_t s = 0.25, TString t = "", TString x = "", TString c = "");

  void AddToLegend(TLegend* leg, Bool_t doyi = 1);
  TH1F* GetVarHistoStatBin(Int_t bin = 0, TString dir = "Up");
};

Histo::Histo(const char *name, const char *title, Int_t nbins, Double_t xlow, Double_t xup)
  : TH1F(name, title, nbins, xlow, xup){
  SetType(0);
  SetStyle();
  SetColor(1);
  tag = TString(name);
  process = TString(title);
}

Histo::Histo(const char *name, const char *title, Int_t nbins, const Float_t* xbins)
 : TH1F(name,title,nbins,xbins){
  SetType(0);
  SetStyle();
  SetColor(1);
  tag = TString(name);
  process = TString(title);
}

#endif
