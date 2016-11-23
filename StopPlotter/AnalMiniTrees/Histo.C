#include "Histo.h"

void Histo::SetType(Int_t tipo){
	type = tipo;
	if(type < 0 || type > 3) type = 3;
	if(type == 0){         // Background
		SetLineStyle(0);
		SetFillStyle(1001);
	}
	else if(type == 1){    // Signal
		SetLineStyle(1);
    SetLineWidth(2);
	}
	else if(type == 2){    // Data
		SetLineStyle(0);
    SetFillStyle(0);
    SetColor(kBlack);
    SetMarkerStyle(20);
    SetMarkerSize(1.1);
	}
  else if(type == 3){   // Other
    SetLineStyle(1);
    SetLineWidth(2);
  }
}

void Histo::StackOverflow(Bool_t doStackOverflow){
	if(!doStackOverflow) return;
	int lastBin = GetNbinsX();
	float lastBinContent = GetBinContent(lastBin);
	float overflow = GetBinContent(lastBin+1);
	SetBinContent(lastBin, lastBinContent + overflow);
  SetBinContent(lastBin+1, 0);
}

void Histo::SetStyle(){
	SetStats(0);
	StackOverflow();
	yield = Integral();
  max = GetMaximum();
  GetXaxis()->SetTitle(xlabel);
}

void Histo::SetTag(TString p, TString t, TString x, TString c){
  if(p != ""){process = p; SetName(p);}
	if(t != "") tag = t;
	if(x != "") xlabel = x;
	if(c != "") cuts = x;
}

void Histo::SetTitles(TString x, TString c){
	xlabel = x; 
	if(cuts != "") cuts = c;
}

void Histo::SetColor(Int_t c){
	color = c;
	SetLineColor(c); SetFillColor(c);
}

void Histo::SetSyst(Float_t s){
  syst = s;
}

void Histo::AddToLegend(TLegend* leg, Bool_t doyi){
  TH1F* h = (TH1F*) Clone();
  TString op = "f";
  if(type == 1) op = "l";
  if(type == 2) op = "pe";
  if(doyi) leg->AddEntry(h, Form(process + ": %1.0f", yield), op);
  else leg->AddEntry(h, process, op);
}

TH1F* Histo::GetVarHistoStatBin(Int_t bin, TString dir){
  Float_t var = GetBinContent(bin);
  Float_t stat = GetBinError(bin);
  TH1F* h = (TH1F*) Clone();
  if      (dir == "up" || dir == "Up" || dir == "UP")  h->SetBinContent(bin, var + stat);
  else if (dir == "down" || dir == "Down" || dir == "DOWN")  h->SetBinContent(bin, var + stat);
  else    cout << " ---> ERROR!!!! No valid direction: " << dir << endl;
  return h;
}
