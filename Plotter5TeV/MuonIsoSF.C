#include "SetPlotter.C"

TFile *f1 = TFile::Open(path + "Tree_" + "DYJetsToLL_M50_aMCatNLO" + ".root");
TFile *f2 = TFile::Open(path + "Tree_" + "Data_SingleMu" + ".root");

const TString v = "Pt";

void MuonIsoSF(){
	TH1F* hDY_PtNoIso; TH1F* hData_PtNoIso; 
	TH1F* hDY_PtIso; TH1F* hData_PtIso; 
  TString hname = TString("fHMuon") + v;
	f1->GetObject(hname + "noIso", hDY_PtNoIso);
	f1->GetObject(hname + "Iso"  , hDY_PtIso);
	f2->GetObject(hname + "noIso", hData_PtNoIso);
	f2->GetObject(hname + "Iso"  , hData_PtIso);

	TCanvas *c = new TCanvas("c","c",10,10,800,600);;
	c->Divide(1,2);
	hDY_PtNoIso->SetLineColor(kRed);	  hDY_PtNoIso->SetStats(0);
	hDY_PtIso->SetLineColor(kBlue);
	hData_PtNoIso->SetLineColor(kRed);	hData_PtNoIso->SetStats(0);
	hData_PtIso->SetLineColor(kBlue);
	hData_PtIso->SetLineWidth(2);
	hData_PtNoIso->SetLineWidth(2);
	hDY_PtIso->SetLineWidth(2);
	hDY_PtNoIso->SetLineWidth(2);

	TLegend* LDY = new TLegend(0.60,0.75,0.93,0.93); TLegend* LData  = new TLegend(0.60,0.75,0.93,0.93);

  float nMCnoIso = hDY_PtNoIso->GetEntries();
  float nMCIso   = hDY_PtIso  ->GetEntries();
  float nDatanoIso = hData_PtNoIso->GetEntries();
  float nDataIso   = hData_PtIso  ->GetEntries();
  float fData = nDataIso/nDatanoIso;
  float fMC   = nMCIso/nMCIso;  

  cout << Form(" - Ratio for MC  : %1.3f", fMC  ) << endl;
  cout << Form(" - Ratio for Data: %1.3f", fData) << endl;
  cout << Form("   Diference:      %1.3f", fabs(fMC-fData)) << endl; 


	c->cd(1);
  hDY_PtNoIso->SetTitle("Comparison for MC");
	hDY_PtNoIso->Draw("le");
	hDY_PtIso->Draw("lesame");
	LDY->AddEntry(hDY_PtNoIso, Form("Not Iso MC Muons: %.0f", nMCnoIso), "le");
	LDY->AddEntry(hDY_PtIso  , Form("    Iso MC Muons: %.0f", nMCIso), "le");
	LDY->Draw();

	c->cd(2);
  hData_PtNoIso->SetTitle("Comparison for Data");
	hData_PtNoIso->Draw("le");
	hData_PtIso->Draw("lesame");

	LData->AddEntry(hData_PtNoIso, Form("Not Iso Data Muons: %.0f", nDatanoIso), "le");
	LData->AddEntry(hData_PtIso  , Form("    Iso Data Muons: %.0f", nDataIso), "le");
	LData->Draw();
}
