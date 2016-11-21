#include "Plot.h"
#include "TSystem.h"
#include <fstream> 

/*Plot::Plot(TString v = "MET", TString ch = "ElMu", Bool_t logy = 1, TString sig = "T2tt_mStop183_mLsp1"){ // Constructor
	var = v;
  chan = ch;
  signal = sig;
	SetAllBkg();
	SetData();
  SetSignal();
}*/
TString Plot::getVar(){  return var;}
TString Plot::getChan(){ return chan;}
TString Plot::getSignal(){ return signal;}
void Plot::setVar(TString v){  var = v;}
void Plot::setChan(TString c){  chan = c;}
void Plot::setSignal(TString s){  signal = s;}

Histo* Plot::GetH(TString sample = "TTbar_Powheg"){
  TH1F* f;
  AnalMiniTree *g = new AnalMiniTree(sample, chan);
  g->Loop(var);
  Histo* h = g->histo;
  if( (sample == "MuonEG") || (sample == "DoubleEG") || (sample == "DoubleMuon") || (sample.Contains("Single")) ){ // Data sample
    h->SetType(2);
  }
  else if(sample.BeginsWith("T2tt")){ // Signal sample
    h->Scale(Lumi*1000);
    h->SetType(1);
  }
  else{                              // Background
    h->Scale(Lumi*1000);
    h->SetType(0);
  }
  h->SetStyle();
  h->SetDirectory(0);
  delete g;
  return h;
}

void Plot::SetTTbar(){ 
  httbar = GetH("TTbar_Powheg");
  httbar->SetTag("ttbar");
  httbar->SetColor(kRed+1);
  httbar->SetStyle();
  httbar->SetSyst(0.25);
}
void Plot::SetDY(){
  Histo* h1 = GetH("DYJetsToLL_M10to50_aMCatNLO");
  hDY = GetH("DYJetsToLL_M50_aMCatNLO");
  hDY -> Add(h1);
  hDY->SetTag("DY");
  hDY->SetColor(kAzure-8);
  hDY->SetStyle();
  hDY->SetSyst(0.25);
}
void Plot::SetWJets(){
  hWJets = GetH("WJetsToLNu_aMCatNLO");
  hWJets->SetTag("WJets");
  hWJets->SetColor(kGreen-3);
  hWJets->SetStyle();
  hWJets->SetSyst(0.25);
}
void Plot::SettW(){
  TH1F* h1 = GetH("TbarW");
  htW = GetH("TW"); htW->Add(h1);
  htW->SetTag("tW");
  htW->SetColor(kMagenta);
  htW->SetStyle();
  htW->SetSyst(0.25);
}
void Plot::SetVV(){
  TH1F* h1 = GetH("WW");
  TH1F* h2 = GetH("WZ");
  hVV = GetH("ZZ");
  hVV->Add(h1); hVV->Add(h2);
  hVV->SetTag("VV");
  hVV->SetColor(kYellow-10);
  hVV->SetStyle();
  hVV->SetSyst(0.25);
}
void Plot::SetttV(){
  TH1F* h1 = GetH("TTWToLNu");
  TH1F* h2 = GetH("TTWToQQ");
  TH1F* h3 = GetH("TTZToQQ");
  httV = GetH("TTZToLLNuNu");
  httV->Add(h1); httV->Add(h2); httV->Add(h3); 
  httV->SetTag("ttV");
  httV->SetColor(kOrange-3);
  httV->SetStyle();
  httV->SetSyst(0.25);
}
void Plot::SetData(){
  if(!doData){ hData = AllBkg; return;} 
  hMuon = GetH("DoubleMuon");
  hElec = GetH("DoubleEG");
  hElMu = GetH("MuonEG");
  if(doSingleLep){
    TH1F* h1 = GetH("SingleElectron");
    TH1F* h2 = GetH("SingleMuon");
    hMuon->Add(h2);
    hElec->Add(h1);
    hElMu->Add(h1); hElMu->Add(h2);
  }
    hMuon->SetLineColor(kBlack); hElec->SetLineColor(kBlack); hElMu->SetLineColor(kBlack);
    hMuon->SetMarkerStyle(20); hElec->SetMarkerStyle(20); hElMu->SetMarkerStyle(20);
    hMuon->SetMarkerSize(1.1); hElec->SetMarkerSize(1.1); hElMu->SetMarkerSize(1.1);
  if(chan == "ElMu") hData = hElMu;
  else if(chan == "Elec") hData = hElec;
  else if(chan == "Muon") hData = hMuon;
  else if(chan == "SF" || chan == "SameF" || chan == "sameF"){
    hData = hElec; 
    hData->Add(hMuon);
  }
  else{
    hData = hElMu; hData->Add(hElec); hData->Add(hMuon);
  }
  hData->SetTag("Data");
  hData->SetStyle();
}

void Plot::SetSignal(){
  hSignal = GetH(signal);
  hSignal->SetTag("signal");
  hSignal->SetColor(kGreen);
  hSignal->SetSyst(0.20);
  VSignals.push_back(hSignal);
}

void Plot::AddToHistos(Histo* p){
  if(p->type == 0) VBkgs.push_back(p);
  else if(p->type == 1) VSignals.push_back(p);
  else return;
}

void Plot::SetAllBkg(){
  SetTTbar(); SetDY(); SetWJets(); SetVV(); SettW(); SetttV();
  //AllBkg = (Histo*) httbar->Clone(); AllBkg->Add(hDY);
  //AllBkg->Add(hWJets); AllBkg->Add(hVV); AllBkg->Add(htW); AllBkg->Add(httV); 
  //AllBkg->SetTag("SumBkgs");
  //AllBkg->SetStyle();
  AddToHistos(httV); AddToHistos(hVV); AddToHistos(hWJets);
  AddToHistos(hDY); AddToHistos(htW); AddToHistos(httbar);
	hStack = new THStack(var, "");
  for(int i = 0; i < VBkgs.size(); i++){
    hStack->Add((TH1F*) VBkgs.at(i));
  }
  AllBkg = new Histo(*(TH1F*) hStack->GetStack()->Last(), 3);
}

void Plot::SetLegend(bool doyi = 1){
  leg = new TLegend(0.70,0.65,0.93,0.93);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  Float_t MinYield = AllBkg->yield/5000;
  for(int i = VBkgs.size()-1; i >= 0; i--){
    if(VBkgs.at(i)->yield < MinYield) continue;
    VBkgs.at(i)->AddToLegend(leg,doyi);
  }
	if(doSignal){
  for(int i = VSignals.size()-1; i >= 0; i--){
			VSignals.at(i)->AddToLegend(leg, doyi);
		}
	}
  hData->AddToLegend(leg,doyi);
}

void Plot::SetTexChan(TString cuts){
  TString t = "";
  if (chan == "ElMu") t += "e^{#pm}#mu^{#mp}";
  else if (chan == "Elec") t += "e^{+}e^{-}";
  else if (chan == "Muon") t += "#mu^{+}#mu^{-}";
  else if (chan == "All") t += "#mu^{+}#mu^{-} + e^{+}e^{-} + e^{#pm}#mu^{#mp}";
  else if (chan == "sameF") t += "#mu^{+}#mu^{-} + e^{+}e^{-}";
  t += cuts;
  texchan = new TLatex(-20, 50, t);
  texchan->SetNDC();
  texchan->SetTextAlign(12);
  texchan->SetX(0.15);
  texchan->SetY(0.81);
  texchan->SetTextFont(42);
  texchan->SetTextSize(0.05);
  texchan->SetTextSizePixels(22);
}

void Plot::SetHRatio(){
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
  }
  TString xvar = "";
    if(var == "DiLepPt") xvar = "P_{T}^{ll} [GeV]";
    else if(var == "Jet0Pt") xvar = "P_{T}^{jet0} [GeV]";
    else if(var == "Jet1Pt") xvar = "P_{T}^{jet1} [GeV]";
    else if(var == "Lep0Pt") xvar = "P_{T}^{lep0} [GeV]";
    else if(var == "Lep1Pt") xvar = "P_{T}^{lep1} [GeV]";
    else if(var == "NJets") xvar = "Jet Multiplicity";
    else if(var == "NBtagJets") xvar = "b-jet multiplicity";
    else if(var == "InvMass") xvar = "M_{ll} [GeV]";
    else if(var == "HT") xvar = "HT [GeV]";
    else if(var == "METHT") xvar = "MET/#sqrt{HT}";
    else if(var == "MinDPhiMetJets") xvar = "Min[#Delta#phi_{jet-MET}] [GeV]";
    else if(var == "MT") xvar = "M_{T} [GeV/c^2]";
    else if(var == "MT2") xvar = "M_{T2} [GeV/c^2]";
    else if(var == "MET") xvar = "E_{T}^{miss} [GeV]";
    else if(var == "Ptllb") xvar = "P_{T}^{llb} [GeV]";
    else if(var == "Meff") xvar = "M_{eff} [GeV]";
    else if(var == "CosDelLepPhi") xvar = "Cos(#Delta#phi_{ll})";
    else if(var == "DelLepPhi") xvar = "#Delta#phi_{ll} [rad]";
    else if(var == "DelPhiLepMet") xvar = "#Delta#phi_{l-MET} [rad]";
    else if(var == "DelPhiJetMet") xvar = "#Delta#phi_{jet-MET} [rad]";
    else if(var == "DelPhiLepJet") xvar = "#Delta#phi_{lep-jet} [rad]";
    else if(var == "DelPhiPllbMet") xvar = "#Delta#phi_{P_{T}^{llb}-MET} [rad]";
    else {xvar = var;}
  hratio->GetXaxis()->SetTitle(xvar);
  hratio->SetMinimum(0.8);
  hratio->SetMaximum(1.2);
}

TCanvas* Plot::SetCanvas(){
    TCanvas* c= new TCanvas("c","c",10,10,800,600);
  c->Divide(1,2);

  plot = (TPad*)c->GetPad(1);
  plot->SetPad(0.0, 0.23, 1.0, 1.0);
  plot->SetTopMargin(0.06); plot->SetRightMargin(0.02);

  pratio = (TPad*)c->GetPad(2);
  pratio->SetPad(0.0, 0.0, 1.0, 0.29);
  pratio->SetGridy();// pratio->SetGridx();
  pratio->SetTopMargin(0.03); pratio->SetBottomMargin(0.4); pratio->SetRightMargin(0.02);

  texlumi = new TLatex(-20.,50., Form("%2.1f fb^{-1}, #sqrt{s} = 13 TeV", Lumi));
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
  texcms->SetTextSize(0.052);
  texcms->SetTextSizePixels(23);
  return c;
}

void Plot::DrawStack(bool sav = save){
  TCanvas* c = SetCanvas();
  plot->cd(); 
  if(!doData) hData = AllBkg;
  float maxData = hData->max;
  float maxMC = AllBkg->max;
  float Max = maxMC > maxData? maxMC : maxData;
  if(doSetLogy){
    hStack->SetMaximum(Max*300);
    hStack->SetMinimum(0.1);
    plot->SetLogy();
  }
  else{
    hStack->SetMaximum(Max*1.15);
    hStack->SetMinimum(0);
  }
  hStack->Draw("hist");
  hStack->GetYaxis()->SetTitle("Number of Events");
  hStack->GetYaxis()->SetTitleSize(0.06);
  hStack->GetYaxis()->SetTitleOffset(0.5);
  hStack->GetYaxis()->SetNdivisions(505);
  hStack->GetXaxis()->SetLabelSize(0.0);
  //if(doData) hData->Draw("pesame");
  hData->Draw("pesame");
  if(doSignal){
    hSignal->Draw("lsame");
  }
  SetLegend(doYields);
  leg->Draw("same");
	texcms->Draw("same");
	texlumi->Draw("same");

	pratio->cd();
//	if(doData){
		hratio = (TH1F*)hData->Clone();
		hratio->Divide(AllBkg);
		SetHRatio();
		hratio->Draw("same");
//	}
	//if(doSys) hratioerror->Draw("same,e2");
	if(sav){
		TString dir = plotfolder + "StopPlots/";
		TString plotname = var + "_" + chan;
		gSystem->mkdir(dir, kTRUE);
		c->Print( dir + plotname + ".pdf", "pdf");
		c->Print( dir + plotname + ".png", "png");
		delete c;
	}
}

void Plot::SaveHistograms(TString tag = "0"){
  if(!doSignal){ std::cout << "No datacards without signal!" << std::endl; return;}
  TFile *f;
  f = new TFile(LimitFolder + var + tag + ".root", "recreate");
  httbar->Write(); hDY->Write(); hWJets->Write(); hVV->Write(); htW->Write(); httV->Write();
  hData->Write(); hSignal->Write(); 
  hStack->Write();  
  f->Close(); delete f;
}

void Plot::MakeDatacard(TString tag = "0"){
  if(!doSignal){ std::cout << "No datacards without signal!" << std::endl; return;}
  hData->SetTag("data_obs");
 // if(!doData){ hData = hAllBkg;}
	SaveHistograms(tag); 
  ofstream outputfile;
  TString filename = "datacard_" + var + tag + ".txt";
  outputfile.open(LimitFolder + filename);

  Int_t nChan = 1;
  Int_t nBkgs = VBkgs.size();
  Int_t nSyst = nBkgs + 1;

  outputfile << Form("imax %i\n", nChan);
	outputfile << Form("jmax %i\n", nBkgs);
	outputfile << Form("kmax %i\n", nSyst);
	outputfile << "##-----------\n";
  outputfile << TString("shapes * ") + chan + " " + LimitFolder + var + tag + TString(".root") + " $PROCESS $PROCESS_$SYSTEMATIC \n";
	outputfile << "##-----------\n";
  outputfile << TString("bin ") + chan + "\n";
  outputfile << Form("observation %1.0f \n",hData->yield);
	outputfile << "##-----------\n";
  TString bin      = "bin ";
  TString process1 = "process ";
  TString process2 = "process ";
  TString rate     = "rate    ";
	for(int i = 0; i < nBkgs+1; i++){ 
		if(i<nBkgs){
      bin += chan + " ";
			process1 += VBkgs[i]->process + " ";
			process2 += Form(" %i ", i+1);
			rate += Form(" %1.2f ", VBkgs[i]->yield);
		}
		else{ 
      bin += chan + " ";
			process1 += hSignal->process + " ";
			process2 += Form(" %i ", -1);
			rate += Form(" %1.2f ", hSignal->yield);
		}    
	}
	outputfile << bin      + "\n";
	outputfile << process1 + "\n";
	outputfile << process2 + "\n";
	outputfile << rate     + "\n";

		outputfile << "##-----------\n";
	TString out = "";

	out = "Lumi  lnN  ";
	for(int i = 0; i < nBkgs+1; i++){ out += Form(" %1.2f ", 1+sys_lumi);}

	for(int i = 0; i < nBkgs+1; i++){
		if(i<nBkgs) out = VBkgs[i]->process + " lnN ";
		else        out = hSignal ->process + " lnN ";
		for(int j = 0; j < nBkgs+1; j++){
			if(j != i) out += " - ";
			else{
				if(i<nBkgs) out += Form(" %1.2f ", 1+VBkgs[i]->syst);
				else        out += Form(" %1.2f ", 1+hSignal->syst); 
			}
		}
		out += "\n";
		outputfile << out;
	}
	outputfile.close();
}



