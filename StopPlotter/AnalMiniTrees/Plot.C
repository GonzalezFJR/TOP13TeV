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
  AnalMiniTree *g = new AnalMiniTree(sample, chan, sys);
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
	if(sys == "0" || sys == "nom"){
		httbar = GetH("TTbar_Powheg");
		httbar->SetTag("ttbar");
		httbar->SetColor(kRed+1);
		httbar->SetStyle();
		httbar->SetSyst(0.15);
	  AddToHistos(httbar);
  }
	else{
		Histo* h = GetH("TTbar_Powheg");
		h->SetTag("ttbar");
		h->SetStyle();
		h->SetTag("ttbar_" + sys);
		h->type = -1;
		AddToHistos(h);
	}
}
void Plot::SetDY(){
	if(sys == "0" || sys == "nom"){
		Histo* h1 = GetH("DYJetsToLL_M10to50_aMCatNLO");
		hDY = GetH("DYJetsToLL_M50_aMCatNLO");
		hDY -> Add(h1);
		hDY->SetTag("DY");
		hDY->SetColor(kAzure-8);
		hDY->SetStyle();
		hDY->SetSyst(0.25);
	  AddToHistos(hDY);
  } 
	else{
		Histo* h  = GetH("DYJetsToLL_M50_aMCatNLO");
		Histo* h1 = GetH("DYJetsToLL_M10to50_aMCatNLO"); h->Add(h1);
		h->SetStyle(); h->SetTag("DY_" + sys); h->type = -1;
		AddToHistos(h);
	}
}
void Plot::SetWJets(){
	if(sys == "0" || sys == "nom"){
		hWJets = GetH("WJetsToLNu_aMCatNLO");
		hWJets->SetTag("WJets");
		hWJets->SetColor(kGreen-3);
		hWJets->SetStyle();
  	hWJets->SetSyst(0.25);
	  AddToHistos(hWJets);
  }
	else{
		Histo* h = GetH("WJetsToLNu_aMCatNLO");
		h->SetStyle(); h->SetTag("WJets_" + sys); h->type = -1;
		AddToHistos(h);
	}
}
void Plot::SettW(){
		Histo* h1 = GetH("TbarW");
	if(sys == "0" || sys == "nom"){
		htW = GetH("TW"); htW->Add(h1);
		htW->SetTag("tW");
		htW->SetColor(kMagenta);
		htW->SetStyle();
		htW->SetSyst(0.25);
    AddToHistos(htW); 
	}
	else{
		Histo* h = GetH("TW"); h->Add(h1);
		h->SetStyle(); h->SetTag("tW_" + sys); h->type = -1;
		AddToHistos(h);
	}
}
void Plot::SetVV(){
	Histo* h1 = GetH("WW");
	Histo* h2 = GetH("WZ");
	if(sys == "0" || sys == "nom"){
		hVV = GetH("ZZ");
		hVV->Add(h1); hVV->Add(h2);
		hVV->SetTag("VV");
		hVV->SetColor(kYellow-10);
		hVV->SetStyle();
		hVV->SetSyst(0.25);
    AddToHistos(hVV); 
	}
	else{
		Histo *h = GetH("ZZ"); h->Add(h1); h->Add(h2);
		h->SetStyle();
		h->SetTag("VV_" + sys);
		h->type = -1;
		AddToHistos(h);
	}
}
void Plot::SetttV(){
	Histo* h1 = GetH("TTWToLNu");
	Histo* h2 = GetH("TTWToQQ");
	Histo* h3 = GetH("TTZToQQ");
	if(sys == "0" || sys == "nom"){
		httV = GetH("TTZToLLNuNu");
		httV->Add(h1); httV->Add(h2); httV->Add(h3); 
		httV->SetTag("ttV");
		httV->SetColor(kOrange-3);
		httV->SetStyle();
		httV->SetSyst(0.25);
		AddToHistos(httV);
	}
	else{
		Histo* h = GetH("TTZToLLNuNu");
		h->Add(h1); h->Add(h2); h->Add(h3); 
		h->SetStyle(); h->SetTag("ttV_" + sys); h->type = -1;
		AddToHistos(h);
	}
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
  hData->type = 2;
}

void Plot::SetSignal(){
  //TString t = signal;
  //t.ReplaceAll("T2tt_mStop", ""); t.ReplaceAll("mLsp", "");
  if(sys == "0" || sys == "nom"){
		hSignal = GetH(signal);
		hSignal->SetTag("Signal");
		hSignal->type = 1;
		hSignal->SetColor(kCyan);
		hSignal->SetSyst(0.20);
		VSignals.push_back(hSignal);
 }
	else{
    Histo* ts = GetH(signal); 
		ts->SetTag("Signal_" + sys);
		ts->type = -3;
		VSignals.push_back(ts);
  }
}

//void Plot::SetOtherSignal(TString signal2 = "T2tt_mStop250_mLsp75", Int_t ccol = kGreen+4){
void Plot::SetOtherSignal(TString signal2, Int_t ccol){
  TString t2 = signal2;
  t2.ReplaceAll("T2tt_mStop", ""); t2.ReplaceAll("mLsp", ""); 
  Histo* sig2 = GetH(signal2);
  sig2->SetTag(t2);
  sig2->SetColor(ccol);
  VSignals.push_back(sig2);
}

void Plot::SetAllProcesses(){
	SetttV(); SettW(); SetVV(); SetWJets(); SetDY(); SetTTbar();
//	AddToHistos(httV); AddToHistos(hVV); AddToHistos(hWJets);
//	AddToHistos(hDY); AddToHistos(htW); AddToHistos(httbar);
}

void Plot::SetAllBkg(){
	if(sys == "0" && nBkgs == 0){
		hStack = new THStack(var, "");
		nBkgs = VBkgs.size();
		for(int i = 0; i < nBkgs; i++){
			hStack->Add((TH1F*) VBkgs.at(i));
		}
		AllBkg = new Histo(*(TH1F*) hStack->GetStack()->Last(), 3);
    AllBkg->SetTag("Total Bkg");
	}
}

void Plot::AddToHistos(Histo* p){
  if(p->type == 0 ) VBkgs.push_back(p);
  else if(p->type == 1) VSignals.push_back(p);
  else if(p->type < 0) VSyst.push_back(p);
  else return;
}


void Plot::SetLegend(bool doyi = 1){
  leg = new TLegend(0.70,0.65,0.93,0.93);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  Float_t MinYield = AllBkg->yield/5000;
  int nVBkgs = VBkgs.size();
  for(int i = nVBkgs-1; i >= 0; i--){
    if(VBkgs.at(i)->yield < MinYield) continue;
    else VBkgs.at(i)->AddToLegend(leg,doyi);
  }
	if(doSignal){
  for(int i = VSignals.size()-1; i >= 0; i--){
      if(VSignals.at(i)->type < 0) continue;
			VSignals[i]->AddToLegend(leg, doyi);
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

void Plot::DrawStack(TString tag = "0", bool sav = save){
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

  // Draw systematics histo
  AllBkg->SetFillStyle(3145);
  AllBkg->SetFillColor(kGray+2);
  AllBkg->SetLineColor(kGray+2);
  AllBkg->SetLineWidth(0);
  AllBkg->SetMarkerSize(0);
  AllBkg->Draw("same,e2");

  hData->Draw("pesame");
  if(doSignal){
		// Draw systematics signal
		TH1F* hSignalerr = (TH1F*) hSignal->Clone();
		hSignalerr->SetFillColor(kGreen-10); 
		hSignalerr->SetMarkerStyle(0);
		hSignalerr->SetFillStyle(3145);
    for(unsigned int  i = 0; i < VSignals.size(); i++){
      if(VSignals.at(i)->type < 0) continue; 
      VSignals.at(i)->Draw("lsame");
    }
    hSignalerr->Draw("same,e2");
  }

  // Draw systematics ratio
  AllBkg->SetFillStyle(3145);
  TH1F* hratioerr =  (TH1F*) AllBkg->Clone();
  Int_t nbins = AllBkg->GetNbinsX(); Float_t binval = 0; Float_t errbin = 0; Float_t totalerror = 0;
  for(int bin = 1; bin <= nbins; bin++){
    totalerror = AllBkg->GetBinError(bin); 
    binval = AllBkg->GetBinContent(bin);
    errbin = binval > 0 ? totalerror/binval : 0.0;
    hratioerr->SetBinContent(bin, 1);
    hratioerr->SetBinError(bin, errbin);
  }
  hratioerr->SetFillColor(kTeal-7);
  hratioerr->SetFillStyle(3144);
  hratioerr->SetMarkerSize(0);


  SetLegend(doYields);
  leg->Draw("same");
	texcms->Draw("same");
	texlumi->Draw("same");

	pratio->cd();
	hratio = (TH1F*)hData->Clone();
	hratio->Divide(AllBkg);
	SetHRatio();
	hratio->Draw("same");
	hratioerr->Draw("same,e2");
	hratio->Draw("same");

	if(sav){
		TString dir = plotfolder + "StopPlots/";
		TString plotname = var + "_" + chan + "_" + tag;
		gSystem->mkdir(dir, kTRUE);
		c->Print( dir + plotname + ".pdf", "pdf");
		c->Print( dir + plotname + ".png", "png");
		delete c;
	}
}

void Plot::SaveHistograms(TString tag = "0"){
	if(!doSignal){ std::cout << "No datacards without signal!" << std::endl; return;}
	TFile *f;
  TString filename =  var + "_" + chan + "_" + tag + ".root";
	f = new TFile(LimitFolder + filename, "recreate");

	TH1F* statup; TH1F* statdown; Histo* nom;
	int nVBkgs = VBkgs.size(); 
	int nSig   = VSignals.size();
	int nbins;
	for(int i = 0; i < nVBkgs; i++){
		nom = VBkgs.at(i);
		nbins = nom->GetNbinsX();
		nom->Write();
		if(nom->type < 0) continue; // no stat for systematics
		for(int j = 1; j <= nbins; j++){
			statup   = nom->GetVarHistoStatBin(j, "up");
			statdown = nom->GetVarHistoStatBin(j, "down");
			statup  ->SetName(nom->process + "_" + nom->process + "_" + chan + Form("_statbin%i", j) + "Up");
			statdown->SetName(nom->process + "_" + nom->process + "_" + chan + Form("_statbin%i", j) + "Down");
			statup  ->Write(); statdown->Write();
			//delete statup; delete statdown;
		}
		//delete nom;
	}
	for(int i = 0; i < nSig; i++){
		nom = VSignals.at(i);
		nbins = nom->GetNbinsX();
		nom->Write();
    //cout << "Saving " << nom->tag << "..." << endl;
		if(nom->type < 0) continue; // no stat for systematics
		for(int j = 1; j <= nbins; j++){
			statup   = nom->GetVarHistoStatBin(j, "up");
			statdown = nom->GetVarHistoStatBin(j, "down");
			statup  ->SetName(nom->process + "_" + nom->process + "_" + chan + Form("_statbin%i", j) + "Up");
			statdown->SetName(nom->process + "_" + nom->process + "_" + chan + Form("_statbin%i", j) + "Down");
			statup  ->Write(); statdown->Write();
			//delete statup; delete statdown;
		}
		//delete nom;
	}
  for(int i = 0; i < (Int_t) VSyst.size(); i++){
    nom = VSyst.at(i);
		nom->Write();
  }
	hData->Write();
	hStack->Write();  
	cout << "-------> Root file created: " << LimitFolder + filename << endl;
	f->Close(); delete f;
}

TString Plot::GetStatUncDatacard(){
	Histo* nom;
	int nbins = 0;
	TString lin = TString("");
	for(int i = 0; i < nBkgs+1; i++){
		if(i<nBkgs) nom = VBkgs.at(i);
		else        nom = hSignal;
		nbins = nom->GetNbinsX();
		for(int j = 1; j <= nbins; j++){
			lin += nom->process + "_" + chan + Form("_statbin%i shape ", j);
			for(int k = 0; k < nBkgs+1; k++){
				if(i==k) lin += TString(" 1 ");
				else lin += TString(" - ");
			}
			lin += TString("\n");
		}
	}
  return lin;
}

void Plot::MakeDatacardBin(Int_t bin, TString tag = "b"){
  if(!doSignal){ std::cout << "No datacards without signal!" << std::endl; return;}
	hData->SetTag("data_obs");

	ofstream outputfile;
	TString filename = TString("datacard_") + var + TString("_") + chan + "_" + tag + Form("%i",bin) + TString(".txt");
	outputfile.open(LimitFolder + filename);

	Int_t nChan = 1;
	Int_t nSyst = nBkgs + 1;

  outputfile << Form("imax %i\n", nChan);
	outputfile << "jmax *\n";
	outputfile << "kmax *\n";
  //outputfile << Form("jmax *", nBkgs);
	//outputfile << Form("kmax %i\n", nSyst);
	outputfile << "##-----------\n";
  outputfile << TString("bin ") + chan + "\n";
  outputfile << Form("observation %1.0f \n",hData->GetBinContent(bin));
	outputfile << "##-----------\n";
  TString bint     = TString("bin ");
  TString process1 = TString("process ");
	TString process2 = TString("process ");
	TString rate     = TString("rate    ");
	for(int i = 0; i < nBkgs; i++){ 
		bint += chan + TString(" ");
		process1 += VBkgs[i]->process + TString(" ");
		process2 += Form(" %i ", i+1);
		rate += Form(" %1.2f ", VBkgs[i]->GetBinContent(bin));
	}
	bint += chan + TString(" ");
	process1 += hSignal->process + TString(" ");
	process2 += Form(" %i ", -1);
	rate += Form(" %1.2f ", hSignal->GetBinContent(bin));

	outputfile << bint      + TString("\n");
	outputfile << process1 + TString("\n");
	outputfile << process2 + TString("\n");
	outputfile << rate     + TString("\n");
		outputfile << "##-----------\n";
	TString out = TString("Lumi  lnN  ");
	for(int i = 0; i < nBkgs+1; i++){ out += Form(" %1.2f ", 1+sys_lumi);}

	for(int i = 0; i < nBkgs+1; i++){
		if(i<nBkgs) 	out = VBkgs[i]->process + " lnN ";
		else	out = hSignal ->process + " lnN ";
		for(int j = 0; j < nBkgs+1; j++){
			if(j != i) out += TString(" - ");
			else{
				if(i<nBkgs) out += Form(" %1.2f ", 1+VBkgs[i]->syst);
				else        out += Form(" %1.2f ", 1+hSignal->syst); 
			}
		}
		out += TString("\n");
		outputfile << out;
	}
	outputfile.close();
  cout << "-------> Datacard created: " << LimitFolder + filename << endl;
}

void Plot::MakeDatacardAllBins(TString tag = "b"){
  Int_t nbins = hSignal->GetNbinsX();
  for(int i = 1; i <= nbins; i++){
    MakeDatacardBin(i, tag);
  } 
}

void Plot::MakeDatacard(TString tag = "0"){
	if(!doSignal){ std::cout << "No datacards without signal!" << std::endl; return;}
	hData->SetTag("data_obs");
	// if(!doData){ hData = hAllBkg;}
	SaveHistograms(tag); 
	ofstream outputfile;
	TString filename = TString("datacard_") + var + TString("_") + chan + "_" + tag + TString(".txt");
	outputfile.open(LimitFolder + filename);

	Int_t nChan = 1;
	Int_t nSyst = nBkgs + 1;

  outputfile << Form("imax %i\n", nChan);
	outputfile << "jmax *\n";
	//outputfile << Form("kmax %i\n", nSyst);
	outputfile << "kmax *\n";
	outputfile << "##-----------\n";
  outputfile << TString("shapes * ") + chan + " " + LimitFolder + var + "_" + chan + "_" + tag + TString(".root") + TString(" $PROCESS $PROCESS_$SYSTEMATIC \n");
	outputfile << "##-----------\n";
  outputfile << TString("bin ") + chan + "\n";
  outputfile << Form("observation %1.0f \n",hData->yield);
	outputfile << "##-----------\n";
  TString bin      = TString("bin ");
  TString process1 = TString("process ");
  TString process2 = TString("process ");
  TString rate     = TString("rate    ");
	for(int i = 0; i < nBkgs; i++){ 
		bin += chan + TString(" ");
		process1 += VBkgs[i]->process + TString(" ");
		process2 += Form(" %i ", i+1);
		rate += Form(" %1.2f ", VBkgs[i]->yield);
	}
	// Add signal
	bin += chan + TString(" ");
	process1 += hSignal->process + TString(" ");
	process2 += Form(" %i ", -1);
	rate += Form(" %1.2f ", hSignal->yield);

	outputfile << bin      + TString("\n");
	outputfile << process1 + TString("\n");
	outputfile << process2 + TString("\n");
	outputfile << rate     + TString("\n");
		outputfile << "##-----------\n";
	TString out = TString("Lumi  lnN  ");
	for(int i = 0; i < nBkgs+1; i++){ out += Form(" %1.2f ", 1+sys_lumi);}

	for(int i = 0; i < nBkgs+1; i++){
		if(i<nBkgs){
			out = VBkgs[i]->process + " lnN ";
    }
		else	out = hSignal ->process + " lnN ";
		for(int j = 0; j < nBkgs+1; j++){
			if(j != i) out += TString(" - ");
			else{
				if(i<nBkgs) out += Form(" %1.2f ", 1+VBkgs[i]->syst);
				else        out += Form(" %1.2f ", 1+hSignal->syst); 
			}
		}
		out += TString("\n");
		outputfile << out;
	}
  TString stat = GetStatUncDatacard();
  outputfile << stat;
  TString systshapes = GetShapeUncDatacard();
  outputfile << systshapes;
	outputfile.close();
  cout << "-------> Datacard created: " << LimitFolder + filename << endl;
}

TString Plot::GetShapeUncDatacard(){
	Histo* nom;
  TString sys;
	TString lin = TString(""); Int_t nsyst = VSystLabel.size();
	for(Int_t gs = 0; gs < nsyst; gs++){
    sys = VSystLabel.at(gs);
	  lin += sys + " shape ";
    if(sys == "FS"){ // FullSim/FastSim -> only signal
			for(int k = 0; k < nBkgs+1; k++){
				if(k < nBkgs) lin += TString(" - ");
				else          lin += TString(" 1 ");
			}
			lin += TString("\n");
		}
		else{
			for(int k = 0; k < nBkgs+1; k++) lin += TString(" 1 ");
			lin += TString("\n");
		}
	}
  return lin;
}

void Plot::AddSystematic(TString s){
	sys = s + "Up"; 
	SetAllProcesses(); 
	SetSignal(); 
	sys = s + "Down"; 
	SetAllProcesses(); 
	SetSignal(); 
	VSystLabel.push_back(s);
	sys = "0";
}

void Plot::AddAllSystematics(){
	AddSystematic("FS");
	AddSystematic("LepEff");
	AddSystematic("Trig");
	AddSystematic("Btag");
	AddSystematic("Misgag");
	AddSystematic("JES");
	//AddSystematic("JER");

	TString pr; TString ps;
	for(int i = 0; i < VBkgs.size(); i++){
		pr = VBkgs.at(i)->process; 
		for(int k = 0; k < VSyst.size(); k++){
			ps = VSyst.at(k)->process;
			if(ps.BeginsWith(pr)) VBkgs.at(i)->AddToSystematics(VSyst.at(k));
		}
	}
	pr = hSignal->process;
	for(int k = 0; k < VSyst.size(); k++){
		ps = VSyst.at(k)->process;
		if(ps.BeginsWith(pr)) hSignal->AddToSystematics(VSyst.at(k));
	}
	Histo* h = VBkgs.at(0);
  Int_t nbins = h->GetNbinsX();
  Float_t valu = 0; Float_t vald = 0; Float_t maxsyst = 0;
  TotalSysUp   = new Float_t[nbins];
  TotalSysDown = new Float_t[nbins];
  for(int i = 0; i < nbins; i++){
    valu = 0; vald = 0;
		for(int k = 0; k < VBkgs.size(); k++){
			valu += VBkgs.at(k)->vsysu[i];
			vald += VBkgs.at(k)->vsysd[i];
		}
   TotalSysUp[i]   = TMath::Sqrt(valu);
   TotalSysDown[i] = TMath::Sqrt(vald);
   maxsyst = TotalSysUp[i] > TotalSysDown[i]? TotalSysUp[i] : TotalSysDown[i];
   AllBkg->SetBinError(i+1, maxsyst);
	 cout << "Bin " << i+1 << ", val = " << AllBkg->GetBinContent(i+1) << ", syst up = " << TotalSysUp[i] << ", syst down = " << TotalSysDown[i] << endl;
  }
}

