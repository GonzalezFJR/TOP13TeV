//////////////////////////////////////////////////////////////////////////
// Use the definitions in SetPlotter.C
// Reads the rootfiles and plot histograms and get yields
// The plots are saved in plotfolder
//
// Usage:
// PlotStackStop(TString var, TString chan, TString level);
// PlotStackStop("MET", "ElMu", "1btag");
// DrawPlots();
// T2ttStackPlots("MT2", "All", "DYVeto", 400, 200);
//
//////////////////////////////////////////////////////////////////////////

#include "SetPlotter.C"

double getYield(TString sample, TString chan, TString level);
double yield(TString process, TString chan, TString level);
void PrintYields(TString chan);
float Get8TeVScale(TString sample);

TH1F* loadHistogram(TString sample, TString var, TString chan, TString level);
TH1F* H_data(TString var, TString chan, TString level);
TH1F* H_ttbar(TString var, TString chan, TString level, TString sys);
TH1F* H_tW(TString var, TString chan, TString level, TString sys);
TH1F* H_VV(TString var, TString chan, TString level, TString sys);
TH1F* H_DY(TString var, TString chan, TString level, TString sys);
TH1F* H_fake(TString var, TString chan, TString level, TString sys);
THStack* TH_bkg(TString var, TString chan, TString level);

void RebinHist(TH1F* hist, TString var);
void StackOverflow(TH1F* hist);

TH1F* BkgSystH(TString var, TString chan, TString level, TString syst);
float* BinVar2(TString var, TString chan, TString level);
TH1F* loadHSyst(TString sample, TString var, TString chan, TString level, TString syst);

void PlotStackStop(TString var, TString chan, TString level);
void DrawPlots();


// Constants:
const TString plotfolder = "/nfs/fanae/user/juanr/StopTOP/plots/";
bool doSys     = false;
bool doSetLogy = true;
bool dodata    = false;
bool do8TeV    = false;

//#######################################//
//################ Yields ###############//

double getYield(TString sample, TString chan, TString level){
  if(chan == "All")   return getYield(sample, "ElMu", level) + getYield(sample, "Muon", level) + getYield(sample, "Elec", level); 
  if(chan == "sameF") return getYield(sample, "Muon", level) + getYield(sample, "Elec", level); 

  int bin = 0;
  if(level == "dilepton") bin = 1;
  if(level == "ZVeto")    bin = 2;
  if(level == "MET")      bin = 3;
  if(level == "2jets")    bin = 4;
  if(level == "1btag")    bin = 5;
  if(level == "DYVeto")    bin = 6;
 
  TH1F* h;
  TFile* inputfile = TFile::Open(path + "/Tree_" + sample + ".root");
  inputfile->GetObject("H_Yields_" + chan, h);
  double yield = h->GetBinContent(bin);
  delete inputfile;
  if(sample != "DoubleMuon" && sample != "DoubleEG" && sample != "MuonEG") return yield*Lumi;
  else return yield;
}

double yield(TString process, TString chan, TString level){
  double y = 0;
//  if       (process == "ttbar"  ) y = getYield("TTbar_Powheg", chan, level) + getYield("TTJetsSemi", chan, level);
  if       (process == "ttbar"  ) y = getYield("TTJets", chan, level);//+ getYield("TTJetsSemi", chan, level);
  else if  (process == "tW"     ) y = getYield("TbarW", chan, level) + getYield("TW", chan, level);
 // else if  (process == "DY"     ) y = getYield("DYJetsToLL_M10to50_aMCatNLO", chan, level) + getYield("DYJetsToLL_M50_aMCatNLO", chan, level);
  else if  (process == "DY"     ) y = getYield("DYJetsToLL_M5to50_MLM", chan, level) + getYield("DYJetsToLL_M50_MLM", chan, level);
  else if  (process == "VV"     ) y = getYield("WW", chan, level) + getYield("WZ", chan, level) + getYield("ZZ", chan, level);
  else if  (process == "ttV"    ) y = getYield("TTWToLNu", chan, level) + getYield("TTWToQQ", chan, level) + getYield("TTZToQQ", chan, level);
  else if  (process == "fake"||process == "NonW/Z"||process == "WJets")  y = getYield("WJetsToLNu_aMCatNLO", chan, level); 
  //else if  (process == "fake"||process == "NonW/Z"||process == "WJets")  y = getYield("WJetsToLNu_MLM", chan, level); 
  else if  (process == "bkg"    ) {y = yield("ttbar", chan, level) + yield("tW", chan, level) + yield("DY", chan, level) + yield("VV", chan, level) + yield("fake", chan, level) + yield("ttV", chan, level); return y;}
  else if  (process == "data" || process == "Data"   ){ 
    if      (chan == "ElMu") y = getYield("MuonEG", chan, level); 
    else if (chan == "Elec") y = getYield("DoubleEG", chan, level);
    else if (chan == "Muon") y = getYield("DoubleMuon", chan, level); 
  }
  else y = getYield(process, chan, level);
  return y;
}

void PrintYields(TString chan){
  TString c;
  cout << "Yields, Lumi = " << Lumi << " pb-1, Channel = " << chan << endl;
  cout << "====================================================================================" << endl;
  cout << "       dilepton  |    ZVeto    |     MET     |    2jets    |   1btag   |    DYVeto" << endl;
  cout << "------------------------------------------------------------------------------------" << endl;
  for(int i = 0; i<5; i++){
    c = "";
    for(int k = 0; k<6; k++){ c += Form("   %5.2f  | ", yield(process[i], chan, level[k])); }
    cout << process[i] << c << endl; 
  } 
  cout << "------------------------------------------------------------------------------------" << endl;
  c = "";
  for(int k = 0; k<6; k++){ c += Form("   %5.0f  | ", yield("data", chan, level[k])); }
  cout << "Data " + c << endl;
  cout << "------------------------------------------------------------------------------------" << endl;
  c = "";
  for(int k = 0; k<6; k++){ c += Form("   %2.4f  | ", yield("T2tt_500_325", chan, level[k])); }
  cout << "S.500_325 " + c << endl;
  c = "";
  for(int k = 0; k<6; k++){ c += Form("   %2.4f  | ", yield("T2tt_850_100", chan, level[k])); }
  cout << "S.850_100 " + c << endl;
  cout << "=====================================================================================" << endl;
}

//################################################//
//################ Load histograms ###############//

TH1F* loadHistogram(TString sample, TString var, TString chan, TString level){
  TH1F* hist;
  if(chan == "All"){
    hist = loadHistogram(sample, var, "ElMu", level);
    hist->Add(loadHistogram(sample, var, "Muon", level));
    hist->Add(loadHistogram(sample, var, "Elec", level));
    hist -> SetDirectory(0);
    return hist;
  }
  if((chan == "SF") || (chan == "sameF")){
    hist = loadHistogram(sample, var, "Elec", level);
    hist->Add(loadHistogram(sample, var, "Muon", level));
    hist -> SetDirectory(0);
    return hist;
  }
	TString thefile;
  if(sample.BeginsWith("T2tt")) thefile = path + "/Susy/Tree_" + sample + ".root";
	else                          thefile = path + "Tree_" + sample + ".root";
  TFile* inputfile = TFile::Open(thefile);
  inputfile->GetObject("H_" + var + "_" + chan + "_" + level,hist);
  hist -> SetStats(0);
  hist->SetLineStyle(0);
  //hist -> SetTitle(var + " " + chan + " " + level);
  hist->SetTitle(" ");
  hist -> GetXaxis() -> SetTitle(var);
  //hist -> GetXaxis() -> CenterTitle();
  if( (sample == "MuonEG") || (sample == "DoubleEG") || (sample == "DoubleMuon") ){
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.1);
  }
  else{
    hist->SetFillStyle(1001);
    hist->Scale(Lumi);
  }
  RebinHist(hist, var);
  if(do8TeV) hist->Scale(Get8TeVScale(sample));
 // if(AddOverflow){
 // }
  hist -> SetDirectory(0);
  delete inputfile;
  return hist;
}

TH1F* H_data(TString var, TString chan, TString level){ 
  TH1F* h; TH1F* h1; TH1F* h2;
  if      (chan == "ElMu") h = loadHistogram("MuonEG", var, chan, level); 
  else if (chan == "Muon") h = loadHistogram("DoubleMuon", var, chan, level);
  else if (chan == "Elec") h = loadHistogram("DoubleEG", var, chan, level);
  else if ( (chan == "sameF") || (chan == "SF") ){ h = loadHistogram("DoubleMuon", var, "Muon", level); h2 = loadHistogram("DoubleEG", var, "Elec", level); h->Add(h2);}
  else{
    h  = loadHistogram("MuonEG", var, "ElMu", level); 
    h1 = loadHistogram("DoubleMuon", var, "Muon", level);  
    h2 = loadHistogram("DoubleEG", var, "Elec", level); 
    h->Add(h1); h->Add(h2);
  }
  return h;
}
TH1F* H_ttbar(TString var, TString chan, TString level, TString sys){
  TH1F* h = loadHSyst("TTJets", var, chan, level, "0"); 
  h->SetFillColor(kRed);
  return h;
}

TH1F* H_ttV(TString var, TString chan, TString level, TString sys){
  TH1F* h  = loadHSyst("TTWToLNu", var, chan, level, "0");
  TH1F* h1 = loadHSyst("TTWToQQ", var, chan, level, "0");
  TH1F* h2 = loadHSyst("TTZToQQ", var, chan, level, "0");
  h->Add(h1); h->Add(h2);
  h->SetFillColor(kOrange+2);
  return h;
}

TH1F* H_VVV(TString var, TString chan, TString level, TString sys){
  TH1F* h  = loadHSyst("WWZ", var, chan, level, "0");
  TH1F* h1 = loadHSyst("WZZ", var, chan, level, "0");
  TH1F* h2 = loadHSyst("ZZZ", var, chan, level, "0");
  h->Add(h1); h->Add(h2);
  h->SetFillColor(kViolet-1);
  return h;
}

TH1F* H_tW(TString var, TString chan, TString level, TString sys){
  TH1F* hist; TH1F* h;
  hist = loadHSyst("TbarW", var, chan, level, "0");
  h    = loadHSyst("TW"   , var, chan, level, "0");
  hist->Add(h);
  hist->SetFillColor(kViolet-3); 
  return hist;
}

TH1F* H_VV(TString var, TString chan, TString level, TString sys){
  TH1F* hist; TH1F* h1; TH1F* h2;
  hist = loadHSyst("WW", var, chan, level, "0");
  h1   = loadHSyst("WZ", var, chan, level, "0");
  h2   = loadHSyst("ZZ", var, chan, level, "0");
  hist->Add(h1); hist->Add(h2);
  hist->SetFillColor(kGreen-7);
  return hist;
}

TH1F* H_DY(TString var, TString chan, TString level, TString sys){
  TH1F* h1; TH1F* h2;
  //h1   = loadHSyst("DYJetsToLL_M10to50_aMCatNLO", var, chan, level, "0");
  //h2   = loadHSyst("DYJetsToLL_M50_aMCatNLO", var, chan, level, "0");
  h1   = loadHSyst("DYJetsToLL_M5to50_MLM", var, chan, level, "0");
  h2   = loadHSyst("DYJetsToLL_M50_MLM", var, chan, level, "0");
  h1->Add(h2);
  h1->SetFillColor(kYellow);
  return h1;
}

TH1F* H_fake(TString var, TString chan, TString level, TString sys){
  TH1F* h; TH1F* h2;
  h = loadHSyst("WJetsToLNu_aMCatNLO", var, chan, level, "0");
  h->SetFillColor(kGreen-9);
  return h;
}

TH1F* loadProcessHistogram(TString process, TString var, TString chan, TString level, TString sys){
	if((process == "tW"))                            return H_tW(var, chan, level, sys);
	if((process == "DY"))                            return H_DY(var, chan, level, sys);
	if((process == "VV") || (process == "dibosons")) return H_VV(var, chan, level, sys);
	if((process == "NonW/Z") || (process == "fake")) return H_fake(var, chan, level, sys);
	if((process == "data") || (process == "Data"))   return H_data(var, chan, level);
	if((process == "ttbar") || (process == "tt"))    return H_ttbar(var, chan, level, sys);
	else    return H_ttbar(var, chan, level, sys);
}

THStack* TH_bkg(TString var, TString chan, TString level){
  THStack* hs = new THStack(var, "");
  TH1F* hDY   = H_DY(var, chan, level, "0");  
  TH1F* htW   = H_tW(var, chan, level, "0");  
  TH1F* hfake = H_fake(var, chan, level, "0"); 
  TH1F* hVV   = H_VV(var, chan, level, "0"); 
  TH1F* ttV   = H_ttV(var, chan, level, "0"); 
  TH1F* httbar= H_ttbar(var, chan, level, "0"); 
  hs->Add(hfake); hs->Add(hVV); hs->Add(ttV); hs->Add(hDY); hs->Add(htW); hs->Add(httbar);
  return hs;
}

// H_ttbar, H_ttV, H_tW, H_VV, H_DY, H_fake
TH1F* loadHSyst(TString sample, TString var, TString chan, TString level, TString syst){
  TH1F* hist;
  if((syst == "0") || (syst == "norm")) return loadHistogram(sample, var, chan, level);
  if(chan == "All"){
    hist = loadHSyst(sample, var, "ElMu", level, syst);
    hist->Add(loadHSyst(sample, var, "Muon", level, syst));
    hist->Add(loadHSyst(sample, var, "Elec", level, syst));
    hist -> SetDirectory(0);
    return hist;
  }
  if((chan == "SF") || (chan == "sameF")){
    hist = loadHSyst(sample, var, "Elec", level, syst);
    hist->Add(loadHSyst(sample, var, "Muon", level, syst));
    hist -> SetDirectory(0);
    return hist;
  }
	TString thefile;
  if(sample.BeginsWith("T2tt")) thefile = path + "/Susy/Tree_" + sample + ".root";
	else                          thefile = path + "Tree_" + sample + ".root";
  TFile* inputfile = TFile::Open(thefile);
  hist -> SetStats(0);

  if( (sample == "MuonEG") || (sample == "DoubleEG") || (sample == "DoubleMuon") ){
    cout << "No systematics for data!" << endl;
    return hist;
  }
  else  hist->Scale(Lumi);
  RebinHist(hist, var);
  hist -> SetDirectory(0);
  delete inputfile;
  return hist;
}

void RebinHist(TH1F* hist, TString var){
  if( (var == "DiLepPt") || (var == "Lep0Pt") || (var == "Lep1Pt") || (var == "Jet0Pt") || (var == "Jet1Pt") || (var == "MET") ) hist -> Rebin(25);
  if( (var == "MT") || (var == "Ptllb") || (var == "MT2bb") || (var == "HT")) hist -> Rebin(30);
  if( (var == "MT2lblb")) hist -> Rebin(20);
  if( (var=="DelLepPhi") || (var == "DelPhiLepMet") || (var == "DelPhiJetMet") || (var == "DelPhiLepJet") || (var == "DelPhiPllbMet")) hist -> Rebin(10);
  if( (var == "Meff") ) hist -> Rebin(100);
  if(var == "InvMass") hist-> Rebin(40);
  if(var == "MT2") hist -> Rebin(10);
  if(var == "MinDPhiMetJets") hist -> Rebin(4);
  if((var == "MT2") || (var == "MT2lblb")) hist->SetBinContent(1, 0);
}

void StackOverflow(TH1F* hist){
  int nbins = hist->GetNbinsX();
  float overflowYield = hist->GetBinContent(nbins+1);
  float lastBinYield = hist->GetBinContent(nbins);
  hist->SetBinContent(nbins+1, lastBinYield+overflowYield);
}

//###########################################//
//############## Systematics  ###############//

TH1F* BkgSystH(TString var, TString chan, TString level, TString syst){
  TH1F* h = loadHSyst(bkg[0], var, chan, level, syst);
  for(int i = 1; i<nBkgs; i++){
    h->Add(loadHSyst(bkg[i], var, chan, level, syst));
  }
  return h;
}

float* BinVar2(TString var, TString chan, TString level){
  float vpos = 0; float vneg = 0; float f; float valnom = 0; float valsys = 0;
  TH1F* hnom = (TH1F*) TH_bkg(var, chan, level)->GetStack()->Last();
  TH1F* hsys[nSyst];
  for(int s = 0; s < nSyst; s++) hsys[s] = BkgSystH(var, chan, level, syst[s]);
  int n = hnom->GetNbinsX()+1;
  float *v; v = new float[2*n+1];
  float staterr = 0;
  for(int bin = 0; bin < n; bin++){
    valnom = hnom->GetBinContent(bin);
    staterr = hnom->GetBinError(bin);
    vpos = 0; vneg = 0;
    for(int s = 0; s < nSyst; s++){
      valsys = hsys[s]->GetBinContent(bin);
      f = (valsys - valnom);
      if (f>0) vpos += f*f;
      else     vneg += f*f;
    }
    vpos = vpos + staterr*staterr;
    vneg = vneg + staterr*staterr;
    v[2*bin] = vpos;
    v[2*bin+1] = vneg;
  }
  return v;
}

//#####################################//
//################ Plot ###############//

void PlotStackStop(TString var, TString chan, TString level){
  TCanvas* c = SetCanvas();
  plot->cd();
  THStack* hs = TH_bkg(var, chan, level);
  TH1F* hdata = H_data(var, chan, level);
  if(doSetLogy){
    plot->SetLogy();  
    hs->SetMinimum(10e-3);
    hs->SetMaximum(hdata->GetMaximum()*500);
  }
  else{
    hs->SetMinimum(0);
    hs->SetMaximum(hdata->GetMaximum()*1.15);
  }
  hs->Draw("hist");
  hs->GetYaxis()->SetTitle("Number of Events");
  hs->GetYaxis()->SetTitleSize(0.06);
  hs->GetYaxis()->SetTitleOffset(0.5);
  hs->GetYaxis()->SetNdivisions(505);
	if (dodata){
		hdata->Draw("pesame");
		hs->GetXaxis()->SetLabelSize(0.0);
	}
  plot->RedrawAxis("same");
 // TH1F* h_stop_1 = loadHistogram("T2tt_500_325", var, chan, level); h_stop_1->SetFillColor(0); h_stop_1->SetLineColor(kBlue+3);    h_stop_1->SetLineWidth(2); h_stop_1->SetLineStyle(1); 
 // TH1F* h_stop_2 = loadHistogram("T2tt_850_100", var, chan, level); h_stop_2->SetFillColor(0); h_stop_2->SetLineColor(kMagenta-4); h_stop_2->SetLineWidth(2); h_stop_2->SetLineStyle(1); 

  SetLegend();
  leg->AddEntry(H_ttbar(var,chan,level, "0"), Form("ttbar: %5.2f", yield("ttbar", chan, level)), "f");
  leg->AddEntry(H_tW(var,chan,level, "0"), Form("tW: %5.2f", yield("tW", chan, level)), "f");
  leg->AddEntry(H_DY(var,chan,level, "0"), Form("DY: %5.2f", yield("DY", chan, level)), "f");
  leg->AddEntry(H_VV(var,chan,level, "0"), Form("VV: %5.2f", yield("VV", chan, level)), "f");
  leg->AddEntry(H_ttV(var,chan,level, "0"), Form("ttV: %5.2f", yield("ttV", chan, level)), "f");
  leg->AddEntry(H_fake(var,chan,level, "0"), Form("W+Jets: %5.2f", yield("fake", chan, level)), "f");
  TH1F* hh = H_ttbar(var,chan,level, "0"); hh->SetFillStyle(1001); hh->SetFillColor(kBlack);
  leg->AddEntry(hh, Form("All bkg: %5.2f", yield("bkg", chan, level)), "");
  if (dodata) leg->AddEntry(hdata, Form("Data: %5.2f", yield("data", chan, level)), "pl");
 // leg->AddEntry(h_stop_1, Form("S.500-325 (x10) : %5.2f", h_stop_1->Integral()), "l");
 // leg->AddEntry(h_stop_2, Form("S.850-100 (x100): %5.2f", h_stop_2->Integral()), "l");

  leg->Draw("same"); texlumi->Draw("same"); texcms->Draw("same");
  SetTexChan(chan, level); texchan->Draw("same");
 // h_stop_1->Scale(20);  h_stop_1->Draw("same,hist");
 // h_stop_2->Scale(100); h_stop_2->Draw("same,hist");


  TH1F* mh = (TH1F*)hs->GetStack()->Last()->Clone();
  
  TH1F* hratioerror;
  if(doSys){
    int n = (mh->GetNbinsX()+1);
    float* v = BinVar2(var, chan, level);
    float staterror; float totalerror; float mhval; float hratioerbin;
    hratioerror = (TH1F*) mh->Clone();
    for(int p = 1; p<n; p++){
      mhval = mh->GetBinContent(p);
      totalerror = TMath::Max( TMath::Sqrt(v[2*p]), TMath::Sqrt(v[2*p+1]));
      mh->SetBinError(p, totalerror);
      hratioerbin = (mhval > 0) ? totalerror/mhval : 0.0;
      hratioerror->SetBinContent(p, 1);
      hratioerror->SetBinError(p, hratioerbin); 
    }
    hratioerror->SetFillColor(kTeal-7);
    hratioerror->SetFillStyle(3144);
    hratioerror->SetMarkerSize(0);
    mh->SetFillColor  (kGray+2);
    mh->SetFillStyle  (   3345);
    mh->SetLineColor  (kGray+2);
    mh->SetMarkerColor(kGray+2);
    mh->SetMarkerSize (      0);
    mh->Draw("same,e2");
    if(dodata)  hdata->Draw("pesame");
  }
  pratio->cd();
  hratio = (TH1F*)hdata->Clone();
  hratio->Divide(mh);
	SetHRatio(var);  
	if(dodata){
		hratio->Draw("same");
		if(doSys) hratioerror->Draw("same,e2");
		hratio->Draw("same");
	}

  c->Print( plotfolder + chan + "/" + var + "_" + chan + "_" + level + ".pdf", "pdf");
  c->Print( plotfolder + chan + "/" + var + "_" + chan + "_" + level + ".png", "png");
  delete c;
}

void DrawPlots(){
 //TString levels[3] = {"dilepton", "2jets", "1btag"};
 TString levels[4] = {"dilepton", "2jets", "1btag","DYVeto"};
 for(int k = 0; k<4; k++){
  for(int i = 0; i<nVars; i++){
    for(int j = 0; j<nChannels; j++){
      PlotStackStop(var[i], chan[j], levels[k]);
    }
  }
 }
}


void T2ttStackPlots(TString var, TString chan, TString level, int mStop, int mLsp){
  TCanvas* c = SetCanvas();
  plot->cd();
  THStack* hs = TH_bkg(var, chan, level);
  TH1F* hdata = H_data(var, chan, level);
  if(doSetLogy){
    plot->SetLogy();  
    hs->SetMinimum(10e-3);
    hs->SetMaximum(hdata->GetMaximum()*500);
  }
  else{
    hs->SetMinimum(0);
    hs->SetMaximum(hdata->GetMaximum()*1.15);
  }
  hs->Draw("hist");
  hs->GetYaxis()->SetTitle("Number of Events");
  hs->GetYaxis()->SetTitleSize(0.06);
  hs->GetYaxis()->SetTitleOffset(0.5);
  hs->GetYaxis()->SetNdivisions(505);
  hs->GetXaxis()->SetLabelSize(0.0);
  if (dodata) hdata->Draw("pesame");
  plot->RedrawAxis("same");
  TString signal = Form("T2tt_mStop%i_mLsp%i", mStop, mLsp);
  TH1F* h_stop_1 = loadHistogram(signal, var, chan, level); h_stop_1->SetFillColor(0); h_stop_1->SetLineColor(kMagenta-4);    h_stop_1->SetLineWidth(2); h_stop_1->SetLineStyle(1); 

  SetLegend();
  leg->AddEntry(H_ttbar(var,chan,level, "0"), Form("ttbar: %5.2f", yield("ttbar", chan, level)), "f");
  leg->AddEntry(H_tW(var,chan,level, "0"), Form("tW: %5.2f", yield("tW", chan, level)), "f");
  leg->AddEntry(H_DY(var,chan,level, "0"), Form("DY: %5.2f", yield("DY", chan, level)), "f");
  leg->AddEntry(H_VV(var,chan,level, "0"), Form("VV: %5.2f", yield("VV", chan, level)), "f");
  leg->AddEntry(H_ttV(var,chan,level, "0"), Form("ttV: %5.2f", yield("ttV", chan, level)), "f");
  leg->AddEntry(H_fake(var,chan,level, "0"), Form("W+Jets: %5.2f", yield("fake", chan, level)), "f");
  TH1F* hh = H_ttbar(var,chan,level, "0"); hh->SetFillStyle(1001); hh->SetFillColor(kBlack);
  leg->AddEntry(hh, Form("All bkg: %5.2f", yield("bkg", chan, level)), "");
  if (dodata) leg->AddEntry(hdata, Form("Data: %5.2f", yield("data", chan, level)), "pl");
  leg->AddEntry(h_stop_1, Form("S.%i-%i (x20) : %5.2f", mStop, mLsp, h_stop_1->Integral()), "l");

  leg->Draw("same"); texlumi->Draw("same"); texcms->Draw("same");
  SetTexChan(chan, level); texchan->Draw("same");
  h_stop_1->Scale(20);  h_stop_1->Draw("same,hist");


  TH1F* mh = (TH1F*)hs->GetStack()->Last()->Clone();
  
  TH1F* hratioerror;
  if(doSys){
    int n = (mh->GetNbinsX()+1);
    float* v = BinVar2(var, chan, level);
    float staterror; float totalerror; float mhval; float hratioerbin;
    hratioerror = (TH1F*) mh->Clone();
    for(int p = 1; p<n; p++){
      mhval = mh->GetBinContent(p);
      totalerror = TMath::Max( TMath::Sqrt(v[2*p]), TMath::Sqrt(v[2*p+1]));
      mh->SetBinError(p, totalerror);
      hratioerbin = (mhval > 0) ? totalerror/mhval : 0.0;
      hratioerror->SetBinContent(p, 1);
      hratioerror->SetBinError(p, hratioerbin); 
    }
    hratioerror->SetFillColor(kTeal-7);
    hratioerror->SetFillStyle(3144);
    hratioerror->SetMarkerSize(0);
    mh->SetFillColor  (kGray+2);
    mh->SetFillStyle  (   3345);
    mh->SetLineColor  (kGray+2);
    mh->SetMarkerColor(kGray+2);
    mh->SetMarkerSize (      0);
    mh->Draw("same,e2");
    if(dodata)  hdata->Draw("pesame");
  }
  pratio->cd();
  hratio = (TH1F*)hdata->Clone();
  hratio->Divide(mh);
  SetHRatio(var);
  if(dodata){
    hratio->Draw("same");
    if(doSys) hratioerror->Draw("same,e2");
    hratio->Draw("same");
  }

  c->Print( plotfolder + chan + Form("/mStop%i_mLsp%i_", mStop, mLsp) + var + "_" + chan + "_" + level + ".pdf", "pdf");
  c->Print( plotfolder + chan + Form("/mStop%i_mLsp%i_", mStop, mLsp) + var + "_" + chan + "_" + level + ".png", "png");
  delete c;
}


float Get8TeVScale(TString sample){
  float fact = 1;
  if      (sample == "TTJets"                      ) fact = 245.8/831.8;
  else if (sample == "TW"                          ) fact = 7.8/35.6;
  else if (sample == "TbarW"                       ) fact = 7.8/35.6;
  else if (sample == "DYJetsToLL_M10to50_aMCatNLO" ) fact = 18610/6025.2*3532;
  else if (sample == "DYJetsToLL_M50_aMCatNLO"     ) fact = 3532/6025.2;
  else if (sample == "WW"                          ) fact = 1/1.96;
  else if (sample == "WZ"                          ) fact = 22.44/47.1;
  else if (sample == "ZZ"                          ) fact = 9.03/16.5;
  else if (sample == "TTWToLNu"                    ) fact = 1/2.8;
  else if (sample == "TTWToQQ"                     ) fact = 1/2.8;
  else if (sample == "TTZToQQ"                     ) fact = 1/3.5;
	else if (sample == "WJetsToLNu_aMCatNLO"         ) fact = 1/1.64;
  else                   fact = 1;
  return fact;
}

