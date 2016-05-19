//////////////////////////////////////////////////////////////////////////
// Use the definitions in SetPlotter.C
// Reads the rootfiles and plot histograms and get yields
// The plots are saved in plotfolder
//
// Usage:
// root -l -q -b 'DrawPlots(500, 200, 850, 100)'
// DrawPlots(500, 100, 850, 100);
// T2ttStackPlots("MET", "ElMu", "1btag", 850, 100, 500, 100);
//
//////////////////////////////////////////////////////////////////////////

#ifndef DrawPlots_h
#define DrawPlots_h
#include "SetPlotter.C"

double getYield(TString sample, TString chan = "ElMu", TString level = "2jets", TString syst = "0", bool isSS = false);
double yield(TString process, TString chan = "ElMu", TString level = "2jets", TString syst = "0", bool isSS = false);
void PrintYields(TString chan);
double getStatError(TString sample, TString chan = "ElMu", TString level = "2jets", TString syst = "0", bool isSS = false);

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
const TString plotfolder = "/nfs/fanae/user/juanr/TOP13TeV/Plotter5TeV/plots/";
bool doSys     = false;
bool doSetLogy = false;
bool dodata    = true;

//#######################################//
//################ Yields ###############//

double getYield(TString sample, TString chan, TString level, TString syst, bool isSS){
  if(chan == "All")   return getYield(sample, "ElMu", level, syst) + getYield(sample, "Muon", level, syst) + getYield(sample, "Elec", level, syst); 
  if(chan == "sameF") return getYield(sample, "Muon", level, syst) + getYield(sample, "Elec", level, syst); 

  int bin = 0;
  if(level == "dilepton") bin = 1;
  if(level == "ZVeto")    bin = 2;
  if(level == "MET")      bin = 3;
  if(level == "2jets")    bin = 4;
  if(level == "1btag")    bin = 5;
  if(level == "DYVeto")    bin = 6;
 
  TH1F* h;
  TFile* inputfile = TFile::Open(path + "/Tree_" + sample + ".root");
  TString prefix = "H_Yields_"; if(isSS) prefix = "H_SSYields_";
  if(syst == "0")  inputfile->GetObject(prefix + chan, h);
  else             inputfile->GetObject(prefix + chan + "_" + syst, h);
  double yield = h->GetBinContent(bin);
  delete inputfile;
  if( (sample.Contains("amcatnlo") || sample.Contains("aMCatNLO")) && ( sample.Contains("TTbar") || sample.Contains("TTJets") ) ) yield *= (1.08*0.9)*(1.08*0.9);
  if(!sample.Contains("Data")) return yield*Lumi;
  else return yield;
}

double yield(TString process, TString chan, TString level, TString syst, bool isSS){
  double y = 0;
  if       (process == "ttbar"  ) y = getYield("TTbar_Powheg", chan, level, syst, isSS);
  else if  (process == "tW"     ) y = getYield("TbarW", chan, level, syst, isSS) + getYield("TW", chan, level, syst, isSS);
  else if  (process == "DY"     ) y = getYield("DYJetsToLL_M10to50_aMCatNLO", chan, level, syst, isSS) + getYield("DYJetsToLL_M50_aMCatNLO", chan, level, syst, isSS);
 // else if  (process == "DY"     ) y = getYield("DYJetsToLL_M50_aMCatNLO", chan, level, syst, isSS);
  else if  (process == "VV"     ) y = getYield("WW", chan, level, syst, isSS) + getYield("WZ", chan, level, syst, isSS);// + getYield("ZZ", chan, level, syst, isSS);
  else if  (process == "fake"||process == "NonW/Z" ||process=="NonW")  y = getYield("WJetsToLNu_aMCatNLO", chan, level, syst, isSS) + getYield("TTbar_PowhegSemi", chan, level, syst, isSS); 
  else if  (process == "bkg"    ) {y = yield("ttbar", chan, level, syst, isSS) + yield("tW", chan, level, syst, isSS) + yield("DY", chan, level, syst, isSS) + yield("VV", chan, level, syst, isSS) + yield("fake", chan, level, syst, isSS); return y;}
  else if  (process == "data" || process == "Data"   ){ 
   // if      (chan == "ElMu") y = getYield("MuonEG", chan, level); 
   // else if (chan == "Elec") y = getYield("DoubleEG", chan, level);
   // else if (chan == "Muon") y = getYield("DoubleMuon", chan, level); 
   y = getYield("Data_SingleMu", chan, level, "0", isSS) + getYield("Data_SingleElec", chan, level, "0", isSS);
  }
  else y = getYield(process, chan, level, syst, isSS);
  return y;
}

double getStatError(TString sample, TString chan, TString level, TString syst, bool isSS){
  if(chan == "All")   return getStatError(sample, "ElMu", level, syst, isSS) + getStatError(sample, "Muon", level, syst, isSS) + getStatError(sample, "Elec", level, syst, isSS);
  if(chan == "sameF") return getStatError(sample, "Muon", level, syst, isSS) + getStatError(sample, "Elec", level, syst, isSS);

  if     (sample == "ttbar"  ) return getStatError("TTbar_Powheg", chan, level, syst, isSS);
  else if(sample == "tW"     ) return getStatError("TbarW", chan, level, syst, isSS) + getStatError("TW", chan, level, syst, isSS);
  else if  (sample == "DY"   ) return getStatError("DYJetsToLL_M50_aMCatNLO", chan, level, syst, isSS) + getStatError("DYJetsToLL_M10to50_aMCatNLO", chan, level, syst, isSS);
  else if  (sample == "VV"   ) return getStatError("WW", chan, level, syst, isSS) + getStatError("WZ", chan, level, syst, isSS);// + getStatError("ZZ", chan, level, syst);
  else if  (sample == "fake"||sample == "NonW/Z" ||sample =="NonW")  return getStatError("WJetsToLNu_aMCatNLO", chan, level, syst, isSS) + getStatError("TTbar_PowhegSemi", chan, level, syst, isSS); 
  else if  (sample == "bkg"  ) return  getStatError("tW", chan, level, syst, isSS) + getStatError("DY", chan, level, syst, isSS) + getStatError("VV", chan, level, syst, isSS) + getStatError("fake", chan, level, syst, isSS);
  else if  (sample == "data" || sample == "Data"   ) return getStatError("Data_SingleMu", chan, level, syst, isSS) + getStatError("Data_SingleElec", chan, level, syst, isSS);


  int bin = 0; if(level == "dilepton") bin = 1; if(level == "ZVeto") bin = 2;  if(level == "MET") bin = 3; if(level == "2jets") bin = 4; if(level == "1btag") bin = 5; if(level == "DYVeto") bin = 6;
  TH1F* h; TFile* inputfile = TFile::Open(path + "/Tree_" + sample + ".root");
  TString prefix = "H_Yields_"; if(isSS) prefix = "H_SSYields_";
  if(syst == "0")  inputfile->GetObject(prefix + chan, h);
  else             inputfile->GetObject(prefix + chan + "_" + syst, h);
  double err = h->GetBinError(bin);
  delete inputfile;
  if(!sample.Contains("Data")) return err*Lumi;
  else return err;
}
  

void PrintYields(){
  TString chan = "ElMu";
  TString c;
  cout << "Yields, Lumi = " << Lumi << " pb-1, Channel = " << chan << endl;
  cout << "===============================" << endl;
  cout << "       dilepton  |    2jets    " << endl;
  cout << "-------------------------------" << endl;
  for(int i = 0; i<5; i++){
    c = "";
    c += Form("   %5.3f  |  %5.3f" , yield(process[i], chan, "dilepton"), yield(process[i], chan, "2jets")); 
    cout << process[i] << c << endl; 
  } 
  cout << "-------------------------------" << endl;
  c = "";
  c += Form("   %5.0f  |  %5.0f ", yield("data", chan, "dilepton"), yield("data", chan, "2jets")); 
  cout << "Data " + c << endl;
  cout << "===============================" << endl;
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
  if(sample.Contains("Data") || sample.Contains("data")){ 
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.1);
  }
  else{
    hist->SetFillStyle(1001);
    hist->Scale(Lumi);
  }
  RebinHist(hist, var);
 // if(AddOverflow){
 // }
  hist -> SetDirectory(0);
  delete inputfile;
  return hist;
}

TH1F* H_data(TString var, TString chan, TString level){ 
  TH1F* h; TH1F* h1; TH1F* h2;
    h  = loadHistogram("Data_SingleMu", var, "ElMu", level); 
    h2 = loadHistogram("Data_SingleElec", var, "ElMu", level);
  return h;
}
TH1F* H_ttbar(TString var, TString chan, TString level, TString sys){
  TH1F* h = loadHSyst("TTbar_Powheg", var, chan, level, "0"); 
  h->SetFillColor(kRed);
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
  //h2   = loadHSyst("ZZ", var, chan, level, "0");
  hist->Add(h1); //hist->Add(h2);
  hist->SetFillColor(kGreen-7);
  return hist;
}

TH1F* H_DY(TString var, TString chan, TString level, TString sys){
  TH1F* h1; TH1F* h2;
  h2   = loadHSyst("DYJetsToLL_M10to50_aMCatNLO", var, chan, level, "0");
  h1   = loadHSyst("DYJetsToLL_M50_aMCatNLO", var, chan, level, "0");
  h1->Add(h2);
  h1->SetFillColor(kYellow);
  return h1;
}

TH1F* H_fake(TString var, TString chan, TString level, TString sys){
  TH1F* h; TH1F* h2;
  h = loadHSyst("WJetsToLNu_aMCatNLO", var, chan, level, "0");
  h2 = loadHSyst("TTbar_PowhegSemi", var, chan, level, "0");
  h->Add(h2);
  h->SetFillColor(kCyan+2);
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
  TH1F* httbar= H_ttbar(var, chan, level, "0"); 
  hs->Add(hVV); hs->Add(hfake); hs->Add(hDY);  hs->Add(htW); hs->Add(httbar);
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
  /*
  if( (var == "DiLepPt") || (var == "Lep0Pt") || (var == "Lep1Pt") || (var == "Jet0Pt") || (var == "Jet1Pt") || (var == "MET") ) hist -> Rebin(100);
  if( (var == "MT") || (var == "Ptllb") || (var == "MT2bb") || (var == "HT")) hist -> Rebin(100);
  if( (var == "MT2lblb")) hist -> Rebin(20);
  if( (var=="DelLepPhi") || (var == "DelPhiLepMet") || (var == "DelPhiJetMet") || (var == "DelPhiLepJet") || (var == "DelPhiPllbMet")) hist -> Rebin(10);
  if( (var == "Meff") ) hist -> Rebin(100);
  if(var == "InvMass") hist-> Rebin(30);
  if(var == "MT2") hist -> Rebin(10);
  if(var == "MinDPhiMetJets") hist -> Rebin(4);
  if((var == "MT2") || (var == "MT2lblb")) hist->SetBinContent(1, 0);
  */
 hist -> Rebin(20); 
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
    float max = hdata->GetMaximum();
    float max2 = ((TH1F*)hs->GetStack()->Last())->GetMaximum();
    if(max2 > max) max = max2 ;
  if(doSetLogy){
    plot->SetLogy();  
    hs->SetMinimum(10e-3);
    hs->SetMaximum(max*500);
  }
  else{
    hs->SetMinimum(0);
    hs->SetMaximum(max*1.16);
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

  SetLegend();
  leg->AddEntry(H_ttbar(var,chan,level, "0"), Form("ttbar: %5.2f", yield("ttbar", chan, level)), "f");
  leg->AddEntry(H_tW(var,chan,level, "0"), Form("tW: %5.2f", yield("tW", chan, level)), "f");
  leg->AddEntry(H_DY(var,chan,level, "0"), Form("DY: %5.2f", yield("DY", chan, level)), "f");
  leg->AddEntry(H_VV(var,chan,level, "0"), Form("WW+WZ: %5.2f", yield("VV", chan, level)), "f");
  leg->AddEntry(H_fake(var,chan,level, "0"), Form("NonW/Z: %5.2f", yield("fake", chan, level)), "f");
  TH1F* hh = H_ttbar(var,chan,level, "0"); hh->SetFillStyle(1001); hh->SetFillColor(kBlack);
  leg->AddEntry(hh, Form("All MC: %5.2f", yield("bkg", chan, level)), "");
  if (dodata) leg->AddEntry(hdata, Form("Data: %5.2f", yield("data", chan, level)), "pl");

  leg->Draw("same"); texlumi->Draw("same"); texcms->Draw("same");
  SetTexChan(chan, level); texchan->Draw("same");

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
 TString levels[2] = {"dilepton", "2jets"};
 TString TheChannel = "ElMu";
 for(int k = 0; k<2; k++){
	 for(int i = 0; i<nUsefulvars; i++){
		 PlotStackStop(usefulvars[i], TheChannel, levels[k]);
	 }
 }
}

#endif 
