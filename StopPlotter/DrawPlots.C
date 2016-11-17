
// Use the definitions in SetPlotter.C
// Reads the rootfiles and plot histograms and get yields
// The plots are saved in plotfolder
//
// Usage:
// PlotStackStop(TString var, TString chan, TString level);
// PlotStackStop("MET", "ElMu", "1btag");
// DrawPlots();
// T2ttStackPlots("MET", "ElMu", "1btag", 850, 100, 500, 100)
//
//////////////////////////////////////////////////////////////////////////

#include "SetPlotter.C"

double getYield(TString sample, TString chan = "ElMu", TString level = "2jets", TString syst = "0", bool isSS = false);
double yield(TString process, TString chan = "ElMu", TString level = "2jets", TString syst = "0", bool isSS = false);
float gettrigeff(TString chan = "ElMu");

//double getYield(TString sample, TString chan, TString level);
//double yield(TString process, TString chan, TString level);
void PrintYields(TString chan = "ElMu");
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
TH1F* RebinVar(TH1F* hist, TString var);
void StackOverflow(TH1F* hist);

TH1F* BkgSystH(TString var, TString chan, TString level, TString syst);
float* BinVar2(TString var, TString chan, TString level);
TH1F* loadHSyst(TString sample, TString var, TString chan, TString level, TString syst);

void PlotStackStop(TString var, TString chan, TString level);
void T2ttStackPlots(TString var, TString chan = "All", TString level = "DYVeto", int mStop1 = 0, int mLsp1 = 0, int mStop2 = 0, int mLsp2 = 0);
void DrawT2ttPlots(int mStop1, int mLsp1, int mStop2 = 0, int mLsp2 = 0);
void DrawPlots();


// Constants:
const TString plotfolder = "/mnt_pool/fanae105/user/juanr/stop80/StopPlotter/plots/";
bool doSys       = false;
bool doSetLogy   = true;
bool dodata      = true;
bool doYields    = false;
bool do8TeV      = false;
bool doSingleLep = false;

//#######################################//
//################ Yields ###############//
double yield(TString process, TString chan, TString level, TString syst, bool isSS){
  double y = 0;
  if       (process == "ttbar"  ) y = getYield("TTbar_Powheg", chan, level);//+ getYield("TTJetsSemi", chan, level);
  else if  (process == "tW"     ) y = getYield("TbarW", chan, level) + getYield("TW", chan, level);
  else if  (process == "DY"     ) y = getYield("DYJetsToLL_M10to50_aMCatNLO", chan, level) + getYield("DYJetsToLL_M50_aMCatNLO", chan, level);
 // else if  (process == "DY"     ) y = getYield("DYJetsToLL_M10to50_aMCatNLO_ext", chan, level) + getYield("DYJetsToLL_M10to50_aMCatNLO_ext2", chan, level) + getYield("DYJetsToLL_M50_aMCatNLO", chan, level) + getYield("DYJetsToLL_M50_aMCatNLO_ext", chan, level);
  else if  (process == "VV"     ) y = getYield("WW", chan, level) + getYield("WZ", chan, level) + getYield("ZZ", chan, level);
  else if  (process == "ttV"    ) y = getYield("TTWToLNu", chan, level) + getYield("TTWToQQ", chan, level) + getYield("TTZToQQ", chan, level) + getYield("TTZToLLNuNu", chan, level);
  else if  (process == "fake"||process == "NonW/Z"||process == "WJets")  y = getYield("WJetsToLNu_aMCatNLO", chan, level); 
  else if  (process == "bkg"    ) {y = yield("ttbar", chan, level) + yield("tW", chan, level) + yield("DY", chan, level) + yield("VV", chan, level) + yield("fake", chan, level) + yield("ttV", chan, level); return y;}
  else if  (process == "data" || process == "Data"   ){ 
    if      (chan == "ElMu"){ y = getYield("MuonEG", chan, level); 
      if(doSingleLep){
         y += getYield("SingleMuon", chan, level);
         y += getYield("SingleElectron", chan, level);
      }
    }
    else if (chan == "Elec"){ y = getYield("DoubleEG", chan, level);
      if(doSingleLep){
         y += getYield("SingleElectron", chan, level);
      }
    }
    else if (chan == "Muon"){ y = getYield("DoubleMuon", chan, level); 
      if(doSingleLep){
         y += getYield("SingleMuon", chan, level);
      }
    }
    else if (chan == "sameF") y = getYield("DoubleMuon", "Muon", level) + getYield("DoubleEG", "Elec", level); 
    else if (chan == "All") y = getYield("MuonEG", chan, level) + yield("data", "sameF", level);
  }
  else y = getYield(process, chan, level);
  return y;
}

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
  if( sample.Contains("DYJetsToLL_M10to50")) yield *= (18610./22635.09);
  if(sample.Contains("Data") || sample.Contains("data") || sample.Contains("Double") || sample.Contains("Muon") || sample.Contains("EG")) return yield;
  else return yield*Lumi*gettrigeff(chan);
}



void PrintYields(TString chan){
  TString c;
  cout << "Yields, Lumi = " << Lumi << " pb-1, Channel = " << chan << endl;
  cout << "====================================================================================" << endl;
  cout << "       dilepton  |    ZVeto    |     MET     |    2jets    |   1btag   |    DYVeto" << endl;
  cout << "------------------------------------------------------------------------------------" << endl;
  for(int i = 0; i<6; i++){
    c = "";
    for(int k = 0; k<6; k++){ c += Form("   %5.2f  | ", yield(process[i], chan, Level[k])); }
    cout << process[i] << c << endl; 
  } 
  cout << "------------------------------------------------------------------------------------" << endl;
  c = "";
    for(int k = 0; k<6; k++){ c += Form("   %5.2f  | ", yield("bkg", chan, Level[k])); }
    cout << "All bkg " << c << endl; 
  cout << "------------------------------------------------------------------------------------" << endl;
  c = "";
  for(int k = 0; k<6; k++){ c += Form("   %5.0f  | ", yield("data", chan, Level[k])); }
  cout << "Data " + c << endl;
//  cout << "------------------------------------------------------------------------------------" << endl;
//  c = "";
//  for(int k = 0; k<6; k++){ c += Form("   %2.4f  | ", yield("T2tt_500_325", chan, level[k])); }
//  cout << "S.500_325 " + c << endl;
//  c = "";
//  for(int k = 0; k<6; k++){ c += Form("   %2.4f  | ", yield("T2tt_850_100", chan, level[k])); }
//  cout << "S.850_100 " + c << endl;
  cout << "=====================================================================================" << endl;
}

float gettrigeff(TString chan){
  float trigeff = 1;
  if(chan == "Elec") trigeff =  0.90;//0.90;//0.97;
  else if(chan == "Muon") trigeff = 0.88;//0.88;//0.94;
  else if(chan == "ElMu") trigeff = 0.90;//0.90;//0.89;
  return 1;
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
  if     (sample.BeginsWith("T2tt"))     thefile = path + "/Susy/Tree_" + sample + ".root";
  else if(sample.BeginsWith("TTDMJets")) thefile = path + "/DM/Tree_"   + sample + ".root";
	else                                   thefile = path + "Tree_" + sample + ".root";
  TFile* inputfile = TFile::Open(thefile);
  TString hname = "H_" + var + "_" + chan + "_" + level;
  if(var == "Yields") hname = "H_Yields_" + chan;
  inputfile->GetObject(hname, hist);
  if(var == "Yields"){
    hist->SetBinContent(7,0);
    hist->SetBinContent(8,0);
  }
  hist -> SetStats(0);
  hist->SetLineStyle(0);
  //hist -> SetTitle(var + " " + chan + " " + level);
  hist->SetTitle(" ");
  hist -> GetXaxis() -> SetTitle(var);
  //hist -> GetXaxis() -> CenterTitle();
  if( (sample == "MuonEG") || (sample == "DoubleEG") || (sample == "DoubleMuon") || (sample.Contains("Single")) ){
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1.1);
  }
  else{
    hist->SetFillStyle(1001);
    hist->Scale(Lumi*gettrigeff(chan));
  }
  if(var == "MET" || var.BeginsWith("MT2") || (var.BeginsWith("Del") && var.Contains("Phi")) || ( (var == "Jet0Pt") || (var == "Jet1Pt")) || ( (var == "Lep0Pt") || (var == "DiLepPt") || (var == "Lep1Pt")) || (var == "InvMass") || (var == "Yields2")){
    hist = RebinVar(hist, var);
  }
  else RebinHist(hist, var);
  if(do8TeV) hist->Scale(Get8TeVScale(sample));
 // if(AddOverflow){
    int nbins = hist->GetNbinsX();
    float fb  = hist->GetBinContent(nbins);
    float of  = hist->GetBinContent(nbins+1);
    hist->SetBinContent(nbins, fb+of);
 // }
  if( sample.Contains("DYJetsToLL_M10to50")) hist->Scale(18610./22635.09);
  hist -> SetDirectory(0);
  delete inputfile;
  return hist;
}

TH1F* H_data(TString var, TString chan, TString level){ 
  TH1F* h; TH1F* h1; TH1F* h2; TH1F* hse; TH1F* hsm;
  if      (chan == "ElMu"){
    h   = loadHistogram("MuonEG", var, chan, level); 
    hse = loadHistogram("SingleElectron", var, chan, level); 
    hsm = loadHistogram("SingleMuon", var, chan, level); 
    if(doSingleLep){ h->Add(hse); h->Add(hsm);}
  }
  else if (chan == "Muon"){
     h   = loadHistogram("DoubleMuon", var, chan, level);
    hsm = loadHistogram("SingleMuon", var, chan, level);
    if(doSingleLep) h->Add(hsm);
  }
  else if (chan == "Elec"){
    h = loadHistogram("DoubleEG", var, chan, level);
    hse = loadHistogram("SingleElectron", var, chan, level); 
		if(doSingleLep) h->Add(hse);
	}
	else if ( (chan == "sameF") || (chan == "SF") ){
		h = loadHistogram("DoubleMuon", var, "Muon", level); 
    hsm = loadHistogram("SingleMuon", var, "Muon", level); 
    hse = loadHistogram("SingleElectron", var, "Elec", level); 
		h2 = loadHistogram("DoubleEG", var, "Elec", level); 
	  h->Add(h2); if(doSingleLep){ h->Add(hse); h->Add(hsm);}
	}
	else{
    h  = loadHistogram("MuonEG", var, "ElMu", level); 
    h1 = loadHistogram("DoubleMuon", var, "Muon", level);  
    h2 = loadHistogram("DoubleEG", var, "Elec", level); 
    hsm = loadHistogram("SingleMuon", var, "Muon", level); 
    hse = loadHistogram("SingleElectron", var, "Elec", level); 
    TH1F* hemsm = loadHistogram("SingleMuon", var, "ElMu", level);
    TH1F* hemse = loadHistogram("SingleElectron", var, "ElMu", level);
    h->Add(h1); h->Add(h2); 
    if(doSingleLep){h->Add(hse); h->Add(hsm); h->Add(hemsm); h->Add(hemse);}
  }
  return h;
}
TH1F* H_ttbar(TString var, TString chan, TString level, TString sys){
  TH1F* h = loadHSyst("TTbar_Powheg", var, chan, level, "0"); 
  h->SetFillColor(kRed);
  return h;
}

TH1F* H_ttV(TString var, TString chan, TString level, TString sys){
  TH1F* h  = loadHSyst("TTWToLNu", var, chan, level, "0");
  TH1F* h1 = loadHSyst("TTWToQQ", var, chan, level, "0");
  TH1F* h2 = loadHSyst("TTZToQQ", var, chan, level, "0");
  TH1F* h3 = loadHSyst("TTZToLLNuNu", var, chan, level, "0");
  h->Add(h1); h->Add(h2); h->Add(h3);
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
  TH1F* h1; TH1F* h2; TH1F* h3; TH1F* h4;
  h1   = loadHSyst("DYJetsToLL_M10to50_aMCatNLO",  var, chan, level, "0"); 
  h3   = loadHSyst("DYJetsToLL_M50_aMCatNLO",          var, chan, level, "0"); 
  h1->Add(h3);
  h1->SetFillColor(kYellow);
  return h1;
}

TH1F* H_fake(TString var, TString chan, TString level, TString sys){
  TH1F* h; TH1F* h2;
  h = loadHSyst("WJetsToLNu_aMCatNLO", var, chan, level, "0");
  h->SetFillColor(kCyan-3);
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
  TH1F* httV   = H_ttV(var, chan, level, "0"); 
  TH1F* httbar= H_ttbar(var, chan, level, "0"); 
  hs->Add(hfake); 
  hs->Add(hVV); hs->Add(httV); hs->Add(hDY); hs->Add(htW); hs->Add(httbar); 
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
  inputfile->GetObject("H_" + var + "_" + chan + "_" + level + "_" + syst,hist);
  hist -> SetStats(0);

  if( (sample == "MuonEG") || (sample == "DoubleEG") || (sample == "DoubleMuon") ){
    cout << "No systematics for data!" << endl;
    return hist;
  }
  else  hist->Scale(Lumi*gettrigeff(chan));
  if(var == "MET" || var.BeginsWith("MT2") || (var.BeginsWith("Del") && var.Contains("Phi")) || ( (var == "Jet0Pt") || (var == "Jet1Pt")) || ( (var == "Lep0Pt") || (var == "DiLepPt") || (var == "Lep1Pt")) || (var == "InvMass") || (var == "Yields")){
    hist = RebinVar(hist, var);
  }
  else RebinHist(hist, var);
  hist -> SetDirectory(0);
  delete inputfile;
  return hist;
}

void RebinHist(TH1F* hist, TString var){
/*  if( (var == "DiLepPt") || (var == "Lep0Pt") || (var == "Lep1Pt") || (var == "Jet0Pt") || (var == "Jet1Pt") || (var == "MET") ) hist -> Rebin(25);
  if( (var == "MT") || (var == "Ptllb") || (var == "MT2bb") || (var == "HT")) hist -> Rebin(30);
  if( (var == "MT2lblb")) hist -> Rebin(20);
  if( (var=="DelLepPhi") || (var == "DelPhiLepMet") || (var == "DelPhiJetMet") || (var == "DelPhiLepJet") || (var == "DelPhiPllbMet")) hist -> Rebin(10);
  if( (var == "Meff") ) hist -> Rebin(100);
  if(var == "InvMass") hist-> Rebin(40);
  if((var == "MT2") || (var.BeginsWith("MT2_")) ) hist -> Rebin(10);
  if(var == "MinDPhiMetJets") hist -> Rebin(4);
*/ // if((var == "MT2") || (var == "MT2lblb")) hist->SetBinContent(1, 0);
  
  //if( (var=="DelLepPhi") || (var == "DelPhiLepMet") || (var == "DelPhiJetMet") || (var == "DelPhiLepJet") || (var == "DelPhiPllbMet")) hist -> Rebin(10);
  if(var == "HT") hist->Rebin(60);
//  if(var == "InvMass") hist->Rebin(50);
  if(var == "Meff") hist->Rebin(100);
//  if( (var == "DiLepPt") || (var == "Lep0Pt") || (var == "Lep1Pt") || (var == "Jet0Pt") || (var == "Jet1Pt") ) hist -> Rebin(40);
  if(var == "CosDelLepPhi") hist->Rebin(2);
}

TH1F* RebinVar(TH1F* hist, TString var){
	TH1F* h; Double_t xbins13[13]; Double_t xbins10 [10]; Double_t xbins11[11]; Double_t xbinsim[18]; Double_t xbinslept[18]; Double_t xbinsjetpt[18]; Double_t xbinscdp[17]; Double_t xbinsmet[17];

  /*if(var == "Yields"){
    Double_t xbinsy[9] ={0,1,2,3,4,5,6,7,8};
    h = (TH1F*) hist->Rebin(8, var + "_rebin", xbinsy);
  }*/

  if( (var.BeginsWith("Del") && var.Contains("Phi")) ){
    xbinscdp[0] = 0; xbinscdp[1] = 0.25; xbinscdp[2] = 0.45; xbinscdp[3] = 0.65; xbinscdp[4] = 0.85; xbinscdp[5] = 1.05; xbinscdp[6] = 1.25; xbinscdp[7] = 1.45; xbinscdp[8] = 1.65; xbinscdp[9] = 1.85; xbinscdp[10] = 2.05; xbinscdp[11] = 2.25; xbinscdp[12] = 2.45; xbinscdp[13] = 2.65; xbinscdp[14] = 2.85; xbinscdp[15] = 3.05; xbinscdp[16] = 3.25; h = (TH1F*) hist->Rebin(16, var + "_rebin", xbinscdp);
  }

  else if( (var == "Jet0Pt") || (var == "Jet1Pt")){
    xbinsjetpt[0] = 0; xbinsjetpt[1] = 30; xbinsjetpt[2] = 60; xbinsjetpt[3] = 100; xbinsjetpt[4] = 150; xbinsjetpt[5] = 200; xbinsjetpt[6] = 250; xbinsjetpt[7] = 300; xbinsjetpt[8] = 350; xbinsjetpt[9] = 400; xbinsjetpt[10] = 450; xbinsjetpt[11] = 500; xbinsjetpt[12] = 550; xbinsjetpt[13] = 600; xbinsjetpt[14] = 650; xbinsjetpt[15] = 700; xbinsjetpt[16] = 750; xbinsjetpt[17] = 800; h = (TH1F*) hist->Rebin(17, var + "_rebin", xbinsjetpt);
  }

  else if( (var == "Lep0Pt") || (var == "DiLepPt") || (var == "Lep1Pt")){
    xbinslept[0] = 0; xbinslept[1] = 20; xbinslept[2] = 60; xbinslept[3] = 100; xbinslept[4] = 150; xbinslept[5] = 200; xbinslept[6] = 250; xbinslept[7] = 300; xbinslept[8] = 350; xbinslept[9] = 400; xbinslept[10] = 450; xbinslept[11] = 500; xbinslept[12] = 550; xbinslept[13] = 600; xbinslept[14] = 650; xbinslept[15] = 700; xbinslept[16] = 750; xbinslept[17] = 800; h = (TH1F*) hist->Rebin(17, var + "_rebin", xbinslept);
  }

	else if(var == "InvMass"){xbinsim[0] = 0; xbinsim[1] = 20; xbinsim[2] = 60; xbinsim[3] = 100; xbinsim[4] = 150; xbinsim[5] = 200; xbinsim[6] = 250; xbinsim[7] = 300; xbinsim[8] = 350; xbinsim[9] = 400; xbinsim[10] = 450; xbinsim[11] = 500; xbinsim[12] = 550; xbinsim[13] = 600; xbinsim[14] = 700; xbinsim[15] = 800; xbinsim[16] = 900; xbinsim[17] = 1000; h = (TH1F*) hist->Rebin(17, var + "_rebin", xbinsim);}

	else if(var == "MT2" || var.BeginsWith("MT2_")){xbins10[0] = 0; xbins10[1] = 20; xbins10[2] = 40; xbins10[3] = 60; xbins10[4] = 80; xbins10[5] = 100; xbins10[6] = 140; xbins10[7] = 200; xbins10[8] = 300; xbins10[9] = 400;   h = (TH1F*) hist->Rebin(9, var + "_rebin", xbins10);}

	else if(var == "MT2bb"){xbins13[0] = 0; xbins13[1] = 50; xbins13[2] = 100; xbins13[3] = 150; xbins13[4] = 200; xbins13[5] = 250; xbins13[6] = 300; xbins13[7] = 350; xbins13[8] = 420; xbins13[9] = 500; xbins13[10] = 600; xbins13[11] = 700; xbins13[12] = 800; h = (TH1F*) hist->Rebin(12, var + "_rebin", xbins13);}

  else if(var == "MT2lblb"){ xbins11[0] = 0; xbins11[1] = 50; xbins11[2] = 100; xbins11[3] = 150; xbins11[4] = 200; xbins11[5] = 250; xbins11[6] = 300; xbins11[7] = 350; xbins11[8] = 420; xbins11[9] = 500; xbins11[10] = 600; h = (TH1F*) hist->Rebin(10, var + "_rebin", xbins11);}

	else if(var == "MET"){xbinsmet[0] = 0; xbinsmet[1] = 50; xbinsmet[2] = 80; xbinsmet[3] = 110; xbinsmet[4] = 140; xbinsmet[5] = 170; xbinsmet[6] = 200; xbinsmet[7] = 230; xbinsmet[8] = 260; xbinsmet[9] = 290; xbinsmet[10] = 320; xbinsmet[11] = 350; xbinsmet[12] = 400; xbinsmet[13] = 450; xbinsmet[14] = 500; xbinsmet[15] = 550; xbinsmet[16] = 600; h = (TH1F*) hist->Rebin(16, var + "_rebin", xbinsmet);}
	//else if(var == "MET"){Double_t xbinsmet[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200}; h = (TH1F*) hist->Rebin(40, var + "_rebin", xbinsmet);}

  //int nbins = sizeof(xbins)/sizeof(Double_t)-1;
  //h = (TH1F*) hist->Rebin(nbins, var + "_rebin", xbins);
  return h;
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
	for(int s = 0; s < nSyst; s++){ 
		hsys[s] = BkgSystH(var, chan, level, syst[s]);
	}
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
	T2ttStackPlots(var, chan, level, 0,0,0,0);
}

void DrawSRs(){
	TString thelevel = "DYVeto"; TString thechan[3] = {"Elec", "ElMu", "Muon"};
	for(int k = 0; k < 3; k++){
		for(int s = 0; s < nSR; s++){
			T2ttStackPlots("MT2_" + SRlabel[s], thechan[k], thelevel, 183, 1, 250, 75);
		}
	}
}


void DrawPlots(){
	TString levels[4] = {"dilepton", "2jets", "1btag","DYVeto"};
	TString Vars[26] = {"Vtx", "NBtagsNJets","NJets", "NBtagJets", "MT2", "MT2bb", "MT2lblb", "InvMass", "HT", "Lep0Pt", "Lep1Pt", "Jet0Pt", "Jet1Pt", "MET", "DelLepPhi", "Meff", "DelPhiLepMet", "DelPhiJetMet", "DelPhiLepJet", "DelPhiPllbMet", "CosDelLepPhi", "METHT", "Lep0Eta", "Lep1Eta", "Jet0Eta", "Jet1Eta"};
	TString Channels[5] = {"ElMu", "Elec", "Muon", "sameF",  "All"};
	for(int k = 0; k<4; k++){
		for(int i = 0; i < 26; i++){ //26
			for(int j = 0; j<5; j++){
        //if(Vars[i] != "MET") continue;
				T2ttStackPlots(Vars[i], Channels[j], levels[k], 183, 1, 250, 75);
        if(levels[k] == "dilepton" && (i == 0))  T2ttStackPlots("Yields", Channels[j], levels[k]);
			}
		}
	}
		//		T2ttStackPlots("MT2", "All", "1btag");

//DrawSRs();
}

//void T2ttStackPlots(TString var, TString chan = "All", TString level = "DYVeto", int mStop1 = 0, int mLsp1 = 0, int mStop2 = 0, int mLsp2 = 0){
void T2ttStackPlots(TString var, TString chan, TString level, int mStop1, int mLsp1, int mStop2, int mLsp2){
	TH1F* h_stop_1;
	TH1F* h_stop_2;
	TCanvas* c = SetCanvas();
  plot->cd();
  THStack* hs = TH_bkg(var, chan, level);
  TH1F* hdata; 
  if(dodata)  hdata = H_data(var, chan, level);
  else hdata = (TH1F*) ((TH1F*) hs->GetStack()->Last())->Clone();
  float maxMC = ((TH1F*) hs->GetStack()->Last())->GetMaximum();
  float maxData = 0;
  if(dodata) maxData = hdata->GetMaximum();
  if(doSetLogy && !var.Contains("Del") && !var.Contains("NBtag") && !var.Contains("NJet") && !var.Contains("Eta") && !(var == "Yields" && chan == "ElMu") && (var != "Vtx")){
    plot->SetLogy();  
    hs->SetMinimum(1.);
    if(maxData > maxMC)  hs->SetMaximum(maxData*500); else hs->SetMaximum(maxMC*500);
  }
  else{
    hs->SetMinimum(0);
    if(maxData > maxMC)  hs->SetMaximum(maxData*1.15); else hs->SetMaximum(maxMC*1.15);
  }
  hs->Draw("hist");
  hs->GetYaxis()->SetTitle("Number of Events");
  hs->GetYaxis()->SetTitleSize(0.06);
  hs->GetYaxis()->SetTitleOffset(0.5);
  hs->GetYaxis()->SetNdivisions(505);
  hs->GetXaxis()->SetLabelSize(0.0);
	if (dodata) hdata->Draw("pesame");
	plot->RedrawAxis("same");
	if(mStop1 != 0){
		TString signal1 = Form("T2tt_mStop%i_mLsp%i", mStop1, mLsp1);
		h_stop_1 = loadHistogram(signal1, var, chan, level); 
		h_stop_1->SetFillColor(0); h_stop_1->SetLineColor(kMagenta-4); h_stop_1->SetLineWidth(2);	h_stop_1->SetLineStyle(1); 
	}
	if(mStop2 != 0){
		TString signal2 = Form("T2tt_mStop%i_mLsp%i", mStop2, mLsp2);
		h_stop_2 = loadHistogram(signal2, var, chan, level);
		h_stop_2->SetFillColor(0); h_stop_2->SetLineColor(kBlue+3); h_stop_2->SetLineWidth(2); h_stop_2->SetLineStyle(1); 
	}
	if(mStop1 == 0 && mLsp1 != 0){ 
		TString signal1 = Form("TTDMJets_EFT_M%i", mLsp1);
		h_stop_1 = loadHistogram(signal1, var, chan, level); 
		h_stop_1->SetFillColor(0); h_stop_1->SetLineColor(kMagenta-4); h_stop_1->SetLineWidth(2);	h_stop_1->SetLineStyle(1); 
	}
	if(mStop2 == 0 && mLsp2 != 0){
		TString signal2 = Form("TTDMJets_EFT_M%i", mLsp2);
		h_stop_2 = loadHistogram(signal2, var, chan, level);
		h_stop_2->SetFillColor(0); h_stop_2->SetLineColor(kBlue+3); h_stop_2->SetLineWidth(2); h_stop_2->SetLineStyle(1); 
	}

	SetLegend();
	if(doYields){
		leg->AddEntry(H_ttbar(var,chan,level, "0"), Form("ttbar: %5.2f", yield("ttbar", chan, level)), "f");
		leg->AddEntry(H_tW(var,chan,level, "0"), Form("tW: %5.2f", yield("tW", chan, level)), "f");
		leg->AddEntry(H_DY(var,chan,level, "0"), Form("DY: %5.2f", yield("DY", chan, level)), "f");
		leg->AddEntry(H_VV(var,chan,level, "0"), Form("VV: %5.2f", yield("VV", chan, level)), "f");
		leg->AddEntry(H_ttV(var,chan,level, "0"), Form("ttV: %5.2f", yield("ttV", chan, level)), "f");
		leg->AddEntry(H_fake(var,chan,level, "0"), Form("W+Jets: %5.2f", yield("fake", chan, level)), "f");
		TH1F* hh = H_ttbar(var,chan,level, "0"); hh->SetFillStyle(1001); hh->SetFillColor(kBlack);
		leg->AddEntry(hh, Form("All bkg: %5.2f", yield("bkg", chan, level)), "");
		if (dodata) leg->AddEntry(hdata, Form("Data: %5.0f",  yield("Data", chan, level)), "pl");
		if(mStop1 != 0)               leg->AddEntry(h_stop_1, Form("S.%i-%i : %5.2f", mStop1, mLsp1, h_stop_1->Integral()), "l");
		if(mStop2 != 0)               leg->AddEntry(h_stop_2, Form("S.%i-%i : %5.2f", mStop2, mLsp2, h_stop_2->Integral()), "l");
    if(mStop1 == 0 && mLsp1 != 0) leg->AddEntry(h_stop_1, Form("DM.%i: %5.2f", mLsp1, h_stop_1->Integral()), "l");
    if(mStop2 == 0 && mLsp2 != 0) leg->AddEntry(h_stop_2, Form("DM.%i: %5.2f", mLsp1, h_stop_2->Integral()), "l");
  }
  else{
		leg->AddEntry(H_ttbar(var,chan,level, "0"), "ttbar", "f");
		leg->AddEntry(H_tW(var,chan,level, "0"), "tW", "f");
		leg->AddEntry(H_DY(var,chan,level, "0"), "DY", "f");
		leg->AddEntry(H_VV(var,chan,level, "0"), "VV", "f");
		leg->AddEntry(H_ttV(var,chan,level, "0"), "ttV", "f");
		leg->AddEntry(H_fake(var,chan,level, "0"), "W+Jets", "f");
		if (dodata) leg->AddEntry(hdata, "Data", "pl");
		if(mStop1 != 0) leg->AddEntry(h_stop_1, Form("S.%i-%i ", mStop1, mLsp1), "l");
		if(mStop2 != 0) leg->AddEntry(h_stop_2, Form("S.%i-%i ", mStop2, mLsp2), "l");
    if(mStop1 == 0 && mLsp1 != 0) leg->AddEntry(h_stop_1, Form("DM.%i ", mLsp1), "l");
    if(mStop2 == 0 && mLsp1 != 0) leg->AddEntry(h_stop_2, Form("DM.%i ", mLsp1), "l");
  }
  leg->Draw("same"); texlumi->Draw("same"); texcms->Draw("same");
  SetTexChan(chan, level); texchan->Draw("same");
  if(mStop1 != 0){ h_stop_1->Scale(1);  h_stop_1->Draw("same,hist");}
  if(mStop2 != 0){ h_stop_2->Scale(1);  h_stop_2->Draw("same,hist");}
	if(mStop1 == 0 && mLsp1 != 0) h_stop_1->Draw("same,hist");
	if(mStop2 == 0 && mLsp2 != 0) h_stop_2->Draw("same,hist");

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
	//if(dodata){
		hratio = (TH1F*)hdata->Clone();
		hratio->Divide(mh);
		SetHRatio(var);
		hratio->Draw("same");
    if(doSys) hratioerror->Draw("same,e2");
    hratio->Draw("same");
  //}
  //TString dir = plotfolder + chan + "/" + level + "/";
  TString dir = plotfolder + level + "/";
  TString plotname = var + "_" + chan + "_" + level;
  gSystem->mkdir(dir, kTRUE);
  c->Print( dir + plotname + ".pdf", "pdf");
  c->Print( dir + plotname + ".png", "png");
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

