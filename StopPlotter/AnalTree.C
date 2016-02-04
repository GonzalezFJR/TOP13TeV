//////////////////////////////////////////////////////////////////////////
// Use the definitions in SetPlotter.C
// Reads the trees, perform a quick analysis and plot the histograms
// The histograms are saved in outputPlots
//
// Ussage:
// AnalAll(TString chan, float metcut, int nbtags);
// DrawAnal(int ivar);
//
//////////////////////////////////////////////////////////////////////////

#include "SetPlotter.C"

TH1F* histo[nSamples][nVars];
TH1F* pttbar;
TH1F* pVV;
TH1F* pData;
TH1F* ptW;
TH1F* pDY;
TH1F* pttV;
TH1F* pWJets;
TH1F* pS500;
TH1F* pS850;
THStack* pStack;

void InitHistos();
void AnalTree(int isample, TString chan, float metcut, int nbtags);
void AnalAll(TString chan, float metcut, int nbtags);
void ConfigHistos(int ivar);
void DrawAnal(int ivar);
float getYield(TH1F* hist);


// Constants:
TString outputPlots = "/nfs/fanae/user/juanr/StopTOP/plots/TreePlots/";

void InitHistos(){
  for(int isample = 0; isample<nSamples; isample ++){
    float nbins = 100; float xmin = 0; float xmax = 100;
    for(int ind = 0; ind < nVars; ind++){
      if( (var[ind] == "MET") || (var[ind] == "DilepPt") || (var[ind] == "Lep0Pt") || (var[ind] == "Lep1Pt") || (var[ind] == "Jet0Pt") || (var[ind] == "Jet1Pt") || (var[ind] == "MT2lblb")) { nbins  = 40; xmin = 0; xmax = 600;}
      if( (var[ind] == "Ptllb") || (var[ind] == "MT2bb") ){ nbins = 40; xmin = 0.; xmax = 800;}
      if(var[ind] == "MT2"){ nbins = 30; xmin = 0.; xmax = 300;}
      if(var[ind] == "InvMass"){ nbins = 40; xmin = 0; xmax = 1000;}
      if(var[ind] == "NJets"){ nbins = 8; xmin = -0.5; xmax = 7.5;}
      if(var[ind] == "NBtagJets"){ nbins = 4; xmin = -0.5; xmax = 3.5;}
      if(var[ind] == "NBtagNJets"){ nbins = 12; xmin = -0.5; xmax = 11.5;}
      if(var[ind] == "Meff"){ nbins = 40; xmin = 0; xmax = 2000;}
      if(var[ind] == "DelLepPhi"){nbins = 20; xmin = -3.2; xmax = 3.2;}
      if( (var[ind] == "DelPhiLepMet") || (var[ind] == "DelPhiJetMet") || (var[ind] == "DelPhiLepJet") || (var[ind] == "DelPhiPllbMet") || (var[ind] == "MinDPhiMetJets") ){nbins = 20; xmin = 0; xmax = 1;}
      if(var[ind] == "METHT"){ nbins = 30, xmin = 0; xmax = 60;}
      histo[isample][ind] = new TH1F(sample[isample] + "_" + var[ind], var[ind], nbins, xmin, xmax);
      histo[isample][ind] -> SetStats(0);
      histo[isample][ind] -> SetLineStyle(0);
    }
  }
}

void AnalTree(int isample, TString chan, float metcut, int nbs){
  if( (chan == "Elec") && ((isample == 13) || (isample == 14)) ) return;
  if( (chan == "ElMu") && ((isample == 12) || (isample == 14)) ) return;
  if( (chan == "Muon") && ((isample == 12) || (isample == 13)) ) return;
  if( ((chan == "sameF") || (chan == "SF")) && (isample == 13) ) return;
  // 12, 13, 14  "DoubleEG", "MuonEG", "DoubleMuon",
  float metCut = metcut; int nBtags = nbs; int ischan = 0;
  float mt2llcut = 0; float mt2bbcut = 0; float mt2lblbcut = 0;
  float ptllbcut = 0; float mindphimetjetscut = 0; float methtcut = 0;

  TTree* t;
  TFile* inputfile = TFile::Open(path + "/Tree_" + sample[isample] + ".root");
  inputfile->GetObject("sTopTree",t);

  float met; int njets; int nbtags;
  int iselel; int ismumu; int iselmu; float weight;
  float mt2ll; float mt2bb; float mt2lblb;
  float meff; float ptllb; float ht; float mll;
  float dphiptllbmet; float dphijetmet; float dphilepmet; float dphilepjet; float dphill;
  float metht; float mindphimetjets; 

  t->SetBranchAddress("TMET", &met);
  t->SetBranchAddress("TWeight", &weight);
  t->SetBranchAddress("TNJets", &njets);
  t->SetBranchAddress("TNJetsBtag", &nbtags);
  t->SetBranchAddress("TIsDoubleMuon", &ismumu);
  t->SetBranchAddress("TIsDoubleElec", &iselel);
  t->SetBranchAddress("TIsElMu", &iselmu);
  t->SetBranchAddress("TMT2ll", &mt2ll);
  t->SetBranchAddress("TMT2bb", &mt2bb);
  t->SetBranchAddress("TMT2lblb", &mt2lblb);
  t->SetBranchAddress("TMll", &mll);
  t->SetBranchAddress("THT", &ht);
  t->SetBranchAddress("TMeff", &meff);
  t->SetBranchAddress("TPtllb", &ptllb);
  t->SetBranchAddress("TdPhiLepMet", &dphilepmet);
  t->SetBranchAddress("TdPhiJetMet", &dphijetmet);
  t->SetBranchAddress("TdPhiLepJet", &dphilepjet);
  t->SetBranchAddress("TdPhiPtllbMET", &dphiptllbmet);
  t->SetBranchAddress("TMinDPhiMetJets", &mindphimetjets);
  t->SetBranchAddress("TdPhill", &dphill);

  metht = met/TMath::Sqrt(ht);
  for(int k = 0; k<t->GetEntries(); k++){
    t->GetEntry(k);
    if      (chan == "ElMu") ischan = iselmu;
    else if (chan == "Elec") ischan = iselel;
    else if (chan == "Muon") ischan = ismumu;
    else if ((chan == "SF") || (chan == "sameF"))  ischan = (iselel + ismumu);
    else if (chan == "All") ischan = 3;
    if( ischan                &&  
        njets >= 2            &&
        met >= metCut          &&
        nbtags >= nBtags       &&
        mt2ll >= mt2llcut      &&
        mt2bb >= mt2bbcut      &&
				(mt2lblb >= mt2lblbcut  || mt2lblb == 0) &&
				metht >= methtcut       &&
				ptllb>= ptllbcut        &&
				mindphimetjets >= mindphimetjetscut
			){
        for(int ind = 0; ind < nVars; ind++){

if(var[ind] == "MET"          )  histo[isample][ind] -> Fill(met, weight);
if(var[ind] == "MT2lblb"      )  histo[isample][ind] -> Fill(mt2lblb, weight);
if(var[ind] == "Ptllb"        )  histo[isample][ind] -> Fill(ptllb, weight);
if(var[ind] == "MT2bb"        )  histo[isample][ind] -> Fill(mt2bb, weight);
if(var[ind] == "MT2"          )  histo[isample][ind] -> Fill(mt2ll, weight);
if(var[ind] == "InvMass"      )  histo[isample][ind] -> Fill(mll, weight);
if(var[ind] == "NJets"        )  histo[isample][ind] -> Fill(njets, weight);
if(var[ind] == "NBtagJets"    )  histo[isample][ind] -> Fill(nbtags, weight);
if(var[ind] == "Meff"         )  histo[isample][ind] -> Fill(meff, weight);
if(var[ind] == "DelLepPhi"    )  histo[isample][ind] -> Fill(dphill, weight);
if(var[ind] == "DelPhiLepMet" )  histo[isample][ind] -> Fill(dphilepmet, weight);
if(var[ind] == "DelPhiJetMet" )  histo[isample][ind] -> Fill(dphijetmet, weight);
if(var[ind] == "DelPhiLepJet" )  histo[isample][ind] -> Fill(dphilepjet, weight);
if(var[ind] == "DelPhiPllbMet")  histo[isample][ind] -> Fill(dphiptllbmet, weight);

/*
if(var[ind] == "NBtagNJets"   )  histo[isample][ind] -> Fill(x, weight);
if(var[ind] == "DilepPt"      )  histo[isample][ind] -> Fill(x, weight);
if(var[ind] == "Lep0Pt"       )  histo[isample][ind] -> Fill(x, weight);
if(var[ind] == "Lep1Pt"       )  histo[isample][ind] -> Fill(x, weight);
if(var[ind] == "Jet0Pt"       )  histo[isample][ind] -> Fill(x, weight);
if(var[ind] == "Jet1Pt"       )  histo[isample][ind] -> Fill(x, weight);
*/
       }
    }
  }
  delete t;
  delete inputfile;
}

void AnalAll(TString chan, float metcut, int nbtags){
  InitHistos();
  for(int isample = 0; isample < nSamples; isample++)  AnalTree(isample, chan, metcut, nbtags);
  cout << "All histograms available in TH1F* histo[isample][ivar] " << endl;
}

void ConfigHistos(TString rvar){
// 0           "TTJets",
// 1, 2        "TW","TbarW",
// 3, 4        "DYJetsToLL_M10to50_aMCatNLO", "DYJetsToLL_M50_aMCatNLO",
// 5, 6, 7     "WW","WZ","ZZ",
// 8, 9, 10    "TTWToLNu", "TTWToQQ", "TTZToQQ",
// 11,         "WJetsToLNu_aMCatNLO",
// 12, 13, 14  "DoubleEG", "MuonEG", "DoubleMuon",
// 15          "T2tt_500_325",
// 16          "T2tt_850_100"
  int ivar = 0; for(int i = 0; i<nVars; i++){ if(rvar == var[i]) ivar = i;}

  pttbar = new TH1F(*histo[0][ivar]);
  pttbar->SetFillColor(kRed);
  pttbar->SetFillStyle(1001); pttbar->Scale(Lumi); 

  ptW = histo[1][ivar]; ptW->Add(histo[2][ivar]);
  ptW->SetFillColor(kViolet-3);
  ptW->SetFillStyle(1001); ptW->Scale(Lumi);

  pDY = histo[3][ivar]; pDY->Add(histo[4][ivar]);
  pDY->SetFillColor(kYellow+1);
  pDY->SetFillStyle(1001); pDY->Scale(Lumi);

  pVV = histo[5][ivar]; pVV->Add(histo[6][ivar]); pVV->Add(histo[7][ivar]);
  pVV->SetFillColor(kGreen-7);
  pVV->SetFillStyle(1001); pVV->Scale(Lumi);

  pttV = histo[8][ivar]; pttV->Add(histo[9][ivar]); pttV->Add(histo[10][ivar]);
  pttV->SetFillColor(kOrange+2);
  pttV->SetFillStyle(1001); pttV->Scale(Lumi);

  pWJets = histo[11][ivar];
  pWJets->SetFillColor(kGreen-9);
  pWJets->SetFillStyle(1001); pWJets->Scale(Lumi);

  pData = histo[12][ivar]; pData->Add(histo[13][ivar]); pData->Add(histo[14][ivar]);
  pData->SetLineColor(kBlack); pData->SetMarkerStyle(20); pData->SetMarkerSize(1.1);

  pS500 = histo[15][ivar]; pS500->SetFillColor(0); pS500->SetLineColor(kBlue+3);    pS500->SetLineWidth(2.5); pS500->SetLineStyle(1);
  pS850 = histo[16][ivar]; pS850->SetFillColor(0); pS850->SetLineColor(kMagenta-4); pS850->SetLineWidth(2.5); pS850->SetLineStyle(1);
  pS500->Scale(Lumi);  pS850->Scale(Lumi);

  pStack = new THStack(rvar, "");
  pStack->Add(pWJets); pStack->Add(pVV); pStack->Add(pttV);  pStack->Add(pDY);  pStack->Add(ptW);  pStack->Add(pttbar);
}

void DrawAnal(TString rvar){
  //AnalAll("All", 80, 1);
  ConfigHistos(rvar);

  TCanvas* c = SetCanvas();
  plot->cd();
  plot->SetLogy();
  THStack* hs = (THStack*)pStack->Clone();
  TH1F* hdata = (TH1F*)pData->Clone();
  hs->SetMinimum(10e-3);
  hs->SetMaximum(hdata->GetMaximum()*500);
  hs->Draw("hist");
  hs->GetYaxis()->SetTitle("Number of Events");
  hs->GetYaxis()->SetTitleSize(0.06);
  hs->GetYaxis()->SetTitleOffset(0.9);
  hs->GetYaxis()->SetNdivisions(505);
  hs->GetXaxis()->SetLabelSize(0.0);
  hdata->Draw("pesame");
  plot->RedrawAxis("same");
  TH1F* h_stop_1 = (TH1F*)pS500->Clone();
  TH1F* h_stop_2 = (TH1F*)pS850->Clone();

  SetLegend();
  leg->AddEntry(pttbar,   Form("ttbar: %5.2f" , getYield(pttbar)), "f");
  leg->AddEntry(ptW,      Form("tW: %5.2f"    , getYield(ptW)),    "f");
  leg->AddEntry(pDY,      Form("DY: %5.2f"    , getYield(pDY)),    "f");
  leg->AddEntry(pVV,      Form("VV: %5.2f"    , getYield(pVV)),    "f");
  leg->AddEntry(pttV,     Form("ttV: %5.2f"   , getYield(pttV)),   "f");
  leg->AddEntry(pWJets,   Form("W+Jets: %5.2f", getYield(pWJets)), "f");
  leg->AddEntry(hdata,    Form("Data: %5.2f"  , getYield(hdata)),  "pl");
  leg->AddEntry(h_stop_1, Form("S.500-325 (x10 ): %5.2f", getYield(h_stop_1)), "l");
  leg->AddEntry(h_stop_2, Form("S.850-100 (x100): %5.2f", getYield(h_stop_2)), "l");

  leg->Draw("same"); texlumi->Draw("same"); texcms->Draw("same");
  SetTexChan("ElMu", "1btag"); texchan->Draw("same");
  h_stop_1->Scale(10);  
  h_stop_1->Draw("same,hist");
  h_stop_2->Scale(100); 
  h_stop_2->Draw("same,hist");

  pratio->cd();
  TH1F* mh = (TH1F*)hs->GetStack()->Last()->Clone();
  hratio = (TH1F*)hdata->Clone();
  hratio->Divide(mh);
  SetHRatio(rvar);
  hratio->Draw("pesame");

  c->Print(outputPlots+ "/" + rvar + ".pdf", "pdf");
  c->Print(outputPlots+ "/" + rvar + ".png", "png");
  delete c;
}

float getYield(TH1F* hist){
  float I = hist->Integral();
  I+= hist->GetBinContent(0);
  I+= hist->GetBinContent(hist->GetNbinsX()+1);
  return I;
}

