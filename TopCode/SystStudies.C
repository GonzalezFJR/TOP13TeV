enum gChannel{
  channels_begin,
  Muon = channels_begin,
  Elec,
  ElMu,
  gNCHANNELS,
};  
TString gChanLabel[gNCHANNELS+1] = {"Muon","Elec","ElMu","EEMM"};
enum iCut{
  iDilepton, 
  iZVeto, 
  iMET, 
  i2jets, 
  i1btag, 
  iNCUTS
  //i2btag,
};
void SystStudies(TString pathtofiles = "TopTrees/ReReco_CSVM_METType0I_No3rdLepV_FullSyst_Mar17/"){
  
  TFile *_file[7];
  _file[0] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin_jet20.root");
  _file[1] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin.root");
  _file[2] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin_jet35.root");
  _file[3] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin_jet40.root");
  _file[4] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin_jet45.root");
  _file[5] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin_jet50.root");
  _file[6] = new TFile(pathtofiles + "/Tree_Legacy_TTJets_MadSpin_jet60.root");
  
  Float_t JES[3][7];
  TH1F *nom;
  TH1F *up;
  TH1F *down;
  Float_t ynom, yup, ydown;
  for (UInt_t i=0; i<7; i++){
    for (UInt_t ch=0; ch<gNCHANNELS; ch++){
      nom  = (TH1F*) _file[i]->Get("H_Yields_"+gChanLabel[ch]);
      up   = (TH1F*) _file[i]->Get("H_Yields_"+gChanLabel[ch]+"_JESUp");
      down = (TH1F*) _file[i]->Get("H_Yields_"+gChanLabel[ch]+"_JESDown");
      
      ynom  = nom ->GetBinContent(i1btag+1);
      yup   = up  ->GetBinContent(i1btag+1);
      ydown = down->GetBinContent(i1btag+1);
      
      cout << " up:   " << yup;
      cout << " nom:  " << ynom;
      cout << " down: " << ydown;
      cout << endl;
      
      JES[ch][i] = 100 * TMath::Max(TMath::Abs(yup - ynom), TMath::Abs(ydown - ynom)) / ynom;
    }
  }
  
  Float_t JetEt[7] = {20., 30., 35., 40., 45., 50., 60.};
  
  TGraph *JESvsJetEt[gNCHANNELS];
  Color_t col[gNCHANNELS] = {kRed+1,kBlue-3,kGreen+3};

  for (UInt_t ch=0; ch<gNCHANNELS; ch++){
    JESvsJetEt[ch] = new TGraph(7,JetEt,JES[ch]);
    JESvsJetEt[ch]->SetTitle("JES vs Jet E_{T}");

    
    JESvsJetEt[ch]->GetYaxis()->SetTitle("JES Unc (%)");
    //    JESvsJetEt[ch]->GetYaxis()->SetTitleOffset(1.1);
    //    JESvsJetEt[ch]->GetYaxis()->SetTitleSize(0.07);
    //    JESvsJetEt[ch]->GetYaxis()->SetLabelSize(0.055);
    JESvsJetEt[ch]->GetYaxis()->SetNdivisions(607);
    JESvsJetEt[ch]->GetYaxis()->SetRangeUser(1.,4.);
    //    JESvsJetEt[ch]->GetXaxis()->SetLabelSize(0.05);
    JESvsJetEt[ch]->GetXaxis()->SetTitle("Jet E_{T} (GeV)");
    
    JESvsJetEt[ch]->SetMarkerStyle(20);
    JESvsJetEt[ch]->SetMarkerColor(col[ch]);
    
  }
  
  
  TCanvas *c1 = new TCanvas("c1","Plot");
  c1->cd();
  
  TLegend *leg = new TLegend(0.73,0.20,0.87,0.50);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);

  leg->AddEntry(JESvsJetEt[Muon], gChanLabel[Muon],"P");
  leg->AddEntry(JESvsJetEt[Elec], gChanLabel[Elec],"P");
  leg->AddEntry(JESvsJetEt[ElMu], gChanLabel[ElMu],"P");
  
  JESvsJetEt[Muon]->Draw("AP");
  JESvsJetEt[Elec]->Draw("P SAME");
  JESvsJetEt[ElMu]->Draw("P SAME");
  
  leg->Draw("SAME");
  
}
