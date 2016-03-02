TString path = "/nfs/fanae/user/palencia/testHeppy/TOP/TopTrees/mar1/";

TString Sample[] = {"TTJets", "TTJets_aMCatNLO", "TTbar_Powheg_Herwig", "TTbar_Powheg_ScaleDown", "TTbar_Powheg_ScaleUp"};
TString Channel[] = {"Elec", "Muon", "ElMu"};
TString level = "1btag";
float Lumi = 2223;
TString Levels[] = {"dilepton", "2jets", "1btag"};
TString level = "1btag";

double getYield(TString sample, TString chan, TString level){
  TFile *f = TFile::Open(path + "Tree_13TeV_EA_" + sample + ".root");
  TH1F *h; f->GetObject("H_Yields_" + chan, h);
  int ilevel = 0;
  if      (level == "dilepton") ilevel = 1;
  else if (level == "2jets")    ilevel = 4;
  else if (level == "1btag")    ilevel = 5;
  double yield = h->GetBinContent(ilevel);
  delete h; f->Close(); delete f;
  if(sample.Contains("aMCatNLO")) return yield*Lumi*(1.08*0.9)*(1.08*0.9);
  return yield*Lumi;
}

double getYieldSys(TString sample, TString chan, TString level, TString sys){
  TFile *f = TFile::Open(path + "Tree_13TeV_EA_" + sample + ".root");
  TH1F *h; f->GetObject("H_Yields_" + chan + "_" + sys, h);
  int ilevel = 0;
  if      (level == "dilepton") ilevel = 1;
  else if (level == "2jets")    ilevel = 4;
  else if (level == "1btag")    ilevel = 5;
  double yield = h->GetBinContent(ilevel);
  delete h; f->Close(); delete f;
  return yield*Lumi;
}

double getNLO(TString chan){
  TString level = "1btag";
  double yPowheg = getYield("TTbar_Powheg", chan, level);	
  double yNLO    = getYield("TTJets_aMCatNLO", chan, level);	
  return TMath::Abs(yPowheg-yNLO); 
}

double getPS(TString chan){
  TString level = "1btag";
  double yPowheg = getYield("TTbar_Powheg", chan, level);	
  double yUp     = getYield("TTbar_Powheg_scaleUp", chan, level);	
  double yDown   = getYield("TTbar_Powheg_scaleDown_ext", chan, level);	
  return TMath::Max(TMath::Abs(yPowheg - yUp), TMath::Abs(yPowheg - yDown));
}

double getHad(TString chan){
  TString level = "1btag";
  double yPythia = getYield("TTbar_Powheg", chan, level);
  double yHerwig = getYield("TTbar_Powheg_Herwig", chan, level);
  return TMath::Abs(yPythia-yHerwig);
}

double getTrig(TString chan){
  TString level = "1btag";
  double yPowheg = getYield("TTbar_Powheg", chan, level);
  double yUp     = getYield("TTbar_Powheg_TrigUp", chan, level);
  double yDown   = getYield("TTbar_Powheg_TrigDown", chan, level);
  return TMath::Max(TMath::Abs(yPowheg - yUp), TMath::Abs(yPowheg - yDown));
}

double getMuon(TString chan){
  double yPowheg = getYield("TTbar_Powheg", chan, level);
  double yUp     = getYield("TTbar_Powheg_MuonUp", chan, level);
  double yDown   = getYield("TTbar_Powheg_MuonDown", chan, level);
  return TMath::Max(TMath::Abs(yPowheg - yUp), TMath::Abs(yPowheg - yDown));
}

double getElec(TString chan){
  double yPowheg = getYield("TTbar_Powheg", chan, level);
  double yUp     = getYield("TTbar_Powheg_ElecUp", chan, level);
  double yDown   = getYield("TTbar_Powheg_ElecDown", chan, level);
  return TMath::Max(TMath::Abs(yPowheg - yUp), TMath::Abs(yPowheg - yDown));
}

double getJES(TString chan){
  double yPowheg = getYield("TTbar_Powheg", chan, level);
  double yUp     = getYieldSys("TTbar_Powheg", chan, level, "JESUp");
  double yDown   = getYieldSys("TTbar_Powheg", chan, level, "JESDown");
  return TMath::Max(TMath::Abs(yPowheg - yUp), TMath::Abs(yPowheg - yDown));
}

void PrintYields(){
	TString chan = "ElMu";
	cout << " Channel = ElMu,  Lumi = 2223 fb-1 " << endl;
  cout << "================================================================================" << endl;
  cout << "                                dilepton     ||      2jets     ||     1btag     ||" << endl; 
  cout << "--------------------------------------------------------------------------------" << endl;
	TString nominal    = " TTbar_Powheg                ";
	TString amcatnlo   = " TTJets_aMCatNLO             ";
	TString amcatnom   = " TTJets_aMCatNLO             ";
	TString scaleup    = " TTbar_Powheg_scaleUp        ";
	TString scaledown  = " TTbar_Powheg_scaleDown      ";
  TString hadron     = " TTbar_Powheg_Herwig       ";
  double nom = 0; double amc = 0; double up = 0; double down = 0; double had = 0; double nom2 = 0;
  for(int i = 0; i<3; i++){
    chan = "ElMu";
    nom  = getYield("TTbar_Powheg", chan, Levels[i]);
    had  = getYield("TTbar_Powheg_Herwig", chan, Levels[i]);
    amc  = getYield("TTJets_aMCatNLO", chan, Levels[i]);
    up   = getYield("TTbar_Powheg_scaleUp", chan, Levels[i]);
    down = getYield("TTbar_Powheg_scaleDown_ext", chan, Levels[i]);
    nominal   += Form("%.1f (%.1f %) || ", nom, 0.0);
    amcatnlo  += Form("%.1f (%.1f %) || ", amc, abs(nom-amc)/nom*100); 
    amcatnom  += Form("%.1f (%.1f %) || ", amc, 0.0); 
    scaleup   += Form("%.1f (%.1f %) || ", up, abs(nom-up)/nom*100); 
    scaledown += Form("%.1f (%.1f %) || ", down, abs(nom-down)/nom*100); 
    hadron    += Form("%.1f (%.1f %) || ", had, abs(nom-had)/nom*100); 
	}
  cout << nominal << endl;
  cout << amcatnlo << endl;
  cout << scaleup << endl;
  cout << scaledown << endl;
  cout << hadron << endl;
  cout << "=================================================================================" << endl;
}

void getOtherSystematics(){
  Double_t normYield; Double_t normYieldamc;
  TString chan = "";
  TString nlo = " NLO Generator     ";
  TString had = " Hadronization     ";
  TString SPS = " Scale PS          ";
  TString muo = " Muon              ";
  TString ele = " Electron          ";
  TString tri = " Trigger           ";
  TString jes = " Jet Energy Scale  ";
  for(int i = 0; i<3; i++){
    chan = Channel[i];
    normYield = getYield("TTbar_Powheg", chan, "1btag"); 
    normYieldamc = getYield("TTJets_aMCatNLO", chan, "1btag"); 
    nlo += Form("%3.1f (%1.2f %)  ||  ", getNLO(chan), getNLO(chan)*100/normYield);
    had += Form("%3.1f (%1.2f %)  ||  ", getHad(chan), getHad(chan)*100/normYield);
    SPS += Form("%3.1f (%1.2f %)  ||  ", getPS(chan) , getPS(chan) *100/normYield);
    muo += Form("%3.1f (%1.2f %)  ||  ", getMuon(chan) , getMuon(chan) *100/normYield);
    ele += Form("%3.1f (%1.2f %)  ||  ", getElec(chan) , getElec(chan) *100/normYield);
    tri += Form("%3.1f (%1.2f %)  ||  ", getTrig(chan) , getTrig(chan) *100/normYield);
    jes += Form("%3.1f (%1.2f %)  ||  ", getJES(chan) , getJES(chan) *100/normYield);
  }


  cout << "=============================================================================" << endl;
  cout << "                       Elec        ||       Muon       ||       ElMu       ||" << endl;
  cout << "-----------------------------------------------------------------------------" << endl;
  cout << nlo << endl;
  cout << had << endl;
  cout << SPS << endl;
  cout << muo << endl;
  cout << ele << endl;
  cout << tri << endl;
  cout << jes << endl;
  cout << "=============================================================================" << endl;

  cout << endl;
  PrintYields();
}

