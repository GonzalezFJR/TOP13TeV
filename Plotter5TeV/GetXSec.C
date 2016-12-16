#include "DrawPlots.C"
#include "GetPDFweights.C"
#include "GetDataDriven.C"

float getSyst(TString tag, TString direc = "Max"); // Returns the variation in the yield due to a systematic uncertainty
float getExpError(TString Process = "ttbar", TString tag = "0", TString direc = "Max"); // Returns the variation in the ttbat yield due to a experimental uncertainty
float getAcc(bool docout = false);
float getXsec(TString tag = "0", TString direc = "Max");

float getNLOsyst(bool isFidu = false);
float getHadSyst(bool isFidu = false); 
float getScalePSSyst(bool isFidu = false, TString direc = "Max"); 
void PrintXsecSyst();
void PrintGoodYields();
void GetXSec();
float getFiducialNorm(TString SampleName = "TTbar_PowhegFiduElMu");
float GetFiducialXSec(bool docout = false);
void PrintXSec();
float getMCstatError();

///////// float DYDD(bool docout = false);
///////// float NonWDD(bool docout = false)

int   doDD = 0;
//const float lumiunc = 0.12;
const float lumiunc = 0.023;
const float xsec_th = 68.9;

const float tWerr = 0.3; const float VVerr = 0.3; const float DYerr = 1.0; const float NonWerr = 0.3;
const TString up = "Up"; const TString down = "Down"; const TString maxi = "Max";

float getSyst(TString tag, TString direc){ // Returns the variation in the yield due to a systematic uncertainty
  // Modeling
  if     (tag == "PDF" || tag == "pdf" || tag == "PDFFidu" || tag == "pdfFidu")          return NNPDFsyst(); // PDF
  else if(tag == "ME"  || tag == "MEFidu" || tag == "scaleME" || tag == "ScalesME")        return MEvariations(); // Scale ME
  else if(tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales")  return getScalePSSyst(0,direc); // Scale PS
	else if(tag == "had" || tag == "Had")                                       return getHadSyst(); // Hadronization
  else if(tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator") return getNLOsyst(); // NLO
  else if(tag == "QCDFidu" || tag == "qcdFidu" || tag == "scaleFidu" || tag == "scalesFidu")  return getScalePSSyst(1, direc); // Scale PS
	else if(tag == "hadFidu" || tag == "HadFidu")                                               return getHadSyst(1); // Hadronization
  else if(tag == "MCstat") return getMCstatError()*yield("ttbar", Chan, Level);
  // Backgrounds
  else if(tag == "tW")                                                          return tWerr*yield("tW", Chan, Level); // tW
  else if(tag == "VV")                                                          return VVerr*yield("VV", Chan, Level); // VV
  else if(doDD && tag == "DY"){
    float staterr = GetDYDD(Chan, Level, 1, 0)/GetDYDD(Chan, Level, 0, 0);
    return TMath::Sqrt(DYerr*DYerr + staterr*staterr)*yield("DY", Chan, Level);
  }
  else if(doDD && (tag == "NonW/Z" || tag == "fake" || tag == "NonW" || tag == "NonWZ")){
    float staterr = GetNonWDD(Chan, Level, 1, 0)/GetNonWDD(Chan, Level, 0, 0);
		return TMath::Sqrt(NonWerr*NonWerr + staterr*staterr)*yield("fake", Chan, Level);
	}
	else if(tag == "DY")                                                          return DYerr*yield("DY", Chan, Level); // DY
  else if(tag == "NonW/Z" || tag == "fake" || tag == "NonW" || tag == "NonWZ")  return NonWerr*yield("fake", Chan, Level); // NonW
  // Experimental
	else if(tag.Contains("trig") || tag.Contains("Trig") || tag.Contains("Elec") || tag.Contains("elec") || tag.Contains("Muon") || tag.Contains("muon") || tag == "JES" || tag == "jes" || tag == "JER" || tag == "jer") return (float )getExpError("ttbar", tag, direc);
	else{ cout << "Systematic " << tag << " not found!" << endl; return 0.;}
}

float getExpError(TString Process, TString tag, TString direc){ // Returns the variation in the ttbat yield due to a experimental uncertainty
  float ynom = yield(Process, Chan, Level, "0"); 
  TString tup = "0"; TString tdown = "0";
  if     (tag.Contains("trig") || tag.Contains("Trig")) { tup = "TrigUp"; tdown = "TrigDown";}
  else if(tag.Contains("elec") || tag.Contains("Elec")) { tup = "ElecUp"; tdown = "ElecDown";}
  else if(tag.Contains("muon") || tag.Contains("Muon")) { tup = "MuonUp"; tdown = "MuonDown";}
  else if(tag == "JES" || tag == "jes")                 { tup = "JESUp" ; tdown = "JESDown" ;}
  else if(tag == "JER" || tag == "jer")                 { tup = "JER"   ; tdown = "JER"     ;}
  float yup   = yield(Process, Chan, Level, tup);  
  float ydown = yield(Process, Chan, Level, tdown);  
  if     (direc == maxi ) return (float )TMath::Max( fabs(ynom - yup), fabs(ynom-ydown) );
  else if(direc == up  ) return (float )(yup   - ynom);
  else if(direc == down) return (float )(ydown - ynom);
  else return 0.;
}

float getMCstatError(){
  TH1F* hist;
  TString thefile = path + "Tree_TTbar_Powheg.root";
  TFile* inputfile = TFile::Open(thefile);
  inputfile->GetObject("H_Lep0Pt_" + Chan + "_" + Level,hist);
  return 1/TMath::Sqrt(hist->GetEntries()); 
}

float getXsec(TString tag, TString direc){
	float N = yield("data", Chan, Level); 
	float VV = yield("VV", Chan, Level); float tW = yield("tW", Chan, Level); float DY = yield("DY", Chan, Level); float NonW  = yield("NonW", Chan, Level); float ttbar = yield("ttbar", Chan, Level);
	if(doDD){
		DY *= GetDYDD();
		NonW = GetNonWDD();
	}
	float Nbkg = NonW + tW + VV + DY;
  float xsec = (N-Nbkg)/ttbar*xsec_th;
	if     (tag == "tW" || tag == "VV" || tag == "DY" || tag == "NonW/Z" || tag == "fake" || tag == "NonW" || tag == "NonWZ") (direc == down)? Nbkg -= getSyst(tag) : Nbkg += getSyst(tag);
  else if(tag.Contains("Fidu") || tag.Contains("fidu")){ float ttbarFidu = getYield("TTbar_PowhegFidu"+Chan, Chan, Level)/getFiducialNorm("TTbar_PowhegFidu"+Chan); return GetFiducialXSec()*fabs(getSyst(tag, direc) - ttbarFidu)/ttbarFidu; }
	else if(tag.Contains("trig") || tag.Contains("Trig") || tag.Contains("Elec") || tag.Contains("elec") || tag.Contains("Muon") || tag.Contains("muon") || tag == "JES" || tag == "jes" || tag == "JER" || tag == "jer" || tag == "PDF" || tag == "pdf" || tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales" || tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator" || tag == "had" || tag == "Had" || tag.Contains("ME") || tag == "MCstat") return xsec*fabs(getSyst(tag, direc)-ttbar)/ttbar; 
  else if(tag == "stat") (direc == down)? N -= TMath::Sqrt(N) : N += TMath::Sqrt(N);
  else if(tag == "lumi" || tag == "Lumi") return xsec*(1+lumiunc);
	return (N-Nbkg)/ttbar*xsec_th;
}

float getAcc(bool docout){
	float y = yield("ttbar", Chan, Level);
	float xsec = getXsec();
	float acc_err = TMath::Sqrt( fabs(xsec-getXsec("MCstat"))*fabs(xsec-getXsec("MCstat")) + fabs(xsec-getXsec("scale"))*fabs(xsec-getXsec("scale")) + fabs(xsec-getXsec("ME"))*fabs(xsec-getXsec("ME")) + fabs(xsec-getXsec("pdf"))*fabs(xsec-getXsec("pdf")) + fabs(xsec-getXsec("had"))*fabs(xsec-getXsec("had")) )/xsec;
	float eff_err = TMath::Sqrt( fabs(xsec-getXsec("Trig"))*fabs(xsec-getXsec("Trig")) + fabs(xsec-getXsec("Elec"))*fabs(xsec-getXsec("Elec")) + fabs(xsec-getXsec("Muon"))*fabs(xsec-getXsec("Muon")) + fabs(xsec-getXsec("JES"))*fabs(xsec-getXsec("JES")) + fabs(xsec-getXsec("JER"))*fabs(xsec-getXsec("JER")) )/xsec;

	int ngenEvt = 498000;

	TFile *ffidu = TFile::Open(path + "Tree_" + "TTbar_PowhegFidu" + Chan + ".root");
	TH1F* hfidu; ffidu->GetObject("YieldFidu", hfidu);
	float nFidu = hfidu->GetBinContent(1);

	float nSelect = y/(xsec_th*Lumi)*ngenEvt;

	double BR = 0.03263;
	double acc = nFidu/ngenEvt/BR;
	double eff = nSelect/nFidu;
	double TotalAcc = y/(xsec_th*Lumi);
  double TotalAcc_err = TMath::Sqrt(acc_err*acc_err + eff_err*eff_err);

	if(docout){
    cout << " Acceptance " << endl;
		cout << " ==================================================== " << endl;
		cout << "  Number of generated events: " << ngenEvt << endl;
		cout << "  Number of fiducial events : " << nFidu << endl;
		cout << "  Number of selected events : " << nSelect << endl;
		cout << endl;
		cout << "  BR         = " << BR << endl;
		cout << "  Acceptance = " << acc << " +/- " << acc_err*acc << " (" << acc_err*100 << " %)" << endl;
		cout << "  Efficiency = " << eff << " +/- " << eff_err*eff << " (" << eff_err*100 << " %)" << endl;
		cout << "  BR·Acc·Eff = " << TotalAcc << " +/- " << TotalAcc_err*TotalAcc << " (" << TotalAcc_err*100 << " %)" << endl;
		cout << " ==================================================== " << endl;
	}

	return BR*acc*eff;
}

float getScalePSSyst(bool isFidu, TString direc){
  //return 1;
	float scaleerror = 0;
	float nom; float eup; float edown;

	nom   = getYield("TTbar_PowhegPart"+Chan, Chan, Level);
	eup   = getYield("TTbar_Powheg_ScaleUpPart"+Chan, Chan, Level);
	edown = getYield("TTbar_Powheg_ScaleDownPart"+Chan, Chan, Level);

  float  fidnorm_n = getFiducialNorm("TTbar_PowhegPart" + Chan);
  float  fidnorm_u = getFiducialNorm("TTbar_Powheg_ScaleUpPart" + Chan);
  float  fidnorm_d = getFiducialNorm("TTbar_Powheg_ScaleDownPart" + Chan);

  nom/=fidnorm_n; eup/=fidnorm_u; edown/=fidnorm_d;

  if      (direc == up)   scaleerror = fabs(eup-nom)/nom;
	else if (direc == down) scaleerror = fabs(edown-nom)/nom;
	else                    scaleerror = (float ) TMath::Max(fabs(nom-eup)/nom, fabs(nom-edown)/nom);
  
  float y = yield("ttbar", Chan, Level);

	return y*scaleerror;
}

float getHadSyst(bool isFidu){ 
	float haderror = 0.;
	float nom; float v;
	if(isFidu){ 
		nom = getYield("TTbar_PowhegFidu" + Chan, Chan, Level)/getFiducialNorm("TTbar_PowhegFidu"+Chan);
		v   = getYield("TTbar_Powheg_HerwigFidu" + Chan, Chan, Level)/getFiducialNorm("TTbar_Powheg_HerwigFidu"+Chan);
	}
	else{
		nom = yield("ttbar", Chan, Level);
		v   = getYield("TTbar_Powheg_Herwig", Chan, Level);
	}
	haderror = fabs(nom - v);
  //cout << "haderror = " << haderror/nom*100 << endl;
  //haderror = nom*1.2/100;
	return haderror;
}

void GetHad(){
	float nom1; float v1; float e1; float w1; float p1; float n1;
	float nom2; float v2; float e2; float w2; float p2; float n2;
	float nom3; float v3; float e3; float w3; float p3; float n3;
	float nom4; float v4; float e4; float w4; float p4; float n4;
    
		nom1 = yield("ttbar", Chan, "dilepton");
		v1   = getYield("TTbar_Powheg_Herwig", Chan, "dilepton");
    e1   = fabs(nom1-v1)/nom1*100;
		nom2 = yield("ttbar", Chan, "ZVeto");
		v2   = getYield("TTbar_Powheg_Herwig", Chan, "ZVeto");
    e2   = fabs(nom2-v2)/nom2*100;
		nom3 = yield("ttbar", Chan, "MET");
		v3   = getYield("TTbar_Powheg_Herwig", Chan, "MET");
    e3   = fabs(nom3-v3)/nom3*100;
		nom4 = yield("ttbar", Chan, "2jets");
		v4   = getYield("TTbar_Powheg_Herwig", Chan, "2jets");
    e4   = fabs(nom4-v4)/nom4*100;
    if(Chan == "ElMu"){
      cout << "Channel: ElMu | Hadronization (%) " << endl; //| Parton Shower (%) " << endl;
      cout << "======================================================" << endl;
      cout << " Dilepton:        " << e1 << endl; //"       |      " << p1 << endl;
      cout << " > 2 jets:        " << e4 << endl; //"       |      " << p4 << endl;
      cout << "======================================================" << endl;
    }
    else if (Chan == "Muon"){
      cout << "Channel: MuMu   | Hadronization (%) " << endl; // | Parton Shower  " << endl;
      cout << "=====================================================" << endl;
      cout << " Dilepton       :    " << e1 << endl; //"       |      " << p1 << endl;
      cout << " Z pm 15 Mll cut:    " << e2 << endl; //"       |      " << p2 << endl;
      cout << " MET cut        :    " << e3 << endl; //"       |      " << p3 << endl;
      cout << " > 2 jets       :    " << e4 << endl; //"       |      " << p4 << endl;
      cout << "=====================================================" << endl;
    } 
  
}


void GetPS(){
	float nom1; float v1; float e1; float w1; float p1; float n1;
	float nom2; float v2; float e2; float w2; float p2; float n2;
	float nom3; float v3; float e3; float w3; float p3; float n3;
	float nom4; float v4; float e4; float w4; float p4; float n4;
    
		nom1 = getYield("TTbar_PowhegPartMuon",           "Muon", "dilepton");
		v1   = getYield("TTbar_Powheg_ScaleUpPartMuon",   "Muon", "dilepton");
		w1   = getYield("TTbar_Powheg_ScaleDownPartMuon", "Muon", "dilepton");

		nom2 = getYield("TTbar_PowhegPartMuon",           "Muon", "ZVeto");
		v2   = getYield("TTbar_Powheg_ScaleUpPartMuon",   "Muon", "ZVeto");
		w2   = getYield("TTbar_Powheg_ScaleDownPartMuon", "Muon", "ZVeto");

		nom3 = getYield("TTbar_PowhegPartMuon",           "Muon", "MET");
		v3   = getYield("TTbar_Powheg_ScaleUpPartMuon",   "Muon", "MET");
		w3   = getYield("TTbar_Powheg_ScaleDownPartMuon", "Muon", "MET");

		nom4 = getYield("TTbar_PowhegPartMuon",           "Muon", "2jets");
		v4   = getYield("TTbar_Powheg_ScaleUpPartMuon",   "Muon", "2jets");
		w4   = getYield("TTbar_Powheg_ScaleDownPartMuon", "Muon", "2jets");

    float fidnorm_n = getFiducialNorm("TTbar_PowhegPartMuon");
    float fidnorm_v = getFiducialNorm("TTbar_Powheg_ScaleUpPartMuon");
    float fidnorm_w = getFiducialNorm("TTbar_Powheg_ScaleDownPartMuon");
/*
    cout << "nom  = " << nom4/Lumi << endl;
    cout << "up   = " << v4  /Lumi << endl;
    cout << "down = " << w4  /Lumi << endl;
    cout << "fid n = " << fidnorm_n << endl;
    cout << "fid u = " << fidnorm_v << endl;
    cout << "fid d = " << fidnorm_w << endl;
  */  
    nom1/=fidnorm_n; nom2/=fidnorm_n; nom3/=fidnorm_n; nom4/=fidnorm_n;
    v1/=fidnorm_v; v2/=fidnorm_v; v3/=fidnorm_v; v4/=fidnorm_v;
    w1/=fidnorm_w; w2/=fidnorm_w; w3/=fidnorm_w; w4/=fidnorm_w;


    e1   = TMath::Max(fabs(nom1-v1)/nom1*100, fabs(nom1-w1)/nom1*100);
    e2   = TMath::Max(fabs(nom2-v2)/nom2*100, fabs(nom2-w2)/nom2*100);
    e3   = TMath::Max(fabs(nom3-v3)/nom3*100, fabs(nom3-w3)/nom3*100);
    e4   = TMath::Max(fabs(nom4-v4)/nom4*100, fabs(nom4-w4)/nom4*100);


		cout << "==================================" << endl;
		cout << "Channel: MuMu   |  PS scale (%) " << endl; // | Parton Shower  " << endl;
		cout << "----------------------------------" << endl;
		cout << " Dilepton       :    " << e1 << endl; //"       |      " << p1 << endl;
		cout << " Z pm 15 Mll cut:    " << e2 << endl; //"       |      " << p2 << endl;
		cout << " MET cut        :    " << e3 << endl; //"       |      " << p3 << endl;
		cout << " > 2 jets       :    " << e4 << endl; //"       |      " << p4 << endl;
		cout << "===================================" << endl;

}

float getNLOsyst(bool isFidu){
	float nloerror = 0.;
	float nom; float v;
	if(isFidu){
		nom = getYield("TTbar_PowhegFidu"+Chan, Chan, Level)*getFiducialNorm("TTbar_PowhegFidu"+Chan);
		v   = getYield("TTJets_aMCatNLOFidu"+Chan, Chan, Level)*getFiducialNorm("TTbar_aMCatNLOFidu"+Chan);
	}
	else{
		nom = yield("ttbar", Chan, Level);
		v   = getYield("TTJets_aMCatNLO", Chan, Level);
	}
	nloerror = fabs(nom - v);
	return nloerror;
} 

void PrintXsecSyst(){
  float xsec = getXsec();
  float totalSys = TMath::Sqrt( (fabs(xsec-getXsec("MCstat"))*fabs(xsec-getXsec("MCstat")) + fabs(xsec-getXsec("Trig"))*fabs(xsec-getXsec("Trig")) + fabs(xsec-getXsec("Elec"))*fabs(xsec-getXsec("Elec")) + fabs(xsec-getXsec("Muon"))*fabs(xsec-getXsec("Muon")) + fabs(xsec-getXsec("JES"))*fabs(xsec-getXsec("JES")) + fabs(xsec-getXsec("JER"))*fabs(xsec-getXsec("JER")) + fabs(xsec-getXsec("scale"))*fabs(xsec-getXsec("scale")) + fabs(xsec-getXsec("ME"))*fabs(xsec-getXsec("ME")) + fabs(xsec-getXsec("pdf"))*fabs(xsec-getXsec("pdf")) + fabs(xsec-getXsec("had"))*fabs(xsec-getXsec("had")) + fabs(xsec-getXsec("tW"))*fabs(xsec-getXsec("tW")) + fabs(xsec-getXsec("VV"))*fabs(xsec-getXsec("VV")) + fabs(xsec-getXsec("DY"))*fabs(xsec-getXsec("DY")) + fabs(xsec-getXsec("NonW"))*fabs(xsec-getXsec("NonW")) )/(xsec*xsec));
  float stat = fabs(xsec-getXsec("stat"))/xsec;
  float lumin = fabs(xsec-getXsec("lumi"))/xsec;
  float total = TMath::Sqrt(totalSys*totalSys + lumin*lumin + stat*stat);
  cout << "==================================================================================" << endl;
  cout << " Source                        #Delta #sigma_{t#bar{t}}(pb)       (   %)          " << endl;
  cout << "----------------------------------------------------------------------------------" << endl;
  cout << Form(" Trigger efficiencies                   %2.3f                      %1.2f          ",fabs(xsec-getXsec("Trig")) , fabs(xsec-getXsec("Trig"))/xsec*100 ) << endl;
  cout << Form(" Electron efficiencies                  %2.3f                      %1.2f          ",fabs(xsec-getXsec("Elec" )), fabs(xsec-getXsec("Elec"))/xsec*100 ) << endl;
  cout << Form(" Muon efficiencies                      %2.3f                      %1.2f          ",fabs(xsec-getXsec("Muon")) , fabs(xsec-getXsec("Muon"))/xsec*100 ) << endl;
  cout << Form(" Jet Energy Scale                       %2.3f                      %1.2f          ",fabs(xsec-getXsec("JES")) , fabs(xsec-getXsec("JES"))/xsec*100 ) << endl;
  cout << Form(" Jer Energy Resolution                  %2.3f                      %1.2f          ",fabs(xsec-getXsec("JER")) , fabs(xsec-getXsec("JER"))/xsec*100 ) << endl;
  cout << Form(" QCD scale PS                           %2.3f                      %1.2f          ",fabs(xsec-getXsec("scale")) , fabs(xsec-getXsec("scale"))/xsec*100 ) << endl;
  cout << Form(" QCD scale ME                           %2.3f                      %1.2f          ",fabs(xsec-getXsec("ME")) , fabs(xsec-getXsec("ME"))/xsec*100 ) << endl;
  cout << Form(" PDF                                    %2.3f                      %1.2f          ",fabs(xsec-getXsec("pdf")) , fabs(xsec-getXsec("pdf"))/xsec*100 ) << endl;
  cout << Form(" Hadronization                          %2.3f                      %1.2f          ",fabs(xsec-getXsec("had")) , fabs(xsec-getXsec("had"))/xsec*100 ) << endl;
  cout << Form(" MC statistics                          %2.3f                      %1.2f          ",fabs(xsec-getXsec("MCstat")) , fabs(xsec-getXsec("MCstat"))/xsec*100 ) << endl;
  cout << "----------------------------------------------------------------------------------" << endl;
  cout << Form(" tW background                          %2.3f                      %1.2f          ",fabs(xsec-getXsec("tW")) , fabs(xsec-getXsec("tW"))/xsec*100 ) << endl;
  cout << Form(" VV background                          %2.3f                      %1.2f          ",fabs(xsec-getXsec("VV")) , fabs(xsec-getXsec("VV"))/xsec*100 ) << endl;
  cout << Form(" DY background                          %2.3f                      %1.2f          ",fabs(xsec-getXsec("DY")) , fabs(xsec-getXsec("DY"))/xsec*100 ) << endl;
  cout << Form(" Non prompt backgrounds                 %2.3f                      %1.2f          ",fabs(xsec-getXsec("NonW")) , fabs(xsec-getXsec("NonW"))/xsec*100 ) << endl;
  cout << "----------------------------------------------------------------------------------" << endl;
  cout << Form(" Total systematic                       %2.3f                     %1.2f          ",xsec*totalSys , totalSys*100) << endl;
  cout << Form(" Stat                                  %2.3f                     %1.2f          ",xsec*stat , stat*100) << endl;
  cout << Form(" Lumi                                   %2.3f                     %1.2f          ",xsec*lumin , lumin*100) << endl;
  cout << "----------------------------------------------------------------------------------" << endl;
  cout << Form(" Total                                 %2.3f                     %1.2f          ",xsec*total , total*100) << endl;
  cout << "==================================================================================" << endl;
}

void PrintGoodYields(){
  TString c;
	float totalbkg_dilep = yield("DY", Chan, "dilepton")+yield("fake", Chan, "dilepton")+yield("tW", Chan, "dilepton")+yield("VV", Chan, "dilepton");
	float totalbkg_zveto = yield("DY", Chan, "ZVeto")+yield("fake", Chan, "ZVeto")+yield("tW", Chan, "ZVeto")+yield("VV", Chan, "ZVeto");
  float totalbkg_met   = yield("DY", Chan, "MET")+yield("fake", Chan, "MET")+yield("tW", Chan, "MET")+yield("VV", Chan, "MET");
  float totalbkg_2jets = yield("DY", Chan, "2jets")+yield("fake", Chan, "2jets")+yield("tW", Chan, "2jets")+yield("VV", Chan, "2jets");

	float totalbkg_dilep_err = getStatError("DY", Chan, "dilepton")+getStatError("fake", Chan, "dilepton")+getStatError("tW", Chan, "dilepton")+getStatError("VV", Chan, "dilepton");
	float totalbkg_zveto_err = getStatError("DY", Chan, "ZVeto")+getStatError("fake", Chan, "ZVeto")+getStatError("tW", Chan, "ZVeto")+getStatError("VV", Chan, "ZVeto");
  float totalbkg_met_err   = getStatError("DY", Chan, "MET")+getStatError("fake", Chan, "MET")+getStatError("tW", Chan, "MET")+getStatError("VV", Chan, "MET");
  float totalbkg_2jets_err = getStatError("DY", Chan, "2jets")+getStatError("fake", Chan, "2jets")+getStatError("tW", Chan, "2jets")+getStatError("VV", Chan, "2jets");
	if(doDD){
		totalbkg_2jets= yield("DY", Chan, "2jets")*GetDYDD()+GetNonWDD()+yield("tW", Chan, "2jets")+yield("VV", Chan, "2jets");
		totalbkg_dilep = yield("DY", Chan, "dilepton")*GetDYDD(Chan, "dilepton")+GetNonWDD(Chan, "dilepton")+yield("tW", Chan, "dilepton")+yield("VV", Chan, "dilepton");
	}
  float yDY_dilep     = yield("DY", Chan, "dilepton");
  float yDY_dilep_err = getStatError("DY", Chan, "dilepton");
  float yDY_2jets     = yield("DY", Chan, "2jets");
  float yDY_2jets_err = getStatError("DY", Chan, "2jets");
	float yNonW_dilep     = yield("fake", Chan, "dilepton");
	float yNonW_dilep_err = getStatError("fake", Chan, "dilepton");
	float yNonW_2jets     = yield("fake", Chan, "2jets");
	float yNonW_2jets_err = getStatError("fake", Chan, "2jets");
  TString DYflag   = " Drell-Yan (MC)    ";
  TString NonWflag = " NonW leptons (MC) ";
	if(doDD){
		DYflag   = " Drell-Yan (DD)    ";
		NonWflag = " NonW leptons (DD) ";
		yDY_dilep*=GetDYDD(Chan, "dilepton");
		yDY_dilep_err = yDY_dilep*GetDYDD(Chan, "dilepton", 1);
		yDY_2jets*=GetDYDD(Chan, "2jets");
		yDY_2jets_err = yDY_dilep*GetDYDD(Chan, "2jets", 1);
		yNonW_dilep     = GetNonWDD(Chan, "dilepton", 0);
		yNonW_dilep_err = GetNonWDD(Chan, "dilepton", 1);
		yNonW_2jets     = GetNonWDD(Chan, "2jets", 0);
		yNonW_2jets_err = GetNonWDD(Chan, "2jets", 1);
	}
	if(Chan == "ElMu"){
		cout <<      "============================================================" << endl;
		cout <<      " Source            |     emu OS pair     |   emu + >= 2jets  " << endl;
		cout <<      "------------------------------------------------------------" << endl;
		cout << DYflag          +  Form("|  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yDY_dilep, yDY_dilep_err, yDY_2jets, yDY_2jets_err)  << endl;
		cout << NonWflag        +  Form("|  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yNonW_dilep, yNonW_dilep_err, yNonW_2jets, yNonW_2jets_err) << endl;
		cout << Form(" tW + tbarW        |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yield("tW", Chan, "dilepton"), getStatError("tW", Chan, "dilepton"), yield("tW", Chan, "2jets"), getStatError("tW", Chan, "2jets"))     << endl;
		cout << Form(" WW + WZ           |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yield("VV", Chan, "dilepton"), getStatError("VV", Chan, "dilepton"), yield("VV", Chan, "2jets"), getStatError("VV", Chan, "2jets"))     << endl;
		cout <<      "------------------------------------------------------------" << endl;
		cout << Form(" Total background  |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , totalbkg_dilep, totalbkg_dilep_err, totalbkg_2jets, totalbkg_2jets_err)   << endl;
		cout <<      "------------------------------------------------------------" << endl;
		cout << Form(" Signal (tt->emu)  |  %5.2f +/- %2.2f     |  %5.2f  +/- %2.2f       " , yield("ttbar", Chan, "dilepton"), getStatError("ttbar", Chan, "dilepton"), yield("ttbar", Chan, "2jets"), getStatError("ttbar", Chan, "2jets")) << endl;
		cout <<      "------------------------------------------------------------" << endl;
		cout << Form(" Data                 |     %5.0f       |   %5.0f         " , yield("data", Chan, "dilepton"), yield("data", Chan, "2jets")) << endl;
		cout <<      "============================================================" << endl;
	}
	else if(Chan == "Muon"){
  float yDY_zveto     = yield("DY", Chan, "ZVeto");
  float yDY_zveto_err = getStatError("DY", Chan, "ZVeto");
  float yDY_met       = yield("DY", Chan, "MET");
  float yDY_met_err = getStatError("DY", Chan, "MET");
	float yNonW_zveto   = yield("fake", Chan, "ZVeto");
	float yNonW_zveto_err = getStatError("fake", Chan, "ZVeto");
	float yNonW_met     = yield("fake", Chan, "MET");
	float yNonW_met_err = getStatError("fake", Chan, "MET");
		cout <<      "========================================================================================================" << endl;
		cout <<      " Source            |     µµ OS pair      |     µµ + ZVeto      |      MET cut      |     µµ, >= 2jets      " << endl;
		cout <<      "--------------------------------------------------------------------------------------------------------" << endl;
		cout << DYflag          +  Form("|  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f       " , yDY_dilep, yDY_dilep_err, yDY_zveto, yDY_zveto_err, yDY_met, yDY_met_err, yDY_2jets, yDY_2jets_err)  << endl;
		cout << NonWflag        +  Form("|  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f       " , yNonW_dilep, yNonW_dilep_err, yNonW_zveto, yNonW_zveto_err, yNonW_met, yNonW_met_err, yNonW_2jets, yNonW_2jets_err) << endl;
		cout << Form(" tW + tbarW        |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f       " , yield("tW", Chan, "dilepton"), getStatError("tW", Chan, "dilepton"),  yield("tW", Chan, "ZVeto"), getStatError("tW", Chan, "ZVeto"), yield("tW", Chan, "MET"), getStatError("tW", Chan, "MET"), yield("tW", Chan, "2jets"), getStatError("tW", Chan, "2jets"))     << endl;
		cout << Form(" WW + WZ           |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f       " , yield("VV", Chan, "dilepton"), getStatError("VV", Chan, "dilepton"),  yield("VV", Chan, "ZVeto"), getStatError("VV", Chan, "ZVeto"), yield("VV", Chan, "MET"), getStatError("VV", Chan, "MET"), yield("VV", Chan, "2jets"), getStatError("VV", Chan, "2jets"))     << endl;
		cout <<      "--------------------------------------------------------------------------------------------------------" << endl;
		cout << Form(" Total background  |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f       " ,totalbkg_dilep, totalbkg_dilep_err, totalbkg_zveto, totalbkg_zveto_err, totalbkg_met, totalbkg_met_err, totalbkg_2jets, totalbkg_2jets_err) << endl; 
		cout <<      "--------------------------------------------------------------------------------------------------------" << endl;
		cout << Form(" Singal (tt->µµ)   |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f    |  %5.2f +/- %2.2f       " , yield("ttbar", Chan, "dilepton"), getStatError("ttbar", Chan, "dilepton"),  yield("ttbar", Chan, "ZVeto"), getStatError("ttbar", Chan, "ZVeto"), yield("ttbar", Chan, "MET"), getStatError("ttbar", Chan, "MET"), yield("ttbar", Chan, "2jets"), getStatError("ttbar", Chan, "2jets"))     << endl;
		cout <<      "--------------------------------------------------------------------------------------------------------" << endl;
		cout << Form(" Data              |      %5.0f          |      %5.0f          |      %5.0f          |      %5.0f              " , yield("data", Chan, "dilepton"), yield("data", Chan, "ZVeto"), yield("data", Chan, "MET"), yield("data", Chan, "2jets")) << endl;
		cout <<      "========================================================================================================" << endl;
	}
}

void PrintXSec(){
	float VV = yield("VV", Chan, Level); float tW = yield("tW", Chan, Level); float DY = yield("DY", Chan, Level); float NonW  = yield("NonW", Chan, Level); float ttbar = yield("ttbar", Chan, Level);
	if(doDD){
		DY *= GetDYDD();
		NonW = GetNonWDD();
	}
	float N_bkg = NonW + tW + VV + DY;
  float xsec = getXsec();
  float totalSys = TMath::Sqrt( fabs(xsec-getXsec("MCstat"))*fabs(xsec-getXsec("MCstat")) + (fabs(xsec-getXsec("Trig"))*fabs(xsec-getXsec("Trig")) + fabs(xsec-getXsec("Elec"))*fabs(xsec-getXsec("Elec")) + fabs(xsec-getXsec("Muon"))*fabs(xsec-getXsec("Muon")) + fabs(xsec-getXsec("JES"))*fabs(xsec-getXsec("JES")) + fabs(xsec-getXsec("JER"))*fabs(xsec-getXsec("JER")) + fabs(xsec-getXsec("scale"))*fabs(xsec-getXsec("scale")) + fabs(xsec-getXsec("ME"))*fabs(xsec-getXsec("ME")) + fabs(xsec-getXsec("pdf"))*fabs(xsec-getXsec("pdf")) + fabs(xsec-getXsec("had"))*fabs(xsec-getXsec("had")) + fabs(xsec-getXsec("tW"))*fabs(xsec-getXsec("tW")) + fabs(xsec-getXsec("VV"))*fabs(xsec-getXsec("VV")) + fabs(xsec-getXsec("DY"))*fabs(xsec-getXsec("DY")) + fabs(xsec-getXsec("NonW"))*fabs(xsec-getXsec("NonW")) ));
  float stat = fabs(xsec-getXsec("stat"));
  float lumin = fabs(xsec-getXsec("lumi"));
  cout <<      "===============================================================" << endl;
  cout << Form("                  N_data = %2.0f", yield("data", Chan, "2jets")) << endl; 
  cout << Form("                  N_bkg  = %2.2f", N_bkg)                        << endl; 
  cout << Form("                  Accept = %2.5f", getAcc())                     << endl;
  cout <<      "---------------------------------------------------------------" << endl;
	cout << Form("  sigma_tt = %2.1f +/- %2.1f (stat) +/- %2.1f (sys) +/- %2.1f (lumi) ", xsec, stat, totalSys, lumin) << endl; 
  cout <<      "===============================================================" << endl;
}

void GetXSec(){
	cout << Form(" ttbar cross section at sqrt(s) = 5 TeV, integrated luminosity = %2.1f pb-1", Lumi) << endl;
  doDD = 0;
  cout << "\n##########################################################" << endl;
  cout << "####### Yields and cross section using MC estimates ######" << endl;
  cout << "##########################################################\n" << endl; 
	PrintGoodYields();
	cout << endl;
	PrintXsecSyst();
	cout << endl;
	PrintXSec();
  cout << "\n\n" << endl;
/*  doDD = 1;
  cout << "##########################################################" << endl;
  cout << "####### Yields and cross section using DD estimates ######" << endl;
  cout << "##########################################################\n" << endl; 
	PrintGoodYields();
	cout << endl;
	PrintXsecSyst();
	cout << endl;
	PrintXSec();
  doDD = 0;
  cout << "\n\n" << endl;
  cout << "#################################" << endl;
  cout << "###### Data Driven methods ######" << endl;
  cout << "#################################\n" << endl; 
  GetDataDriven();
  cout << "\n\n" << endl;
  cout << "################################" << endl;
  cout << "########## Acceptance ##########" << endl;
  cout << "################################\n" << endl;
	getAcc(1);
*/ //  GetFiducialXSec(1);
}

float getFiducialNorm(TString SampleName){
  float TotalEvents = 0.;
  if     (SampleName == "TTbar_PowhegFidu"          +Chan || SampleName == "TTbar_PowhegPart"          + Chan)  TotalEvents = 498000;
  else if(SampleName == "TTbar_Powheg_ScaleUpFidu"  +Chan || SampleName == "TTbar_Powheg_ScaleUpPart"  + Chan)  TotalEvents = 462000;
  else if(SampleName == "TTbar_Powheg_ScaleDownFidu"+Chan || SampleName == "TTbar_Powheg_ScaleDownPart"+ Chan)  TotalEvents = 378500;
  else if(SampleName == "TTbar_Powheg_HerwigFidu"   +Chan)    TotalEvents = 494088;
  else return 1.;
  TFile *fidufile = TFile::Open(path + "Tree_" + SampleName + ".root");
  TH1F* h; fidufile->GetObject("YieldFidu", h);
  float FiduEvents = h->GetBinContent(1);
  float fidunorm = FiduEvents/TotalEvents;
  return fidunorm;
}


float GetFiducialXSec(bool docout){
  IsFidu = 1;
  float TotalEvents = 498000;
  TFile *fidufile = TFile::Open(path + "Tree_" + "TTbar_PowhegFidu"+ Chan + ".root");
  TH1F* h; fidufile->GetObject("YieldFidu", h);
  float FiduEvents = h->GetBinContent(1);
  float fidunorm = FiduEvents/TotalEvents; 

  float ttbar = getYield("TTbar_PowhegFidu"+Chan, Chan, "2jets");
  float N     = yield("data", Chan, Level);
  float Nbkg  = yield("bkg", Chan, Level);
  float xsec  = (N-Nbkg)/(ttbar/fidunorm)*xsec_th;
  if(docout){
    cout << " ##################################" << endl;
    cout << " ##### FIDUCIAL CROSS SECTION #####" << endl;
    cout << " ##################################" << endl;
    cout << endl;
    cout << Form(" # Total events in the sample        = %8.0f", TotalEvents) << endl;
    cout << Form(" # Events in the fiducial region     = %8.0f",  FiduEvents) << endl;
    cout << Form(" # Signal events after the selection = %8.2f",       ttbar) << endl;
    cout << Form(" # Fiducial cross section (pb)       = %8.4f",        xsec) << endl;

		cout << endl;
		cout <<      "================================================================================" << endl;
    cout <<      " Modiling uncertainties     #Delta #sigma_{t#bar{t}}(pb)       (   %)          " << endl;
		cout <<      "-------------------------------------------------------------------------------" << endl;
		cout << Form(" QCD scale PS                           %2.4f                      %1.2f       ",fabs(xsec-getXsec("scaleFidu")), fabs(xsec-getXsec("scaleFidu"))/xsec*100) << endl;
		cout << Form(" QCD scale ME                           %2.4f                      %1.2f       ",fabs(xsec-getXsec("MEFidu")) , fabs(xsec-getXsec("MEFidu"))/xsec*100 ) << endl;
		cout << Form(" PDF                                    %2.4f                      %1.2f       ",fabs(xsec-getXsec("pdfFidu")) , fabs(xsec-getXsec("pdfFidu"))/xsec*100 ) << endl;
		cout << Form(" Hadronization                          %2.4f                      %1.2f       ",fabs(xsec-getXsec("hadFidu")) , fabs(xsec-getXsec("hadFidu"))/xsec*100 ) << endl;
		cout <<      "===============================================================================" << endl;
		cout << endl;
  }
  return xsec;
}
