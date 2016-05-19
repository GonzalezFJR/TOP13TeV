#include "DrawPlots.C"
#include "GetPDFweights.C"
#include "GetDataDriven.C"

float getSyst(TString tag, TString direc = "Max"); // Returns the variation in the yield due to a systematic uncertainty
float getExpError(TString Process = "ttbar", TString tag = "0", TString direc = "Max"); // Returns the variation in the ttbat yield due to a experimental uncertainty
float getAcc(TString tag = "0", TString direc = "Max");
float getXsec(TString tag = "0", TString direc = "Max");

float getNLOsyst(bool isFidu = false);
float getHadSyst(bool isFidu = false); 
float getScalePSSyst(bool isFidu = false, TString direc = "Max"); 
void PrintXsecSyst();
void PrintGoodYields();
void GetXSec();
float getFiducialNorm(TString SampleName = "TTbar_PowhegFidu");
float GetFiducialXSec(bool docout = false);

///////// float DYDD(bool docout = false);
///////// float NonWDD(bool docout = false)

const float lumiunc = 0.12;
const float xsec_th = 67.8;

const float tWerr = 0.3; const float VVerr = 0.3; const float DYerr = 0.15; const float NonWerr = 0.3;
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
  // Backgrounds
  else if(tag == "tW")                                                          return tWerr*yield("tW", Chan, Level); // tW
  else if(tag == "VV")                                                          return VVerr*yield("VV", Chan, Level); // VV
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

float getXsec(TString tag, TString direc){
	float N = yield("data", Chan, Level); 
	float VV = yield("VV", Chan, Level); float tW = yield("tW", Chan, Level); float DY = yield("DY", Chan, Level); float NonW  = yield("NonW", Chan, Level); float ttbar = yield("ttbar", Chan, Level);
	float Nbkg = NonW + tW + VV + DY;
  float xsec = (N-Nbkg)/ttbar*xsec_th;
	if     (tag == "tW" || tag == "VV" || tag == "DY" || tag == "NonW/Z" || tag == "fake" || tag == "NonW" || tag == "NonWZ") (direc == down)? Nbkg -= getSyst(tag) : Nbkg += getSyst(tag);
  else if(tag.Contains("Fidu") || tag.Contains("fidu")){ float ttbarFidu = getYield("TTbar_PowhegFidu", Chan, Level)/getFiducialNorm("TTbar_PowhegFidu"); return GetFiducialXSec()*fabs(getSyst(tag, direc) - ttbarFidu)/ttbarFidu; }
	else if(tag.Contains("trig") || tag.Contains("Trig") || tag.Contains("Elec") || tag.Contains("elec") || tag.Contains("Muon") || tag.Contains("muon") || tag == "JES" || tag == "jes" || tag == "JER" || tag == "jer" || tag == "PDF" || tag == "pdf" || tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales" || tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator" || tag == "had" || tag == "Had" || tag.Contains("ME")) return xsec*fabs(getSyst(tag, direc)-ttbar)/ttbar; 
  else if(tag == "stat") (direc == down)? N -= TMath::Sqrt(N) : N += TMath::Sqrt(N);
  else if(tag == "lumi" || tag == "Lumi") return xsec*(1+lumiunc);
	return (N-Nbkg)/ttbar*xsec_th;
}

float getAcc(TString tag, TString direc){
	float y = yield("ttbar", Chan, Level);
	if(tag.Contains("trig") || tag.Contains("Trig") || tag.Contains("Elec") || tag.Contains("elec") || tag.Contains("Muon") || tag.Contains("muon") || tag == "JES" || tag == "jes" || tag == "JER" || tag == "jer" || tag == "PDF" || tag == "pdf" || tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales" || tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator" || tag == "had" || tag == "Had") y += getSyst(tag, direc);
	return y/(xsec_th*Lumi);
}

float getScalePSSyst(bool isFidu, TString direc){
	float scaleerror = 0;
	float nom; float eup; float edown;
/*	if(isFidu){
		nom = getYield("TTbar_PowhegFidu", Chan, Level)/getFiducialNorm("TTbar_PowhegFidu");
		eup   = getYield("TTbar_Powheg_ScaleUpFidu", Chan, Level)/getFiducialNorm("TTbar_Powheg_ScaleUpFidu");
		edown = getYield("TTbar_Powheg_ScaleDownFidu", Chan, Level)/getFiducialNorm("TTbar_Powheg_ScaleDownFidu");
	}
	else{
		nom = getYield("TTbar_Powheg", Chan, Level);
		eup   = getYield("TTbar_Powheg_ScaleUp", Chan, Level);
		edown = getYield("TTbar_Powheg_ScaleDown", Chan, Level);
	}*/
		nom = getYield("TTbar_PowhegPart", Chan, Level);
		eup   = getYield("TTbar_Powheg_ScaleUpPart", Chan, Level);
		edown = getYield("TTbar_Powheg_ScaleDownPart", Chan, Level);
	if     (direc == maxi) scaleerror = (float ) TMath::Max(fabs(nom-eup), fabs(nom-edown));
	else if(direc == up)   scaleerror = (eup-nom);
	else if(direc == down) scaleerror = (edown-nom);
	return scaleerror;
}

float getHadSyst(bool isFidu){ 
	float haderror = 0.;
	float nom; float v;
	if(isFidu){ 
		nom = getYield("TTbar_PowhegFidu", Chan, Level)/getFiducialNorm("TTbar_PowhegFidu");
		v   = getYield("TTbar_Powheg_HerwigFidu", Chan, Level)/getFiducialNorm("TTbar_Powheg_HerwigFidu");
	}
	else{
		nom = yield("ttbar", Chan, Level);
		v   = getYield("TTbar_Powheg_Herwig", Chan, Level);
	}
	haderror = fabs(nom - v);
	return haderror;
}

float getNLOsyst(bool isFidu){
	float nloerror = 0.;
	float nom; float v;
	if(isFidu){
		nom = getYield("TTbar_PowhegFidu", Chan, Level)*getFiducialNorm("TTbar_PowhegFidu");
		v   = getYield("TTJets_aMCatNLOFidu", Chan, Level)*getFiducialNorm("TTbar_aMCatNLOFidu");
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
  float totalSys = TMath::Sqrt( (fabs(xsec-getXsec("Trig"))*fabs(xsec-getXsec("Trig")) + fabs(xsec-getXsec("Elec"))*fabs(xsec-getXsec("Elec")) + fabs(xsec-getXsec("Muon"))*fabs(xsec-getXsec("Muon")) + fabs(xsec-getXsec("JES"))*fabs(xsec-getXsec("JES")) + fabs(xsec-getXsec("JER"))*fabs(xsec-getXsec("JER")) + fabs(xsec-getXsec("scale"))*fabs(xsec-getXsec("scale")) + fabs(xsec-getXsec("ME"))*fabs(xsec-getXsec("ME")) + fabs(xsec-getXsec("pdf"))*fabs(xsec-getXsec("pdf")) + fabs(xsec-getXsec("had"))*fabs(xsec-getXsec("had")) + fabs(xsec-getXsec("tW"))*fabs(xsec-getXsec("tW")) + fabs(xsec-getXsec("VV"))*fabs(xsec-getXsec("VV")) + fabs(xsec-getXsec("DY"))*fabs(xsec-getXsec("DY")) + fabs(xsec-getXsec("NonW"))*fabs(xsec-getXsec("NonW")) )/(xsec*xsec));
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
  TString Chan = "ElMu";
  TString c;
  float totalbkg_2jets= yield("DY", Chan, "2jets")+yield("fake", Chan, "2jets")+yield("tW", Chan, "2jets")+yield("VV", Chan, "2jets");
  float totalbkg_dilep = yield("DY", Chan, "dilepton")+yield("fake", Chan, "dilepton")+yield("tW", Chan, "dilepton")+yield("VV", Chan, "dilepton");
	cout <<      "============================================================" << endl;
  cout <<      " Source            |     emu OS pair     |   emu + >= 2jets  " << endl;
  cout <<      "------------------------------------------------------------" << endl;
  cout << Form(" Drell-Yan (MC)    |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yield("DY", Chan, "dilepton"), getStatError("DY", Chan, "dilepton"), yield("DY", Chan, "2jets"), getStatError("DY", Chan, "2jets"))     << endl;
  cout << Form(" NonW leptons (MC) |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yield("fake", Chan, "dilepton"), getStatError("fake", Chan, "dilepton"), yield("fake", Chan, "2jets"), getStatError("fake", Chan, "2jets")) << endl;
  cout << Form(" tW + tbarW        |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yield("tW", Chan, "dilepton"), getStatError("tW", Chan, "dilepton"), yield("tW", Chan, "2jets"), getStatError("tW", Chan, "2jets"))     << endl;
  cout << Form(" WW + WZ           |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , yield("VV", Chan, "dilepton"), getStatError("VV", Chan, "dilepton"), yield("VV", Chan, "2jets"), getStatError("VV", Chan, "2jets"))     << endl;
  cout <<      "------------------------------------------------------------" << endl;
  cout << Form(" Total background  |  %5.2f +/- %2.2f     |  %5.2f +/- %2.2f        " , totalbkg_dilep, getStatError("bkg", Chan, "dilepton"), totalbkg_2jets, getStatError("bkg", Chan, "2jets"))   << endl;
  cout <<      "------------------------------------------------------------" << endl;
  cout << Form(" Signal (tt->emu)  |  %5.2f +/- %2.2f     |  %5.2f  +/- %2.2f       " , yield("ttbar", Chan, "dilepton"), getStatError("ttbar", Chan, "dilepton"), yield("ttbar", Chan, "2jets"), getStatError("ttbar", Chan, "2jets")) << endl;
  cout <<      "------------------------------------------------------------" << endl;
  cout << Form(" Data                 |     %5.0f       |   %5.0f         " , yield("data", Chan, "dilepton"), yield("data", Chan, "2jets")) << endl;
	cout <<      "============================================================" << endl;
}

void GetXSec(){
  TString Chan = "ElMu";
  float N_bkg = yield("DY", Chan, "2jets")+yield("fake", Chan, "2jets")+yield("tW", Chan, "2jets")+yield("VV", Chan, "2jets");
  float xsec = getXsec();
  float totalSys = TMath::Sqrt( (fabs(xsec-getXsec("Trig"))*fabs(xsec-getXsec("Trig")) + fabs(xsec-getXsec("Elec"))*fabs(xsec-getXsec("Elec")) + fabs(xsec-getXsec("Muon"))*fabs(xsec-getXsec("Muon")) + fabs(xsec-getXsec("JES"))*fabs(xsec-getXsec("JES")) + fabs(xsec-getXsec("JER"))*fabs(xsec-getXsec("JER")) + fabs(xsec-getXsec("scale"))*fabs(xsec-getXsec("scale")) + fabs(xsec-getXsec("ME"))*fabs(xsec-getXsec("ME")) + fabs(xsec-getXsec("pdf"))*fabs(xsec-getXsec("pdf")) + fabs(xsec-getXsec("had"))*fabs(xsec-getXsec("had")) + fabs(xsec-getXsec("tW"))*fabs(xsec-getXsec("tW")) + fabs(xsec-getXsec("VV"))*fabs(xsec-getXsec("VV")) + fabs(xsec-getXsec("DY"))*fabs(xsec-getXsec("DY")) + fabs(xsec-getXsec("NonW"))*fabs(xsec-getXsec("NonW")) ));
  float stat = fabs(xsec-getXsec("stat"));
  float lumin = fabs(xsec-getXsec("lumi"));
	cout << Form(" ttbar cross section at sqrt(s) = 5 TeV, integrated luminosity = %2.1f pb-1", Lumi) << endl;
	cout << endl;
	PrintGoodYields();
	cout << endl;
	PrintXsecSyst();
	cout << endl;
  cout <<      "===============================================================" << endl;
  cout << Form("                  N_data = %2.0f", yield("data", Chan, "2jets")) << endl; 
  cout << Form("                  N_bkg  = %2.2f", N_bkg)                        << endl; 
	cout << Form("                  Acceptance = %2.4f", getAcc())                 << endl;
  cout <<      "---------------------------------------------------------------" << endl;
	cout << Form("  sigma_tt = %2.1f +/- %2.1f (stat) +/- %2.1f (sys) +/- %2.1f (lumi) ", xsec, stat, totalSys, lumin) << endl; 
  cout <<      "===============================================================" << endl;
  cout << endl;
  GetDataDriven();
  cout << endl;
  GetFiducialXSec(1);
}

float getFiducialNorm(TString SampleName){
  float TotalEvents = 0.;
  if     (SampleName == "TTbar_PowhegFidu")           TotalEvents = 498000;
  else if(SampleName == "TTbar_Powheg_ScaleUpFidu")   TotalEvents = 462000;
  else if(SampleName == "TTbar_Powheg_ScaleDownFidu") TotalEvents = 378500;
  else if(SampleName == "TTbar_Powheg_HerwigFidu")    TotalEvents = 494088;
  else return 0.;
  TFile *fidufile = TFile::Open(path + "Tree_" + SampleName + ".root");
  TH1F* h; fidufile->GetObject("YieldFidu", h);
  float FiduEvents = h->GetBinContent(1);
  float fidunorm = FiduEvents/TotalEvents;
  return fidunorm;
}


float GetFiducialXSec(bool docout){
  IsFidu = 1;
  float TotalEvents = 498000;
  TFile *fidufile = TFile::Open(path + "Tree_" + "TTbar_PowhegFidu" + ".root");
  TH1F* h; fidufile->GetObject("YieldFidu", h);
  float FiduEvents = h->GetBinContent(1);
  float fidunorm = FiduEvents/TotalEvents; 

  float ttbar = getYield("TTbar_PowhegFidu", "ElMu", "2jets");
  float N     = yield("data", Chan, Level);
  float Nbkg  = yield("bkg", Chan, Level);
  float xsec  = (N-Nbkg)/(ttbar/fidunorm)*xsec_th;
  if(docout){
    cout << " #### FIDUCIAL CROSS SECTION ####\n" << endl;
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
