#include "DrawPlots.C"

//double getYield(TString sample, TString Chan, TString Level);
//double yield(TString Process, TString Chan, TString Level);

float getSyst(TString tag, TString direc = "Max"); // Returns the variation in the yield due to a systematic uncertainty
float getExpError(TString Process = "ttbar", TString tag = "0", TString direc = "Max"); // Returns the variation in the ttbat yield due to a experimental uncertainty
float getAcc(TString tag = "0", TString direc = "Max");
float getXsec(TString tag = "0", TString direc = "Max");

float getPDFsyst(); 
float getScaleMESyst();
float getNLOsyst();
float getHadSyst(); 
float getScalePSSyst(TString direc = "Max"); 
void PrintXsecSyst();
void PrintGoodYields();
void GetXSec();

const TString Level = "2jets";
const TString Chan = "ElMu"; 
const float xsec_th = 67.8;
const float lumiunc = 0.12;

const float tWerr = 0.3; const float VVerr = 0.3; const float DYerr = 0.15; const float NonWerr = 0.3;
const TString up = "Up"; const TString down = "Down"; const TString maxi = "Max";


float getSyst(TString tag, TString direc ){ // Returns the variation in the yield due to a systematic uncertainty
  // Modeling
  if     (tag == "PDF" || tag == "pdf")                                       return getPDFsyst(); // PDF
  else if(tag == "ME" || tag == "scaleME" || tag == "ScalesME")               return getScaleMESyst(); // Scale ME
  else if(tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales")  return getScalePSSyst(direc); // Scale PS
	else if(tag == "had" || tag == "Had")                                       return getHadSyst(); // Hadronization
  else if(tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator") return getNLOsyst(); // NLO
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
  else if(tag == "JES" || tag == "jes")                 { tup = "JESUp"; tdown = "JESDown";}
  else if(tag == "JER" || tag == "jer")                 { tup = "JER"; tdown = "JER";}
  float yup   = yield(Process, Chan, Level, tup);  
  float ydown = yield(Process, Chan, Level, tdown);  
  if     (direc == maxi ) return (float )TMath::Max( abs(ynom - yup), abs(ynom-ydown) );
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
	else if(tag.Contains("trig") || tag.Contains("Trig") || tag.Contains("Elec") || tag.Contains("elec") || tag.Contains("Muon") || tag.Contains("muon") || tag == "JES" || tag == "jes" || tag == "JER" || tag == "jer" || tag == "PDF" || tag == "pdf" || tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales" || tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator" || tag == "had" || tag == "Had") return xsec*abs(getSyst(tag, direc)-ttbar)/ttbar; 
  else if(tag == "stat") (direc == down)? N -= TMath::Sqrt(N) : N += TMath::Sqrt(N);
  else if(tag == "lumi" || tag == "Lumi") return xsec*(1+lumiunc);
	return (N-Nbkg)/ttbar*xsec_th;
}

float getAcc(TString tag, TString direc){
	float y = yield("ttbar", Chan, Level);
	if(tag.Contains("trig") || tag.Contains("Trig") || tag.Contains("Elec") || tag.Contains("elec") || tag.Contains("Muon") || tag.Contains("muon") || tag == "JES" || tag == "jes" || tag == "JER" || tag == "jer" || tag == "PDF" || tag == "pdf" || tag == "QCD" || tag == "qcd" || tag == "scale" || tag == "scales" || tag == "NLO" || tag == "nlo" || tag == "gen" || tag == "generator" || tag == "had" || tag == "Had") y += getSyst(tag, direc);
	return y/(xsec_th*Lumi);
}

float getPDFsyst(){
	float pdferror = 0.;
	return pdferror;
} 

float getScaleMESyst(){
	float MEerror = 0.;
	return MEerror;
} 
 
float getScalePSSyst(TString direc){
	float scaleerror = 0;
  float nom = yield("ttbar", Chan, Level);
  float eup   = getYield("TTbar_Powheg_ScaleUp", Chan, Level);
  float edown = getYield("TTbar_Powheg_ScaleDown", Chan, Level);
  if     (direc == maxi) scaleerror = (float ) TMath::Max(abs(nom-eup), abs(nom-edown));
  else if(direc == up)   scaleerror = (eup-nom);
  else if(direc == down) scaleerror = (edown-nom);
	return scaleerror;
}

float getHadSyst(){ 
	float haderror = 0.;
	return haderror;
  float nom = yield("ttbar", Chan, Level);
  float var = getYield("TTbar_Powheg_Herwig", Chan, Level);
  haderror = abs(nom - var);
	return haderror;
}

float getNLOsyst(){
	float nloerror = 0.;
  float nom = yield("ttbar", Chan, Level);
  float var = getYield("TTJets_aMCatNLO", Chan, Level);
  nloerror = abs(nom - var);
	return nloerror;
} 

void PrintXsecSyst(){
  float xsec = getXsec();
  float totalSys = TMath::Sqrt( (abs(xsec-getXsec("Trig"))*abs(xsec-getXsec("Trig")) + abs(xsec-getXsec("Elec"))*abs(xsec-getXsec("Elec")) + abs(xsec-getXsec("Muon"))*abs(xsec-getXsec("Muon")) + abs(xsec-getXsec("JES"))*abs(xsec-getXsec("JES")) + abs(xsec-getXsec("JER"))*abs(xsec-getXsec("JER")) + abs(xsec-getXsec("scale"))*abs(xsec-getXsec("scale")) + abs(xsec-getXsec("ME"))*abs(xsec-getXsec("ME")) + abs(xsec-getXsec("pdf"))*abs(xsec-getXsec("pdf")) + abs(xsec-getXsec("had"))*abs(xsec-getXsec("had")) + abs(xsec-getXsec("tW"))*abs(xsec-getXsec("tW")) + abs(xsec-getXsec("VV"))*abs(xsec-getXsec("VV")) + abs(xsec-getXsec("DY"))*abs(xsec-getXsec("DY")) + abs(xsec-getXsec("NonW"))*abs(xsec-getXsec("NonW")) )/(xsec*xsec));
  float stat = abs(xsec-getXsec("stat"))/xsec;
  float lumin = abs(xsec-getXsec("lumi"))/xsec;
  float total = TMath::Sqrt(totalSys*totalSys + lumin*lumin + stat*stat);
  cout << "==================================================================================" << endl;
  cout << " Source                        #Delta #sigma_{t#bar{t}}(pb)       (   %)          " << endl;
  cout << "----------------------------------------------------------------------------------" << endl;
  cout << Form(" Trigger efficiencies                   %2.3f                      %1.2f          ",abs(xsec-getXsec("Trig")) , abs(xsec-getXsec("Trig"))/xsec*100 ) << endl;
  cout << Form(" Electron efficiencies                  %2.3f                      %1.2f          ",abs(xsec-getXsec("Elec" )), abs(xsec-getXsec("Elec"))/xsec*100 ) << endl;
  cout << Form(" Muon efficiencies                      %2.3f                      %1.2f          ",abs(xsec-getXsec("Muon")) , abs(xsec-getXsec("Muon"))/xsec*100 ) << endl;
  cout << Form(" Jet Energy Scale                       %2.3f                      %1.2f          ",abs(xsec-getXsec("JES")) , abs(xsec-getXsec("JES"))/xsec*100 ) << endl;
  cout << Form(" Jer Energy Resolution                  %2.3f                      %1.2f          ",abs(xsec-getXsec("JER")) , abs(xsec-getXsec("JER"))/xsec*100 ) << endl;
  cout << Form(" QCD scale PS                           %2.3f                      %1.2f          ",abs(xsec-getXsec("scale")) , abs(xsec-getXsec("scale"))/xsec*100 ) << endl;
  cout << Form(" QCD scale ME                           %2.3f                      %1.2f          ",abs(xsec-getXsec("ME")) , abs(xsec-getXsec("ME"))/xsec*100 ) << endl;
  cout << Form(" PDF                                    %2.3f                      %1.2f          ",abs(xsec-getXsec("pdf")) , abs(xsec-getXsec("pdf"))/xsec*100 ) << endl;
  cout << Form(" Hadronization                          %2.3f                      %1.2f          ",abs(xsec-getXsec("had")) , abs(xsec-getXsec("had"))/xsec*100 ) << endl;
  cout << "----------------------------------------------------------------------------------" << endl;
  cout << Form(" tW background                          %2.3f                      %1.2f          ",abs(xsec-getXsec("tW")) , abs(xsec-getXsec("tW"))/xsec*100 ) << endl;
  cout << Form(" VV background                          %2.3f                      %1.2f          ",abs(xsec-getXsec("VV")) , abs(xsec-getXsec("VV"))/xsec*100 ) << endl;
  cout << Form(" DY background                          %2.3f                      %1.2f          ",abs(xsec-getXsec("DY")) , abs(xsec-getXsec("DY"))/xsec*100 ) << endl;
  cout << Form(" Non prompt backgrounds                 %2.3f                      %1.2f          ",abs(xsec-getXsec("NonW")) , abs(xsec-getXsec("NonW"))/xsec*100 ) << endl;
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
  cout << Form(" Data              |     %5.0f       |   %5.0f         " , yield("data", Chan, "dilepton"), yield("data", Chan, "2jets")) << endl;
	cout <<      "============================================================" << endl;
}

void GetXSec(){
  TString Chan = "ElMu";
  float N_bkg = yield("DY", Chan, "2jets")+yield("fake", Chan, "2jets")+yield("tW", Chan, "2jets")+yield("VV", Chan, "2jets");
  float xsec = getXsec();
  float totalSys = TMath::Sqrt( (abs(xsec-getXsec("Trig"))*abs(xsec-getXsec("Trig")) + abs(xsec-getXsec("Elec"))*abs(xsec-getXsec("Elec")) + abs(xsec-getXsec("Muon"))*abs(xsec-getXsec("Muon")) + abs(xsec-getXsec("JES"))*abs(xsec-getXsec("JES")) + abs(xsec-getXsec("JER"))*abs(xsec-getXsec("JER")) + abs(xsec-getXsec("scale"))*abs(xsec-getXsec("scale")) + abs(xsec-getXsec("ME"))*abs(xsec-getXsec("ME")) + abs(xsec-getXsec("pdf"))*abs(xsec-getXsec("pdf")) + abs(xsec-getXsec("had"))*abs(xsec-getXsec("had")) + abs(xsec-getXsec("tW"))*abs(xsec-getXsec("tW")) + abs(xsec-getXsec("VV"))*abs(xsec-getXsec("VV")) + abs(xsec-getXsec("DY"))*abs(xsec-getXsec("DY")) + abs(xsec-getXsec("NonW"))*abs(xsec-getXsec("NonW")) ));
  float stat = abs(xsec-getXsec("stat"));
  float lumin = abs(xsec-getXsec("lumi"));
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
}
