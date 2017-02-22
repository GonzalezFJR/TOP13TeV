///////////////////////////////////////////////////////////////////////
//
//    FILE: LeptonSF.h
//   CLASS: LeptonSF
// AUTHORS: S. Folgueras, I. Gonzalez
//    DATE: February, 2014
//
// CONTENT: An utility class to handle Lepton and Trigger Scale Factors.
//
///////////////////////////////////////////////////////////////////////
#ifndef LEPTONSF_H
#define LEPTONSF 1


#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"

const Float_t lumiBCDEF = 19.706;
const Float_t lumiGH    = 16.1454;

class LeptonSF {
 public:
  LeptonSF(bool loadhistos = true);
  ~LeptonSF() {}

  
  
  // Muon SFs
  double GetTightMuonIDSF (double pt, double eta) { 
    return fTightMuonIDSF->GetBinContent(fTightMuonIDSF->FindBin(TMath::Abs(eta), pt));
  }
  double GetTightMuonIsoSF(double pt, double eta) { 
    return fTightMuonIsoSF->GetBinContent(fTightMuonIsoSF->FindBin(TMath::Abs(eta), pt));
  }
  double GetTrackerMuonSF(double eta){
    float val = 0; float xlow = 0; float xhigh = 0;
    for(int i = 0; i < fTrackerMuonSF->GetN(); i++){
      xlow  = fTrackerMuonSF->GetX()[i] - fTrackerMuonSF->GetErrorXlow(i);
      xhigh = fTrackerMuonSF->GetX()[i] + fTrackerMuonSF->GetErrorXhigh(i);
      if(xlow <= eta && xhigh > eta) val = fTrackerMuonSF->GetY()[i];
    }
    return val;
  }

  
  // Electron SFs
  float GetTightElectronIDSF  (float pt, float eta) const { // binned in eta, pt
    return fTightElectronIDSF->GetBinContent(fTightElectronIDSF->FindBin(TMath::Abs(eta), pt));
  }
  float GetTrackerElectronSF(float eta){
    return fTrackerElectronSF->GetBinContent(fTrackerElectronSF->FindBin(eta, 50.));
  }
  
  
/*  double GetTightMuonSF   (double pt, double eta) { 
    if(pt > 120) pt = 119.9;
    float trackerSF = GetTrackerMuonSF(eta);
    eta = TMath::Abs(eta);
    return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt))*trackerSF; 
    
    //double statErr = fTightMuonSF->GetBinError(fTightMuonSF->FindBin(eta,pt));
    //double systErr = 0.015 * fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt));
    //double totErr = sqrt(statErr*statErr + systErr*systErr);
     
    return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt)) + fTightMuonSF->GetBinError(fTightMuonSF->FindBin(eta,pt));  
    //return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt)) - totErr;  
  }*/
  double GetTightMuonSF   (double pt, double eta) { 
    if(pt > 120) pt = 119.9;
    float SF_BCDEF = 1; float SF_GH = 1;
    float SF_iso_BCDEF = 1; float SF_iso_GH = 1;
    //float trackerSF = GetTrackerMuonSF(eta);
    eta = TMath::Abs(eta);
    //return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(pt,eta))*trackerSF; 
    SF_BCDEF = fTightMuonSF_BCDEF->GetBinContent(fTightMuonSF_BCDEF->FindBin(pt,eta));
    SF_GH    = fTightMuonSF_GH   ->GetBinContent(fTightMuonSF_GH   ->FindBin(pt,eta));
    SF_iso_BCDEF = fTightIsoMuonSF_BCDEF->GetBinContent(fTightIsoMuonSF_BCDEF->FindBin(pt,eta));
    SF_iso_GH    = fTightIsoMuonSF_GH   ->GetBinContent(fTightIsoMuonSF_GH   ->FindBin(pt,eta));
    SF_BCDEF *= SF_iso_BCDEF; SF_GH *= SF_iso_GH;
    return (SF_BCDEF*lumiBCDEF + SF_GH*lumiGH)/(lumiBCDEF+lumiGH);
  }
  double GetTightMuonSF_err(double pt, double eta) {
    if(pt > 120) pt = 119.9; eta = TMath::Abs(eta);
    float SF_BCDEF = 1; float SF_GH = 1;
    float SF_iso_BCDEF = 1; float SF_iso_GH = 1;
    //return fTightMuonSF->GetBinError(fTightMuonSF->FindBin(pt, TMath::Abs(eta)));
    SF_BCDEF = fTightMuonSF_BCDEF->GetBinError(fTightMuonSF_BCDEF->FindBin(pt,eta)); 
    SF_GH    = fTightMuonSF_GH   ->GetBinError(fTightMuonSF_GH   ->FindBin(pt,eta));
    SF_iso_BCDEF = fTightIsoMuonSF_BCDEF->GetBinError(fTightIsoMuonSF_BCDEF->FindBin(pt,eta));
    SF_iso_GH    = fTightIsoMuonSF_GH   ->GetBinError(fTightIsoMuonSF_GH   ->FindBin(pt,eta));
    SF_BCDEF = TMath::Sqrt(SF_BCDEF*SF_BCDEF + SF_iso_BCDEF*SF_iso_BCDEF);
    SF_GH = TMath::Sqrt(SF_GH*SF_GH + SF_iso_GH*SF_iso_GH);
    return (SF_BCDEF*lumiBCDEF + SF_GH*lumiGH)/(lumiBCDEF+lumiGH);
  }

  float GetTightElectronSF  (float pt, float eta) { // binned in eta, pt
    //if(pt > 500) pt = 499.9;
    if(pt > 200) pt = 199.9;
    float trackerSF = GetTrackerElectronSF(eta);
    //trackerSF = 1;
    return fTightElectronSF->GetBinContent(fTightElectronSF->FindBin(eta,pt))*trackerSF;    
    //return fTightElectronSF->GetBinContent(fTightElectronSF->FindBin(eta, pt)) - fTightElectronSF->GetBinError(fTightElectronSF->FindBin(eta, pt));
  }

  
  // Trigger SFs
  float GetDoubleMuSF(float eta1, float eta2) const { // binned in eta1, eta2
    eta1 = TMath::Abs(eta1);
    eta2 = TMath::Abs(eta2);
    return fDoubleMuSF->GetBinContent(fDoubleMuSF->FindBin(eta1, eta2));
    //return fDoubleMuSF->GetBinContent(fDoubleMuSF->FindBin(eta1, eta2)) - fDoubleMuSF->GetBinError(fDoubleMuSF->FindBin(eta1, eta2));
  }
  float GetDoubleElSF(float eta1, float eta2) const { // binned in eta1, eta2
    eta1 = TMath::Abs(eta1);
    eta2 = TMath::Abs(eta2);
    return fDoubleElSF->GetBinContent(fDoubleElSF->FindBin(eta1, eta2));
    //return fDoubleElSF->GetBinContent(fDoubleElSF->FindBin(eta1, eta2)) - fDoubleElSF->GetBinError(fDoubleElSF->FindBin(eta1, eta2));
  }
  float GetMuEGSF    (float eta1, float eta2) const { // binned in eta1, eta2
    eta1 = TMath::Abs(eta1);
    eta2 = TMath::Abs(eta2);
    return fMuEGSF->GetBinContent(fMuEGSF->FindBin(eta1, eta2) );
    //return fMuEGSF->GetBinContent(fMuEGSF->FindBin(eta1, eta2)) - fMuEGSF->GetBinError(fMuEGSF->FindBin(eta1, eta2));
  }
  
  
  // Trigger SF errors
  float GetDoubleMuSF_err(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleMuSF->GetBinError(fDoubleMuSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetDoubleElSF_err(float eta1, float eta2) const { // binned in eta1, eta2
    return fDoubleElSF->GetBinError(fDoubleElSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }
  float GetMuEGSF_err(float eta1, float eta2) const { // binned in eta1, eta2
    return fMuEGSF->GetBinError(fMuEGSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
  }

  

  
  /// Methods to get the ERROR 
/*  double GetTightMuonSF_err(double pt, double eta) {
    return fTightMuonSF->GetBinError(fTightMuonSF->FindBin(pt,TMath::Abs(eta)));
  }*/
  double GetTightElectronSF_err(double pt, double eta) {
    return fTightElectronSF->GetBinError(fTightElectronSF->FindBin(eta, pt));
  }
  
  // Methods to load Scale Factors
  // + Muon SFs https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2#Results_for_CMSSW_7_6_X_dataset
  // MuonIso_Z_RunCD_Reco76X_Feb15.root --> MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin/abseta_pt_ratio
  // MuonID_Z_RunCD_Reco76X_Feb15.root  --> MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_spliteta_bin/abseta_pt_ratio
  TH2D* LoadTightMuonIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/MuonIDSF.root",
			   const char* histo = "DATA_over_MC_Tight");
  TH2D* LoadTightMuonIsoSF(const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/MuonISOSF.root",
			   const char* histo = "DATA_over_MC_Tight");
  //TH2D* LoadTightMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Muons/GlobalMuonSFs_T.root",  const char* histo = "GlobalSF_T");
  //TGraphAsymmErrors* LoadTrackerMuonSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Muons/ratios.root",	   const char* histo = "ratio_eta");
  TGraphAsymmErrors* LoadTrackerMuonSF (const char* file = "/nfs/fanae/user/palencia/sfs13TeV/2016/muonTrackerEffRatios.root",	   const char* histo = "ratio_eta");
					TH2D* LoadTightMuonSF(const char* file = "/nfs/fanae/user/palencia/sfs13TeV/2016/GlobalMuonSFs.root",		const char* histo = "GlobalSF");
 

  TH2D* LoadTightMuonSF_BCDEF(const char* file = "/nfs/fanae/user/juanr/SFs/topM17/MuonSFId_BCDEF.root",  const char* histo = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
  TH2D* LoadTightMuonSF_GH(const char* file = "/nfs/fanae/user/juanr/SFs/topM17/MuonSFId_GH.root",  const char* histo = "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
  TH2D* LoadTightIsoMuonSF_BCDEF(const char* file = "/nfs/fanae/user/juanr/SFs/topM17/MuonSFTightIso_BCDEF.root",  const char* histo = "TightISO_TightID_pt_eta/pt_abseta_ratio");
  TH2D* LoadTightIsoMuonSF_GH(const char* file = "/nfs/fanae/user/juanr/SFs/topM17/MuonSFTightIso_GH.root",  const char* histo = "TightISO_TightID_pt_eta/pt_abseta_ratio");



 
  // + Elec SFs
  // http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/moriond2016_76X/eleCutBasedID/
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
  //TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronSF.root",
  //	  		        const char* histo = "sfTIGHT");
  TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronMVAIDIsoSF.root",
			       const char* histo = "electronsDATAMCratio_FO_ID_ISO");
 //TH2D* LoadTrackerElectronSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Electrons/trakingElectronSF.root",
 TH2D* LoadTrackerElectronSF(const char* file = "/nfs/fanae/user/juanr/SFs/topM17/RecoElec.root",
         const char* histo = "EGamma_SF2D");
  
  //TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/elec_tight_sf2D_13TeV_RunD.root",
  //TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Electrons/GlobalElectronSFs_T.root",                          		     const char* histo = "GlobalSF_T");
TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/juanr/SFs/topM17/TightCBElec.root", const char* histo = "EGamma_SF2D");

  // + Trigger SFs
  TH2F* LoadDoubleElectronSF(const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/76X/triggerSummary_ee.root",
 			     const char* histo = "scalefactor_eta2d_with_syst");
  TH2F* LoadDoubleMuonSF(const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/76X/triggerSummary_mumu.root",
 			 const char* histo = "scalefactor_eta2d_with_syst");
  TH2F* LoadElMuSF(const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/76X/triggerSummary_emu.root",
 		   const char* histo = "scalefactor_eta2d_with_syst");


  // Plot all histograms
  TCanvas* Draw();

 protected:
  TH2D* GetHistogramFromFileD(const char* file, const char* histo, 
			      const char* hname);
  TH2F* GetHistogramFromFileF(const char* file, const char* histo, 
			      const char* hname) const;
  
 private:
  // Lepton SFs
  TH2D*  fTightMuonIDSF;     //Muon ID Scale Factors
  TGraphAsymmErrors*  fTrackerMuonSF;     //Muon ID Scale Factors
  TH2D*  fTightMuonIsoSF;    //Muon Isolation Scale Factors
  TH2D*  fTightMuonSF;    //Muon Isolation Scale Factors
  TH2D*  fTightElectronIDSF; //Electron ID Scale Factors
  TH2D*  fTightElectronSF; //Electron ID Scale Factors
  TH2D*  fTrackerElectronSF;
  TH2D*  fTightMuonSF_BCDEF;  
  TH2D*  fTightMuonSF_GH;   
  TH2D*  fTightIsoMuonSF_BCDEF;  
  TH2D*  fTightIsoMuonSF_GH;   

  // Trigger SFs
  TH2F *fDoubleMuSF;
  TH2F *fDoubleElSF;
  TH2F *fMuEGSF;
};
#endif
