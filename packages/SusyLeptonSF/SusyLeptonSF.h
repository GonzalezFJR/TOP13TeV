///////////////////////////////////////////////////////////////////////
//
//    FILE: SusyLeptonSF.h
//   CLASS: SusyLeptonSF
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


class SusyLeptonSF {
	public:
		SusyLeptonSF(bool loadhistos = true);
		~SusyLeptonSF() {}

  //----------------------------------------------------------------------  
  //--- Get Muon SFs
  //----------------------------------------------------------------------  

  double GetTightMuonIDSF (double pt, double eta) { return fTightMuonIDSF ->GetBinContent(fTightMuonIDSF->FindBin(TMath::Abs(eta), pt)); }
  double GetTightMuonIsoSF(double pt, double eta) { return fTightMuonIsoSF->GetBinContent(fTightMuonIsoSF->FindBin(TMath::Abs(eta), pt));}

  double GetTrackerMuonSF(double eta){
    float val = 0; float xlow = 0; float xhigh = 0;
    for(int i = 0; i < fTrackerMuonSF->GetN(); i++){
      xlow  = fTrackerMuonSF->GetX()[i] - fTrackerMuonSF->GetErrorXlow(i);
      xhigh = fTrackerMuonSF->GetX()[i] + fTrackerMuonSF->GetErrorXhigh(i);
      if(xlow <= eta && xhigh > eta) val = fTrackerMuonSF->GetY()[i];
    }
    return val;
  }
  double GetTightMuonSF(double pt, double eta){
    if(pt >= 120) pt = 119.9;
    float trackerSF = GetTrackerMuonSF(eta);
    eta = TMath::Abs(eta);
    float SFid  = fMediumMuonSF  ->GetBinContent(fMediumMuonSF  ->FindBin(pt,eta));
    float SFiso = fTightIsoMuonSF->GetBinContent(fTightIsoMuonSF->FindBin(pt,eta));
    return SFid*SFiso*trackerSF;
  }
  double GetTightMuonSF_err(double pt, double eta) {
    // Instructions: 0.3% for each muon...
    eta = TMath::Abs(eta);
    float SF = GetTightMuonSF(pt, eta);
    return SF*0.03;
  }

  double GetTightMuonSFPOG   (double pt, double eta) { 
    if(pt > 120) pt = 119.9;
    float SF_BCDEF = 1; float SF_GH = 1;
    float SF_iso_BCDEF = 1; float SF_iso_GH = 1;
    float trackerSF = GetTrackerMuonSF(eta);
    eta = TMath::Abs(eta);
    //return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(pt,eta))*trackerSF; 
    SF_BCDEF = fMediumMuonSF_BCDEF->GetBinContent(fMediumMuonSF_BCDEF->FindBin(pt,eta));
    SF_GH    = fMediumMuonSF_GH   ->GetBinContent(fMediumMuonSF_GH   ->FindBin(pt,eta));
    SF_iso_BCDEF = fTightIsoMuonSF_BCDEF->GetBinContent(fTightIsoMuonSF_BCDEF->FindBin(pt,eta));
    SF_iso_GH    = fTightIsoMuonSF_GH   ->GetBinContent(fTightIsoMuonSF_GH   ->FindBin(pt,eta));
    SF_BCDEF *= SF_iso_BCDEF; SF_GH *= SF_iso_GH;
    return (SF_BCDEF*lumiBCDEF + SF_GH*lumiGH)/(lumiBCDEF+lumiGH)*trackerSF;
  }
  double GetTightMuonSFPOG_err(double pt, double eta) {
    if(pt > 120) pt = 119.9; eta = TMath::Abs(eta);
    float SF_BCDEF = 1; float SF_GH = 1;
    float SF_iso_BCDEF = 1; float SF_iso_GH = 1;
    //return fTightMuonSF->GetBinError(fTightMuonSF->FindBin(pt, TMath::Abs(eta)));
    SF_BCDEF = fMediumMuonSF_BCDEF->GetBinError(fMediumMuonSF_BCDEF->FindBin(pt,eta)); 
    SF_GH    = fMediumMuonSF_GH   ->GetBinError(fMediumMuonSF_GH   ->FindBin(pt,eta));
    SF_iso_BCDEF = fTightIsoMuonSF_BCDEF->GetBinError(fTightIsoMuonSF_BCDEF->FindBin(pt,eta));
    SF_iso_GH    = fTightIsoMuonSF_GH   ->GetBinError(fTightIsoMuonSF_GH   ->FindBin(pt,eta));
    SF_BCDEF = TMath::Sqrt(SF_BCDEF*SF_BCDEF + SF_iso_BCDEF*SF_iso_BCDEF);
    SF_GH = TMath::Sqrt(SF_GH*SF_GH + SF_iso_GH*SF_iso_GH);
    return (SF_BCDEF*lumiBCDEF + SF_GH*lumiGH)/(lumiBCDEF+lumiGH);
  }

	float GetFastSimMuonSF(float pt, float eta){
		if(pt > 200) pt = 199.9;
		return fFastSimMuonSF->GetBinContent(fFastSimMuonSF->FindBin(pt, eta));	
	}
	float GetFastSimMuonSF_err(float pt, float eta){
		if(pt > 200) pt = 199.9;
		return fFastSimMuonSF->GetBinError(fFastSimMuonSF->FindBin(pt, eta));	
	}

  
  //----------------------------------------------------------------------  
  //--- Get Electron SFs
  //----------------------------------------------------------------------  

  float GetTightElectronIDSF  (float pt, float eta) const{ return fTightElectronIDSF->GetBinContent(fTightElectronIDSF->FindBin(TMath::Abs(eta), pt)); }

  float GetTrackerElectronSF(float eta){ return fTrackerElectronSF->GetBinContent(fTrackerElectronSF->FindBin(eta, 50.));	}
  float GetTrackerElectronSF_err(float eta){ return fTrackerElectronSF->GetBinError(fTrackerElectronSF->FindBin(eta, 50.));	}

	float GetFastSimElectronSF(float pt, float eta){
		if(pt > 200) pt = 199.9;
		return fFastSimElectronSF->GetBinContent(fFastSimElectronSF->FindBin(pt, eta));
	}
	float GetFastSimElectronSF_err(float pt, float eta){
		if(pt > 200) pt = 199.9;
		return fFastSimElectronSF->GetBinError(fFastSimElectronSF->FindBin(pt, eta));
	}

  float GetTightElectronSF  (float pt, float eta) { // binned in eta, pt
    if(pt > 200) pt = 199.9;
    float trackerSF = GetTrackerElectronSF(eta);
    float idsf = fTightElectronSF->GetBinContent(fTightElectronSF->FindBin(pt, eta));
    float isosf = fTightElectronIsoSF->GetBinContent(fTightElectronIsoSF->FindBin(pt, eta));
    return idsf*isosf*trackerSF;    
  }
  double GetTightElectronSF_err(double pt, double eta) {
    if(pt > 200) pt = 199.9;
    float trackerSF_err = GetTrackerElectronSF_err(eta);
    float iso_err = fTightElectronIsoSF->GetBinError(fTightElectronIsoSF->FindBin(pt, eta));
    float id_err = fTightElectronSF->GetBinError(fTightElectronSF->FindBin(eta, pt));
    return TMath::Sqrt(trackerSF_err*trackerSF_err+id_err*id_err+iso_err*iso_err);
  }

  /// Methods to get the ERROR 

  
  //----------------------------------------------------------------------  
	//--- Get Trigger SFs               (binned eta1, eta2)
	//----------------------------------------------------------------------  
	float GetDoubleMuSF(float pt, float eta) const {
		if(pt > 200) pt = 199.9;
		return fDoubleMuSF->GetBinContent(fDoubleMuSF->FindBin(TMath::Abs(pt), TMath::Abs(eta)) ); }
	float GetDoubleElSF(float pt, float eta) const { 
		if(pt > 200) pt = 199.9;
		return fDoubleElSF->GetBinContent(fDoubleElSF->FindBin(TMath::Abs(pt), TMath::Abs(eta)) ); }
	float GetMuEGSF    (float pt, float eta) const { 
		if(pt > 200) pt = 199.9;
		return fMuEGSF    ->GetBinContent(fMuEGSF    ->FindBin(TMath::Abs(pt), TMath::Abs(eta)) ); }
	// Trigger SF Errors (binned eta1, eta2)
	float GetDoubleMuSF_err(float pt, float eta) const { 
		if(pt > 200) pt = 199.9;
		return fDoubleMuSF->GetBinError(fDoubleMuSF->FindBin(TMath::Abs(pt),TMath::Abs(eta)));}
	float GetDoubleElSF_err(float pt, float eta) const { 
		if(pt > 200) pt = 199.9;
		return fDoubleElSF->GetBinError(fDoubleElSF->FindBin(TMath::Abs(pt),TMath::Abs(eta)));}
	float GetMuEGSF_err(float pt, float eta)     const {
		if(pt > 200) pt = 199.9;
		return fMuEGSF    ->GetBinError(fMuEGSF    ->FindBin(TMath::Abs(pt),TMath::Abs(eta)));}



	//######################################################################
	//### Methods to load histograms with Scale Factors
  //######################################################################
  // Muon SFs https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2#Results_for_CMSSW_7_6_X_dataset
  TH2D* LoadTightMuonIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/MuonIDSF.root", const char* histo = "DATA_over_MC_Tight");
  TH2D* LoadTightMuonIsoSF(const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/MuonISOSF.root", const char* histo = "DATA_over_MC_Tight");
  //TH2D* LoadTightMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Muons/GlobalMuonSFs_T.root",  const char* histo = "GlobalSF_T");
  TH2D* LoadTightMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/SUS_MuonMediumIdM17.root",  const char* histo = "SF");
  TH2D* LoadMediumMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/SUS_MuonMediumIdM17.root",  const char* histo = "SF");
  TH2D* LoadTightIsoMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/SUS_MuonRelIso03b012.root",  const char* histo = "pt_abseta_ratio");
  TH2D* LoadMediumMuonSF_BCDEF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/MuonSFMediumId_BCDEF.root",  const char* histo = "MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
  TH2D* LoadMediumMuonSF_GH(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/MuonSFMediumId_GH.root",  const char* histo = "MC_NUM_MediumID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio");
  TH2D* LoadTightIsoMuonSF_BCDEF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/MuonSFTightIso_BCDEF.root",  const char* histo = "TightISO_MediumID_pt_eta/pt_abseta_ratio");
  TH2D* LoadTightIsoMuonSF_GH(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/MuonSFTightIso_GH.root",  const char* histo = "TightISO_MediumID_pt_eta/pt_abseta_ratio");
 // TGraphAsymmErrors* LoadTrackerMuonSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Muons/ratios.root",	   const char* histo = "ratio_eta");
  TGraphAsymmErrors* LoadTrackerMuonSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/Tracking_EfficienciesAndSF_BCDEFGH.root",	   const char* histo = "ratio_eff_eta3_dr030e030_corr");
	//"/nfs/fanae/user/palencia/sfs13TeV/2016/GlobalMuonSFs.root"
  
  // Elec SFs
  // http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/moriond2016_76X/eleCutBasedID/
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
  TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronMVAIDIsoSF.root",                 const char* histo = "electronsDATAMCratio_FO_ID_ISO");
  //TH2D* LoadTrackerElectronSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Electrons/trakingElectronSF.root",             const char* histo = "EGamma_SF2D");
  TH2D* LoadTrackerElectronSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/ElecRecoM17.root",             const char* histo = "EGamma_SF2D");
  TH2D* LoadFastSimElectronSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Electrons_FastSim/sf_el_tightCB_MultiVT.root", const char* histo = "histo2D");
  TH2D* LoadFastSimMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Muons_FastSim/sf_mu_mediumID_multiVT.root",         const char* histo = "histo2D");
  //TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/Electrons/GlobalElectronSFs_T.root",            const char* histo = "GlobalSF_T");
  TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/SUS_electonSFs.root",const char* histo = "GsfElectronToCutBasedSpring15T");
  TH2D* LoadTightElectronIsoSF (const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/M17/SUS_electonSFs.root",            const char* histo = "CutBasedStopsDileptonToRelIso012");
  //"/nfs/fanae/user/palencia/sfs13TeV/25ns/elec_tight_sf2D_13TeV_RunD.root",
  //"/nfs/fanae/user/palencia/sfs13TeV/2016/egammaEffiMedium.txt_SF2D_runBCD_12p9fb.root"

  // + Trigger SFs
  TH2F* LoadDoubleElectronSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/trigger/triggEff.root", const char* histo = "eff_ee");
  TH2F* LoadDoubleMuonSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/trigger/triggEff.root",   const char* histo = "eff_mm");
  TH2F* LoadElMuSF(const char* file = "/nfs/fanae/user/juanr/SFs/StopSFs/trigger/triggEff.root",          const char* histo = "eff_OF");

  // Plot all histograms
  TCanvas* Draw();

 //protected:
  TH2D* GetHistogramFromFileD(const char* file, const char* histo, const char* hname);
  TH2F* GetHistogramFromFileF(const char* file, const char* histo, const char* hname) const;
  
 private:
  // Muon SFs
  TH2D*  fTightMuonIDSF;     //Muon ID Scale Factors
  TGraphAsymmErrors*  fTrackerMuonSF;     //Muon ID Scale Factors
  TH2D*  fTightMuonIsoSF;    //Muon Isolation Scale Factors
  TH2D*  fTightMuonSF;    //Muon Isolation Scale Factors
  TH2D*  fMediumMuonSF;  
  TH2D*  fTightIsoMuonSF;  
  TH2D*  fMediumMuonSF_BCDEF;  
  TH2D*  fMediumMuonSF_GH;   
  TH2D*  fTightIsoMuonSF_BCDEF;  
  TH2D*  fTightIsoMuonSF_GH;   

  // Electron SFs
  TH2D*  fTightElectronIDSF; //Electron ID Scale Factors
  TH2D*  fTightElectronSF; //Electron ID Scale Factors
  TH2D*  fTightElectronIsoSF;
  TH2D*  fTrackerElectronSF;
	TH2D*  fFastSimElectronSF;
	TH2D*  fFastSimMuonSF;

  // Trigger SFs
  TH2F *fDoubleMuSF;
  TH2F *fDoubleElSF;
  TH2F *fMuEGSF;
};
#endif
