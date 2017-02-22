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

#include "TH2D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"

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
  
  // Electron SFs
  float GetTightElectronIDSF  (float pt, float eta) const { // binned in eta, pt
    return fTightElectronIDSF->GetBinContent(fTightElectronIDSF->FindBin(TMath::Abs(eta), pt));
  }
  
  
  double GetTightMuonSF   (double pt, double eta) { 
    if(pt > 120) pt = 119.9;
    eta = TMath::Abs(eta);
    return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt)); 
    
    double statErr = fTightMuonSF->GetBinError(fTightMuonSF->FindBin(eta,pt));
    double systErr = 0.015 * fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt));
    double totErr = sqrt(statErr*statErr + systErr*systErr);
     
    //return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt)) + fTightMuonSF->GetBinError(fTightMuonSF->FindBin(eta,pt));  
    //return fTightMuonSF->GetBinContent(fTightMuonSF->FindBin(eta, pt)) - totErr;  
  }

  float GetTightElectronSF  (float pt, float eta) const { // binned in eta, pt
    //if(pt > 500) pt = 499.9;
    if(pt > 200) pt = 199.9;
    eta = TMath::Abs(eta);
    return fTightElectronSF->GetBinContent(fTightElectronSF->FindBin(eta, pt));    
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
  double GetTightMuonSF_err(double pt, double eta) {
    return fTightMuonSF->GetBinError(fTightMuonSF->FindBin(pt,TMath::Abs(eta)));
  }
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
  TH2D* LoadTightMuonSF(const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/76X/MuonSF_IdIso_Z_RunCD_Reco76X_Feb15.root",
  			  const char* histo = "GlobalSF");
  
  // + Elec SFs
  // http://fcouderc.web.cern.ch/fcouderc/EGamma/scaleFactors/moriond2016_76X/eleCutBasedID/
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Electron_efficiencies_and_scale
  //TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronSF.root",
  //	  		        const char* histo = "sfTIGHT");
  TH2D* LoadTightElectronIDSF (const char* file = "http://www.hep.uniovi.es/iglez/CMS/WZ/ElectronMVAIDIsoSF.root",
			       const char* histo = "electronsDATAMCratio_FO_ID_ISO");
  
  //TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/elec_tight_sf2D_13TeV_RunD.root",
 // TH2D* LoadTightElectronSF (const char* file = "/nfs/fanae/user/palencia/sfs13TeV/25ns/76X/CutBasedID_MediumWP_76X_18Feb.txt_SF2D.root", const char* histo = "EGamma_SF2D");
//  TH2D* LoadTightElectronSF (const char* file = "/mnt_pool/fanae105/user/juanr/TOP13TeV/CutBasedID_LooseWP_76X_18Feb.txt_SF2D.root", const char* histo = "EGamma_SF2D");
  TH2D* LoadTightElectronSF (const char* file = "CutBasedID_LooseWP_fromTemplates_withSyst_Final.txt_SF2D.root", const char* histo = "EGamma_SF2D");

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
  TH2D*  fTightMuonIsoSF;    //Muon Isolation Scale Factors
  TH2D*  fTightMuonSF;    //Muon Isolation Scale Factors
  TH2D*  fTightElectronIDSF; //Electron ID Scale Factors
  TH2D*  fTightElectronSF; //Electron ID Scale Factors

  // Trigger SFs
  TH2F *fDoubleMuSF;
  TH2F *fDoubleElSF;
  TH2F *fMuEGSF;
};
#endif
