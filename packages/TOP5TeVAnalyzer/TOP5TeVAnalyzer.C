//========================================================
//  TOP5TeVAnalyzer selector
//========================================================

#include "TOP5TeVAnalyzer.h"
#include <iostream>
#include <math.h>

ClassImp(TOP5TeVAnalyzer);
const float gJetEtCut = 25.;

//#define DEBUG

//------------------------------------------------------------------------------
// GetParameters
//------------------------------------------------------------------------------
void TOP5TeVAnalyzer::GetParameters(){
  gSampleName    = GetParam<TString>("sampleName");
  gIsData        = GetParam<bool>("IsData");
  gWeight        = GetParam<float>("weight"); // cross section / events in the sample
  gLumiForPU     = GetParam<float>("LumiForPU");
  gTotalLumi     = GetParam<float>("TotalLumi");
  gDoSystStudies = GetParam<bool>("DoSystStudies");
  gUseCSVM       = GetParam<bool>("UseCSVM");
  gSelection     = GetParam<Int_t>("Selection");
  gIsMCatNLO     = GetParam<bool>("IsMCatNLO");
  gDoSF = 1; //if(gSelection == 2 || gSelection == 4 || gSelection == 5) gDoSF = 0; else gDoSF = 1;
  gDoDF = 1; //if(gSelection == 3) gDoDF = 0;


  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gSampleName = %s",gSampleName.Data()));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gIsData = %d",gIsData ));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gWeight = %e", gWeight));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gLumiForPU = %f", gLumiForPU));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gTotalLumi = %f", gTotalLumi));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gDoSystStudies = %d", gDoSystStudies));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gUseCSVM = %d",gUseCSVM ));
  PAF_INFO("TOP5TeVAnalyzer::GetParameters()", Form("gSelection = %i", gSelection));
}

//-----------------------------------------------------------------------------------
// GetTreeVariables
//-----------------------------------------------------------------------------------
void TOP5TeVAnalyzer::GetTreeVariables(){
  evt = Get<Int_t>("EventNumber");
	if(gIsData){
		run = Get<Int_t>("RunNumber");
		lum = Get<Int_t>("LumiSection");
	}
  nElec                = Get<Int_t>("nElec");
  nMuon                = Get<Int_t>("nMuon");
  nJet                 = Get<Int_t>("nJet");
  if(!gIsData){
    ngenLep              = Get<Int_t>("ngenLep");
    ngenJet              = Get<Int_t>("ngenJet");
    genWeight            = Get<Float_t>("genWeight");
  }
  for(int k = 0; k<nMuon; k++){
    MuonPx[k]      = Get<Float_t>("MuonPx", k);
    MuonPy[k]      = Get<Float_t>("MuonPy", k);
    MuonPz[k]      = Get<Float_t>("MuonPz", k);
    MuonEnergy[k]  = Get<Float_t>("MuonEnergy", k);
    MuonPt[k]      = Get<Float_t>("MuonPt", k);
    MuonEta[k]     = Get<Float_t>("MuonEta", k);
    MuonCharge[k]  = Get<Int_t>("MuonCharge", k);
    MuonDxy[k]     = Get<Float_t>("MuonDxy", k);
    MuonDz[k]      = Get<Float_t>("MuonDz", k);
    MuonChi2NDF[k]    = Get<Float_t>("MuonChi2NDF", k);
    MuonPixelHits[k]  = Get<Int_t>("MuonPixelHits", k);
    MuonHits[k]  = Get<Int_t>("MuonHits", k);
    MuonStations[k]   = Get<Int_t>("MuonStations", k);
    MuonTrkLayers[k]  = Get<Int_t>("MuonTrkLayers", k);
  }
  for(int k = 0; k<nElec; k++){
    ElecPx[k]      = Get<Float_t>("ElecPx", k);
    ElecPy[k]      = Get<Float_t>("ElecPy", k);
    ElecPz[k]      = Get<Float_t>("ElecPz", k);
    ElecEnergy[k]  = Get<Float_t>("ElecEnergy", k);
    ElecPt[k]      = Get<Float_t>("ElecPt", k);
    ElecEta[k]     = Get<Float_t>("ElecEta", k);
    ElecCharge[k]  = Get<Int_t>("ElecCharge", k);
    ElecIDVeto[k]  = Get<Int_t>("ElecIDVeto", k);
    ElecIDLoose[k] = Get<Int_t>("ElecIDLoose", k);
    ElecIDMedium[k]= Get<Int_t>("ElecIDMedium", k);
    ElecIDTight[k] = Get<Int_t>("ElecIDTight", k);
	}
  for(int k = 0; k<nJet; k++){
    Jet_px[k]          = Get<Float_t>("Jet_px", k);
    Jet_py[k]          = Get<Float_t>("Jet_py", k);
    Jet_pz[k]          = Get<Float_t>("Jet_pz", k);
    Jet_energy[k]      = Get<Float_t>("Jet_energy", k);
    Jet_eta[k]         = Get<Float_t>("Jet_eta", k);
    Jet_btagCSV[k]     = Get<Float_t>("Jet_btagCSV", k);
    Jet_mcFlavour[k]   = Get<Int_t>("Jet_mcFlavour", k);
    rawpt[k]           = Get<Float_t>("RawPt", k);
    neutralSum[k]      = Get<Float_t>("NeutralSum", k);
    eSum[k]            = Get<Float_t>("ESum", k);
    chargedHardSum[k]  = Get<Float_t>("ChargedHardSum", k);
    chargedSum[k]      = Get<Float_t>("ChargedSum", k);
    chargedMult[k]     = Get<Int_t>("chargedMult", k);
    totalMult[k]       = Get<Int_t>("totalMult", k);

    CHF[k] = Get<Float_t>("JetPfCHF", k);
    NHF[k] = Get<Float_t>("JetPfNHF", k);
    CEF[k] = Get<Float_t>("JetPfCEF", k);
    NEF[k] = Get<Float_t>("JetPfNEF", k);
    MUF[k] = Get<Float_t>("JetPfMUF", k);
    CHM[k] = Get<Int_t>("JetPfCHM", k);
    NHM[k] = Get<Int_t>("JetPfNHM", k);
    CEM[k] = Get<Int_t>("JetPfCEM", k);
    NEM[k] = Get<Int_t>("JetPfNEM", k);
    MUM[k] = Get<Int_t>("JetPfMUM", k);
  }

  if(!gIsData){
    for(int k = 0; k<ngenLep; k++){
      genLep_pdgId[k]    = TMath::Abs(Get<Int_t>("genLep_pdgId", k));
      genLep_pt[k]       = Get<Float_t>("genLep_pt", k);
      genLep_eta[k]      = Get<Float_t>("genLep_eta", k);
      genLep_phi[k]      = Get<Float_t>("genLep_phi", k);
      genLep_energy[k]   = Get<Float_t>("genLep_energy", k);
      genLep_status[k]   = Get<Int_t>("genLep_status", k);
      genLep_MomPID[k]   = Get<Int_t>("genLep_MomPID", k);
      genLep_GMomPID[k]   = Get<Int_t>("genLep_GMomPID", k);
    }
    for(int k = 0; k < ngenJet; k++){
  genJet_pt[k]       = Get<Float_t>("genJet_pt", k);
  genJet_eta[k]      = Get<Float_t>("genJet_eta", k);;
  genJet_phi[k]      = Get<Float_t>("genJet_phi", k);;
  genJet_m[k]        = Get<Float_t>("genJet_m", k);;
  genJet_matchId[k]  = Get<Int_t>("genJet_matchId", k);;

    }

  }
  //met_pt      = Get<Float_t>("pfTypeIMET");
  //met_phi     = Get<Float_t>("pfTypeIMETPhi");
  met_pt      = Get<Float_t>("met_pt");
  met_phi     = Get<Float_t>("met_phi");
  met_trk      = Get<Float_t>("met_trk");
  met_phi_trk  = Get<Float_t>("met_phi_trk");
  TrueMet_pt  = Get<Float_t>("TrueMet_pt");
  TrueMet_phi = Get<Float_t>("TrueMet_phi");
  TrueMet_pt_noHF  = Get<Float_t>("TrueMet_pt_noHF");
  TrueMet_phi_noHF = Get<Float_t>("TrueMet_phi_noHF");
  HLT_HIL2Mu15_v1  = Get<Int_t>("HLT_HIL2Mu15_v1");
  HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1 = Get<Int_t>("HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1");
}

int getNFPtBins(gChannel chan){ // fake ratios
  if(chan == Muon || chan == ElMu) return gNMuFPtBins;
  if(chan == Elec)                 return gNElFPtBins;
  else return -99;
};

const double *getFPtBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuFPtBins;
  else                             return gElFPtBins;
  //  if(chan == Elec)                 return gElFPtBins;
  //  else return *-999;
};

int getNPPtBins(gChannel chan){
  if(chan == Muon || chan == ElMu) return gNMuPPtbins;
  if(chan == Elec)                 return gNElPPtbins;
  else return -99;
};

const double *getPPtBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuPPtbins;
  else                             return gElPPtbins;
  //  if(chan == Elec)                 return gElPPtbins;
  //  else return -99;
};

int getNEtaBins(gChannel chan){
  if(chan == Muon || chan == ElMu) return gNMuEtabins;
  if(chan == Elec)                 return gNElEtabins;
  else return -99;
};

const double *getEtaBins (gChannel chan){
  if(chan == Muon || chan == ElMu) return gMuEtabins;
  else                             return gElEtabins;
  //  if(chan == Elec)                 return gElEtabins;
  //  else return *-99.;
};

//------------------------------------------------------------------------------------
// TOP5TeVAnalyzer class constructor (make sure the pointers are initialized to zero)
//------------------------------------------------------------------------------------
TOP5TeVAnalyzer::TOP5TeVAnalyzer() : PAFChainItemSelector() {
	fHDummy = 0;
	fHFidu = 0;
	fHBR = 0;
	fHTopPtWeight = 0;
	fHnGenEle = 0;
	fHnGenMuo = 0;
	fHGenElePt = 0;
	fHGenMuoPt = 0;
	fHWeightsFidu  = 0;
	fHnGenLeptons = 0;
  fHMuonPtnoIso = 0;
  fHMuonPtIso   = 0;
  fHMuonEtanoIso = 0;
  fHMuonEtaIso   = 0;
  fJetDataPt     = 0;
  fJetDataL2L3Pt = 0;
  hmeteff = 0;
  trkhmeteff = 0;

	for (unsigned int ichan = 0; ichan < gNCHANNELS; ichan++) {
		for (unsigned int isyst = 0; isyst < gNSYST; isyst++) {
			fHyields     [ichan][isyst] = 0;
			fHSSyields   [ichan][isyst] = 0;
		}
		for (unsigned int iweight = 0; iweight < gNWEIGHT; iweight++) {
			fHWeightyield[ichan][iweight] = 0;
		}
		for (unsigned int icut = 0; icut < iNCUTS; icut++) {
			fHLepSys[ichan][icut] = 0;
			fHTrigSys[ichan][icut] = 0;

			fHDY_InvMassVsNPV   [ichan][icut] = 0;
			fHDY_InvMassVsMET   [ichan][icut] = 0;
			fHDY_InvMassVsNjets [ichan][icut] = 0;
			fHDY_InvMassVsNbtags[ichan][icut] = 0;
			fHDY_InvMass        [ichan][icut] = 0;

			//++ Origin Histos
			//  fHSSOrigins[ichan][icut] = 0;
			//  fHOrigins[ichan][icut] = 0;

			//++ Kinematic  
			fHWeights[ichan][icut] = 0;
			fHMET[ichan][icut] = 0;       
			fHtrkMET[ichan][icut] = 0;       
			fHMET_noHF[ichan][icut] = 0;       
			fHMETPhi[ichan][icut] = 0;       
			fHTrueMET[ichan][icut] = 0;       
			fHTrueMETPhi[ichan][icut] = 0;       
			fHTrueMET_noHF[ichan][icut] = 0;       
			fHTrueMETPhi_noHF[ichan][icut] = 0;       
			fHLep0Eta[ichan][icut] = 0;    
			fHLep1Eta[ichan][icut] = 0;    
			fHDelLepPhi[ichan][icut] = 0; 
			fHHT[ichan][icut] = 0;        
			fHHT2[ichan][icut] = 0;        
			fHHT3[ichan][icut] = 0;        
			fHHT4[ichan][icut] = 0;        
			fHHT5[ichan][icut] = 0;        
			fHJet0Eta[ichan][icut] = 0;    
			fHJet1Eta[ichan][icut] = 0;    
			fHBtagJet0Pt[ichan][icut] = 0;

			fHMT[ichan][icut] = 0;       
			
			fHJetCSV[ichan][icut]  = 0;
			fHCSVTag[ichan][icut] = 0; 
			fHDelPhillJet[ichan][icut] = 0;

			fHDRLep[ichan][icut] = 0;
			fHDRLep0Jet[ichan][icut] = 0;
			fHDPhiLep0Jet[ichan][icut] = 0;
			fHDRLep1Jet[ichan][icut] = 0;
			fHDPhiLep1Jet[ichan][icut] = 0;

			fHStopMass[ichan][icut] = 0;
			fHChi0Mass[ichan][icut] = 0;
			fHChi0StopMass[ichan][icut] = 0;
			fHvertices[ichan][icut] = 0;
			fHgoodvertices[ichan][icut] = 0;
			fnGenLep[ichan][icut] = 0;

      fMuonIsoCharged[ichan][icut]   = 0;
      fMuonIsoNeutral[ichan][icut]   = 0; 
      fMuonIsoPhotons[ichan][icut]   = 0; 
      fMuonIsoPU[ichan][icut]        = 0; 
      fMuonIso[ichan][icut]          = 0; 
      fSSMuonIso[ichan][icut]        = 0; 

			for (unsigned int isyst = 0; isyst < gNSYST; isyst++) {
				fHInvMass[ichan][icut][isyst] = 0;   
				fHInvMass2[ichan][icut][isyst] = 0;   
				fHSSInvMass[ichan][icut][isyst] = 0;   
				fHNBtagsNJets[ichan][icut][isyst] = 0; 
				fHSSNBtagsNJets[ichan][icut][isyst] = 0; 

				/// STOP
				//fHAbsDelPhiLep[ichan][icut] = 0;
				fHminDelRJetsLeps[ichan][icut][isyst] = 0;
				fHSSminDelRJetsLeps[ichan][icut][isyst] = 0;
				fHdelPhi2LeadJets[ichan][icut][isyst] = 0;
				fHSSdelPhi2LeadJets[ichan][icut][isyst] = 0;
				fHAbsDelPhiLeps[ichan][icut][isyst] = 0;
				fHSSAbsDelPhiLeps[ichan][icut][isyst] = 0;	

				fHNJets[ichan][icut][isyst]       = 0;
				fHNBtagJets[ichan][icut][isyst]   = 0;
				fHJet0Pt[ichan][icut][isyst]      = 0;
				fHJet1Pt[ichan][icut][isyst]      = 0;
				fHDiLepPt[ichan][icut][isyst]     = 0;
				fHLep0Pt[ichan][icut][isyst]      = 0;
				fHLep1Pt[ichan][icut][isyst]      = 0; 
			}
		}
	}

	//++ Gen Info
	for (unsigned int ichan = 0; ichan < gNCHANNELS - 1; ichan++) {
		fHDeltaRLepJet[ichan] = 0;
	}
}

//-------------------------------------------------------------------
// Initialise
//-------------------------------------------------------------------
void TOP5TeVAnalyzer::Initialise() {
	PAF_INFO("TOP5TeVAnalyzer", "+ Initializing...");
	//PAF_INFO("TOP5TeVAnalyzer", "+ Initializing paramenters...");
	GetParameters();
	//PAF_INFO("TOP5TeVAnalyzer", "+ Sumw2 set for all histograms...");
	TH1::SetDefaultSumw2();
	fHDummy = CreateH1F("fHDummy","",1,0,1);
  hmeteff = CreateH1F("hmeteff", "hmeteff", 17, -0.5, 16.5);
  trkhmeteff = CreateH1F("trkhmeteff", "trkhmeteff", 17, -0.5, 16.5);
	fHFidu = CreateH1F("YieldFidu","",1,0,1);
	fHBR = CreateH1F("fHBR","",1,0,1);;
	fHWeightsFidu  = CreateH1F("hPDFweightsFidu","hPDFweightsFidu", nWeights, -0.5, nWeights - 0.5);
	//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise Yield histograms...");
	InitialiseYieldsHistos();
	//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise Kinematic histograms...");
	InitialiseKinematicHistos();
	if (!gIsData) {
		//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise Gen histograms...");
		InitialiseGenHistos();
	}
	//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise other histograms...");
	fHTopPtWeight  = CreateH1F("H_TopPtWeight" ,"TopPt Weight",100, 0, 2);

	fHnGenEle  = CreateH1F("fHnGenEle" , "nGenPromptElecs"  , 11, -1.5, 9.5);
	fHnGenMuo  = CreateH1F("fHnGenMuo" , "nGenPromptMuons"  , 11, -1.5, 9.5);
	fHGenElePt = CreateH1F("fHGenElePt", "GenPromptElecs Pt", 500, 0, 500);
	fHGenMuoPt = CreateH1F("fHGenMuoPt", "GenPromptMuons Pt", 500, 0, 500);
	fHnGenLeptons  = CreateH1F("fHnGenLeptons" , "nGenLeptons"  , 9, -1.5, 7.5);
  fHMuonPtnoIso = CreateH1F("fHMuonPtnoIso", "MuonPt_noIso", 24, 0, 120);
  fHMuonPtIso   = CreateH1F("fHMuonPtIso", "MuonPt_Iso", 24, 0, 120);
  fHMuonEtanoIso = CreateH1F("fHMuonEtanoIso", "MuonEta_noIso", 40, 0, 2.2);
  fHMuonEtaIso   = CreateH1F("fHMuonEtaIso", "MuonEta_Iso", 40, 0, 2.2);
  fJetDataPt     = CreateH1F("fJetDataPt", "fJetDataPt", 40, 0, 120);
  fJetDataL2L3Pt = CreateH1F("fJetDataL2L3Pt", "fJetDataL2L3Pt", 40, 0, 120);

	if (gSampleName.Contains("Data") ||     
			gSampleName == "DoubleEG"        || 
			gSampleName == "SingleMu"        || 
			gSampleName == "SingleEle"  || 
			gSampleName == "MuonEG"	       ||	    
			gSampleName == "TTJets    "      ||  
			gSampleName == "DYJetsToLL_M50_aMCatNLO"     ||
			gSampleName == "DYJetsToLL_M10to50_aMCatNLO" ||
			gSampleName == "DYJetsToLL_M10to50_aMCatNLO_ext" ||
			gSampleName == "ZJets_MLM"       ||
			gSampleName == "DYJets_MLM"      ||
			gSampleName == "TTbar_Powheg"    ||  
			gSampleName == "TTJets_aMCatNLO"  ||  
			gSampleName == "TTbar_Powheg_Pythia6")  {  
		//	PAF_INFO("TOP5TeVAnalyzer", "+ Initialise Drell-Yan histograms...");
		InitialiseDYHistos();
	}
	PAF_INFO("TOP5TeVAnalyzer", "+ Initialise histograms for systematics studies...");
	InitialiseSystematicHistos();

//----xxxxxxxxxx------>	//	PU Reweight
	//--------------------------------------
	//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise Pile-Up reweighting tool...");
/*	fPUWeight     = new PUWeight(gLumiForPU, Fall2015_25ns_matchData_poisson,"2015_25ns_76"); 
	if (!gIsData) {
		fPUWeightUp   = new PUWeight(18494.9,    Fall2015_25ns_matchData_poisson,"2015_25ns_76"); //  18494.9 
		fPUWeightDown = new PUWeight(20441.7,    Fall2015_25ns_matchData_poisson,"2015_25ns_76"); //  20441.7 
	}*/

	//if (gUseCSVM) fBTagSF   = new BTagSFUtil("CSVM","ABCD");//ReReco
	//else          fBTagSF   = new BTagSFUtil("CSVT","ABCD");//ReReco 

	//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise b-tag scale factors...");
	//-------xxxxxxxxxxxx----------
	/*
	if (gUseCSVM){
		fBTagSFnom = new BTagSFUtil("mujets", "CSVv2", "Medium",  0);
		fBTagSFbUp = new BTagSFUtil("mujets", "CSVv2", "Medium",  1);
		fBTagSFbDo = new BTagSFUtil("mujets", "CSVv2", "Medium", -1);
		fBTagSFlUp = new BTagSFUtil("mujets", "CSVv2", "Medium",  3);
		fBTagSFlDo = new BTagSFUtil("mujets", "CSVv2", "Medium", -3);
	}
	else{
		fBTagSFnom = new BTagSFUtil("mujets", "CSVv2", "Tight",  0);
		fBTagSFbUp = new BTagSFUtil("mujets", "CSVv2", "Tight",  1);
		fBTagSFbDo = new BTagSFUtil("mujets", "CSVv2", "Tight", -1);
		fBTagSFlUp = new BTagSFUtil("mujets", "CSVv2", "Tight",  3);
		fBTagSFlDo = new BTagSFUtil("mujets", "CSVv2", "Tight", -3);
	}*/

	PAF_INFO("TOP5TeVAnalyzer", "+ Initialise lepton scale factors...");
	fLeptonSF = new LeptonSF();

	//PAF_INFO("TOP5TeVAnalyzer", "+ Initialise random 3...");
	//fRand3 = new TRandom3(50);

	// No systematics activaded...
	gSysSource = Norm;
	PAF_INFO("TOP5TeVAnalyzer", "+ Initialisation DONE.");
}

void TOP5TeVAnalyzer::InitialiseGenHistos(){
	fHDeltaRLepJet[Muon] = CreateH1F("H_DeltaRLepJet_"+gChanLabel[Muon],"",1000,0.,5.);
	fHDeltaRLepJet[Elec] = CreateH1F("H_DeltaRLepJet_"+gChanLabel[Elec],"",1000,0.,5.);  
	fHnGenLep0 = CreateH1F("H_nGenLep_ElMu_raw","",8, -0.5, 7.5);  
}

void TOP5TeVAnalyzer::InitialiseDYHistos(){
	for (size_t ch=0; ch<gNCHANNELS; ch++){
		for (size_t cut=0; cut<iNCUTS; cut++){
			TString name = "_"+gChanLabel[ch]+"_"+sCut[cut];
			fHDY_InvMassVsNPV   [ch][cut] = CreateH2F("H_DY_InvMassVsNPV"   +name,"",50, -0.5, 49.5, 200, 0, 200);
			fHDY_InvMassVsMET   [ch][cut] = CreateH2F("H_DY_InvMassVsMET"   +name,"",200, 0  , 200 , 200, 0, 200);
			fHDY_InvMassVsNjets [ch][cut] = CreateH2F("H_DY_InvMassVsNjets" +name,"",10, -0.5,  9.5, 200, 0, 200);
			fHDY_InvMassVsNbtags[ch][cut] = CreateH2F("H_DY_InvMassVsNbtags"+name,"",10, -0.5,  9.5, 200, 0, 200);
			fHDY_InvMass        [ch][cut] = CreateH1F("H_DY_InvMass"        +name,"",                200, 0, 200);
		}
	}
}

void TOP5TeVAnalyzer::InitialiseYieldsHistos(){
	//++ Yields histograms
	if (gDoSF) {
		fHyields[Muon][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5); 
		fHyields[Elec][Norm]   = CreateH1F("H_Yields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHSSyields[Muon][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Muon],"", iNCUTS, -0.5, iNCUTS-0.5); 
		fHSSyields[Elec][Norm] = CreateH1F("H_SSYields_"+gChanLabel[Elec],"", iNCUTS, -0.5, iNCUTS-0.5);
	}
	if (gDoDF) {
		fHyields[ElMu][Norm]   = CreateH1F("H_Yields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
		fHSSyields[ElMu][Norm] = CreateH1F("H_SSYields_"+gChanLabel[ElMu],"", iNCUTS, -0.5, iNCUTS-0.5);
	}

	if (gDoSystStudies){
		for (size_t chan=0; chan<gNCHANNELS; chan++){
			if (!gDoSF && chan==Muon) continue;
			if (!gDoSF && chan==Elec) continue;
			if (!gDoDF && chan==ElMu) continue;
			for (size_t sys=1; sys<gNSYST; sys++){
				fHyields[chan][sys]   = CreateH1F("H_Yields_"+gChanLabel[chan]+"_"+SystName[sys],"",iNCUTS,-0.5,iNCUTS-0.5);
				fHSSyields[chan][sys] = CreateH1F("H_SSYields_"+gChanLabel[chan]+"_"+SystName[sys],"", iNCUTS, -0.5, iNCUTS-0.5);
			}
		}
	}
	for (size_t chan=0; chan<gNCHANNELS; chan++){
		//hEventsWeight[chan] = CreateH1F("hEventsWeight"+gChanLabel[chan],"",2*iNCUTS,-0.5,2*iNCUTS-0.);
		for (int wei = 0; wei < gNWEIGHT; ++wei){
			fHWeightyield[chan][wei] = CreateH1F("H_Yields_Wei_"+gChanLabel[chan]+"_"+WeiName[wei],"",iNCUTS,-0.5,iNCUTS-0.5);
		}
	}
	for (size_t chan=0; chan<gNCHANNELS; chan++){
		if (!gDoSF && chan==Muon) continue;
		if (!gDoSF && chan==Elec) continue;
		if (!gDoDF && chan==ElMu) continue;
		for (size_t cut=0; cut<iNCUTS; cut++){
			fHLepSys [chan][cut] = CreateH1F("H_LepSys_" +gChanLabel[chan]+"_"+sCut[cut],"LepSys" , 400, 0, 0.04);
			fHTrigSys[chan][cut] = CreateH1F("H_TrigSys_"+gChanLabel[chan]+"_"+sCut[cut],"TrigSys", 400, 0, 0.04);
		}
	}
}

void TOP5TeVAnalyzer::InitialiseKinematicHistos(){
	//  PAF_DEBUG("TOP5TeVAnalyzer::InitialiseKinematicHistos()",Form("nWeights = %i", nWeights));
	//++ Kinematic histograms
	for (size_t ch=0; ch<gNCHANNELS; ch++){
		if (!gDoSF && ch==Muon) continue;
		if (!gDoSF && ch==Elec) continue;
		if (!gDoDF && ch==ElMu) continue;

		for (size_t cut=0; cut<iNCUTS; cut++){
			//PAF_DEBUG("TOP5TeVAnalyzer::InitialiseKinematicHistos()",Form("cut = %i", cut));
			fHWeights[ch][cut]  = CreateH1F("hPDFweights_"  +gChanLabel[ch]+"_"+sCut[cut],"hPDFweights", nWeights, -0.5, nWeights - 0.5);
			fHMET[ch][cut]         = CreateH1F("H_MET_"        +gChanLabel[ch]+"_"+sCut[cut],"MET"                       , 150, 0,150);
			fHtrkMET[ch][cut]         = CreateH1F("H_trkMET_"        +gChanLabel[ch]+"_"+sCut[cut],"trkMET"              , 150, 0,150);
			fHMET_noHF[ch][cut]  = CreateH1F("H_MET_noHF_"   +gChanLabel[ch]+"_"+sCut[cut],"MET_noHF"                    , 150, 0,150);
			fHTrueMET[ch][cut]     = CreateH1F("H_TrueMET_"        +gChanLabel[ch]+"_"+sCut[cut],"TrueMET"               , 150, 0,150);
			fHTrueMET_noHF[ch][cut]     = CreateH1F("H_TrueMET_noHF_"        +gChanLabel[ch]+"_"+sCut[cut],"TrueMET_noHF", 150, 0,150);
			fHMETPhi[ch][cut]      = CreateH1F("H_METPhi_"        +gChanLabel[ch]+"_"+sCut[cut],"METPhi"       , 32, 0,3.2);
			fHTrueMETPhi[ch][cut]  = CreateH1F("H_TrueMETPhi_"        +gChanLabel[ch]+"_"+sCut[cut],"TrueMETPhi"       , 32, 0,3.2);
			fHTrueMETPhi_noHF[ch][cut]  = CreateH1F("H_TrueMETPhi_noHF_"        +gChanLabel[ch]+"_"+sCut[cut],"TrueMETPhi_noHF"       , 32, 0,3.2);
			fHLep0Eta[ch][cut]     = CreateH1F("H_Lep0Eta_"    +gChanLabel[ch]+"_"+sCut[cut],"Lep0Eta"   , 50  ,0 ,2.5);
			fHLep1Eta[ch][cut]     = CreateH1F("H_Lep1Eta_"    +gChanLabel[ch]+"_"+sCut[cut],"Lep1Eta"   , 50  ,0 ,2.5);
			fHDelLepPhi[ch][cut]   = CreateH1F("H_DelLepPhi_"  +gChanLabel[ch]+"_"+sCut[cut],"DelLepPhi" , 100,-3.2, 3.2);
			fHHT[ch][cut]          = CreateH1F("H_HT_"         +gChanLabel[ch]+"_"+sCut[cut],"HT"        , 300,0,300);
			fHHT2[ch][cut]         = CreateH1F("H_HT2_"        +gChanLabel[ch]+"_"+sCut[cut],"HT2"       , 5000,0,500);
			fHHT3[ch][cut]         = CreateH1F("H_HT3_"        +gChanLabel[ch]+"_"+sCut[cut],"HT3"       , 3000,0,300);
			fHHT4[ch][cut]         = CreateH1F("H_HT4_"        +gChanLabel[ch]+"_"+sCut[cut],"HT4"       , 1000,0,1000);
			fHHT5[ch][cut]         = CreateH1F("H_HT5_"        +gChanLabel[ch]+"_"+sCut[cut],"HT5"       , 1200,0,1200);
			fHBtagJet0Pt[ch][cut]  = CreateH1F("H_BtagJet0Pt_" +gChanLabel[ch]+"_"+sCut[cut],"BtagJet0Pt", 2700,30,300);
			fHJet0Eta[ch][cut]     = CreateH1F("H_Jet0Eta_"	 +gChanLabel[ch]+"_"+sCut[cut],"Jet0Eta"   , 60,0,3.0);
			fHJet1Eta[ch][cut]     = CreateH1F("H_Jet1Eta_"	 +gChanLabel[ch]+"_"+sCut[cut],"Jet1Eta"   , 60,0,3.0);

			fHMT[ch][cut]            = CreateH1F("H_MT_"           +gChanLabel[ch]+"_"+sCut[cut],"MT"           , 800,0.,800);

			fHInvMass[ch][cut][0]       = CreateH1F("H_InvMass_"    +gChanLabel[ch]+"_"+sCut[cut],"InvMass"   ,  300,  0., 300.);
			fHInvMass2[ch][cut][0]      = CreateH1F("H_InvMass2_"   +gChanLabel[ch]+"_"+sCut[cut],"InvMass2"  ,  400, 70., 110.);
			fHSSInvMass[ch][cut][0]     = CreateH1F("HSS_InvMass_"  +gChanLabel[ch]+"_"+sCut[cut],"InvMass"   ,  300,  0., 300.);

			fHNBtagsNJets[ch][cut][0]   = CreateH1F("H_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]  ,"NBtagsNJets"   ,15 , -0.5, 14.5);
			fHSSNBtagsNJets[ch][cut][0] = CreateH1F("HSS_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut],"SS_NBtagsNJets",15 , -0.5, 14.5);

			fHAbsDelPhiLeps[ch][cut][0]   = CreateH1F("H_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut]  ,"AbsDelPhiLeps"   , 20 ,0., 1.);
			fHSSAbsDelPhiLeps[ch][cut][0] = CreateH1F("HSS_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut],"SS_AbsDelPhiLeps", 20 ,0., 1.);

			fHdelPhi2LeadJets[ch][cut][0]   = CreateH1F("H_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut]  ,"delPhi2LeadJets"   , 28,-0.2, 1.2);
			fHSSdelPhi2LeadJets[ch][cut][0] = CreateH1F("HSS_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut],"SS_delPhi2LeadJets", 28,-0.2, 1.2);

			fHminDelRJetsLeps[ch][cut][0]   = CreateH1F("H_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut]  ,"minDelRJetsLeps"   ,   500,0.0, 5.0);
			fHSSminDelRJetsLeps[ch][cut][0] = CreateH1F("HSS_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut],"SS_minDelRJetsLeps",   500,0.0, 5.0);

			fHNJets[ch][cut][0]       = CreateH1F("H_NJets_"      +gChanLabel[ch]+"_"+sCut[cut],"NJets"     , 6 ,-0.5, 5.5);
			fHNBtagJets[ch][cut][0]   = CreateH1F("H_NBtagJets_"  +gChanLabel[ch]+"_"+sCut[cut],"NBtagJets" , 4 ,-0.5, 3.5);
			fHJet0Pt[ch][cut][0]      = CreateH1F("H_Jet0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet0Pt"    , 200,0,200);
			fHJet1Pt[ch][cut][0]      = CreateH1F("H_Jet1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet1Pt"    , 150,0,150);
			fHDiLepPt[ch][cut][0]     = CreateH1F("H_DiLepPt_"    +gChanLabel[ch]+"_"+sCut[cut],"DiLepPt"   , 120,0,120); 
			fHLep0Pt[ch][cut][0]      = CreateH1F("H_Lep0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep0Pt"    , 120,0,120);
			fHLep1Pt[ch][cut][0]      = CreateH1F("H_Lep1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep1Pt"    , 120,0,120);
			
			fHJetCSV[ch][cut]  = CreateH1F("H_JetCSV_" +gChanLabel[ch]+"_"+sCut[cut],"CSV" , 100,0, 1.3);
			// other variables 
			fHCSVTag[ch][cut]      = CreateH1F("H_CSVTag_"     +gChanLabel[ch]+"_"+sCut[cut], "NBtagsNJets"     , 1000, 0.0, 1.0);
			fHDelPhillJet[ch][cut] = CreateH1F("H_DelPhillJet_"+gChanLabel[ch]+"_"+sCut[cut], "DeltaPhi"        , 1000,0.0, TMath::Pi());

			// Different Top / Z topologies
			fHDRLep[ch][cut]       = CreateH1F("H_DRLep_"       +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep",       1000,0.0, 5.0);
			fHDRLep0Jet[ch][cut]   = CreateH1F("H_DRLep0Jet_"   +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep0Jet",   1000,0.0, 5.0);
			fHDPhiLep0Jet[ch][cut] = CreateH1F("H_DPhiLep0Jet_" +gChanLabel[ch]+"_"+sCut[cut], "DeltaPhiLep0Jet", 1000,0.0, TMath::Pi());
			fHDRLep1Jet[ch][cut]   = CreateH1F("H_DRLep1Jet_"   +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep1Jet",   1000,0.0, 5.0);
			fHDPhiLep1Jet[ch][cut] = CreateH1F("H_DPhiLep1Jet_" +gChanLabel[ch]+"_"+sCut[cut], "DeltaPhiLep1Jet", 1000,0.0, TMath::Pi());

			fHvertices[ch][cut]     = CreateH1F("H_Vtx_"+gChanLabel[ch]+"_"+sCut[cut],"", 51, -0.5, 50.5); 
			fHgoodvertices[ch][cut] = CreateH1F("H_goodVtx_"+gChanLabel[ch]+"_"+sCut[cut],"", 51, -0.5, 50.5); 
      fnGenLep[ch][cut] = CreateH1F("H_nGenLep_" +gChanLabel[ch]+"_"+sCut[cut],"NGenLeps", 8, -0.5, 7.5);

      fMuonIsoCharged[ch][cut]   = CreateH1F("H_MuonIsoCharged_" + gChanLabel[ch]+"_"+sCut[cut],"", 150, 0, 150);
      fMuonIsoNeutral[ch][cut]   = CreateH1F("H_MuonIsoNeutral_" + gChanLabel[ch]+"_"+sCut[cut],"",  50, 0,  50);
      fMuonIsoPhotons[ch][cut]   = CreateH1F("H_MuonIsoPhotons_" + gChanLabel[ch]+"_"+sCut[cut],"", 100, 0, 100);
      fMuonIsoPU[ch][cut]        = CreateH1F("H_MuonIsoPU_"      + gChanLabel[ch]+"_"+sCut[cut],"",  20, 0,  20);
      fMuonIso[ch][cut]          = CreateH1F("H_MuonIso_"        + gChanLabel[ch]+"_"+sCut[cut],"",  30, 0,  10);
      fSSMuonIso[ch][cut]        = CreateH1F("HSS_MuonIso_"      + gChanLabel[ch]+"_"+sCut[cut],"",  30, 0,  10);

		}
	}
}

void TOP5TeVAnalyzer::InitialiseSystematicHistos(){
	TString histoname = "";
	for (size_t ch=0; ch<gNCHANNELS; ch++){
		if (!gDoSF && ch==Muon) continue;
		if (!gDoSF && ch==Elec) continue;
		if (!gDoDF && ch==ElMu) continue;

		for (size_t cut=0; cut<iNCUTS; cut++){
			for (size_t sys=1; sys<gNSYST; sys++){
				histoname = "H_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHNBtagsNJets[ch][cut][sys] = CreateH1F(histoname,"NBtagsNJets", 15 , -0.5, 14.5);

				histoname = "H_InvMass_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHInvMass[ch][cut][sys]     = CreateH1F(histoname, "InvMass"   , 300, 0., 300.);

				histoname = "H_InvMass2_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHInvMass2[ch][cut][sys]     = CreateH1F(histoname, "InvMass2"   , 400, 70., 110.);

				histoname = "H_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHAbsDelPhiLeps[ch][cut][sys] = CreateH1F(histoname,"AbsDelPhiLeps", 28,-0.2, 1.2);

				histoname = "H_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHdelPhi2LeadJets[ch][cut][sys] = CreateH1F(histoname,"delPhi2LeadJets", 28,-0.2, 1.2);

				histoname = "H_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHminDelRJetsLeps[ch][cut][sys] = CreateH1F(histoname,"minDelRJetsLeps", 500, 0., 5.0);

				histoname = "HSS_NBtagsNJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHSSNBtagsNJets[ch][cut][sys] = CreateH1F(histoname,"SS_NBtagsNJets", 15 , -0.5, 14.5);

				histoname = "HSS_InvMass_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHSSInvMass[ch][cut][sys]     = CreateH1F(histoname,"SS_InvMass"   , 300, 0., 300.);

				histoname = "HSS_AbsDelPhiLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHSSAbsDelPhiLeps[ch][cut][sys] = CreateH1F(histoname,"SS_AbsDelPhiLeps", 28,-0.2, 1.2);

				histoname = "HSS_delPhi2LeadJets_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHSSdelPhi2LeadJets[ch][cut][sys] = CreateH1F(histoname,"SS_delPhi2LeadJets", 28,-0.2, 1.2);

				histoname = "HSS_minDelRJetsLeps_"+gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys];
				fHSSminDelRJetsLeps[ch][cut][sys] = CreateH1F(histoname,"SS_minDelRJetsLeps", 500, 0., 5.0);

				fHNJets[ch][cut][sys]       = CreateH1F("H_NJets_"      +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"NJets"     , 6 ,-0.5, 5.5);
				fHNBtagJets[ch][cut][sys]   = CreateH1F("H_NBtagJets_"  +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"NBtagJets" , 4 ,-0.5, 3.5);
				fHJet0Pt[ch][cut][sys]      = CreateH1F("H_Jet0Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Jet0Pt"    , 200,0,200);
				fHJet1Pt[ch][cut][sys]      = CreateH1F("H_Jet1Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Jet1Pt"    , 150,0,150);
				fHDiLepPt[ch][cut][sys]     = CreateH1F("H_DiLepPt_"    +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"DiLepPt"   , 120,0,120); 
				fHLep0Pt[ch][cut][sys]      = CreateH1F("H_Lep0Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Lep0Pt"    , 120,0,120);
				fHLep1Pt[ch][cut][sys]      = CreateH1F("H_Lep1Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Lep1Pt"    , 120,0,120);
			}
		}
	}
}

//---------------------------------------------------------------------------------------------------
// Set objets, to be called once per event, saving information in tmp vectors for systematic studies.
//---------------------------------------------------------------------------------------------------
void TOP5TeVAnalyzer::SetOriginalObjects(){
	ResetHypLeptons();
	gSysSource = Norm;

	// SAVING ORIGINAL VALUES FOR MET, JET, LEPTONS for SYST
	JetEt.clear();
	JetPt.clear();
	JetPhi.clear();
	MuPx.clear();
	MuPy.clear();
	MuPz.clear();
	MuEnergy.clear();
	ElPx.clear();
	ElPy.clear();
	ElPz.clear();
	ElEnergy.clear();
	MET  = 0.;
	MET_noHF  = 0.;
	MET_Phi = 0.;
	TrueMET  = 0.;
	TrueMET_Phi = 0.;
	TrueMET_noHF  = 0.;
	TrueMET_Phi_noHF = 0.;
	int k = 0;
	// Save original values for MET, Jets and Leptons
	TLorentzVector j;
	for (Int_t i=0; i < nJet; i++){    
		j.SetPxPyPzE(Jet_px[i], Jet_py[i], Jet_pz[i], Jet_energy[i]);
		JetEt.push_back(j.Et());
		JetPt.push_back(j.Pt());
		JetPhi.push_back(j.Phi());
    if(!gIsData) JetPt.at(i) *= getJetSF(i);
    else{
			fJetDataPt    -> Fill(JetPt.at(i));
			JetPt.at(i) = L2L3->get_corrected_pt(JetPt.at(i), Jet_eta[i]);
			fJetDataL2L3Pt -> Fill(JetPt.at(i));
    }
	}
	for (Int_t i=0; i < nElec; i++){     
		ElPx.push_back(ElecPx[i]); 
		ElPy.push_back(ElecPy[i]);
		ElPz.push_back(ElecPz[i]);
		ElEnergy.push_back(ElecEnergy[i]);
	}
	for (Int_t i=0; i<nMuon; i++){
		MuPx.push_back(MuonPx[i]);
		MuPy.push_back(MuonPy[i]); 
		MuPz.push_back(MuonPz[i]); 
		MuEnergy.push_back(MuonEnergy[i]); 
	}
	MET     = met_pt; //met
	MET_noHF  = 0;//Get<Float_t>("pfTypeIMETnoHF"); //met
  MET_Phi = met_phi; //met
	TrueMET     = TrueMet_pt; //met
	TrueMET_Phi = TrueMet_phi; //met
	TrueMET_noHF     = TrueMet_pt_noHF; //met
	TrueMET_Phi_noHF = TrueMet_phi_noHF; //met
}

void TOP5TeVAnalyzer::SetEventObjects(){
	ResetHypLeptons();

	fChargeSwitch = false;
	EventWeight = 1.;

	// USEFUL COUNTERS
	nGenLepton = 0;
  isEMu      = 0;
	nGenElec   = 0;
	nGenMuon   = 0;
	nGenTau    = 0;
	nTauElec   = 0;
	nTauMuon   = 0;
	nGoodVertex = 0;
	nVertex     = 0;
	nBtags      = 0;
	nJets       = 0;
	nMuons      = 0;
	nElecs      = 0;
	nLeptons    = 0;

	//// READ AND SAVE OBJETS...
	Jet.clear();
	Lepton.clear();

	nLeptons = getSelectedLeptons();
	nJets    = getSelectedJets();
	nBtags   = getNBTags();
}

void TOP5TeVAnalyzer::ResetOriginalObjects(){
	// Save original values for MET, Jets and Leptons
	TLorentzVector j;
	for (Int_t i=0; i < nJet; i++){    
		j.SetPxPyPzE(Jet_px[i],	Jet_py[i], Jet_pz[i],	Jet_energy[i]);
		JetEt[i]  = j.Et();
		JetPt[i]  = j.Pt();
		JetPhi[i] = j.Phi();
    JetPt.at(i)*=getJetSF(i);
	}
	for (Int_t i=0; i<nElec; i++){
		ElPx[i] = ElecPx[i]; 
		ElPy[i] = ElecPy[i];
		ElPz[i] = ElecPz[i];
		ElEnergy[i] = ElecEnergy[i];
	}
	for (Int_t i=0; i<nMuon; i++){ 
		MuPx[i] = MuonPx[i]; 
		MuPy[i] = MuonPy[i]; 
		MuPz[i] = MuonPz[i];
		MuEnergy[i] = MuonEnergy[i];
	}
	setMET(met_pt); //met
}

void TOP5TeVAnalyzer::ResetHypLeptons(){
  TLorentzVector vec(0., 0., 0., 0.);
  fHypLepton1 = lepton(vec, 0, -1, -1);
  fHypLepton2 = lepton(vec, 0, -1, -1);
}


void TOP5TeVAnalyzer::CoutEvent(long unsigned int en, TString t){
  //if(en == 59521 || en == 436079){
  //if(en == 377858 || en == 154320 || en == 39567 || en == 39691){
  if(0){
    cout << t << endl;
  }
  else return;
}

//-----------------------------------------------------------------------
// InsideLoop
//-----------------------------------------------------------------------
void TOP5TeVAnalyzer::InsideLoop() {
	fHDummy->Fill(0.5);

	CoutEvent(evt, Form("Event number = %i", evt));

	// Calculate PU Weight
	PUSF = 1.;
	if (!gIsData)
		PUSF = 1.; //fPUWeight->GetWeight(Get<Float_t>("nTrueInt")); //True       //nTruePU

	// Init data members ........................................................
	GetTreeVariables();
	SetOriginalObjects();
	SetEventObjects();

	// Get number of generated leptons ........................................................
	if (!gIsData) {
		SelectedGenLepton();

// ------------xxxxxxxxxx-------------------
/*			TLorentzVector top;
			for (Int_t t=0; t<Get<Int_t>("nGenTop"); t++){ 
				top.SetPtEtaPhiM(Get<Float_t>("GenTop_pt",   t), 
						Get<Float_t>("GenTop_eta",  t), 
						Get<Float_t>("GenTop_phi",  t),
						Get<Float_t>("GenTop_mass", t));   
				Float_t pt    = TMath::Min(top.Pt(), 400.);
				Float_t topSF = TMath::Exp(0.156 - 0.00137 * pt);
				Weight *= topSF;
			}
			fHTopPtWeight->Fill(TMath::Sqrt(Weight));*/

    bool isEmuDilepton = 1;
    if(gSelection == 0                   ) isEmuDilepton = (nGenElec+nGenMuon) >= 2;    // for dileptonic selection
    if(gSelection == 1                   ) isEmuDilepton = (nGenElec+nGenMuon) <  2;    // for semileptonic selection
    if(gSelection == 2 || gSelection == 4) isEmuDilepton = nGenElec>=1 && nGenMuon>=1;  // for emu selection
    if(gSelection == 5                   ) isEmuDilepton = nGenElec>=1 && nGenMuon>=1;  // for emu selection
    if(gSelection == 3                   ) isEmuDilepton = nGenElec>=2 || nGenMuon>=2;  // for SF selection
    if(gSelection == 6                   ) isEmuDilepton = nGenMuon >= 2;
    if(gSelection == 7                   ) isEmuDilepton = nGenMuon >= 2;
		if(!isEmuDilepton){
			CoutEvent(evt, " >>> IS NOT DILEPTON EVENT <<<");
			if (gSampleName.Contains("TTbar")  ) return;
			if (gSampleName.Contains("TTJets") ) return;
		}
		// Fill Gen Info 
		//----------------------------------------------------------------------------
		TLorentzVector lep,jet;
// -----------------xxxxxxxxxxxxxxxxx---------------
/*
		Float_t minDRmu(999.),minDRel(999.);
		for (Int_t jt = 0; jt < Get<Int_t>("nGenPart"); jt++){ 
			if (abs(Get<Int_t>("GenPart_pdgId",jt)) != 5) continue; // b-quarks!!
			jet.SetPtEtaPhiM(Get<Float_t>("GenPart_pt",   jt), 
					Get<Float_t>("GenPart_eta",  jt), 
					Get<Float_t>("GenPart_phi",  jt),
					Get<Float_t>("GenPart_mass", jt));

			for (Int_t mu=0; mu<ngenLep; mu++){
				if ((TMath::Abs(LepGood_pdgId[mu])) != 13) continue;
				lep.SetPtEtaPhiM(genLep_pt[mu], genLep_eta[mu], genLep_phi[mu], genLep_mass[mu]);

				if (minDRmu > lep.DeltaR(jet))  minDRmu = lep.DeltaR(jet);
			}
			for (Int_t el=0; el < ngenLep; el++){ 
				if ((TMath::Abs(LepGood_pdgId[el])) != 11) continue;
				lep.SetPtEtaPhiM(genLep_pt[el],	genLep_eta[el],	genLep_phi[el],	genLep_mass[el]);
				if (minDRel > lep.DeltaR(jet))  minDRel = lep.DeltaR(jet);
			}
		}
		fHDeltaRLepJet[Muon] -> Fill(minDRmu);
		fHDeltaRLepJet[Elec] -> Fill(minDRel);
*/

		if(gSelection == 2) fHBR->Fill(0.5);

		if((gSelection >= 4) && !gIsData){
/*			TLorentzVector l1; l1.SetPtEtaPhiE(genLep_pt[0], genLep_eta[0], genLep_phi[0], genLep_energy[0]); int l1Id = genLep_pdgId[0];
			for(int k = 1; k< ngenLep; k++){ // Search for the emu pair!
				if(genLep_pdgId[0] == genLep_pdgId[k]) continue;
				else{
					TLorentzVector l2; l2.SetPtEtaPhiE(genLep_pt[k], genLep_eta[k], genLep_phi[k], genLep_energy[k]);
					if      ( (l1+l2).M() < 20) return;
          break;
				}
			}*/
			TLorentzVector l1; TLorentzVector l2; int ngenlepmass = 0;
			for(int k = 0; k < ngenLep; k++){
				if((abs(genLep_MomPID[k]) == 24 && abs(genLep_GMomPID[k]) == 6) || (abs(genLep_MomPID[k]) != 15 && abs(genLep_GMomPID[k]) == 24)){
          //if(gSelection == 7 && abs(genLep_pdgId[k]) != 13) continue; 
					if (ngenlepmass == 0){ l1.SetPtEtaPhiE(genLep_pt[k], genLep_eta[k], genLep_phi[k], genLep_energy[k]); ngenlepmass++;}
					else if (ngenlepmass == 1){ l2.SetPtEtaPhiE(genLep_pt[k], genLep_eta[k], genLep_phi[k], genLep_energy[k]); ngenlepmass++;}
					else break;
					}
				}
			if      ( (l1+l2).M() < 20) return;
      if((l1+l2).M() < 76 || (l1+l2).M() > 106) passGenZVetoCut = true;
      else passGenZVetoCut = false;


		}
		fHFidu->Fill(0.5);
		for(int i = 0; i < nWeights; i++){
			fHWeightsFidu->Fill(i, EventWeight*Get<Float_t>("hWeight", i));
		}
		fHnGenLeptons->Fill(nGenLepton);
	}
	// Fill Yields ...............................................................
#ifdef DEBUG
	cout << "N Leptons: " << Lepton.size() << endl;
	cout << "PassTriggerEMu/EE/MuMu= " 
		<< PassTriggerEMu() << "/"
		<< PassTriggerEE() << "/"
		<< PassTriggerMuMu() << endl;
	cout << "Is ElMu/ElEl/MuMu Event= " 
		<< IsElMuEvent() << "/" 
		<< IsElElEvent() << "/" 
		<< IsMuMuEvent() << endl;
#endif
	FillYields();

	// Get SS Yields...
	fChargeSwitch = true;
	FillYields(); /// Get SS yields....
	fChargeSwitch = false;

	// Fill DY DD histograms
	if (gSampleName.Contains("Data") || gSampleName.Contains("data") || gSampleName.Contains("DY")) FillDYHistograms();

	// Fill Yields for syst. studies (only for MC) ..............................
	if (gIsData)         return;
	if (!gDoSystStudies) return;

	// B-tagging systematics .................................................................
	ResetOriginalObjects();
	gSysSource = BtagUp;
	SetEventObjects();
	FillYields(BtagUp);
	fChargeSwitch = true;  
	FillYields(BtagUp); /* Get SS yields....*/  
	fChargeSwitch = false;

	ResetOriginalObjects();
	gSysSource = BtagDown;
	SetEventObjects();  
	FillYields(BtagDown);
	fChargeSwitch = true;
	FillYields(BtagDown); /// Get SS yields....
	fChargeSwitch = false;

	ResetOriginalObjects();
	gSysSource = MisTagUp;
	SetEventObjects();
	FillYields(MisTagUp);
	fChargeSwitch = true;
	FillYields(MisTagUp); /// Get SS yields....
	fChargeSwitch = false;

	ResetOriginalObjects();
  gSysSource = MisTagDown;
  SetEventObjects();
  FillYields(MisTagDown);
  fChargeSwitch = true;
  FillYields(MisTagDown); /// Get SS yields....
  fChargeSwitch = false;
  
	// JES/JER sytematics ....................................................................
  ResetOriginalObjects();
  SmearJetPts(1);
  gSysSource = JESUp;
  SetEventObjects();
  FillYields(JESUp);
  fChargeSwitch = true;
  FillYields(JESUp); /// Get SS yields....
  fChargeSwitch = false;
  
  ResetOriginalObjects();
  SmearJetPts(2);
  gSysSource = JESDown;
  SetEventObjects();
  FillYields(JESDown);
  fChargeSwitch = true;
  FillYields(JESDown); /// Get SS yields....
  fChargeSwitch = false;
 
  ResetOriginalObjects();
  SmearJetPts(3);
  gSysSource = JER;
  SetEventObjects();
  FillYields(JER);
  fChargeSwitch = true;
  FillYields(JER); /// Get SS yields....
  fChargeSwitch = false;

  // Lepton Scale  sytematics ....................................................................
  ResetOriginalObjects();
  ScaleLeptons(1); //up
  gSysSource = LESUp;
  SetEventObjects();
  FillYields(LESUp);
  fChargeSwitch = true;
  FillYields(LESUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  ScaleLeptons(2); //down
  gSysSource = LESDown;
  SetEventObjects();
  FillYields(LESDown);
  fChargeSwitch = true;
  FillYields(LESDown); /// Get SS yields....
  fChargeSwitch = false;

  // Electon  sytematics ....................................................................
  ResetOriginalObjects();
  gSysSource = ElecUp;
  SetEventObjects();
  FillYields(ElecUp);
  fChargeSwitch = true;
  FillYields(ElecUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = ElecDown;
  SetEventObjects();
  FillYields(ElecDown);
  fChargeSwitch = true;
  FillYields(ElecDown); /// Get SS yields....
  fChargeSwitch = false;

  // Trigger  sytematics ....................................................................
  ResetOriginalObjects();
  gSysSource = TrigUp;
  SetEventObjects();
  FillYields(TrigUp);
  fChargeSwitch = true;
  FillYields(TrigUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = TrigDown;
  SetEventObjects();
  FillYields(TrigDown);
  fChargeSwitch = true;
  FillYields(TrigDown); /// Get SS yields....
  fChargeSwitch = false;

  // Muon   sytematics ....................................................................
  ResetOriginalObjects();
  gSysSource = MuonUp;
  SetEventObjects();
  FillYields(MuonUp);
  fChargeSwitch = true;
  FillYields(MuonUp); /// Get SS yields....
  fChargeSwitch = false;

  ResetOriginalObjects();
  gSysSource = MuonDown;
  SetEventObjects();
  FillYields(MuonDown);
  fChargeSwitch = true;
  FillYields(MuonDown); /// Get SS yields....
  fChargeSwitch = false;

  // Pile Up sytematics ....................................................................
  ResetOriginalObjects();
  if (!gIsData)
    PUSF = 1.; //fPUWeightUp->GetWeight(Get<Float_t>("nTrueInt")); //nTruePU
  gSysSource = PUUp;
  SetEventObjects();
  FillYields(PUUp);
  
  ResetOriginalObjects();
  if (!gIsData)
    PUSF = 1.; //fPUWeightDown->GetWeight(Get<Float_t>("nTrueInt")); //nTruePU
  gSysSource = PUDown;
  SetEventObjects();
  FillYields(PUDown); 
  
  // Top PT ...............................................................................
  ResetOriginalObjects();
  gSysSource = TopPt;
  SetEventObjects();
  FillYields(TopPt);
  fChargeSwitch = true;
  FillYields(TopPt); /// Get SS yields....
  fChargeSwitch = false;
}

void TOP5TeVAnalyzer::Summary(){}

//------------------------------------------------------------------------------
// TRIGGER INFORMATION
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
bool TOP5TeVAnalyzer::PassTriggerMuMu() {
  Bool_t pass = false;
 // if      (gSampleName == "Data_SingleMu") pass = HLT_HIL2Mu15_v1;
 // else if (gSampleName == "Data_SingleElec") pass = (HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1 && !HLT_HIL2Mu15_v1);
 // else    pass = (HLT_HIL2Mu15_v1);
  pass = HLT_HIL2Mu15_v1;
	return pass;
}

bool TOP5TeVAnalyzer::PassTriggerEE(){ 
  Bool_t pass = false;
  //if      (gSampleName == "Data_SingleMu") pass = HLT_HIL2Mu15_v1;
  //else if (gSampleName == "Data_SingleElec") pass = (HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1 && !HLT_HIL2Mu15_v1);
  //else    pass = (HLT_HIL2Mu15_v1);
  pass =  HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1;
	return pass;
}

bool TOP5TeVAnalyzer::PassTriggerEMu(){ 
  Bool_t pass = false;
  if      (gSampleName == "Data_SingleMu") pass = HLT_HIL2Mu15_v1;
  else if (gSampleName == "Data_SingleElec") pass = (HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1 && !HLT_HIL2Mu15_v1);
  else    pass = (HLT_HIL2Mu15_v1);// || HLT_HIDoublePhoton15_Eta2p5_Mass50_1000_R9SigmaHECut_v1;
	return pass;
}

//------------------------------------------------------------------------------
// Get METHODS
//------------------------------------------------------------------------------
float TOP5TeVAnalyzer::getHT(){
	float ht(0);
	for (unsigned int i=0; i<Jet.size(); i++) ht+=Jet[i].p.Pt();
	return ht;
}
float TOP5TeVAnalyzer::getJetPtIndex(unsigned int ind){
	if (Jet.size() <= ind) return -999.;
	return Jet[ind].p.Pt();
}
float TOP5TeVAnalyzer::getJetEtaIndex(unsigned int ind){
	if (Jet.size() <= ind) return -999.;
	return TMath::Abs(Jet[ind].p.Eta());
}
float TOP5TeVAnalyzer::getBtagJetPtIndex(unsigned int ind){
	if (Jet.size() <= ind) return -999.;
	Int_t btagInd = 0;
	if (ind==0) btagInd = getLeadingJetbTag();
	else  return -999.;
	return Jet[btagInd].p.Pt();
}

float TOP5TeVAnalyzer::getMT(gChannel chan){
	float ptl1 = fHypLepton1.p.Pt();
	float ptl2 = fHypLepton2.p.Pt();
	float dphi = getDelPhill();
	return TMath::Sqrt(2*ptl1*ptl2*(1-TMath::Cos(dphi)));
}

float TOP5TeVAnalyzer::getDelPhill(){ return fHypLepton1.p.DeltaPhi(fHypLepton2.p);}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TOP5TeVAnalyzer::getJERScale(int jet){
	float eta = Jet_eta[ jet];
	// 8 TeV
	if(     TMath::Abs(eta) < 0.5) return 1.079;
	else if(TMath::Abs(eta) < 1.1) return 1.099;
	else if(TMath::Abs(eta) < 1.7) return 1.121;
	else if(TMath::Abs(eta) < 2.3) return 1.208;
	else if(TMath::Abs(eta) < 2.8) return 1.254;
	else if(TMath::Abs(eta) < 3.2) return 1.395;
	else                           return 1.056;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TOP5TeVAnalyzer::getJERScaleUp(int jet){
	float eta = Jet_eta[ jet];
	// up, 8 TeV
	if(	  TMath::Abs(eta) < 0.5) return 1.053;
	else if(TMath::Abs(eta) < 1.1) return 1.071;
	else if(TMath::Abs(eta) < 1.7) return 1.092;
	else if(TMath::Abs(eta) < 2.3) return 1.162;
	else if(TMath::Abs(eta) < 2.8) return 1.192;
	else if(TMath::Abs(eta) < 3.2) return 1.332;
	else  			 return 0.865;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TOP5TeVAnalyzer::getJERScaleDown(int jet){
	float eta = Jet_eta[ jet];
	// down, 8 TeV
	if(	  TMath::Abs(eta) < 0.5) return 1.105;
	else if(TMath::Abs(eta) < 1.1) return 1.127;
	else if(TMath::Abs(eta) < 1.7) return 1.150;
	else if(TMath::Abs(eta) < 2.3) return 1.254;
	else if(TMath::Abs(eta) < 2.8) return 1.316;
	else if(TMath::Abs(eta) < 3.2) return 1.458;
	else  			 return 1.247;
}

float TOP5TeVAnalyzer::getErrPt(float Pt, float Eta) {
	float InvPerr2;
	float N(0.), S(0.), C(0.), m(0.);

	if(TMath::Abs(Eta) < 0.5 ) {
		N = 3.96859;
		S = 0.18348;
		C = 0.;
		m = 0.62627;
	} else if( TMath::Abs(Eta) < 1.  ) {
		N = 3.55226;
		S = 0.24026;
		C = 0.;
		m = 0.52571;
	} else if( TMath::Abs(Eta) < 1.5  ) {
		N = 4.54826;
		S = 0.22652;
		C = 0.;
		m = 0.58963;
	} else if( TMath::Abs(Eta) < 2.  ) {
		N = 4.62622;
		S = 0.23664;
		C = 0.;
		m = 0.48738;
	} else if( TMath::Abs(Eta) < 2.5  ) {
		N = 2.53324;
		S = 0.34306;
		C = 0.;
		m = 0.28662;
	} else if( TMath::Abs(Eta) < 3.  ) {
		N = -3.33814;
		S = 0.73360;
		C = 0.;
		m = 0.08264;
	} else if( TMath::Abs(Eta) < 5.  ) {
		N = 2.95397;
		S = 0.11619;
		C = 0.;
		m = 0.96086;
	}
	// this is the absolute resolution (squared), not sigma(pt)/pt
	// so have to multiply by pt^2, thats why m+1 instead of m-1
	InvPerr2 =  (N * TMath::Abs(N) ) + (S * S) * pow(Pt, m+1) + (C * C) * Pt * Pt ;
	return sqrt(InvPerr2);
}

/*
float TOP5TeVAnalyzer::getLeptonError(gChannel chan){
	float err1(0.), err2(0.);
	if (chan==Muon){
		err1 = fLeptonSF->GetTightMuonSF_err(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		err2 = fLeptonSF->GetTightMuonSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
	if (chan==ElMu){
		err1 = fLeptonSF->GetTightMuonSF_err    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		err2 = fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
	if (chan==Elec){
		err1 = fLeptonSF->GetTightElectronSF_err(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		err2 = fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	}
	return TMath::Sqrt(err1*err1+err2*err2);
}

float TOP5TeVAnalyzer::getTriggerError(gChannel chan){
	float trig(0.);
	if (chan==Muon) trig = fLeptonSF->GetDoubleMuSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	if (chan==ElMu) trig = fLeptonSF->GetMuEGSF_err    (fHypLepton2.p.Eta(),fHypLepton1.p.Eta());
	if (chan==Elec) trig = fLeptonSF->GetDoubleElSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	return trig;
}*/

float TOP5TeVAnalyzer::getSF(gChannel chan, gSystFlag sys ) {
	if (gIsData)              return 1.; //Don't scale data
	float id1(1.),id2(1.), trig(1.);
	float err1(0.), err2(0.), err_trg(0.);
  float muoSF(1.), eleSF(1.03), trigSF(1.0);//trigSF(0.985);
  float muoSFerr(0.03), eleSFerr(0.02), trigSFerr(0.0);//trigSFerr(0.015);
  if(sys == TrigUp  ) trigSF += trigSFerr;
  if(sys == TrigDown) trigSF -= trigSFerr;
  if(sys == MuonUp  ) muoSF += muoSFerr;
  if(sys == MuonDown) muoSF -= muoSFerr;
  //if(sys == ElecUp  ) eleSF  += eleSFerr;
  //if(sys == ElecDown) eleSF  -= eleSFerr;
	if (chan == Muon){
		id1  = muoSF;
		id2  = muoSF;
		trig = trigSF;
	} 
	else if (chan == Elec){
		id1  = eleSF;
		id2  = eleSF;
		trig = trigSF;
	}
	else if (chan == ElMu){
		id1  = muoSF;
		//id2  = eleSF;
		id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
		trig = trigSF;
	}
  if(sys == ElecUp  ) id2  += fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
  if(sys == ElecDown) id2  -= fLeptonSF->GetTightElectronSF_err(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
	if     (chan == ElMu) return (id1*id2*trig);
	else if(chan == Elec) return (id2*id2*trig);
	else if(chan == Muon) return (id1*id1*trig);
}

//--------------------------------------------------------------------------
// Fill histograms      
//------------------------------------------------------------------------
void TOP5TeVAnalyzer::FillDYHistograms(){
	float Mll = 0.;
	if (PassTriggerEMu()  && IsElMuEvent()){
		// Define Hypothesis Leptons...
		EventWeight = gWeight * getSF(ElMu);
		if(gIsMCatNLO)  EventWeight = EventWeight * genWeight; 
		Mll = (fHypLepton1.p+fHypLepton2.p).M();
		if (PassesMllVeto()){
			fHDY_InvMassVsNPV   [ElMu][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
 			fHDY_InvMassVsMET   [ElMu][iDilepton]->Fill(getMET()   , Mll, EventWeight);
			fHDY_InvMassVsNjets [ElMu][iDilepton]->Fill(getNJets() , Mll, EventWeight);
			fHDY_InvMassVsNbtags[ElMu][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
			fHDY_InvMass        [ElMu][iDilepton]->Fill(             Mll, EventWeight);
			//if (1)   { // No cut in MET for e-mu channel
			if (PassesMETCut())   {  
				fHDY_InvMassVsNPV   [ElMu][iMET]->Fill(nGoodVertex, Mll, EventWeight);
				fHDY_InvMassVsMET   [ElMu][iMET]->Fill(getMET()   , Mll, EventWeight);
				fHDY_InvMassVsNjets [ElMu][iMET]->Fill(getNJets() , Mll, EventWeight);
				fHDY_InvMassVsNbtags[ElMu][iMET]->Fill(getNBTags(), Mll, EventWeight);
				fHDY_InvMass	    [ElMu][iMET]->Fill( 	    Mll, EventWeight);

				if (PassesNJetsCut()) {
					fHDY_InvMassVsNPV   [ElMu][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [ElMu][i2jets]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [ElMu][i2jets]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[ElMu][i2jets]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [ElMu][i2jets]->Fill(             Mll, EventWeight);

					if (PassesNBtagCut()) {
						fHDY_InvMassVsNPV   [ElMu][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
						fHDY_InvMassVsMET   [ElMu][i1btag]->Fill(getMET()   , Mll, EventWeight);
						fHDY_InvMassVsNjets [ElMu][i1btag]->Fill(getNJets() , Mll, EventWeight);
						fHDY_InvMassVsNbtags[ElMu][i1btag]->Fill(getNBTags(), Mll, EventWeight);
						fHDY_InvMass        [ElMu][i1btag]->Fill(             Mll, EventWeight);
					}
				}
			}
			if (getNBTags() == 1){
				fHDY_InvMassVsNPV   [ElMu][iExact1btag]->Fill(nGoodVertex, Mll, EventWeight);
				fHDY_InvMassVsMET   [ElMu][iExact1btag]->Fill(getMET()   , Mll, EventWeight);
				fHDY_InvMassVsNjets [ElMu][iExact1btag]->Fill(getNJets() , Mll, EventWeight);
				fHDY_InvMassVsNbtags[ElMu][iExact1btag]->Fill(getNBTags(), Mll, EventWeight);
				fHDY_InvMass        [ElMu][iExact1btag]->Fill(             Mll, EventWeight);
			}
			if (getNBTags() == 2){
				fHDY_InvMassVsNPV   [ElMu][iExact2btag]->Fill(nGoodVertex, Mll, EventWeight);
				fHDY_InvMassVsMET   [ElMu][iExact2btag]->Fill(getMET()   , Mll, EventWeight);
				fHDY_InvMassVsNjets [ElMu][iExact2btag]->Fill(getNJets() , Mll, EventWeight);
				fHDY_InvMassVsNbtags[ElMu][iExact2btag]->Fill(getNBTags(), Mll, EventWeight);
				fHDY_InvMass        [ElMu][iExact2btag]->Fill(             Mll, EventWeight);
			}
		}
	}
	ResetHypLeptons(); 
	if (PassTriggerMuMu() && IsMuMuEvent()){
		EventWeight = gWeight * getSF(Muon);
		if(gIsMCatNLO)   EventWeight = EventWeight * genWeight; 
		Mll = (fHypLepton1.p+fHypLepton2.p).M();

		if (PassesMllVeto()){
			fHDY_InvMassVsNPV   [Muon][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
			fHDY_InvMassVsMET   [Muon][iDilepton]->Fill(getMET()   , Mll, EventWeight);
			fHDY_InvMassVsNjets [Muon][iDilepton]->Fill(getNJets() , Mll, EventWeight);
			fHDY_InvMassVsNbtags[Muon][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
			fHDY_InvMass        [Muon][iDilepton]->Fill(             Mll, EventWeight);

			if (PassesMETCut())   {
			//if (1)   {  // No cut in MET for e-mu channel
				fHDY_InvMassVsNPV   [Muon][iMET]->Fill(nGoodVertex, Mll, EventWeight);
				fHDY_InvMassVsMET   [Muon][iMET]->Fill(getMET()   , Mll, EventWeight);
				fHDY_InvMassVsNjets [Muon][iMET]->Fill(getNJets() , Mll, EventWeight);
				fHDY_InvMassVsNbtags[Muon][iMET]->Fill(getNBTags(), Mll, EventWeight);
				fHDY_InvMass	    [Muon][iMET]->Fill( 	    Mll, EventWeight);
				if (getNBTags() == 1){
					fHDY_InvMassVsNPV   [Muon][iExact1btag]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [Muon][iExact1btag]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [Muon][iExact1btag]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[Muon][iExact1btag]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [Muon][iExact1btag]->Fill(             Mll, EventWeight);
				}
				if (getNBTags() == 2){
					fHDY_InvMassVsNPV   [Muon][iExact2btag]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [Muon][iExact2btag]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [Muon][iExact2btag]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[Muon][iExact2btag]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [Muon][iExact2btag]->Fill(             Mll, EventWeight);
				}
				if (PassesNJetsCut()) {
					fHDY_InvMassVsNPV   [Muon][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [Muon][i2jets]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [Muon][i2jets]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[Muon][i2jets]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [Muon][i2jets]->Fill(             Mll, EventWeight);

					if (PassesNBtagCut()) {
						fHDY_InvMassVsNPV   [Muon][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
						fHDY_InvMassVsMET   [Muon][i1btag]->Fill(getMET()   , Mll, EventWeight);
						fHDY_InvMassVsNjets [Muon][i1btag]->Fill(getNJets() , Mll, EventWeight);
						fHDY_InvMassVsNbtags[Muon][i1btag]->Fill(getNBTags(), Mll, EventWeight);
						fHDY_InvMass        [Muon][i1btag]->Fill(             Mll, EventWeight);
					}
				}
			}
		}
	}
	ResetHypLeptons(); 
	if (PassTriggerEE()   && IsElElEvent()){
		EventWeight = gWeight * getSF(Elec);
		if(gIsMCatNLO)    EventWeight = EventWeight * genWeight; 
		Mll = (fHypLepton1.p+fHypLepton2.p).M();

		if (PassesMllVeto()){
			fHDY_InvMassVsNPV   [Elec][iDilepton]->Fill(nGoodVertex, Mll, EventWeight);
			fHDY_InvMassVsMET   [Elec][iDilepton]->Fill(getMET()   , Mll, EventWeight);
			fHDY_InvMassVsNjets [Elec][iDilepton]->Fill(getNJets() , Mll, EventWeight);
			fHDY_InvMassVsNbtags[Elec][iDilepton]->Fill(getNBTags(), Mll, EventWeight);
			fHDY_InvMass        [Elec][iDilepton]->Fill(             Mll, EventWeight);

				if (PassesMETCut())   {
			//if (1)   {  // No vut in MET for e-mu channel
				fHDY_InvMassVsNPV   [Elec][iMET]->Fill(nGoodVertex, Mll, EventWeight);
				fHDY_InvMassVsMET   [Elec][iMET]->Fill(getMET()   , Mll, EventWeight);
				fHDY_InvMassVsNjets [Elec][iMET]->Fill(getNJets() , Mll, EventWeight);
				fHDY_InvMassVsNbtags[Elec][iMET]->Fill(getNBTags(), Mll, EventWeight);
				fHDY_InvMass	    [Elec][iMET]->Fill( 	    Mll, EventWeight);

				if (getNBTags() == 1){
					fHDY_InvMassVsNPV   [Elec][iExact1btag]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [Elec][iExact1btag]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [Elec][iExact1btag]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[Elec][iExact1btag]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [Elec][iExact1btag]->Fill(             Mll, EventWeight);
				}
 				if (getNBTags() == 2){
					fHDY_InvMassVsNPV   [Elec][iExact2btag]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [Elec][iExact2btag]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [Elec][iExact2btag]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[Elec][iExact2btag]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [Elec][iExact2btag]->Fill(             Mll, EventWeight);
				}
				if (PassesNJetsCut()) {
					fHDY_InvMassVsNPV   [Elec][i2jets]->Fill(nGoodVertex, Mll, EventWeight);
					fHDY_InvMassVsMET   [Elec][i2jets]->Fill(getMET()   , Mll, EventWeight);
					fHDY_InvMassVsNjets [Elec][i2jets]->Fill(getNJets() , Mll, EventWeight);
					fHDY_InvMassVsNbtags[Elec][i2jets]->Fill(getNBTags(), Mll, EventWeight);
					fHDY_InvMass        [Elec][i2jets]->Fill(             Mll, EventWeight);
					if (PassesNBtagCut()) {
						fHDY_InvMassVsNPV   [Elec][i1btag]->Fill(nGoodVertex, Mll, EventWeight);
						fHDY_InvMassVsMET   [Elec][i1btag]->Fill(getMET()   , Mll, EventWeight);
						fHDY_InvMassVsNjets [Elec][i1btag]->Fill(getNJets() , Mll, EventWeight);
						fHDY_InvMassVsNbtags[Elec][i1btag]->Fill(getNBTags(), Mll, EventWeight);
						fHDY_InvMass        [Elec][i1btag]->Fill(             Mll, EventWeight);
					}
				}
			}
		}
	}
	ResetHypLeptons();
}
void TOP5TeVAnalyzer::FillKinematicHistos(gChannel chan, iCut cut){
#ifdef DEBUG
	cout << "Filling KinematicHistos("<<chan<<","<<cut<<")... ";
	cout << fHypLepton1.index << " , " << fHypLepton2.index << endl;
#endif

	if (gSysSource != Norm)      return;  //only fill histograms for nominal distributions...
	if (fChargeSwitch == true  ){
		for(int i = 0; i < nMuon; i++){
			fSSMuonIso[chan][cut]     ->Fill(getMuonIso(i), EventWeight);
		}
			//if(chan==0 && cut==0)	cout << evt << "\t" << endl;
		return;
	}

	if (!gIsData) {
		Int_t nWeightsTree = Get<Int_t>("nWeights");
		for(int i = 0; i < nWeightsTree; i++){
			fHWeights[chan][cut]->Fill(i, EventWeight*Get<Float_t>("hWeight", i));
		}
	}
	//++ met info
	fHMET[chan][cut]           ->Fill(getMET(),             EventWeight);
	fHtrkMET[chan][cut]  ->Fill(met_trk, EventWeight);       
	fHMET_noHF[chan][cut]           ->Fill(Get<Float_t>("met_pt_noHF"),             EventWeight);
	fHMETPhi[chan][cut]           ->Fill(getMETPhi(),             EventWeight);
	fHTrueMET[chan][cut]           ->Fill(getTrueMET(),             EventWeight);
	fHTrueMETPhi[chan][cut]           ->Fill(getTrueMETPhi(),             EventWeight);
	fHTrueMET_noHF[chan][cut]           ->Fill(getTrueMET_noHF(),             EventWeight);
	fHTrueMETPhi_noHF[chan][cut]           ->Fill(getTrueMETPhi_noHF(),             EventWeight);
	fHMT[chan][cut]            ->Fill(getMT(chan),          EventWeight);

	//++ lepton info
	// fHDiLepPt[chan][cut]    ->Fill((fHypLepton1.p+fHypLepton2.p).Pt(),      EventWeight);
	// fHLep0Pt[chan][cut]     ->Fill(fHypLepton1.p.Pt(),                      EventWeight);
	// fHLep1Pt[chan][cut]     ->Fill(fHypLepton2.p.Pt(),                      EventWeight);
	fHLep0Eta[chan][cut]    ->Fill(TMath::Abs(fHypLepton1.p.Eta()),         EventWeight);
	fHLep1Eta[chan][cut]    ->Fill(TMath::Abs(fHypLepton2.p.Eta()),         EventWeight);
	fHDelLepPhi[chan][cut]  ->Fill(fHypLepton1.p.DeltaPhi((fHypLepton2.p)), EventWeight);
	//++ jet info		  
	int njets = getNJets(); 
	//fHNJets[chan][cut]      ->Fill(njets,                                   EventWeight);
	if(getHT()>0){
		fHHT[chan][cut]      ->Fill(getHT(),                                 EventWeight);
		fHHT2[chan][cut]     ->Fill(getHT(),                                 EventWeight);
		fHHT3[chan][cut]     ->Fill(getHT(),                                 EventWeight);
		fHHT4[chan][cut]     ->Fill(getHT(),                                 EventWeight);
		fHHT5[chan][cut]     ->Fill(getHT(),                                 EventWeight);
	}
	// fHNBtagJets[chan][cut]  ->Fill(getNBTags(),                             EventWeight);
	// fHJet0Pt[chan][cut]     ->Fill(getJetPtIndex(0),                        EventWeight);
	// fHJet1Pt[chan][cut]     ->Fill(getJetPtIndex(1),                        EventWeight);
	fHJet0Eta[chan][cut]    ->Fill(getJetEtaIndex(0),			  EventWeight);
	fHJet1Eta[chan][cut]    ->Fill(getJetEtaIndex(1),			  EventWeight);
  fHBtagJet0Pt[chan][cut] ->Fill(getBtagJetPtIndex(0),                    EventWeight);

	//if(chan==0 && cut==0){// && (fHypLepton1.p+fHypLepton2.p).M() > 175){ 
	//cout << evt << endl; //"\t" << run << "\t" << lum << endl; //<< "     InvMass = " << (fHypLepton1.p+fHypLepton2.p).M() << endl; //" : " << T_Event_LuminosityBlock <<  endl;	
	//cout << " ## genWeight                = " << genWeight << endl;
	//	cout << " ## norm (xsec/SumOfWeights) = " << gWeight    << endl;
	//}
  

	for(int k = 0; k<nJets; k++){
		fHJetCSV[chan][cut] -> Fill(Jet_btagCSV[k], EventWeight);
	}
	int ib = getLeadingJetbTag();
	if (ib>=0)  fHCSVTag[chan][cut] ->Fill(Jet_btagCSV[ib], EventWeight);

	fHDelPhillJet[chan][cut]->Fill(getDeltaPhillJet(), EventWeight);

	//// Top/Z diff topology.
	fHDRLep[chan][cut]        ->Fill(fHypLepton1.p.DeltaR(fHypLepton2.p),     EventWeight);

	for(int i = 0; i < nMuon; i++){
		fMuonIsoCharged[chan][cut]->Fill( Get<Float_t>("MuonPFChIso", i), EventWeight);
		fMuonIsoNeutral[chan][cut]->Fill( Get<Float_t>("MuonPFNeuIso", i), EventWeight);
		fMuonIsoPhotons[chan][cut]->Fill( Get<Float_t>("MuonPFPhoIso", i), EventWeight);
		fMuonIsoPU[chan][cut]     ->Fill( Get<Float_t>("MuonPFPUIso", i), EventWeight);
		fMuonIso[chan][cut]       ->Fill(getMuonIso(i), EventWeight);
	}

	if (njets > 0) {
		fHDRLep0Jet[chan][cut]    ->Fill(getDRClosestJet(fHypLepton1.p),   EventWeight);
		fHDRLep1Jet[chan][cut]    ->Fill(getDRClosestJet(fHypLepton2.p),   EventWeight);
		fHDPhiLep0Jet[chan][cut]  ->Fill(getDPhiClosestJet(fHypLepton1.p), EventWeight);
		fHDPhiLep1Jet[chan][cut]  ->Fill(getDPhiClosestJet(fHypLepton2.p), EventWeight);
	}
	//-----------xxxxxxxxxxxxxx---------
	/*
		 Int_t nVert =  Get<Int_t>("nPu");
		 fHvertices[chan][cut]     ->Fill(nVert, EventWeight); // for now, using the same
		 fHgoodvertices[chan][cut] ->Fill(nVert, EventWeight); // for now, using the same
		 */
#ifdef DEBUG
	cout << " DONE!" << endl;
#endif

}

void TOP5TeVAnalyzer::FillYieldsHistograms(gChannel chan, iCut cut, gSystFlag sys){
#ifdef DEBUG
	cout << "FillYieldsHistograms("<<chan<<","<<cut<<","<<sys<<")...";
#endif
	if (fChargeSwitch){   fHSSyields[chan][sys]->Fill(cut, EventWeight);  }
	else                  fHyields[chan][sys]  ->Fill(cut, EventWeight);  
	

	/// FOR SYSTEMATIC STUDIES
	int njets  = 0; njets  = getNJets();
	int nbtags = 0; nbtags = getNBTags();

	if (fChargeSwitch) { 
		//cout << "!!!!!!!"<< endl; 
		fHSSInvMass[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
		if (njets == 0) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags,        EventWeight);
		if (njets == 1) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+1,      EventWeight);
		if (njets == 2) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+3,      EventWeight);
		if (njets == 3) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+6,      EventWeight);
		if (njets >= 4) fHSSNBtagsNJets[chan][cut][sys]->Fill(nbtags+10,     EventWeight);
		fHSSAbsDelPhiLeps[chan][cut][sys]->Fill( TMath::Abs(fHypLepton1.p.DeltaPhi(fHypLepton2.p))/TMath::Pi(), EventWeight);
		if(njets >= 2) fHSSdelPhi2LeadJets[chan][cut][sys]->Fill( TMath::Abs(Jet[0].p.DeltaPhi(Jet[1].p))/TMath::Pi(), EventWeight);
		if (njets > 1) {
			double deltaR_temp1 = 999., deltaR_temp2 = 999.;       
			if( fHypLepton1.p.DeltaR(Jet[0].p) < fHypLepton2.p.DeltaR(Jet[0].p) ) {
				deltaR_temp1 = fHypLepton1.p.DeltaR(Jet[0].p);
				deltaR_temp2 = fHypLepton2.p.DeltaR(Jet[1].p);  
			}else{
				deltaR_temp1 = fHypLepton2.p.DeltaR(Jet[0].p);
				deltaR_temp2 = fHypLepton1.p.DeltaR(Jet[1].p);  
			}
			fHSSminDelRJetsLeps[chan][cut][sys]->Fill(TMath::Min(deltaR_temp1, deltaR_temp2), EventWeight);
		}
	}
	else {
		fHInvMass[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
		fHInvMass2[chan][cut][sys]->Fill((fHypLepton1.p+fHypLepton2.p).M(), EventWeight);
		if (njets == 0) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags,        EventWeight);
		if (njets == 1) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+1,      EventWeight);
		if (njets == 2) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+3,      EventWeight);
		if (njets == 3) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+6,      EventWeight);
		if (njets >= 4) fHNBtagsNJets[chan][cut][sys]->Fill(nbtags+10,     EventWeight);
		fHAbsDelPhiLeps[chan][cut][sys]->Fill( TMath::Abs(fHypLepton1.p.DeltaPhi(fHypLepton2.p))/TMath::Pi(), EventWeight);
		if( njets >= 2) fHdelPhi2LeadJets[chan][cut][sys]->Fill( TMath::Abs(Jet[0].p.DeltaPhi(Jet[1].p))/TMath::Pi(), EventWeight);
		if (njets > 1) {
			double deltaR_temp1 = 999., deltaR_temp2 = 999.;       
			if( fHypLepton1.p.DeltaR(Jet[0].p) < fHypLepton2.p.DeltaR(Jet[0].p) ) {
				deltaR_temp1 = fHypLepton1.p.DeltaR(Jet[0].p);
				deltaR_temp2 = fHypLepton2.p.DeltaR(Jet[1].p);  
			}else{
				deltaR_temp1 = fHypLepton2.p.DeltaR(Jet[0].p);
				deltaR_temp2 = fHypLepton1.p.DeltaR(Jet[1].p);  
			}
			fHminDelRJetsLeps[chan][cut][sys]->Fill(TMath::Min(deltaR_temp1, deltaR_temp2), EventWeight);
		}
		fHNJets[chan][cut][sys]       ->Fill(getNJets(), EventWeight);
		fHNBtagJets[chan][cut][sys]   ->Fill(getNBTags(), EventWeight);
		fHJet0Pt[chan][cut][sys]      ->Fill(getJetPtIndex(0), EventWeight);
		fHJet1Pt[chan][cut][sys]      ->Fill(getJetPtIndex(1), EventWeight);
		fHDiLepPt[chan][cut][sys]     ->Fill((fHypLepton1.p+fHypLepton2.p).Pt(), EventWeight);
		fHLep0Pt[chan][cut][sys]      ->Fill(fHypLepton1.p.Pt(), EventWeight);
		fHLep1Pt[chan][cut][sys]      ->Fill(fHypLepton2.p.Pt(), EventWeight); 
	}

	if (!gIsData){
		//fHLepSys[chan][cut] ->Fill(getLeptonError(chan), EventWeight);
		//fHTrigSys[chan][cut]->Fill(getTriggerError(chan),EventWeight);
	}
#ifdef DEBUG
	cout << " DONE! " << endl;
#endif
	return;
}
void TOP5TeVAnalyzer::FillYields(gSystFlag sys){
	ResetHypLeptons();  

#ifdef DEBUG
	cout << "gDoDF= " << gDoDF << endl;
	cout << "PassTriggerEMu= " << PassTriggerEMu() << endl;
	cout << "Is ElMu/ElEl/MuMu Event= " 
		<< IsElMuEvent() << "/" 
		<< IsElElEvent() << "/" 
		<< IsMuMuEvent() << endl;
#endif

	if (gSelection == 5){
		EventWeight = gWeight;
		FillYieldsHistograms(ElMu, iDilepton, sys);
		if(PassesNGenJetsCut()) FillYieldsHistograms(ElMu, i2jets, sys);
		return;
	}
	if (gSelection == 7){
		EventWeight = gWeight;
    //if(nGenMuon <2) return;
		FillYieldsHistograms(Muon, iDilepton, sys);
		if(PassesGenZVetoCut()){
			FillYieldsHistograms(Muon, iZVeto, sys);
			if(PassesGenMetCut()){
				FillYieldsHistograms(Muon, iMET, sys);
				if(PassesNGenJetsCut()) FillYieldsHistograms(Muon, i2jets, sys);
			}
		}
		return;
	}


	if (gDoDF && PassTriggerEMu()  && IsElMuEvent()){
		// Define Hypothesis Leptons...
		EventWeight = gWeight * getSF(ElMu, sys);// * getTopPtSF();
#ifdef DEBUG
		cout << " pass trigger + emu, ";
#endif
		// 0.115 = Fraction events with negative weight
		if(gIsMCatNLO) EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115)); 
		if (PassesMllVeto()){
#ifdef DEBUG
			cout << " pass mll, ";
#endif
			FillYieldsHistograms(ElMu, iDilepton, sys);
			if(sys==Norm) FillKinematicHistos(ElMu,iDilepton);

			FillYieldsHistograms(ElMu, iZVeto, sys);      
			if(sys==Norm) FillKinematicHistos(ElMu, iZVeto);

			FillYieldsHistograms(ElMu, iMET, sys);      
			if(sys==Norm) FillKinematicHistos(ElMu,iMET);

			if (PassesNJetsCut()) {
#ifdef DEBUG
				cout << " pass njets with njets = "<<getNJets()<<", ";
#endif
				FillYieldsHistograms(ElMu, i2jets, sys);      
				if(sys==Norm) FillKinematicHistos(ElMu,i2jets);
				if (PassesNBtagCut()) {
#ifdef DEBUG
					cout << " pass nbjets with nbtags = "<<getNBTags()<<", ";
#endif
					//if (sys == LESUp) cout << evt<< endl;  //LESup    8 //EventNumber
					FillYieldsHistograms(ElMu, i1btag, sys);
					if(sys==Norm) FillKinematicHistos(ElMu,i1btag);
				}
			}
			if (getNBTags() == 1){
#ifdef DEBUG
				cout << " pass nbjets=1";
#endif
				FillYieldsHistograms(ElMu, iExact1btag, sys);      
				if(sys==Norm) FillKinematicHistos(ElMu,iExact1btag);
			}
			if (getNBTags() == 2){
#ifdef DEBUG
				cout << " pass nbjets=2";
#endif
				FillYieldsHistograms(ElMu, iExact2btag, sys);      
				if(sys==Norm) FillKinematicHistos(ElMu,iExact2btag);
			}
		}
	}

	ResetHypLeptons(); 
	if (gDoSF && PassTriggerMuMu() && IsMuMuEvent()){
		EventWeight = gWeight * getSF(Muon, sys); //  * getTopPtSF();
		//EventWeight = 1.;
#ifdef DEBUG
		cout << " pass trigger + mumu, ";
#endif
		// 0.115 = Fraction events with negative weight
		if(gIsMCatNLO) EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115)); 
		if (PassesMllVeto()){
#ifdef DEBUG
			cout << " pass mll, ";
#endif
			FillYieldsHistograms(Muon,iDilepton, sys);
			if(sys==Norm) FillKinematicHistos(Muon,iDilepton);
			if (PassesZVeto())    {
			//if (1){
				FillYieldsHistograms(Muon,iZVeto, sys);      
				if(sys==Norm) FillKinematicHistos(Muon,iZVeto);
          if(sys == Norm){
						hmeteff->Fill(0., EventWeight);
						if (getMET_noHF() > 5.)  hmeteff->Fill(1, EventWeight); 
						if (getMET_noHF() > 10.) hmeteff->Fill(2, EventWeight); 
						if (getMET_noHF() > 15.) hmeteff->Fill(3, EventWeight); 
						if (getMET_noHF() > 20.) hmeteff->Fill(4, EventWeight); 
						if (getMET_noHF() > 25.) hmeteff->Fill(5, EventWeight); 
						if (getMET_noHF() > 30.) hmeteff->Fill(6, EventWeight); 
						if (getMET_noHF() > 35.) hmeteff->Fill(7, EventWeight); 
						if (getMET_noHF() > 40.) hmeteff->Fill(8, EventWeight); 
						if (getMET_noHF() > 45.) hmeteff->Fill(9, EventWeight); 
						if (getMET_noHF() > 50.) hmeteff->Fill(10, EventWeight); 
						if (getMET_noHF() > 55.) hmeteff->Fill(11, EventWeight); 
						if (getMET_noHF() > 60.) hmeteff->Fill(12, EventWeight); 
						if (getMET_noHF() > 65.) hmeteff->Fill(13, EventWeight); 
						if (getMET_noHF() > 70.) hmeteff->Fill(14, EventWeight); 
						if (getMET_noHF() > 75.) hmeteff->Fill(15, EventWeight); 
						if (getMET_noHF() > 80.) hmeteff->Fill(16, EventWeight); 

						trkhmeteff->Fill(0., EventWeight);
						if (getMET() > 5.)  trkhmeteff->Fill(1, EventWeight); 
						if (getMET() > 10.) trkhmeteff->Fill(2, EventWeight); 
						if (getMET() > 15.) trkhmeteff->Fill(3, EventWeight); 
						if (getMET() > 20.) trkhmeteff->Fill(4, EventWeight); 
						if (getMET() > 25.) trkhmeteff->Fill(5, EventWeight); 
						if (getMET() > 30.) trkhmeteff->Fill(6, EventWeight); 
						if (getMET() > 35.) trkhmeteff->Fill(7, EventWeight); 
						if (getMET() > 40.) trkhmeteff->Fill(8, EventWeight); 
						if (getMET() > 45.) trkhmeteff->Fill(9, EventWeight); 
						if (getMET() > 50.) trkhmeteff->Fill(10, EventWeight); 
						if (getMET() > 55.) trkhmeteff->Fill(11, EventWeight); 
						if (getMET() > 60.) trkhmeteff->Fill(12, EventWeight); 
						if (getMET() > 65.) trkhmeteff->Fill(13, EventWeight); 
						if (getMET() > 70.) trkhmeteff->Fill(14, EventWeight); 
						if (getMET() > 75.) trkhmeteff->Fill(15, EventWeight); 
						if (getMET() > 80.) trkhmeteff->Fill(16, EventWeight); 

					}

				if (PassesMETCut())   {
					//if (PassesNJetsCut()) {
					FillYieldsHistograms(Muon,iMET, sys);      
					if(sys==Norm) FillKinematicHistos(Muon,iMET);

					if (getNBTags() == 1){
						FillYieldsHistograms(Muon, iExact1btag, sys);      
						if(sys==Norm) FillKinematicHistos(Muon,iExact1btag);
					}
					if (getNBTags() == 2){
						FillYieldsHistograms(Muon, iExact2btag, sys);      
						if(sys==Norm) FillKinematicHistos(Muon,iExact2btag);
					}
					if (PassesNJetsCut()) {
						FillYieldsHistograms(Muon,i2jets, sys);      
						if(sys==Norm) FillKinematicHistos(Muon,i2jets);
						if (PassesNBtagCut()) {
							FillYieldsHistograms(Muon,i1btag, sys);      
							if(sys==Norm) FillKinematicHistos(Muon,i1btag);
						}
					}
				}
      }
    }
  }

  ResetHypLeptons(); 
  if (gDoSF && PassTriggerEE() && IsElElEvent()){
		EventWeight = gWeight * getSF(Elec, sys);// * getTopPtSF();     
		//EventWeight = 1.;     
#ifdef DEBUG
		cout << " pass trigger + ee, ";
#endif
		// 0.115 = Fraction events with negative weight
		if(gIsMCatNLO) EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115)); 
		if (PassesMllVeto()){
			FillYieldsHistograms(Elec,iDilepton, sys);
			if(sys==Norm) FillKinematicHistos(Elec,iDilepton);
			if (PassesZVeto())    {
				FillYieldsHistograms(Elec,iZVeto, sys);      
				if(sys==Norm) FillKinematicHistos(Elec,iZVeto);
				if (PassesMETCut())   {
					//if (PassesNJetsCut()) {
					FillYieldsHistograms(Elec,iMET, sys);      
					if(sys==Norm) FillKinematicHistos(Elec,iMET);

					if (getNBTags() == 1){
						FillYieldsHistograms(Elec, iExact1btag, sys);      
						if(sys==Norm) FillKinematicHistos(Elec,iExact1btag);
					}
					if (getNBTags() == 2){
						FillYieldsHistograms(Elec, iExact2btag, sys);      
						if(sys==Norm) FillKinematicHistos(Elec,iExact2btag);
					}
					if (PassesNJetsCut()) {
						//if (PassesMETCut())   {
						FillYieldsHistograms(Elec,i2jets, sys);      
						if(sys==Norm) FillKinematicHistos(Elec,i2jets);
						if (PassesNBtagCut()) {
							FillYieldsHistograms(Elec,i1btag, sys);      
							if(sys==Norm) FillKinematicHistos(Elec,i1btag);
						}
					}
				}	  
			}
    }
  }
  ResetHypLeptons();
#ifdef DEBUG
  cout << " DONE!"<<endl;
#endif
}

//----------------------------------------------------------------------
// Passes
//----------------------------------------------------------------------
bool TOP5TeVAnalyzer::PassesMuonEta2p1(gChannel chan){
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;

	if (chan == Muon){  
		if (TMath::Abs(fHypLepton1.p.Eta()) < 2.1) return true; 
		if (TMath::Abs(fHypLepton2.p.Eta()) < 2.1) return true;
	}
	else if (chan == ElMu){
		if (TMath::Abs(fHypLepton1.p.Eta()) < 2.1) return true;
	}
	else if (chan == Elec){    
		return true;  
	}
	return false;
}

bool TOP5TeVAnalyzer::Passes3rdLeptonVeto(){
	return true; // don't apply third lepton veto...
	// Return false if there are not 2 signal leptons
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;  

	//  Int_t nvetoleptons = 0;
	for(Int_t i = 0; i < nMuon; ++i){ //elec
		if (fHypLepton1.index == i && fHypLepton1.type == 0) continue;
		if (fHypLepton2.index == i && fHypLepton2.type == 0) continue;
		//    if (IsVetoMuon(i)) return false;
		//     nvetoleptons++;
	}

	for(Int_t i = 0; i < nElec; ++i){ // muon
		if (fHypLepton1.index == i && fHypLepton1.type == 1) continue;
		if (fHypLepton2.index == i && fHypLepton2.type == 1) continue;
		//     if (IsVetoElectron(i)) return false;
		//    nvetoleptons++;
	}
	//  if (nvetoleptons > 0) return false;
	return true;
}

bool TOP5TeVAnalyzer::PassesMllVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
  CoutEvent(evt, Form(" # Event: InvMass = %2.2f", InvMass));
	if (InvMass < 20.)            return false; 
	return true;
}

bool TOP5TeVAnalyzer::PassesZVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
  if (InvMass > 76. && InvMass < 106.) return false;
  return true;
}

bool TOP5TeVAnalyzer::PassesNGenJetsCut(){
	TLorentzVector genjet;
  int nGenJets = 0;
	for(int k = 0; k < ngenJet; k++){ 
		genjet.SetPtEtaPhiM(genJet_pt[k], genJet_eta[k], genJet_phi[k], genJet_m[k]);
		if(genjet.Pt()>gJetEtCut && abs(genjet.Eta()) < 2.4) nGenJets++;
	}
	//cout << "[nGenJets ngenJet] = [" << nGenJets << " " << ngenJet << "]" << endl;
  if(nGenJets > 1) return true;
  return false;
}

bool TOP5TeVAnalyzer::PassesGenMetCut(){
  if(getTrueMET() <= 35) return false;
  return true;
}
bool TOP5TeVAnalyzer::PassesGenZVetoCut(){
  return passGenZVetoCut;
}

bool TOP5TeVAnalyzer::PassesNJetsCut(){
	//if (getNJets() <= 1) return false;
	CoutEvent(evt, Form(" # Event: njets = %i", getNJets()) );
	if (getNJets() <= 1) return false;
	return true;
}

bool TOP5TeVAnalyzer::PassesMETCut(){
	if (getMET() < 35.) return false;
	return true;
}

bool TOP5TeVAnalyzer::PassesNBtagCut(){
	if (getNBTags() < 1) return false;
	return true;
}

bool TOP5TeVAnalyzer::IsElMuEvent(){
	if (fChargeSwitch){      return (IsDileptonEvent()  == 3);   }
	return (IsDileptonEvent() == -3);
}

bool TOP5TeVAnalyzer::IsMuMuEvent(){
	if (fChargeSwitch){  return (IsDileptonEvent()  == 1); }
	return (IsDileptonEvent() == -1);
}

bool TOP5TeVAnalyzer::IsElElEvent(){
  if (fChargeSwitch){    return (IsDileptonEvent()  == 2); }
  return (IsDileptonEvent() == -2);
}

int TOP5TeVAnalyzer::IsDileptonEvent(){
#ifdef DEBUG
	cout << "IsDileptonEvent(): NLeptons =" << Lepton.size()<< endl;
#endif
  int nSelectedLeptons = Lepton.size();
	int result = 0;
	if     (nSelectedLeptons <  2) return 0;
	else if(nSelectedLeptons == 2){
		int select = Lepton[0].charge*Lepton[1].charge;
		if      (Lepton[0].type == 0 && Lepton[1].type == 0) result = 1; // mu/mu
		else if (Lepton[0].type == 1 && Lepton[1].type == 1) result = 2; // el/el
		else result = 3; // mu/el
		fHypLepton1 = lepton(Lepton[0]);
		fHypLepton2 = lepton(Lepton[1]);

		if(Lepton[0].type == 1 && Lepton[1].type == 0){
			fHypLepton1 = lepton(Lepton[1]);
			fHypLepton2 = lepton(Lepton[0]);
		}
		result *= select; // Add charge to result
	}
	else if(nSelectedLeptons >  2){
		Double_t dileppt = 0; Double_t dilepptmax = 0; int index1 = 0; int index2 = 0;
		for(int k = 1; k < nSelectedLeptons; k++){
			for(int h = 0; h < k; h++){
				dileppt = (Lepton[k].p + Lepton[h].p).Pt();
				if(dileppt > dilepptmax){
					index1 = h;
					index2 = k;
					dilepptmax = dileppt;
				}
			}
		}
    int select = Lepton[index1].charge*Lepton[index2].charge;
    if      (Lepton[index1].type == 0 && Lepton[index2].type == 0) result = 1; // mu/mu
    else if (Lepton[index1].type == 1 && Lepton[index2].type == 1) result = 2; // el/el
    else result = 3; // mu/el
    fHypLepton1 = lepton(Lepton[index1]);
    fHypLepton2 = lepton(Lepton[index2]);

    if(Lepton[index1].type == 1 && Lepton[index2].type == 0){
      fHypLepton1 = lepton(Lepton[index2]);
      fHypLepton2 = lepton(Lepton[index1]);
    }
    result *= select; // Add charge to result
	} 
#ifdef DEBUG
	cout << result;
	cout << " DONE!" << endl;
#endif
	return result;
}
//------------------------------------------------------------------------------
// LEPTON SELECTORS
//------------------------------------------------------------------------------
bool momentumComparator(lepton i, lepton j){ return (i.p.Pt()>j.p.Pt()); }

vector<lepton> TOP5TeVAnalyzer::SortLeptonsByPt(vector<lepton>& leptons){
  vector<lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;
}

int TOP5TeVAnalyzer::getSelectedLeptons(){
  // Loops over the total number of Muons and Electrons and returns the Number of Leptons.
  if (Lepton.size() > 0) {
    cout << "[WARNING]: you have called this function previously... RESETTING..."<<endl;
    Lepton.clear();
  }
  vector<lepton> tmp_lepton;
	nMuons = 0;
  nElecs = 0;
	TLorentzVector lep;
	Int_t thetype = 0;
  CoutEvent(evt, "-- Leptons --");
	for (Int_t i=0; i<nMuon;i++){
  CoutEvent(evt,Form(" >> # %i (MUON), pt: %2.2f, eta: %2.2f, charge: %i", i, MuonPt[i], MuonEta[i], MuonCharge[i]));
		if(IsTightMuon(i)){
			CoutEvent(evt, "   --> Is GOOD Muon!!!");
			thetype = 0;
			nMuons ++; 
			lep.SetPxPyPzE(MuonPx[i], MuonPy[i], MuonPz[i], MuonEnergy[i]); 
			lepton tmpLepton(lep, MuonCharge[i], thetype, i); 
			tmp_lepton.push_back(tmpLepton);
		}
	}
	for (Int_t i=0; i<nElec;i++){
  CoutEvent(evt,Form(" >> # %i (ELECTRON), pt: %2.2f, eta: %2.2f, charge: %i", i+nMuon, ElecEnergy[i], ElecEta[i], ElecCharge[i]));
		if(IsTightElectron(i)){
			CoutEvent(evt, "   --> Is GOOD Electon!!!");
			thetype = 1;
			nElecs++;
			lep.SetPxPyPzE(ElecPx[i], ElecPy[i], ElecPz[i], ElecEnergy[i]); 
			lepton tmpLepton(lep, ElecCharge[i], thetype, i); 
			tmp_lepton.push_back(tmpLepton);
		}
	}
  Lepton = SortLeptonsByPt(tmp_lepton);
  return Lepton.size();
}

//------------------------------------------------------------------------------
// Muon Selectors  
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
bool TOP5TeVAnalyzer::IsTightMuon(unsigned int iMuon,float ptcut){
   if (MuonPt[iMuon]              < ptcut) return false;
   if (TMath::Abs(MuonEta[iMuon]) > 2.1)   return false; CoutEvent(evt, Form(" -----------> Dxy            = %1.2f, ", MuonDxy[iMuon]));
   if (TMath::Abs(MuonDxy[iMuon]) >= 0.2 || MuonDxy[iMuon] == -99)  return false; CoutEvent(evt, Form(" -----------> Dz             = %1.2f, ", MuonDz[iMuon]));
   if (TMath::Abs(MuonDz[iMuon])  >= 0.5 || MuonDz[iMuon] == -99)  return false; CoutEvent(evt, Form(" -----------> MuonChi2NDF    = %1.2f, ", MuonChi2NDF[iMuon]));
   if (MuonChi2NDF[iMuon] > 10 || MuonChi2NDF[iMuon] == -99) return false;            CoutEvent(evt, Form(" -----------> MuonPixelHits  = %i, ", MuonPixelHits[iMuon]));
   if(MuonPixelHits[iMuon] < 1) return false;            CoutEvent(evt, Form(" -----------> Muon Hits = %i, ", MuonHits[iMuon]));
   if(MuonHits[iMuon] < 1) return false;            CoutEvent(evt, Form(" -----------> MuonStations   = %i, ", MuonStations[iMuon]));
	 if(MuonStations[iMuon] < 2) return false;             CoutEvent(evt, Form(" -----------> MuonTrkLayers  = %i, ", MuonTrkLayers[iMuon]));
	 if(MuonTrkLayers[iMuon] <= 5) return false;           CoutEvent(evt, Form(" -----------> MuonIso        = %1.2f, ", getMuonIso(iMuon))); 
	 if(gSysSource == Norm) fHMuonPtnoIso->Fill(MuonPt[iMuon]);
	 if(gSysSource == Norm) fHMuonEtanoIso->Fill(fabs(MuonEta[iMuon]));
   if(getMuonIso(iMuon) > 0.15) return false;
	 if(gSysSource == Norm) fHMuonPtIso->Fill(MuonPt[iMuon]);
	 if(gSysSource == Norm) fHMuonEtaIso->Fill(fabs(MuonEta[iMuon]));
	 //if(getMuonIso(iMuon) > 5) return false; // relaxin this cut to perform SS data driven estimate
   return true;
}

//------------------------------------------------------------------------------
// Electron Selectors
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/view/CMS/TopEGM#Spring15_selection_25ns
bool TOP5TeVAnalyzer::IsTightElectron(unsigned int iElec, float ptcut){
 	if (ElecPt[iElec] < ptcut) return false;
	if (TMath::Abs(ElecEta[iElec]) > 2.4) return false;
  if (ElecIDLoose[iElec] < 1) return false;
	return true;
}

void  TOP5TeVAnalyzer::setMET(float newmet){ MET = newmet;}
float TOP5TeVAnalyzer::getMET(){ return MET; }
float TOP5TeVAnalyzer::getMET_noHF(){ return MET_noHF; }
float TOP5TeVAnalyzer::getMETPhi(){ return MET_Phi;}
float TOP5TeVAnalyzer::getTrueMET(){ return TrueMET; }
float TOP5TeVAnalyzer::getTrueMETPhi(){ return TrueMET_Phi;}
float TOP5TeVAnalyzer::getTrueMET_noHF(){ return TrueMET_noHF; }
float TOP5TeVAnalyzer::getTrueMETPhi_noHF(){ return TrueMET_Phi_noHF;}
int   TOP5TeVAnalyzer::getNJets(){ return nJets;}
float TOP5TeVAnalyzer::getDRClosestJet(TLorentzVector lep){
	float minDR = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDR > lep.DeltaR(Jet[i].p)) minDR = lep.DeltaR(Jet[i].p);
	}
	return minDR;
}

float TOP5TeVAnalyzer::getMuonIso(unsigned int iMuon){
  Float_t charged = Get<Float_t>("MuonPFChIso", iMuon);
  Float_t neutral = Get<Float_t>("MuonPFNeuIso", iMuon);
  Float_t photons = Get<Float_t>("MuonPFPhoIso", iMuon);
  Float_t PUchar  = Get<Float_t>("MuonPFPUIso", iMuon);
  return (charged + TMath::Max(0., neutral + photons - 0.5*PUchar))/MuonPt[iMuon];
}

float TOP5TeVAnalyzer::getDPhiClosestJet(TLorentzVector lep){
	float minDphi = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDphi > TMath::Abs(lep.DeltaPhi(Jet[i].p))) minDphi = TMath::Abs(lep.DeltaPhi(Jet[i].p));
	}
	return minDphi;
}

int TOP5TeVAnalyzer::getLeadingJetbTag(){
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (Jet[i].isbtag) return i;
	}
	return  -1;
}

int TOP5TeVAnalyzer::getNBTags(){
	int ntags(0);
	for(UInt_t i = 0; i <Jet.size(); i++){
		if (Jet[i].isbtag) ntags++;
	}
	return ntags;
}

float TOP5TeVAnalyzer::getDeltaPhillJet(){
	if (fHypLepton1.index == -1) return -999.;
	if (fHypLepton2.index == -1) return -999.;
	Int_t ij = getLeadingJetbTag();
	if (ij < 0) return -999.; 
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	TLorentzVector jet = Jet[ij].p; 
	return TMath::Abs(dilep.DeltaPhi(jet));
}

int TOP5TeVAnalyzer::getSelectedJets(){
	int nj(0);
	if (Jet.size() > 0) {
		cout << "[WARNING]: you have called this function previously, RESETTING..."<<endl;
		Jet.clear();
	}
  TLorentzVector jt;
  for (Int_t i=0; i<nJet; i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;
	   CoutEvent(evt,"   >>> Is Good Jet! ");

    Float_t jetbtagi      = Jet_btagCSV[i];
    Float_t jetetai       = Jet_eta[i];
    Float_t jetenergyi    = Jet_energy[i];
    
    jt.SetPtEtaPhiE(JetPt.at(i), jetetai, JetPhi.at(i), jetenergyi);
    bool isbtag = false;
    if (gIsData) {
      //isbtag = fBTagSFnom->IsTagged(Jet_btagCSV[i], -999999, JetPt.at(i), jetetai);
      isbtag = 0;
    }
    else {
      Int_t   jetmcflavouri = Jet_mcFlavour[i];
      // official b-tag recommendation: use JetHadronFlavour instead of JetPartonFlavor
      /*
	if(TMath::Abs(Jet_mcFlavour[i]) == 5 || TMath::Abs(Jet_mcFlavour[i]) == 4){ 
	if (gSysSource == BtagUp)     btagSys =  1;
	if (gSysSource == BtagDown)   btagSys = -1;
	if (gSysSource == MisTagUp)   btagSys =  0;
	if (gSysSource == MisTagDown) btagSys =  0;
	}
	else {
	if (gSysSource == BtagUp)     btagSys =  0;
	if (gSysSource == BtagDown)   btagSys =  0;
	if (gSysSource == MisTagUp)   btagSys =  1;
	if (gSysSource == MisTagDown) btagSys = -1;
	}*/

//------------xxxxxxxxxxxxxxxxxx----------------
/*
      if      (gSysSource == BtagUp)      isbtag = fBTagSFbUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai); 
      else if (gSysSource == BtagDown)    isbtag = fBTagSFbDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      else if (gSysSource == MisTagUp)    isbtag = fBTagSFlUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      else if (gSysSource == MisTagDown)  isbtag = fBTagSFlDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      else                                isbtag = fBTagSFnom->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
*/
      // Use this line only to get raw numbers for syncronization
      //if(Get<Float_t>("Jet_btagCSV", i) > 0.89) isbtag=true;  // WP for 74
    }
    jet tmpjet(jt, isbtag, i);
    Jet.push_back(tmpjet);
    nj++;
  }
  return nj;
}

bool TOP5TeVAnalyzer::IsGoodJet(unsigned int ijet, float ptcut){
	float minDR = 0.3;
	// https://twiki.cern.ch/twiki/bin/view/CMS/TopJME
	if(ijet == 0) CoutEvent(evt,"-- JETS --");
	TLorentzVector jet;
  float SF = getJerSF(Jet_eta[ijet]);
	jet.SetPtEtaPhiE(JetPt.at(ijet), Jet_eta[ijet], JetPhi.at(ijet), Jet_energy[ijet]);
	CoutEvent(evt,Form(" >> # %i, pt: %2.2f, eta: %2.2f", ijet, jet.Pt(), jet.Eta()));
	if (jet.Pt() < ptcut)     return false;
	CoutEvent(evt, "   --> Pass Jet pt cut");
	if (abs(jet.Eta()) > 3.0) return false;
	CoutEvent(evt,"   --> Pass Eta cut ");
	if (getJet_id(ijet) <= 0)     return false;
	CoutEvent(evt,"   --> Pass Jet ID ");

	// https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_4_12/PhysicsTools/Heppy/python/physicsobjects/Jet.py#L66
	// Remove jets close to all selected leptons... 
	CoutEvent(evt,"   * Cleaning..."); 
	for(unsigned int i = 0; i < Lepton.size(); ++i){
		TString label = "Muon"; if(Lepton[i].type) label = "Electron"; 
		CoutEvent(evt, Form("   * Lep %i (", i) + label + Form("), LepPt: %2.2f, LepEta: %2.2f, DR: %1.3f", Lepton[i].p.Pt(), Lepton[i].p.Eta(), jet.DeltaR(Lepton[i].p)));
		if(jet.DeltaR(Lepton[i].p) < minDR) return false;
	}
	CoutEvent(evt,"   --> Pass Cleaning ");
  return true;
  //return false;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------
void TOP5TeVAnalyzer::SelectedGenLepton() {
	if (!gIsData) {
		nGenElec = 0; nGenMuon = 0; bool isOS = 0;
		for(int n = 0; n<ngenLep; n++){
			Int_t id = TMath::Abs(genLep_pdgId[n]);
			if(abs(genLep_MomPID[n]) != 24 && abs(genLep_MomPID[n]) != 15) continue; // lepton comming from W or tau
			if(abs(genLep_MomPID[n]) == 15){
				if(abs(genLep_GMomPID[n]) != 24) continue; // tau comming from W
			}
			if(gSampleName.Contains("Herwig") && abs(genLep_MomPID[n]) == 24){
		  	if(abs(genLep_GMomPID[n]) != 6) continue; // W comming from top
			} 
			if(gSelection >= 4){ //FIUDIAL
				if((id == 13) && ((genLep_pt[n] < 18) || (abs(genLep_eta[n]) > 2.1)) ) continue;
				if((id == 11) && ((genLep_pt[n] < 20) || (abs(genLep_eta[n]) > 2.4)) ) continue;
			}
			if (id == 11)  nGenElec++;
			if (id == 13)  nGenMuon++;
		}

		fHnGenEle->Fill(nGenElec);
		fHnGenMuo->Fill(nGenMuon);
		//nGenTau  = ngenTau;
		nGenLepton = nGenElec + nGenMuon;
		for (Int_t i=0; i<ngenLep; i++) {
			if (TMath::Abs(genLep_pdgId[i]) == 11){
				fHGenElePt->Fill(genLep_pt[i]);
			}
		}
		for (Int_t i=0; i<ngenLep; i++) {
			if (TMath::Abs(genLep_pdgId[i]) == 13){
				fHGenMuoPt->Fill(genLep_pt[i]);
			}
		}
  }
}


void TOP5TeVAnalyzer::propagateMET(TLorentzVector nVec, TLorentzVector oVec){
	TLorentzVector met;
	met.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	// set the pfMET to the old MET minus original vector plus new vector 
	setMET( (met+oVec-nVec).Pt() );
}

std::vector<int> TOP5TeVAnalyzer::CleanedJetIndices(float pt){
	std::vector<int> cleanJetsInd;
	for(Int_t i = 0; i <nJet; i++){
		if (IsGoodJet(i,pt)) cleanJetsInd.push_back(i);
	}
	return cleanJetsInd;
}

void TOP5TeVAnalyzer::SmearJetPts(int flag){
	// Modify the jet pt for systematics studies. Either shifted or smeared propagate to the MET!!
	TLorentzVector genjet;
	if(gIsData)   return; // don't smear data
	if(flag == 0) return; // 0 makes no sense
	for (Int_t i=0; i < nJet; i++){
		if (!gIsData) {
			float JESunc = 2.8 /100;  // 2.8 %
			if(flag == 1) JetPt.at(i) += JetPt.at(i)*JESunc;   // vary up   for flag 1 
			if(flag == 2) JetPt.at(i) -= JetPt.at(i)*JESunc;  // vary down for flag 2
			if(flag == 3){   // smear for flag 3 
        JetPt.at(i) = JetPt.at(i)*getJetSF(i, 1)/getJetSF(i, 0);
			}
		}
	}
}

float TOP5TeVAnalyzer::getJetSF(int index, bool up){
  if(gIsData) return 1;
	float JerSF = getJerSF(Jet_eta[index]);
  if(up) JerSF += getJerSFerror(Jet_eta[index]); 
	float factor = 0.;
	int genIndex = -1;
	TLorentzVector genjet;
	for(int k = 0; k < ngenJet; k++){ // Guess the index for the matched genJet
		if(genJet_matchId[k] == index){
			genIndex = k;
			genjet.SetPtEtaPhiM(genJet_pt[k], genJet_eta[k], genJet_phi[k], genJet_m[k]);
			break;
		}
	}
  TString report = "   ::: Not matched ";
	if (genIndex >= 0){
    report = "   ::: Matched jet ";
		factor = max(0.0, (genjet.Pt() + JerSF*(JetPt.at(index) - genjet.Pt()))/JetPt.at(index) );
		if(factor == 0) factor = 1;
	}
	else                factor = 1.; // Not matched..
  report += Form(" (%i) genJet_Pt-Eta, JetPt_eta = [%2.2f, %2.2f], [%2.2f, %2.2f]... SF = %2.2f", index, genjet.Pt(), genjet.Eta(), JetPt.at(index), Jet_eta[index], factor);

  CoutEvent(evt, report);
	return factor;
}

float TOP5TeVAnalyzer::getJerSF(float eta){
//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
	float aeta = abs(eta);
	if     (aeta < 0.5) return 1.095;
	else if(aeta < 0.8) return 1.120;
	else if(aeta < 1.1) return 1.097;
	else if(aeta < 1.3) return 1.103;
	else if(aeta < 1.7) return 1.118;
	else if(aeta < 1.9) return 1.100;
	else if(aeta < 2.1) return 1.162;
	else if(aeta < 2.3) return 1.160;
	else if(aeta < 2.5) return 1.161;
	else if(aeta < 2.8) return 1.209;
	else if(aeta < 3.0) return 1.567;
	else if(aeta < 3.2) return 1.384;
	else return 1;
}

float TOP5TeVAnalyzer::getJerSFerror(float eta){
	float aeta = abs(eta);
	if     (aeta < 0.5) return 0.018;
	else if(aeta < 0.8) return 0.028;
	else if(aeta < 1.1) return 0.017;
	else if(aeta < 1.3) return 0.033;
	else if(aeta < 1.7) return 0.014;
	else if(aeta < 1.9) return 0.033;
	else if(aeta < 2.1) return 0.044;
	else if(aeta < 2.3) return 0.048;
	else if(aeta < 2.5) return 0.060;
	else if(aeta < 2.8) return 0.059;
	else if(aeta < 3.0) return 0.321;
	else if(aeta < 3.2) return 0.033;
	else return 0.;
}

void TOP5TeVAnalyzer::ScaleLeptons(int flag){
	// Shift the lepton pts for systematics studies
	if(gIsData) return; // don't smear data
	if(flag == 0) return;
	//float scale = 0.003; // 0.3% for muons
	float scale = 0.005; 
	TLorentzVector oleps, leps, tmp;
	for(Int_t k = 0; k < nMuons; ++k){ //xxx Muon
		//if(TMath::Abs(LepGood_pdgId[i]) != 13) continue; 
		tmp.SetPxPyPzE(MuPx.at(k), MuPy.at(k), MuPz.at(k), MuEnergy.at(k)); 
		oleps += tmp;
		if(flag == 1) { MuPx.at(k) += scale*MuPx.at(k); MuPy.at(k) += scale*MuPy.at(k); }
		if(flag == 2) { MuPx.at(k) -= scale*MuPx.at(k); MuPy.at(k) -= scale*MuPy.at(k); }
		tmp.SetPxPyPzE(MuPx.at(k), MuPy.at(k), MuPz.at(k), MuEnergy.at(k)); 
		leps += tmp;
	}
	//scale = 0.0015; // 0.15% for electrons
	scale = 0.01; 
	for(Int_t i = 0; i < nElecs; ++i){ //xxx Elec
		tmp.SetPxPyPzE(ElPx.at(i), ElPy.at(i), ElPz.at(i), ElEnergy.at(i)); 
		oleps += tmp;
		if(flag == 1) { ElPx.at(i) += scale*ElPx.at(i); ElPy.at(i) += scale*ElPy.at(i); }
		if(flag == 2) { ElPx.at(i) -= scale*ElPx.at(i); ElPy.at(i) -= scale*ElPy.at(i); }
		tmp.SetPxPyPzE(ElPx.at(i), ElPy.at(i), ElPz.at(i), ElEnergy.at(i)); 
		leps += tmp;
	}
	propagateMET(leps, oleps);
	return;
}

Int_t TOP5TeVAnalyzer::getJet_id(Int_t ind){
/*	if(neutralSum[ind]/rawpt[ind] >= 0.99){CoutEvent(evt, "        - Does Not Pass  [neutralSum/rawjet >= 0.99]"); return 0;}
	if(eSum[ind]/rawpt[ind] >= 0.99){CoutEvent(evt, "        - Does Not Pass  [eSum/rawpt >= 0.99]"); return 0;}              
	if(totalMult[ind] <= 1){CoutEvent(evt, "        - Does Not Pass  [totalMult <= 1]"); return 0.; }                     

	if (abs(Jet_eta[ind]) < 2.4){
		if(chargedHardSum[ind]/rawpt[ind] <= 0.){CoutEvent(evt, "        - Does Not Pass  [chargedHardSum/rawjet <= 0.]"); return 0; }
		if(chargedSum[ind]/rawpt[ind] >= 0.99){CoutEvent(evt, "        - Does Not Pass  [chargedSum/rawpt >= 0.99]"); return 0;     } 
		if(chargedMult[ind] <= 0){CoutEvent(evt, "        - Does Not Pass  [chargedMult <= 0]"); return 0;    }               CoutEvent(evt,"       - Pass neutralSum/rawjet  ");
		}*/
  int NumberOfConstituents = CEM[ind] + CHM[ind] + MUM[ind] + NHM[ind] + NEM[ind];
	if(NumberOfConstituents <= 1) return 0;
	if(NHF[ind]>0.99) return 0;
	if(NEF[ind]>0.99) return 0;
	if (abs(Jet_eta[ind]) < 2.4){
		if(CHF[ind] <= 0) return 0;
		if(CHM[ind] <= 0) return 0; // Charged multiplicity
		if(CEF[ind] >= 0.99) return 0;
	}


	return 1;
}

void TOP5TeVAnalyzer::ScaleMET(int flag){
	// first try on MET uncertainty
	if(gIsData) return; // don't scale data
	TLorentzVector umet, jets, leps, tmp;
	jets.SetPtEtaPhiM(0., 0., 0., 0.); // init
	leps.SetPtEtaPhiM(0., 0., 0., 0.); // init
	tmp.SetPtEtaPhiM(0., 0., 0., 0.);  // init
	umet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.); // add met
	// subtract uncleaned jets
	for (Int_t i=0; i<nJet; i++) {
		if (!IsGoodJet(i, 15.)) continue; // do this on all jets in the event, not only the good jets with pT > 40
		tmp.SetPxPyPzE(Jet_px[i], Jet_py[i], Jet_pz[i], Jet_energy[i]);
		umet += tmp;
		jets += tmp;
		tmp.SetPtEtaPhiE(0., 0., 0., 0.);
	}
	// subtract muons
	for (Int_t i=0; i<nMuon; i++) { 
		if (!IsTightMuon(i)) continue;
		tmp.SetPxPyPzE(MuonPx[i], MuonPy[i], MuonPz[i], MuonEnergy[i]); 
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// subtract electrons
	for (Int_t i=0; i<nElec; i++) { 
		if (!IsTightElectron(i)) continue;
		tmp.SetPxPyPzE(ElecPx[i], ElecPy[i], ElecPz[i], ElecEnergy[i]); 
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// scale the unclustered energy by 10%
	if (flag == 0) tmp.SetPtEtaPhiE(1.1 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
	if (flag == 1) tmp.SetPtEtaPhiE(0.9 * umet.Pt(), umet.Eta(), umet.Phi(), umet.E());
	// subtract the leptons and jets again
	tmp -= leps;
	tmp -= jets;
	setMET(tmp.Pt());
	return;
}
