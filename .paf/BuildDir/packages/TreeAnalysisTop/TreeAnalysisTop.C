//========================================================
//  TreeAnalysisTop selector
//========================================================

#include "TreeAnalysisTop.h"
#include <iostream>
#include <math.h>

ClassImp(TreeAnalysisTop);
const float gJetEtCut = 30.;

//#define DEBUG

//------------------------------------------------------------------------------
// GetParameters
//------------------------------------------------------------------------------
void TreeAnalysisTop::GetParameters(){
  gSampleName    = GetParam<TString>("sampleName");
  gIsData        = GetParam<bool>("IsData");
  gWeight        = GetParam<float>("weight"); // cross section / events in the sample
  gLumiForPU     = GetParam<float>("LumiForPU");
  gTotalLumi     = GetParam<float>("TotalLumi");
  gDoSystStudies = GetParam<bool>("DoSystStudies");
  gUseCSVM       = GetParam<bool>("UseCSVM");
  gDoSF          = GetParam<bool>("DoSF");
  gDoDF          = GetParam<bool>("DoDF");
  gStopMass      = GetParam<float>("stopMass");
  gLspMass       = GetParam<float>("lspMass");

  gIsMCatNLO     = GetParam<bool>("IsMCatNLO");

  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gSampleName = %s",gSampleName.Data()));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gIsData = %d",gIsData ));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gWeight = %e", gWeight));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gLumiForPU = %f", gLumiForPU));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gTotalLumi = %f", gTotalLumi));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gDoSystStudies = %d", gDoSystStudies));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gUseCSVM = %d",gUseCSVM ));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gDoSF = %d", gDoSF));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gDoDF = %d",gDoDF ));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gStopMass = %f", gStopMass));
  PAF_INFO("TreeAnalysisTop::GetParameters()", Form("gLspMass = %f",gLspMass ));
}

//-----------------------------------------------------------------------------------
// GetTreeVariables
//-----------------------------------------------------------------------------------
void TreeAnalysisTop::GetTreeVariables(){
  nLepGood             = Get<Int_t>("nLepGood");
  nJet                 = Get<Int_t>("nJet");
  if(!gIsData){
    ngenLep              = Get<Int_t>("ngenLep");
    genWeight            = Get<Float_t>("genWeight");
  }
  for(int k = 0; k<nLepGood; k++){
    LepGood_px[k]      = Get<Float_t>("LepGood_px", k);
    LepGood_py[k]      = Get<Float_t>("LepGood_py", k);
    LepGood_pz[k]      = Get<Float_t>("LepGood_pz", k);
    LepGood_energy[k]  = Get<Float_t>("LepGood_energy", k);
    LepGood_pt[k]      = Get<Float_t>("LepGood_pt", k);
    LepGood_etaSc[k]   = Get<Float_t>("LepGood_etaSc", k);
    LepGood_eta[k]     = Get<Float_t>("LepGood_eta", k);
    LepGood_dxy[k]     = Get<Float_t>("LepGood_dxy", k);
    LepGood_dz[k]      = Get<Float_t>("LepGood_dz", k);
    LepGood_relIso03[k]= Get<Float_t>("LepGood_relIso03", k);
    LepGood_relIso04[k]= Get<Float_t>("LepGood_relIso04", k);
    LepGood_pdgId[k]   = Get<Int_t>("LepGood_pdgId", k);
    LepGood_charge[k]  = Get<Int_t>("LepGood_charge", k);
  }
  for(int k = 0; k<nJet; k++){
    Jet_px[k]          = Get<Float_t>("Jet_px", k);
    Jet_py[k]          = Get<Float_t>("Jet_py", k);
    Jet_pz[k]          = Get<Float_t>("Jet_pz", k);
    Jet_energy[k]      = Get<Float_t>("Jet_energy", k);
    Jet_eta[k]         = Get<Float_t>("Jet_eta", k);
    Jet_btagCSV[k]     = Get<Float_t>("Jet_btagCSV", k);
  }
  if(!gIsData){
    for(int k = 0; k<ngenLep; k++){
      genLep_pdgId[k]    = Get<Int_t>("genLep_pdgId", k);
      genLep_pt[k]       = Get<Float_t>("genLep_pt", k);
      genLep_eta[k]      = Get<Float_t>("genLep_eta", k);
      genLep_phi[k]      = Get<Float_t>("genLep_phi", k);
      genLep_mass[k]     = Get<Float_t>("genLep_mass", k);
    }
  }
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
// TreeAnalysisTop class constructor (make sure the pointers are initialized to zero)
//------------------------------------------------------------------------------------
TreeAnalysisTop::TreeAnalysisTop() : PAFChainItemSelector() {
	fHDummy = 0;
	hWeight = 0;
	fHTopPtWeight = 0;
	fHnGenEle = 0;
	fHnGenMuo = 0;
	fHGenElePt = 0;
	fHGenMuoPt = 0;

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
			fHLHEweights[ichan][icut] = 0;
			fHMET[ichan][icut] = 0;       
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
			fHMT2[ichan][icut] = 0;       
			fHPtllb[ichan][icut] = 0;       
			fHMeff[ichan][icut] = 0;       
			fHDelPhiLepMet[ichan][icut] = 0; 
			fHDelPhiJetMet[ichan][icut] = 0; 
			fHDelPhiPllbMet[ichan][icut] = 0; 
			fHDelPhiLepJet[ichan][icut] = 0; 

			fHCSVTag[ichan][icut] = 0; 
			fHTopD[ichan][icut] = 0;
			fHDelPhillJet[ichan][icut] = 0;

			fHDRLep[ichan][icut] = 0;
			fHDRLep0Jet[ichan][icut] = 0;
			fHDPhiLep0Jet[ichan][icut] = 0;
			fHLep0Iso[ichan][icut] = 0;
			fHDRLep1Jet[ichan][icut] = 0;
			fHDPhiLep1Jet[ichan][icut] = 0;
			fHLep1Iso[ichan][icut] = 0;

			fHStopMass[ichan][icut] = 0;
			fHChi0Mass[ichan][icut] = 0;
			fHChi0StopMass[ichan][icut] = 0;
			fHvertices[ichan][icut] = 0;
			fHgoodvertices[ichan][icut] = 0;

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
void TreeAnalysisTop::Initialise() {
	PAF_INFO("TreeAnalysisTop", "+ Initializing...");
	//PAF_INFO("TreeAnalysisTop", "+ Initializing paramenters...");
	GetParameters();
	//PAF_INFO("TreeAnalysisTop", "+ Sumw2 set for all histograms...");
	TH1::SetDefaultSumw2();
	fHDummy = CreateH1F("fHDummy","",1,0,1);
	//PAF_INFO("TreeAnalysisTop", "+ Initialise Yield histograms...");
	InitialiseYieldsHistos();
	//PAF_INFO("TreeAnalysisTop", "+ Initialise Kinematic histograms...");
	InitialiseKinematicHistos();
	if (!gIsData) {
		//PAF_INFO("TreeAnalysisTop", "+ Initialise Gen histograms...");
		InitialiseGenHistos();
	}
	//PAF_INFO("TreeAnalysisTop", "+ Initialise other histograms...");
	fHTopPtWeight  = CreateH1F("H_TopPtWeight" ,"TopPt Weight",100, 0, 2);

	fHnGenEle  = CreateH1F("fHnGenEle" , "nGenPromptElecs"  , 11, -1.5, 9.5);
	fHnGenMuo  = CreateH1F("fHnGenMuo" , "nGenPromptMuons"  , 11, -1.5, 9.5);
	fHGenElePt = CreateH1F("fHGenElePt", "GenPromptElecs Pt", 500, 0, 500);
	fHGenMuoPt = CreateH1F("fHGenMuoPt", "GenPromptMuons Pt", 500, 0, 500);

	if (gSampleName == "DoubleMuon"      ||     
			gSampleName == "DoubleEG"        || 
			gSampleName == "SingleMu"        || 
			gSampleName == "SingleElectron"  || 
			gSampleName == "MuonEG"	       ||	    
			gSampleName == "TTJets    "      ||  
			gSampleName == "DYJetsToLL_M50_aMCatNLO"     ||
			gSampleName == "DYJetsToLL_M10to50_aMCatNLO" ||
			gSampleName == "ZJets_MLM"       ||
			gSampleName == "DYJets_MLM"      ||
			gSampleName == "TTbar_Powheg"    ||  
			gSampleName == "TTbar_aMCatNLO"  ||  
			gSampleName == "TTbar_Powheg_Pythia6")  {  
		//	PAF_INFO("TreeAnalysisTop", "+ Initialise Drell-Yan histograms...");
		InitialiseDYHistos();
	}
	PAF_INFO("TreeAnalysisTop", "+ Initialise histograms for systematics studies...");
	InitialiseSystematicHistos();

	//	PU Reweight
	//--------------------------------------
	//PAF_INFO("TreeAnalysisTop", "+ Initialise Pile-Up reweighting tool...");
	fPUWeight     = new PUWeight(gLumiForPU,Summer2015_25ns_poisson,"2015_25ns");
	if (!gIsData) {
		fPUWeightUp   = new PUWeight(18494.9,   Summer2015_25ns_poisson,"2015_25ns"); //  18494.9  (5% down)
		fPUWeightDown = new PUWeight(20441.7,   Summer2015_25ns_poisson,"2015_25ns"); //  20441.7  (5% up  )
	}

	//if (gUseCSVM) fBTagSF   = new BTagSFUtil("CSVM","ABCD");//ReReco
	//else          fBTagSF   = new BTagSFUtil("CSVT","ABCD");//ReReco 

	//PAF_INFO("TreeAnalysisTop", "+ Initialise b-tag scale factors...");
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
	}

	//PAF_INFO("TreeAnalysisTop", "+ Initialise lepton scale factors...");
	fLeptonSF = new LeptonSF();

	//PAF_INFO("TreeAnalysisTop", "+ Initialise random 3...");
	fRand3 = new TRandom3(50);

	// No systematics activaded...
	gSysSource = Norm;
	PAF_INFO("TreeAnalysisTop", "+ Initialisation DONE.");
}

void TreeAnalysisTop::InitialiseGenHistos(){
	fHDeltaRLepJet[Muon] = CreateH1F("H_DeltaRLepJet_"+gChanLabel[Muon],"",1000,0.,5.);
	fHDeltaRLepJet[Elec] = CreateH1F("H_DeltaRLepJet_"+gChanLabel[Elec],"",1000,0.,5.);  
}

void TreeAnalysisTop::InitialiseDYHistos(){
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

void TreeAnalysisTop::InitialiseYieldsHistos(){
	hWeight = CreateH1F("hWeight","",200,0,1);
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

void TreeAnalysisTop::InitialiseKinematicHistos(){
	//  PAF_DEBUG("TreeAnalysisTop::InitialiseKinematicHistos()",Form("nWeights = %i", nWeights));
	//++ Kinematic histograms
	for (size_t ch=0; ch<gNCHANNELS; ch++){
		if (!gDoSF && ch==Muon) continue;
		if (!gDoSF && ch==Elec) continue;
		if (!gDoDF && ch==ElMu) continue;

		for (size_t cut=0; cut<iNCUTS; cut++){
			//PAF_DEBUG("TreeAnalysisTop::InitialiseKinematicHistos()",Form("cut = %i", cut));
			fHLHEweights[ch][cut]  = CreateH1F("H_LHEweights"  +gChanLabel[ch]+"_"+sCut[cut],"LHEweights", nWeights, -0.5, nWeights - 0.5);

			fHMET[ch][cut]         = CreateH1F("H_MET_"        +gChanLabel[ch]+"_"+sCut[cut],"MET"       , 3000, 0,300);
			fHLep0Eta[ch][cut]     = CreateH1F("H_Lep0Eta_"    +gChanLabel[ch]+"_"+sCut[cut],"Lep0Eta"   , 50  ,0 ,2.5);
			fHLep1Eta[ch][cut]     = CreateH1F("H_Lep1Eta_"    +gChanLabel[ch]+"_"+sCut[cut],"Lep1Eta"   , 50  ,0 ,2.5);
			fHDelLepPhi[ch][cut]   = CreateH1F("H_DelLepPhi_"  +gChanLabel[ch]+"_"+sCut[cut],"DelLepPhi" , 100,-3.2, 3.2);
			fHHT[ch][cut]          = CreateH1F("H_HT_"         +gChanLabel[ch]+"_"+sCut[cut],"HT"        , 4700,30,500);
			fHHT2[ch][cut]         = CreateH1F("H_HT2_"        +gChanLabel[ch]+"_"+sCut[cut],"HT2"       , 5000,0,500);
			fHHT3[ch][cut]         = CreateH1F("H_HT3_"        +gChanLabel[ch]+"_"+sCut[cut],"HT3"       , 3000,0,300);
			fHHT4[ch][cut]         = CreateH1F("H_HT4_"        +gChanLabel[ch]+"_"+sCut[cut],"HT4"       , 1000,0,1000);
			fHHT5[ch][cut]         = CreateH1F("H_HT5_"        +gChanLabel[ch]+"_"+sCut[cut],"HT5"       , 1200,0,1200);
			fHBtagJet0Pt[ch][cut]  = CreateH1F("H_BtagJet0Pt_" +gChanLabel[ch]+"_"+sCut[cut],"BtagJet0Pt", 2700,30,300);
			fHJet0Eta[ch][cut]     = CreateH1F("H_Jet0Eta_"	 +gChanLabel[ch]+"_"+sCut[cut],"Jet0Eta"   , 50,0,2.5);
			fHJet1Eta[ch][cut]     = CreateH1F("H_Jet1Eta_"	 +gChanLabel[ch]+"_"+sCut[cut],"Jet1Eta"   , 50,0,2.5);

			// Susy useful variables: 
			//fHMET2[ch][cut]          = CreateH1F("H_MET_"          +gChanLabel[ch]+"_"+sCut[cut],"MET"          , 800,0.,800);
			//fHInvMass2[ch][cut][0]   = CreateH1F("H_InvMass_"      +gChanLabel[ch]+"_"+sCut[cut],"InvMass"      , 800,0.,800);
			//fHDelLepPhi[ch][cut]     = CreateH1F("H_DelLepPhi_"    +gChanLabel[ch]+"_"+sCut[cut],"DelLepPhi"    , 100,0, 1);
			fHMT[ch][cut]            = CreateH1F("H_MT_"           +gChanLabel[ch]+"_"+sCut[cut],"MT"           , 800,0.,800);
			fHMT2[ch][cut]           = CreateH1F("H_MT2_"          +gChanLabel[ch]+"_"+sCut[cut],"MT2"          , 800,0.,800);
			fHPtllb[ch][cut]         = CreateH1F("H_Ptllb_"        +gChanLabel[ch]+"_"+sCut[cut],"PTllb"        , 800,0.,800);
			fHMeff[ch][cut]          = CreateH1F("H_Meff_"         +gChanLabel[ch]+"_"+sCut[cut],"Meff"         , 800,0.,800);
			fHDelPhiJetMet[ch][cut]  = CreateH1F("H_DelPhiJetMet_" +gChanLabel[ch]+"_"+sCut[cut],"DelPhiJetMet" , 100,0, 1);
			fHDelPhiLepMet[ch][cut]  = CreateH1F("H_DelPhiLepMet_" +gChanLabel[ch]+"_"+sCut[cut],"DelPhiLepMet" , 100,0, 1);
			fHDelPhiPllbMet[ch][cut] = CreateH1F("H_DelPhiPllbMet_"+gChanLabel[ch]+"_"+sCut[cut],"DelPhiPllbMet", 100,0, 1);
			fHDelPhiLepJet[ch][cut]  = CreateH1F("H_DelPhiLepJet_" +gChanLabel[ch]+"_"+sCut[cut],"DelPhiLepJet" , 100,0, 1);

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

			fHNJets[ch][cut][0]       = CreateH1F("H_NJets_"      +gChanLabel[ch]+"_"+sCut[cut],"NJets"     , 8 ,-0.5, 7.5);
			fHNBtagJets[ch][cut][0]   = CreateH1F("H_NBtagJets_"  +gChanLabel[ch]+"_"+sCut[cut],"NBtagJets" , 4 ,-0.5, 3.5);
			fHJet0Pt[ch][cut][0]      = CreateH1F("H_Jet0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet0Pt"    , 2700,30,300);
			fHJet1Pt[ch][cut][0]      = CreateH1F("H_Jet1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Jet1Pt"    , 2200,30,250);
			fHDiLepPt[ch][cut][0]     = CreateH1F("H_DiLepPt_"    +gChanLabel[ch]+"_"+sCut[cut],"DiLepPt"   , 1600,20,180); 
			fHLep0Pt[ch][cut][0]      = CreateH1F("H_Lep0Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep0Pt"    , 1800,20,200);
			fHLep1Pt[ch][cut][0]      = CreateH1F("H_Lep1Pt_"     +gChanLabel[ch]+"_"+sCut[cut],"Lep1Pt"    , 1800,20,200);

			// other variables 
			fHCSVTag[ch][cut]      = CreateH1F("H_CSVTag_"     +gChanLabel[ch]+"_"+sCut[cut], "NBtagsNJets"     , 1000, 0.0, 1.0);
			fHTopD[ch][cut]        = CreateH1F("H_TopD_"       +gChanLabel[ch]+"_"+sCut[cut], "TopDiscriminator", 1000,0.0,1.0);
			fHDelPhillJet[ch][cut] = CreateH1F("H_DelPhillJet_"+gChanLabel[ch]+"_"+sCut[cut], "DeltaPhi"        , 1000,0.0, TMath::Pi());

			// Different Top / Z topologies
			fHDRLep[ch][cut]       = CreateH1F("H_DRLep_"       +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep",       1000,0.0, 5.0);
			fHDRLep0Jet[ch][cut]   = CreateH1F("H_DRLep0Jet_"   +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep0Jet",   1000,0.0, 5.0);
			fHDPhiLep0Jet[ch][cut] = CreateH1F("H_DPhiLep0Jet_" +gChanLabel[ch]+"_"+sCut[cut], "DeltaPhiLep0Jet", 1000,0.0, TMath::Pi());
			fHLep0Iso[ch][cut]     = CreateH1F("H_Lep0Iso_"     +gChanLabel[ch]+"_"+sCut[cut], "Lep0Iso",         1000,0.0, 0.5);
			fHDRLep1Jet[ch][cut]   = CreateH1F("H_DRLep1Jet_"   +gChanLabel[ch]+"_"+sCut[cut], "DeltaRLep1Jet",   1000,0.0, 5.0);
			fHDPhiLep1Jet[ch][cut] = CreateH1F("H_DPhiLep1Jet_" +gChanLabel[ch]+"_"+sCut[cut], "DeltaPhiLep1Jet", 1000,0.0, TMath::Pi());
			fHLep1Iso[ch][cut]     = CreateH1F("H_Lep1Iso_"     +gChanLabel[ch]+"_"+sCut[cut], "Lep1Iso",         1000,0.0, 0.5);

			// STOP HISTOGRAMS:
			fHvertices[ch][cut]     = CreateH1F("H_Vtx_"+gChanLabel[ch]+"_"+sCut[cut],"", 51, -0.5, 50.5); 
			fHgoodvertices[ch][cut] = CreateH1F("H_goodVtx_"+gChanLabel[ch]+"_"+sCut[cut],"", 51, -0.5, 50.5); 
		}
	}
}

void TreeAnalysisTop::InitialiseSystematicHistos(){
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

				fHNJets[ch][cut][sys]       = CreateH1F("H_NJets_"      +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"NJets"     , 8 ,-0.5, 7.5);
				fHNBtagJets[ch][cut][sys]   = CreateH1F("H_NBtagJets_"  +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"NBtagJets" , 4 ,-0.5, 3.5);
				fHJet0Pt[ch][cut][sys]      = CreateH1F("H_Jet0Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Jet0Pt"    , 2700,30,300);
				fHJet1Pt[ch][cut][sys]      = CreateH1F("H_Jet1Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Jet1Pt"    , 2200,30,250);
				fHDiLepPt[ch][cut][sys]     = CreateH1F("H_DiLepPt_"    +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"DiLepPt"   , 1600,20,180); 
				fHLep0Pt[ch][cut][sys]      = CreateH1F("H_Lep0Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Lep0Pt"    , 1800,20,200);
				fHLep1Pt[ch][cut][sys]      = CreateH1F("H_Lep1Pt_"     +gChanLabel[ch]+"_"+sCut[cut]+"_"+SystName[sys],"Lep1Pt"    , 1800,20,200);
			}
		}
	}
}

//---------------------------------------------------------------------------------------------------
// Set objets, to be called once per event, saving information in tmp vectors for systematic studies.
//---------------------------------------------------------------------------------------------------
void TreeAnalysisTop::SetOriginalObjects(){
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
	MET_Phi = 0.;
	int k = 0;
	// Save original values for MET, Jets and Leptons
	TLorentzVector j;
	for (Int_t i=0; i < nJet; i++){    
		j.SetPxPyPzE(Jet_px[i], Jet_py[i], Jet_pz[i], Jet_energy[i]);
		JetEt.push_back(j.Et());
		JetPt.push_back(j.Pt());
		JetPhi.push_back(j.Phi());
	}
	for (Int_t i=0; i < nLepGood; i++){     
		if(TMath::Abs(LepGood_pdgId[i]) == 11){
			ElPx.push_back(LepGood_px[i]); 
			ElPy.push_back(LepGood_py[i]);
			ElPz.push_back(LepGood_pz[i]);
			ElEnergy.push_back(LepGood_energy[i]);
		}
	}
	for (Int_t i=0; i<nLepGood; i++){
		if(TMath::Abs(LepGood_pdgId[i]) == 13){
			MuPx.push_back(LepGood_px[i]);
			MuPy.push_back(LepGood_py[i]); 
			MuPz.push_back(LepGood_pz[i]); 
			MuEnergy.push_back(LepGood_energy[i]); 
		}
	}
	MET     = Get<Float_t>("met_pt"); //met
	MET_Phi = Get<Float_t>("met_phi"); //met
}

void TreeAnalysisTop::SetEventObjects(){
	ResetHypLeptons();

	fChargeSwitch = false;
	EventWeight = 1.;

	// USEFUL COUNTERS
	nGenLepton = 0;
	nGenElec   = 0;
	nGenMuon   = 0;
	nGenTau    = 0;
	nTauElec   = 0;
	nTauMuon   = 0;
	nGoodVertex = 0;
	nVertex     = 0;
	nBtags      = 0;
	nJets       = 0;
	nMuon       = 0;
	nElec       = 0;
	nLeptons    = 0;

	//// READ AND SAVE OBJETS...
	Jet.clear();
	Lepton.clear();

	nLeptons = getSelectedLeptons();
	nJets    = getSelectedJets();
	nBtags   = getNBTags();
}

void TreeAnalysisTop::ResetOriginalObjects(){
	// Save original values for MET, Jets and Leptons
	TLorentzVector j;
	for (Int_t i=0; i < nJet; i++){    
		j.SetPxPyPzE(Jet_px[i],	Jet_py[i], Jet_pz[i],	Jet_energy[i]);
		JetEt[i]  = j.Et();
		JetPt[i]  = j.Pt();
		JetPhi[i] = j.Phi();
	}
	int k = 0;
	for (Int_t i=0; i<nLepGood; i++){
		if(TMath::Abs(LepGood_pdgId[i]) == 11){
			ElPx[k] = LepGood_px[i]; 
			ElPy[k] = LepGood_py[i];
			ElPz[k] = LepGood_pz[i];
			ElEnergy[k] = LepGood_energy[i];
			k++;
		}
	}
	k = 0;
	for (Int_t i=0; i<nLepGood; i++){ 
		if(TMath::Abs(LepGood_pdgId[i]) == 13){
			MuPx[k] = LepGood_px[i]; 
			MuPy[k] = LepGood_py[i]; 
			MuPz[k] = LepGood_pz[i];
			MuEnergy[k] = LepGood_energy[i];
			k++;
		}
	}
	setMET(Get<Float_t>("met_pt")); //met
}

void TreeAnalysisTop::ResetHypLeptons(){
  TLorentzVector vec(0., 0., 0., 0.);
  fHypLepton1 = lepton(vec, 0, -1, -1);
  fHypLepton2 = lepton(vec, 0, -1, -1);
}

//-----------------------------------------------------------------------
// InsideLoop
//-----------------------------------------------------------------------
void TreeAnalysisTop::InsideLoop() {
	fHDummy->Fill(0.5);
	//if (METFilter() == false) return;

	// Calculate PU Weight
	PUSF = 1.;
	if (!gIsData)
		PUSF = fPUWeight->GetWeight(Get<Int_t>("nTrueInt")); //True       //nTruePU

	// Init data members ........................................................
	GetTreeVariables();
	SetOriginalObjects();
	SetEventObjects();


	// Get number of generated leptons ........................................................
	if (!gIsData) {
		SelectedGenLepton();
		if (gSampleName == "TTJets_MadSpin"       ||
				gSampleName == "TTbar_Powheg"         ||        
				gSampleName == "TTbar_Powheg_Pythia6" ||        
				gSampleName == "TTJets_aMCatNLO"      ||        
				gSampleName == "TTJets"               ){
			Float_t Weight = 1.; 
			TLorentzVector top;
			for (Int_t t=0; t<Get<Int_t>("nGenTop"); t++){ 
				top.SetPtEtaPhiM(Get<Float_t>("GenTop_pt",   t), 
						Get<Float_t>("GenTop_eta",  t), 
						Get<Float_t>("GenTop_phi",  t),
						Get<Float_t>("GenTop_mass", t));   
				Float_t pt    = TMath::Min(top.Pt(), 400.);
				Float_t topSF = TMath::Exp(0.156 - 0.00137 * pt);
				Weight *= topSF;
			}
			fHTopPtWeight->Fill(TMath::Sqrt(Weight));
		}
		bool isEmuDilepton = (nGenElec+nGenMuon)>=2;  // for dileptonic selection
		//bool isEmuDilepton = (nGenElec+nGenMuon)< 2;  // for semileptonic selection
		//bool isEmuDilepton = nGenElec>=1 && nGenMuon>=1;  // for dileptonic e-mu selection
		//bool isEmuDilepton = ( (nGenElec==1 && nGenMuon==0) || (nGenElec==0 && nGenMuon==1) ); // for semileptonic selection
		if(!isEmuDilepton){
			if (gSampleName == "TTbar_Powheg"              ) return;
			if (gSampleName == "TTJets"                    ) return;
			if (gSampleName == "TTJets_MLM"                ) return;
			if (gSampleName == "TTbar_Powheg_ScaleUp"      ) return;
			if (gSampleName == "TTbar_Powheg_ScaleDown"    ) return;
			if (gSampleName == "TTbar_Powheg_mtop1695"     ) return;
			if (gSampleName == "TTbar_Powheg_mtop1755"     ) return;
			if (gSampleName == "TTbar_Powheg_2L"	         ) return;
			if (gSampleName == "TTbar_Powheg_Pythia6"      ) return;
			if (gSampleName == "TTbar_Powheg_Herwig"       ) return;
			if (gSampleName == "TTJets_aMCatNLO"	         ) return;
			if (gSampleName == "TTJets_aMCatNLO_ScaleUp"   ) return;
			if (gSampleName == "TTJets_aMCatNLO_ScaleDown" ) return;
			if (gSampleName == "TTJets_aMCatNLO_mtop1695"  ) return;
			if (gSampleName == "TTJets_aMCatNLO_mtop1755"  ) return;
		}
		// Fill Gen Info 
		//----------------------------------------------------------------------------
		TLorentzVector lep,jet;
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
	if (gSampleName == "DoubleMuon"      ||     
			gSampleName == "DoubleEG"        || 
			gSampleName == "SingleMu"        || 
			gSampleName == "SingleElectron"  || 
			gSampleName == "MuonEG"	       ||	  
			gSampleName == "ZJets_MLM"       ||
			gSampleName == "DYJets_MLM"      ||
			gSampleName == "DYJetsToLL_M50_aMCatNLO"     ||
			gSampleName == "DYJetsToLL_M10to50_aMCatNLO" 
		 )  {
		FillDYHistograms();
	}

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

  // Pile Up sytematics ....................................................................
  ResetOriginalObjects();
  if (!gIsData)
    PUSF = fPUWeightUp->GetWeight(Get<Int_t>("nTrueInt")); //nTruePU
  gSysSource = PUUp;
  SetEventObjects();
  FillYields(PUUp);
  
  ResetOriginalObjects();
  if (!gIsData)
    PUSF = fPUWeightDown->GetWeight(Get<Int_t>("nTrueInt")); //nTruePU
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

void TreeAnalysisTop::Summary(){}

//------------------------------------------------------------------------------
// TRIGGER INFORMATION
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
bool TreeAnalysisTop::PassTriggerMuMu() {
	Bool_t pass = Get<Float_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
	return pass;
}

bool TreeAnalysisTop::PassTriggerEE(){ 
	Bool_t pass = Get<Float_t>("HLT_BIT_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
	return pass;
}

bool TreeAnalysisTop::PassTriggerEMu(){ 
	Bool_t pass = (Get<Float_t>("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v") || 
			Get<Float_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v")); 
	return pass;
}

//------------------------------------------------------------------------------
// Get METHODS
//------------------------------------------------------------------------------
float TreeAnalysisTop::getHT(){
  float ht(0);
  for (unsigned int i=0; i<Jet.size(); i++) ht+=Jet[i].p.Pt();
  return ht;
}
float TreeAnalysisTop::getJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return Jet[ind].p.Pt();
}
float TreeAnalysisTop::getJetEtaIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return TMath::Abs(Jet[ind].p.Eta());
}
float TreeAnalysisTop::getBtagJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  Int_t btagInd = 0;
  if (ind==0) btagInd = getLeadingJetbTag();
  else  return -999.;
  return Jet[btagInd].p.Pt();
}

float TreeAnalysisTop::getMT(gChannel chan){
  float ptl1 = fHypLepton1.p.Pt();
	float ptl2 = fHypLepton2.p.Pt();
	float dphi = getDelPhill();
	return TMath::Sqrt(2*ptl1*ptl2*(1-TMath::Cos(dphi)));
}

float TreeAnalysisTop::getDPhiLepJet(){
	if (fHypLepton1.index == -1) return -999.; if (fHypLepton2.index == -1) return -999.;
	// Int_t ij = getLeadingJetbTag(); if (ij < 0) return -999.; 
	if(Jet.size()<1) return -999.;
	TLorentzVector jet = Jet[0].p;
	TLorentzVector plep = fHypLepton1.p;
	if (plep.Pt() < fHypLepton2.p.Pt()) plep = fHypLepton2.p;
	return TMath::Abs(plep.DeltaPhi(jet));
}

float TreeAnalysisTop::getDelPhill(){ return fHypLepton1.p.DeltaPhi(fHypLepton2.p);}

float TreeAnalysisTop::getDPhiJetMet(){
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	return getDPhiClosestJet(pmet);
}

float TreeAnalysisTop::getDPhiLepMet(){
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	TLorentzVector plep = fHypLepton1.p;
	if (plep.Pt() < fHypLepton2.p.Pt()) plep = fHypLepton2.p;
	return TMath::Abs(plep.DeltaPhi(pmet));
}

float TreeAnalysisTop::getMT2(gChannel chan){
	TLorentzVector plep1, plep2;
	if (chan == Muon) {
		plep1.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
		plep2.SetPxPyPzE(MuPx.at(1), MuPy.at(1), MuPz.at(1), MuEnergy.at(1));
	}
	if (chan == Elec) {
		plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
		plep2.SetPxPyPzE(ElPx.at(1), ElPy.at(1), ElPz.at(1), ElEnergy.at(1));
	}
	if (chan == ElMu) {
		plep1.SetPxPyPzE(ElPx.at(0), ElPy.at(0), ElPz.at(0), ElEnergy.at(0));
		plep2.SetPxPyPzE(MuPx.at(0), MuPy.at(0), MuPz.at(0), MuEnergy.at(0));
	}
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	double pa[3]; double pb[3]; double pmiss[3];
	pmiss[0] = 0.; // irrelevant
	pmiss[1] = pmet.Px(); pmiss[2] = pmet.Py();
	pa[0] = 0.; pa[1] = plep1.Px(); pa[2] = plep1.Py();
	pb[0] = 0.; pb[1] = plep2.Px(); pb[2] = plep2.Py();
	mt2 MT2bisect;
	MT2bisect.set_momenta(pa, pb, pmiss);
	MT2bisect.set_mn(0.); // testmass
	double MT2 = MT2bisect.get_mt2();
	return MT2;
}

float TreeAnalysisTop::getMeff(){
	if(Jet.size()<2) return -999.;
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	return fHypLepton1.p.Pt() + fHypLepton2.p.Pt() + Jet[0].p.Pt() + Jet[1].p.Pt() + pmet.Pt();
}

TLorentzVector TreeAnalysisTop::getPtllb(){
	TLorentzVector pmet; pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	TLorentzVector pl1 = fHypLepton1.p; TLorentzVector pl2 = fHypLepton2.p;
	return pl1 + pl2 + pmet;
}

float TreeAnalysisTop::getDPhibMet(){
	TLorentzVector pmet; pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	TLorentzVector Ptllb = getPtllb();
	return pmet.DeltaPhi(Ptllb); 
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float TreeAnalysisTop::getJERScale(int jet){
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
float TreeAnalysisTop::getJERScaleUp(int jet){
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
float TreeAnalysisTop::getJERScaleDown(int jet){
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

float TreeAnalysisTop::getErrPt(float Pt, float Eta) {
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

float TreeAnalysisTop::getLeptonError(gChannel chan){
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

float TreeAnalysisTop::getTriggerError(gChannel chan){
	float trig(0.);
	if (chan==Muon) trig = fLeptonSF->GetDoubleMuSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	if (chan==ElMu) trig = fLeptonSF->GetMuEGSF_err    (fHypLepton2.p.Eta(),fHypLepton1.p.Eta());
	if (chan==Elec) trig = fLeptonSF->GetDoubleElSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	return trig;
}

float TreeAnalysisTop::getSF(gChannel chan) {
	if (gIsData)              return 1.; //Don't scale data
	float id1(1.),id2(1.), trig(1.);
	float err1(0.), err2(0.), err_trg(0.);
	if (chan == Muon){
		id1  = fLeptonSF->GetTightMuonSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		id2  = fLeptonSF->GetTightMuonSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
		trig = fLeptonSF->GetDoubleMuSF (fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
  } 
  else if (chan == Elec){
    id1  = fLeptonSF->GetTightElectronSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta()); 
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta()); 
    trig = fLeptonSF->GetDoubleElSF     (fHypLepton1.p.Eta(),fHypLepton2.p.Eta()); 
  }
  else if (chan == ElMu){
    id1  = fLeptonSF->GetTightMuonSF    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta()); 
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
    trig = fLeptonSF->GetMuEGSF         (fHypLepton2.p.Eta(),fHypLepton1.p.Eta());
  }
  return (PUSF*id1*id2*trig);
}

float TreeAnalysisTop::getTopPtSF(){
	// Return SF of the pt pt of the top 
	// Only apply SF if the process is ttbar...
	if(!gSampleName.Contains("TTJets")) return 1.;
	if (gSysSource==TopPt) {
		TLorentzVector top;
		Float_t topSF = 0.;
		Float_t Weight = 1.; 
		if (!gIsData) {
			Int_t nGenTop = Get<Int_t>("nGenTop");
			if (nGenTop != 2) return 1.; 
			for (Int_t t=0; t<nGenTop; t++){ 
				top.SetPtEtaPhiM(Get<Float_t>("GenTop_pt",   t), 
						Get<Float_t>("GenTop_eta",  t), 
						Get<Float_t>("GenTop_phi",  t), 
						Get<Float_t>("GenTop_mass", t));   
				Float_t pt = TMath::Min(top.Pt(), 400.);
				topSF = TMath::Exp(0.156 - 0.00137 * pt);
				Weight *= topSF;
			}
			Weight = TMath::Sqrt(Weight);
		}
		return Weight;
	}
	return 1.;
}

//--------------------------------------------------------------------------
// Fill histograms      
//------------------------------------------------------------------------
void TreeAnalysisTop::FillDYHistograms(){
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
			if (1)   {  // No vut in MET for e-mu channel
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

			//if (PassesMETCut())   {
			if (1)   {  // No cut in MET for e-mu channel
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

			//if (PassesMETCut())   {
			if (1)   {  // No vut in MET for e-mu channel
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
void TreeAnalysisTop::FillKinematicHistos(gChannel chan, iCut cut){
#ifdef DEBUG
  cout << "Filling KinematicHistos("<<chan<<","<<cut<<")... ";
  cout << fHypLepton1.index << " , " << fHypLepton2.index << endl;
#endif

	if (gSysSource != Norm)      return;  //only fill histograms for nominal distributions...
	if (fChargeSwitch == true  ) return;

	if (!gIsData) {
		for(int i = 0; i<Get<Int_t>("nLHEweight"); i++){
			fHLHEweights[chan][cut]->Fill(i, EventWeight*Get<Float_t>("LHEweight_wgt", i));
		}
	}

	//++ met info
  fHMET[chan][cut]           ->Fill(getMET(),             EventWeight);
  fHMT[chan][cut]            ->Fill(getMT(chan),          EventWeight);
  fHMT2[chan][cut]           ->Fill(getMT2(chan),         EventWeight);
  fHMeff[chan][cut]          ->Fill(getMeff(),            EventWeight);
  fHPtllb[chan][cut]         ->Fill(getPtllb().Pt(),      EventWeight);
  fHDelPhiLepJet[chan][cut]  ->Fill(TMath::Abs(getDPhiLepJet())/pi, EventWeight);
  fHDelPhiJetMet[chan][cut]  ->Fill(TMath::Abs(getDPhiJetMet())/pi, EventWeight); 
  fHDelPhiLepMet[chan][cut]  ->Fill(TMath::Abs(getDPhiLepMet())/pi, EventWeight); 
  fHDelPhiPllbMet[chan][cut] ->Fill(TMath::Abs(getDPhibMet()  )/pi, EventWeight); 
  
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

	int ib = getLeadingJetbTag();
	if (ib>=0)  fHCSVTag[chan][cut] ->Fill(Jet_btagCSV[ib], EventWeight);

	fHTopD[chan][cut] ->Fill(getTopD(), EventWeight);
	fHDelPhillJet[chan][cut]->Fill(getDeltaPhillJet(), EventWeight);

	//// Top/Z diff topology.
	fHDRLep[chan][cut]        ->Fill(fHypLepton1.p.DeltaR(fHypLepton2.p),     EventWeight);
	if (chan == Muon){
		fHLep0Iso[chan][cut]      ->Fill(getMuonIso(fHypLepton1.index) ,EventWeight);
		fHLep1Iso[chan][cut]      ->Fill(getMuonIso(fHypLepton2.index), EventWeight);
	}
	else if (chan==Elec){
		fHLep0Iso[chan][cut]      ->Fill(getElecIso(fHypLepton1.index) ,EventWeight);
		fHLep1Iso[chan][cut]      ->Fill(getElecIso(fHypLepton2.index), EventWeight);
	}
	else if (chan == ElMu){
		fHLep0Iso[chan][cut]      ->Fill(getMuonIso(fHypLepton1.index) ,EventWeight);
		fHLep1Iso[chan][cut]      ->Fill(getElecIso(fHypLepton2.index), EventWeight);
	}

	if (njets > 0) {
		fHDRLep0Jet[chan][cut]    ->Fill(getDRClosestJet(fHypLepton1.p),   EventWeight);
		fHDRLep1Jet[chan][cut]    ->Fill(getDRClosestJet(fHypLepton2.p),   EventWeight);
		fHDPhiLep0Jet[chan][cut]  ->Fill(getDPhiClosestJet(fHypLepton1.p), EventWeight);
		fHDPhiLep1Jet[chan][cut]  ->Fill(getDPhiClosestJet(fHypLepton2.p), EventWeight);
	}
	Int_t nVert =  Get<Int_t>("nVert");
	fHvertices[chan][cut]     ->Fill(nVert, EventWeight); // for now, using the same
	fHgoodvertices[chan][cut] ->Fill(nVert, EventWeight); // for now, using the same

#ifdef DEBUG
	cout << " DONE!" << endl;
#endif

}

void TreeAnalysisTop::FillYieldsHistograms(gChannel chan, iCut cut, gSystFlag sys){
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
		fHLepSys[chan][cut] ->Fill(getLeptonError(chan), EventWeight);
		fHTrigSys[chan][cut]->Fill(getTriggerError(chan),EventWeight);
	}
#ifdef DEBUG
	cout << " DONE! " << endl;
#endif
	return;
}
void TreeAnalysisTop::FillYields(gSystFlag sys){
	ResetHypLeptons();  

#ifdef DEBUG
	cout << "gDoDF= " << gDoDF << endl;
	cout << "PassTriggerEMu= " << PassTriggerEMu() << endl;
	cout << "Is ElMu/ElEl/MuMu Event= " 
		<< IsElMuEvent() << "/" 
		<< IsElElEvent() << "/" 
		<< IsMuMuEvent() << endl;
#endif

	if (gDoDF && PassTriggerEMu()  && IsElMuEvent()){
		// Define Hypothesis Leptons...
		EventWeight = gWeight * getSF(ElMu);// * getTopPtSF();
		hWeight -> Fill(EventWeight,1.);
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
		EventWeight = gWeight * getSF(Muon); //  * getTopPtSF();
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
				FillYieldsHistograms(Muon,iZVeto, sys);      
				if(sys==Norm) FillKinematicHistos(Muon,iZVeto);
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
						//if (PassesMETCut())   {
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
		EventWeight = gWeight * getSF(Elec);// * getTopPtSF();     
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
bool TreeAnalysisTop::PassesMuonEta2p1(gChannel chan){
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

bool TreeAnalysisTop::Passes3rdLeptonVeto(){
	return true; // don't apply third lepton veto...
	// Return false if there are not 2 signal leptons
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;  

	//  Int_t nvetoleptons = 0;
	for(Int_t i = 0; i < nLepGood; ++i){ //elec
		if (fHypLepton1.index == i && fHypLepton1.type == 0) continue;
		if (fHypLepton2.index == i && fHypLepton2.type == 0) continue;
		//    if (IsVetoMuon(i)) return false;
		//     nvetoleptons++;
	}

	for(Int_t i = 0; i < nLepGood; ++i){ // muon
		if (fHypLepton1.index == i && fHypLepton1.type == 1) continue;
		if (fHypLepton2.index == i && fHypLepton2.type == 1) continue;
		//     if (IsVetoElectron(i)) return false;
		//    nvetoleptons++;
	}
	//  if (nvetoleptons > 0) return false;
	return true;
}

bool TreeAnalysisTop::PassesMllVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
	if (InvMass < 20.)            return false; 
	return true;
}

bool TreeAnalysisTop::PassesZVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
  if (InvMass > 76. && InvMass < 106.) return false;
  return true;
}

bool TreeAnalysisTop::PassesTopDCut(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	if (getTopD() < 0.8) return false;
	return true;
}

bool TreeAnalysisTop::PassesNJetsCut(){
	if (getNJets() <= 1) return false;
	return true;
}

bool TreeAnalysisTop::PassesMETCut(){
	if (getMET() < 40.) return false;
	return true;
}

bool TreeAnalysisTop::PassesNBtagCut(){
	if (getNBTags() < 1) return false;
	return true;
}

bool TreeAnalysisTop::IsElMuEvent(){
	if (fChargeSwitch){      return (IsDileptonEvent()  == 3);   }
	return (IsDileptonEvent() == -3);
}

bool TreeAnalysisTop::IsMuMuEvent(){
	if (fChargeSwitch){  return (IsDileptonEvent()  == 1); }
	return (IsDileptonEvent() == -1);
}

bool TreeAnalysisTop::IsElElEvent(){
  if (fChargeSwitch){    return (IsDileptonEvent()  == 2); }
  return (IsDileptonEvent() == -2);
}

int TreeAnalysisTop::IsDileptonEvent(){
#ifdef DEBUG
	cout << "IsDileptonEvent(): NLeptons =" << Lepton.size()<< endl;
#endif
	if(Lepton.size() < 2) return 0;
	int select = Lepton[0].charge*Lepton[1].charge;
	int result = 0;
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

vector<lepton> TreeAnalysisTop::SortLeptonsByPt(vector<lepton>& leptons){
  vector<lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;
}

int TreeAnalysisTop::getSelectedLeptons(){
  // Loops over the total number of Muons and Electrons and returns the Number of Leptons.
  if (Lepton.size() > 0) {
    cout << "[WARNING]: you have called this function previously... RESETTING..."<<endl;
    Lepton.clear();
  }
  vector<lepton> tmp_lepton;
  nMuon = 0;
  TLorentzVector lep;
  Int_t thetype = 0;
  for (Int_t i=0; i<nLepGood;i++){
    if(IsTightMuon(i)){
      thetype = 0;
      nMuon ++; 
    }
    else if(IsTightElectron(i)){
      thetype = 1;
      nElec++;
    }
    else continue;
    lep.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]); 
    lepton tmpLepton(lep, LepGood_charge[i], thetype, i); 
    tmp_lepton.push_back(tmpLepton);
  }
  Lepton = SortLeptonsByPt(tmp_lepton);
  return Lepton.size();
}

bool TreeAnalysisTop::METFilter(){
  if (Get<Float_t>("Flag_HBHENoiseFilter") && 
      Get<Float_t>("Flag_HBHENoiseIsoFilter") && 
      Get<Float_t>("Flag_CSCTightHaloFilter") && 
      Get<Float_t>("Flag_goodVertices") && 
      Get<Float_t>("Flag_eeBadScFilter")) 
    return true;
  return false;
}

//------------------------------------------------------------------------------
// Muon Selectors  
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
bool TreeAnalysisTop::IsTightMuon(unsigned int iMuon,float ptcut){
   if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return false;
   if (LepGood_pt[iMuon]              < ptcut) return false;
   if (TMath::Abs(LepGood_eta[iMuon]) > 2.4)   return false;
   if (LepGood_relIso04[iMuon]        > 0.15)  return false; 
   if (TMath::Abs(LepGood_dxy[iMuon]) >= 0.2)  return false; 
   if (TMath::Abs(LepGood_dz[iMuon])  >= 0.5)  return false;
   if (Get<Int_t>("LepGood_tightId",iMuon)           < 1   )  return false;
   return true;
}

float TreeAnalysisTop::getMuonIso(int iMuon){
	if (iMuon < 0) return 9999.;
	if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return 9999.;
	return LepGood_relIso04[iMuon];
}

//------------------------------------------------------------------------------
// Electron Selectors
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/view/CMS/TopEGM#Spring15_selection_25ns
bool TreeAnalysisTop::IsTightElectron(unsigned int iElec, float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return false; 
	if (LepGood_pt[iElec] < ptcut) return false;
	if (TMath::Abs(LepGood_eta[iElec]) > 2.4) return false;
	if (TMath::Abs(LepGood_etaSc[iElec])  <= 1.479){ // Barrel
		if (LepGood_relIso03[iElec]        > 0.0766 ) return false;
		if (TMath::Abs(LepGood_dxy[iElec]) >= 0.0118) return false; 
		if (TMath::Abs(LepGood_dz[ iElec]) >= 0.373 ) return false;
		if (Get<Int_t>("LepGood_lostHits", iElec)         > 2      ) return false; 
	}
	if (TMath::Abs(LepGood_etaSc[iElec])  > 1.479){ // Endcaps
		if (LepGood_relIso03[iElec]        >  0.0678) return false;
		if (TMath::Abs(LepGood_dxy[iElec]) >= 0.0739) return false; 
		if (TMath::Abs(LepGood_dz[ iElec]) >= 0.602 ) return false;
		if (Get<Int_t>("LepGood_lostHits", iElec)         > 1      ) return false; 
	}
	if (TMath::Abs(LepGood_etaSc[iElec]) > 1.4442 && 
			TMath::Abs(LepGood_etaSc[iElec]) < 1.566) return false;
	if (Get<Int_t>("LepGood_convVeto", iElec) <  1) return false;
	if(TMath::Abs(Get<Int_t>("LepGood_tightId", iElec)) < 2) return false;
	return true;
}

float TreeAnalysisTop::getElecIso(int iElec){
	if (iElec < 0) return 9999.;
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return 9999.;
	return LepGood_relIso03[iElec];
}

float TreeAnalysisTop::getEACorrection(float eta){  // for a 0.3 CONE
	float abseta = TMath::Abs(eta);
	// numbers from https://indico.cern.ch/event/370494/contribution/2/material/slides/0.pdf
	if      (abseta < 0.8)                  return 0.1013; 
	else if (abseta >= 0.8 && abseta < 1.3) return 0.0988; 
	else if (abseta >= 1.3 && abseta < 2.0) return 0.0572; 
	else if (abseta >= 2.0 && abseta < 2.2) return 0.0842; 
	else if (abseta >= 2.2 && abseta < 5.0) return 0.1530; 
	else                                    return 0.1530; 

	cout << "[ERROR] getEACorrection(): No correction factor found!!" << endl;
	return -999999.;
}

void TreeAnalysisTop::setMET(float newmet){ MET = newmet;}

float TreeAnalysisTop::getMET(){ return MET; }

float TreeAnalysisTop::getMETPhi(){ return MET_Phi;}

int TreeAnalysisTop::getNJets(){ return nJets;}

float TreeAnalysisTop::getDRClosestJet(TLorentzVector lep){
	float minDR = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDR > lep.DeltaR(Jet[i].p)) minDR = lep.DeltaR(Jet[i].p);
	}
	return minDR;
}

float TreeAnalysisTop::getDPhiClosestJet(TLorentzVector lep){
	float minDphi = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDphi > TMath::Abs(lep.DeltaPhi(Jet[i].p))) minDphi = TMath::Abs(lep.DeltaPhi(Jet[i].p));
	}
	return minDphi;
}

int TreeAnalysisTop::getLeadingJetbTag(){
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (Jet[i].isbtag) return i;
	}
	return  -1;
}

int TreeAnalysisTop::getNBTags(){
	int ntags(0);
	for(UInt_t i = 0; i <Jet.size(); i++){
		if (Jet[i].isbtag) ntags++;
	}
	return ntags;
}

float TreeAnalysisTop::getDeltaPhillJet(){
	if (fHypLepton1.index == -1) return -999.;
	if (fHypLepton2.index == -1) return -999.;
	Int_t ij = getLeadingJetbTag();
	if (ij < 0) return -999.; 
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	TLorentzVector jet = Jet[ij].p; 
	return TMath::Abs(dilep.DeltaPhi(jet));
}

float TreeAnalysisTop::getTopD(){
	if (fHypLepton1.index == -1) return -999;
	if (fHypLepton2.index == -1) return -999;
	// Make Dilepton system
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	Float_t DeltaPhi(0.),TopD(0.);
	TLorentzVector jet;
	if (nJets == 0) return -999.;
	DeltaPhi = TMath::Abs(dilep.DeltaPhi(Jet[0].p));
	TopD     = 1 - (DeltaPhi/TMath::Pi()) * (1 - Jet_btagCSV[(Jet[0].index)]);
	return TopD;
}

int TreeAnalysisTop::getSelectedJets(){
	int nj(0);
	if (Jet.size() > 0) {
		cout << "[WARNING]: you have called this function previously, RESETTING..."<<endl;
		Jet.clear();
	}
	//int btagSys = 0;
	//TLorentzVector jtDisc;
	//for (UInt_t i=0; i<nDiscJet; i++) {
	//  jtDisc.SetPtEtaPhiE(DiscJet_pt[i], DiscJet_eta[i], DiscJet_phi[i], DiscJet_energy[i]);
	//}
  TLorentzVector jt;
  for (Int_t i=0; i<nJet; i++) {
    if(!IsGoodJet(i,gJetEtCut)) continue;

    Float_t jetbtagi      = Jet_btagCSV[i];
    Float_t jetetai       = Jet_eta[i];
    Float_t jetenergyi    = Jet_energy[i];
    
    jt.SetPtEtaPhiE(JetPt.at(i), jetetai, JetPhi.at(i), jetenergyi);
    bool isbtag = false;
    if (gIsData) {
      isbtag = fBTagSFnom->IsTagged(Jet_btagCSV[i], -999999, JetPt.at(i), jetetai);
    }
    else {
      Int_t   jetmcflavouri = Get<Int_t>  ("Jet_mcFlavour", i);
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
      if      (gSysSource == BtagUp)      isbtag = fBTagSFbUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai); 
      else if (gSysSource == BtagDown)    isbtag = fBTagSFbDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      else if (gSysSource == MisTagUp)    isbtag = fBTagSFlUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      else if (gSysSource == MisTagDown)  isbtag = fBTagSFlDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      else                                isbtag = fBTagSFnom->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
      // Use this line only to get raw numbers for syncronization
      //if(Get<Float_t>("Jet_btagCSV", i) > 0.89) isbtag=true;  // WP for 74
    }
    jet tmpjet(jt, isbtag, i);
    Jet.push_back(tmpjet);
    nj++;
  }
  return nj;
}

bool TreeAnalysisTop::IsGoodJet(unsigned int ijet, float ptcut){
  float minDR = 0.4;
// https://twiki.cern.ch/twiki/bin/view/CMS/TopJME
  TLorentzVector jet;
  jet.SetPtEtaPhiE(JetPt.at(ijet), Jet_eta[ijet], JetPhi.at(ijet), Jet_energy[ijet]);
  if (jet.Pt() < ptcut)     return false;
  if (abs(jet.Eta()) > 2.4) return false;
  if (Get<Int_t>("Jet_id",ijet) <= 0)    return false;
  //if (Jet_id[ijet] > 0) return true;
  // https://github.com/CERN-PH-CMG/cmg-cmssw/blob/CMGTools-from-CMSSW_7_4_12/PhysicsTools/Heppy/python/physicsobjects/Jet.py#L66
  // Remove jets close to all selected leptons... 
  for(unsigned int i = 0; i < Lepton.size(); ++i){
    if(jet.DeltaR(Lepton[i].p) < minDR) return false;
  }
  return true;
  //return false;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------
void TreeAnalysisTop::SelectedGenLepton() {
	if (!gIsData) {

		for(int n = 0; n<ngenLep; n++){
			Int_t id = TMath::Abs(genLep_pdgId[n]);
			if (id == 11)  nGenElec++;
			if (id == 13)  nGenMuon++;
		}
		// add generated leptons (e/mu) from decays of taus from W/Z/h decays
		for(int n = 0; n<Get<Int_t>("ngenLepFromTau"); n++){
			Int_t id = TMath::Abs(Get<Int_t>("genLepFromTau_pdgId", n));
			if (id == 11)  nGenElec++;
			if (id == 13)  nGenMuon++;
		}
		fHnGenEle->Fill(nGenElec);
		fHnGenMuo->Fill(nGenMuon);
		//nGenTau  = ngenTau;
		nGenLepton = nGenElec + nGenMuon + nGenTau;
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

void TreeAnalysisTop::propagateMET(TLorentzVector nVec, TLorentzVector oVec){
	TLorentzVector met;
	met.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	// set the pfMET to the old MET minus original vector plus new vector 
	setMET( (met+oVec-nVec).Pt() );
}

std::vector<int> TreeAnalysisTop::CleanedJetIndices(float pt){
	std::vector<int> cleanJetsInd;
	for(Int_t i = 0; i <nJet; i++){
		if (IsGoodJet(i,pt)) cleanJetsInd.push_back(i);
	}
	return cleanJetsInd;
}

void TreeAnalysisTop::SmearJetPts(int flag){
	// Modify the jet pt for systematics studies. Either shifted or smeared propagate to the MET!!
	if(gIsData)   return; // don't smear data
	if(flag == 0) return; // 0 makes no sense
	// select the jets you want to have... 
	if (!gIsData) {
		std::vector<int> cleanJets = CleanedJetIndices(15.);
		TLorentzVector ojets, jets, tmp, genJet;                            // 4-vec of old jets, newjets and a tmp-vector

		std::vector<int>::const_iterator it = cleanJets.begin();

		for( it = cleanJets.begin(); it != cleanJets.end(); ++it) {
			tmp.SetPtEtaPhiE(JetPt.at(*it), Jet_eta[*it], JetPhi.at(*it), Jet_energy[*it]);         // set temp to the jet
			// Gen info for jets... 
			genJet.SetPxPyPzE(Get<Float_t>("Jet_mcPx",*it), 
					Get<Float_t>("Jet_mcPy",*it), 
					Get<Float_t>("Jet_mcPz",*it),  
					Get<Float_t>("Jet_mcEnergy",*it));        
			ojets += tmp; 
			if(flag == 1) JetPt.at(*it) *= Get<Float_t>("Jet_corr_JECUp",*it);   // vary up   for flag 1 
			if(flag == 2) JetPt.at(*it) *= Get<Float_t>("Jet_corr_JECDown",*it); // vary down for flag 2;
			if(flag == 3){
				float jerScaleUp   = getJERScaleUp(*it);	  
				float jerScale     = getJERScale(*it);	  
				float factor1 = 0.;
				if (genJet.DeltaR(tmp) < 0.5) factor1 = max(0.0, genJet.Pt() + jerScale*(tmp.Pt() - genJet.Pt()) );
				else                          factor1 = tmp.Pt();
				float sigmaMC  = getErrPt(JetPt.at(*it), Jet_eta[*it]) / JetPt.at(*it);
				float factor   = fRand3->Gaus(1., sqrt(jerScale*jerScale -1.) * sigmaMC );
				JetPt.at(*it) = JetPt.at(*it) * factor;		// smear for flag 3
			}
			// set tmp to the scaled/smeared jet
			tmp.SetPtEtaPhiE(JetPt.at(*it), Jet_eta[*it], JetPhi.at(*it), Jet_energy[*it]);
			jets += tmp;  // add scaled/smeared jet to the new jets
		}
		propagateMET(jets, ojets);  // propagate this change to the MET
	}
}

void TreeAnalysisTop::ScaleLeptons(int flag){
	// Shift the lepton pts for systematics studies
	if(gIsData) return; // don't smear data
	if(flag == 0) return;
	//float scale = 0.003; // 0.3% for muons
	float scale = 0.005; 
	TLorentzVector oleps, leps, tmp;
	for(Int_t k = 0; k < nMuon; ++k){ //xxx Muon
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
	for(Int_t i = 0; i < nElec; ++i){ //xxx Elec
		if(TMath::Abs(LepGood_pdgId[i]) != 11) continue; 
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

float TreeAnalysisTop::weightNvtx(int nvtx){
	float weight = 1.0;
	if(gIsData) return weight;
	// weights from single lepton region based on nVertex
	if(nvtx== 1)      weight = 3.55238;
	else if(nvtx== 2) weight = 3.55238;
	else if(nvtx== 3) weight = 3.55238;
	else if(nvtx== 4) weight = 3.55238;
	else if(nvtx== 5) weight = 3.55238;
	else if(nvtx== 6) weight = 1.14672;
	else if(nvtx== 7) weight = 1.23742;
	else if(nvtx== 8) weight = 1.18588;
	else if(nvtx== 9) weight = 1.18755;
	else if(nvtx==10) weight = 1.18737;
	else if(nvtx==11) weight = 1.17951;
	else if(nvtx==12) weight = 1.19148;
	else if(nvtx==13) weight = 1.23285;
	else if(nvtx==14) weight = 1.18385;
	else if(nvtx==15) weight = 1.19163;
	else if(nvtx==16) weight = 1.18142;
	else if(nvtx==17) weight = 1.12507;
	else if(nvtx==18) weight = 1.09837;
	else if(nvtx==19) weight = 1.0351;
	else if(nvtx==20) weight = 0.976103;
	else if(nvtx==21) weight = 0.876084;
	else if(nvtx==22) weight = 0.807204;
	else if(nvtx==23) weight = 0.71911;
	else if(nvtx==24) weight = 0.624918;
	else if(nvtx==25) weight = 0.512945;
	else if(nvtx==26) weight = 0.449579;
	else if(nvtx==27) weight = 0.356727;
	else if(nvtx==28) weight = 0.288227;
	else if(nvtx==29) weight = 0.237089;
	else if(nvtx==30) weight = 0.187185;
	else if(nvtx==31) weight = 0.155465;
	else if(nvtx==32) weight = 0.106072;
	else if(nvtx==33) weight = 0.0863709;
	else if(nvtx==34) weight = 0.0784446;
	else if(nvtx==35) weight = 0.0731524;
	else if(nvtx==36) weight = 0.0533728;
	else if(nvtx==37) weight = 0.0478281;
	else if(nvtx==38) weight = 0.0174096;
	else if(nvtx==39) weight = 0.0205018;
	else if(nvtx==40) weight = 0.0164379;
	else if(nvtx==42) weight = 0.0156376;
	else if(nvtx==43) weight = 0.0239829;
	return weight;
}

void TreeAnalysisTop::ScaleMET(int flag){
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
	for (Int_t i=0; i<nLepGood; i++) { 
		if (!IsTightMuon(i)) continue;
		tmp.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]); 
		umet += tmp;
		leps += tmp;
		tmp.SetPtEtaPhiM(0., 0., 0., 0.);
	}
	// subtract electrons
	for (Int_t i=0; i<nLepGood; i++) { 
		if (!IsTightElectron(i)) continue;
		tmp.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]); 
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
