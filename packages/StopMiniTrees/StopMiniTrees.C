//========================================================
//  StopMiniTrees selector
//========================================================

#include "StopMiniTrees.h"
#include <iostream>
#include <math.h>

ClassImp(StopMiniTrees);
const float gJetEtCut = 30.;

//#define DEBUG

//------------------------------------------------------------------------------
// GetParameters
//------------------------------------------------------------------------------
void StopMiniTrees::GetParameters(){
  gSampleName    = GetParam<TString>("sampleName");
  gIsData        = GetParam<bool>("IsData");
  gWeight        = GetParam<float>("weight"); // cross section / events in the sample
  gLumiForPU     = GetParam<float>("LumiForPU");
  gDoSystStudies = GetParam<bool>("DoSystStudies");
  gTotalLumi     = GetParam<float>("TotalLumi");
  gUseCSVM       = GetParam<bool>("UseCSVM");
  gStopMass      = GetParam<Int_t>("stopMass");
  gLspMass       = GetParam<Int_t>("lspMass");
  gIsT2tt        = false;
  if(gSampleName.BeginsWith("T2tt")) gIsT2tt = true;

  gIsMCatNLO     = GetParam<bool>("IsMCatNLO");
  gCreateTree    = GetParam<bool>("CreateTree");

  PAF_INFO("StopMiniTrees::GetParameters()", Form("gSampleName = %s",gSampleName.Data()));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gIsData = %d",gIsData ));
  PAF_INFO("StopAnalyzer::GetParameters()", Form("gDoSystStudies = %d", gDoSystStudies));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gWeight = %e", gWeight));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gLumiForPU = %f", gLumiForPU));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gTotalLumi = %f", gTotalLumi));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gUseCSVM = %d",gUseCSVM ));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gStopMass = %i", gStopMass));
  PAF_INFO("StopMiniTrees::GetParameters()", Form("gLspMass = %i",gLspMass ));
}

//-----------------------------------------------------------------------------------
// GetTreeVariables
//-----------------------------------------------------------------------------------
void StopMiniTrees::GetTreeVariables(){
  nLepGood             = Get<Int_t>("nLepGood");
	nJet                 = Get<Int_t>("nJet");
  evt                  = Get<ULong64_t>("evt");
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

//------------------------------------------------------------------------------------
// StopMiniTrees class constructor (make sure the pointers are initialized to zero)
//------------------------------------------------------------------------------------
StopMiniTrees::StopMiniTrees() : PAFChainItemSelector() {
	fHDummy = 0;
	hWeight = 0;
	for (unsigned int ichan = 0; ichan < gNCHANNELS; ichan++) {
		for (unsigned int icut = 0; icut < iNCUTS; icut++) {
			fHDY_InvMass        [ichan][icut] = 0;
		}
	}
}

//-------------------------------------------------------------------
// Initialise
//-------------------------------------------------------------------
void StopMiniTrees::Initialise() {
	GetParameters();
	TH1::SetDefaultSumw2();
	fHDummy = CreateH1F("fHDummy","",1,0,1);
	hWeight = CreateH1F("hWeight","",200,0,1);

	InitialiseTree();

	//	PU Reweight
  fPUWeight     = new PUWeight(gLumiForPU, Spring2016_25ns_poisson_OOTPU, "2016_ichep");
  if (!gIsData) {
    fPUWeightUp   = new PUWeight(18494.9,  Spring2016_25ns_poisson_OOTPU, "2016_ichep"); //  18494.9 
    fPUWeightDown = new PUWeight(20441.7,  Spring2016_25ns_poisson_OOTPU, "2016_ichep"); //  20441.7 
  }

	// Initialise b-tag scale factors...
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

	fLeptonSF = new SusyLeptonSF(); // Initialise lepton scale factors
	fRand3 = new TRandom3(50);	// Initialise random 3
	gSysSource = Norm;
}

void StopMiniTrees::InitialiseTree(){
    fTree = CreateTree("sTopTree","Optimization tree");

    fTree->Branch("TWeight",      &TWeight,      "TWeight/F");
    fTree->Branch("TIsDoubleMuon",&TIsDoubleMuon,"TIsDoubleMuon/I");
    fTree->Branch("TIsDoubleElec",&TIsDoubleElec,"TIsDoubleElec/I");
    fTree->Branch("TIsElMu",      &TIsElMu,      "TIsElMu/I");
    fTree->Branch("TNJets",           &TNJets,         "TNJets/I");
    fTree->Branch("TNISRJets",        &TNISRJets,      "TNISRJets/I");
    fTree->Branch("TNJetsJESUp",           &TNJetsJESUp,         "TNJetsJESUp/I");
    fTree->Branch("TNJetsJESDown",           &TNJetsJESDown,         "TNJetsJESDown/I");
    fTree->Branch("TNJetsJER",           &TNJetsJER,         "TNJetsJER/I");
    fTree->Branch("TNJetsBtag",       &TNJetsBtag,     "TNJetsBtag/I");
    fTree->Branch("TNJetsBtagUp",     &TNJetsBtagUp,   "TNJetsBtagUp/I");
    fTree->Branch("TNJetsBtagDown",   &TNJetsBtagDown, "TNJetsBtagDown/I");
    fTree->Branch("TNJetsBtagMisTagUp",     &TNJetsBtagMisTagUp,   "TNJetsBtagMisTagUp/I");
    fTree->Branch("TNJetsBtagMisTagDown",   &TNJetsBtagMisTagDown, "TNJetsBtagMisTagDown/I");

    fTree->Branch("TMET",         &TMET,         "TMET/F");
    fTree->Branch("TGenMET",         &TGenMET,   "TGenMET/F");
    fTree->Branch("TMETJESUp",    &TMETJESUp,    "TMETJESUp/F");
    fTree->Branch("TMETJESDown",  &TMETJESDown,  "TMETJESDown/F");
    fTree->Branch("TMT2ll",       &TMT2ll,       "TMT2ll/F");
    fTree->Branch("TMT2bb",       &TMT2bb,       "TMT2bb/F");
    fTree->Branch("TMT2lblb",     &TMT2lblb,     "TMT2lblb/F");
    fTree->Branch("TMll",         &TMll,         "TMll/F");
    fTree->Branch("TMeff",        &TMeff,        "TMeff/F");
    fTree->Branch("TPtllb",       &TPtllb,       "TPtllb/F");
    fTree->Branch("THT",          &THT,          "THT/F");
    fTree->Branch("THTJESUp",     &THTJESUp,     "THTJESUp/F");
    fTree->Branch("THTJESDown",   &THTJESDown,   "THTJESDown/F");

    fTree->Branch("TMET_Phi",     &TMET_Phi,     "TMET_Phi/F");
    fTree->Branch("TdPhiPtllbMET",&TdPhiPtllbMET,"TdPhiPtllbMET/F");
    fTree->Branch("TMinDPhiMetJets",&TMinDPhiMetJets,"TMinDPhiMetJets/F");
    fTree->Branch("TdPhiJetMet",  &TdPhiJetMet,  "TdPhiJetMet/F");
    fTree->Branch("TdPhiLepMet",  &TdPhiLepMet,  "TdPhiLepMet/F");
    fTree->Branch("TdPhiLepJet",  &TdPhiLepJet,  "TdPhiLepJet/F");
    fTree->Branch("TdPhill",      &TdPhill,      "TdPhill/F");

    fTree->Branch("TLep1_Px",     &TLep1_Px,     "TLep1_Px/F");
    fTree->Branch("TLep1_Py",     &TLep1_Py,     "TLep1_Py/F");
    fTree->Branch("TLep1_Pz",     &TLep1_Pz,     "TLep1_Pz/F");
    fTree->Branch("TLep1_E" ,     &TLep1_E ,     "TLep1_E/F");
    fTree->Branch("TLep1_Carge",  &TLep1_Charge, "TLep1_Charge/F");
    fTree->Branch("TLep2_Px",     &TLep2_Px,     "TLep2_Px/F");
    fTree->Branch("TLep2_Py",     &TLep2_Py,     "TLep2_Py/F");
    fTree->Branch("TLep2_Pz",     &TLep2_Pz,     "TLep2_Pz/F");
    fTree->Branch("TLep2_E" ,     &TLep2_E ,     "TLep2_E/F");
    fTree->Branch("TLep2_Carge",  &TLep2_Charge, "TLep2_Charge/F");

    fTree->Branch("TJet_isBJet",       TJet_isBJet,       "TJet_isBJet[TNJets]/I");
    fTree->Branch("TJet_Px",           TJet_Px,           "TJet_Px[TNJets]/F");
    fTree->Branch("TJet_Py",           TJet_Py,           "TJet_Py[TNJets]/F");
    fTree->Branch("TJet_Pz",           TJet_Pz,           "TJet_Pz[TNJets]/F");
    fTree->Branch("TJet_E",            TJet_E,            "TJet_E[TNJets]/F");
    fTree->Branch("TJetJESUp_Pt",      TJetJESUp_Pt,      "TJetJESUp_Pt[TNJetsJESUp]/F");
    fTree->Branch("TJetJESDown_Pt",    TJetJESDown_Pt,    "TJetJESDown_Pt[TNJetsJESDown]/F");
    fTree->Branch("TJetJER_Pt",        TJetJER_Pt,        "TJetJER_Pt[TNJetsJER]/F");

    fTree->Branch("TWeight_LepEffUp",      &TWeight_LepEffUp,      "TWeight_LepEffUp/F");
    fTree->Branch("TWeight_LepEffDown",    &TWeight_LepEffDown,    "TWeight_LepEffDown/F");
    fTree->Branch("TWeight_TrigUp",        &TWeight_TrigUp,        "TWeight_TrigUp/F");
    fTree->Branch("TWeight_TrigDown",      &TWeight_TrigDown,      "TWeight_TrigDown/F");
    fTree->Branch("TWeight_FSUp",          &TWeight_FSUp,          "TWeight_FSUp/F");
    fTree->Branch("TWeight_FSDown",        &TWeight_FSDown,        "TWeight_FSDown/F");
    fTree->Branch("TWeight_PUUp",        &TWeight_PUUp,        "TWeight_PUUp/F");
    fTree->Branch("TWeight_PUDown",        &TWeight_PUDown,        "TWeight_PUDown/F");
}

void StopMiniTrees::SetOriginalObjects(){
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

void StopMiniTrees::SetEventObjects(){
	ResetHypLeptons();
	fChargeSwitch = false;
	EventWeight = 1.;

	// USEFUL COUNTERS
	nGenLepton = 0; nGenElec = 0; nGenMuon = 0;
	nGenTau = 0; nTauElec = 0; nTauMuon = 0;
	nGoodVertex = 0; nVertex = 0;
	nBtags = 0; nJets = 0;
	nMuon = 0; nElec = 0; nLeptons = 0;

	//// READ AND SAVE OBJETS...
	Jet.clear();
	Lepton.clear();
	nLeptons = getSelectedLeptons();
	nJets    = getSelectedJets();
	nBtags   = getNBTags();
}

void StopMiniTrees::ResetOriginalObjects(){
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

void StopMiniTrees::ResetHypLeptons(){
  TLorentzVector vec(0., 0., 0., 0.);
  fHypLepton1 = lepton(vec, 0, -1, -1);
  fHypLepton2 = lepton(vec, 0, -1, -1);
}

void StopMiniTrees::SetTreeVariables(gChannel chan){ 
  TWeight     = EventWeight;
  TNJets      = getNJets();
  TNISRJets   = gIsT2tt ? 0 : Get<Float_t>("nISRJet30");
  TNJetsBtag  = getNBTags(); 

  TIsDoubleElec = 0; TIsDoubleMuon = 0; TIsElMu = 0;
  if(chan == Muon) TIsDoubleMuon = 1;
  if(chan == Elec) TIsDoubleElec = 1;
  if(chan == ElMu) TIsElMu       = 1;
  TMET          = getMET();
  TGenMET       = gIsData? 0 : Get<Float_t>("met_genPt");
	TMT2ll        = getMT2ll(chan);
	TMT2bb        = getMT2b(chan);
	TMT2lblb      = getMT2lb(chan);
  TMll          = (fHypLepton1.p+fHypLepton2.p).M();
  TPtllb        = getPtllb().Pt();
  TMeff         = getMeff();
  THT           = getHT();
  THTJESDown    = getHT();
  THTJESUp      = getHT();
  TdPhiPtllbMET = getDPhibMet();
  TMinDPhiMetJets = getMinDPhiMetJets();
  TdPhiJetMet   = getDPhiJetMet();
  TdPhiLepMet   = getDPhiLepMet();
  TdPhiLepJet   = getDPhiLepJet();
  TdPhill       = getDelPhill();
  TMET_Phi      = getMETPhi();

  TLep1_Px      = fHypLepton1.p.Px();
  TLep1_Py      = fHypLepton1.p.Py();
  TLep1_Pz      = fHypLepton1.p.Pz();
  TLep1_E       = fHypLepton1.p.E();
  TLep1_Charge  = fHypLepton1.charge;
  TLep2_Px      = fHypLepton2.p.Px();
  TLep2_Py      = fHypLepton2.p.Py();
  TLep2_Pz      = fHypLepton2.p.Pz();
  TLep2_E       = fHypLepton2.p.E();
  TLep2_Charge  = fHypLepton2.charge;

	TWeight_LepEffUp   = EventWeight_LepEffUp;
	TWeight_LepEffDown = EventWeight_LepEffDown;
	TWeight_TrigUp     = EventWeight_TrigUp;
	TWeight_TrigDown   = EventWeight_TrigDown;
	TWeight_FSUp       = EventWeight_FSUp;
	TWeight_FSDown     = EventWeight_FSDown;
	TWeight_PUUp       = EventWeight_PUUp;
	TWeight_PUDown     = EventWeight_PUDown;

  for(int k = 0; k<40; k++){
    if(k<TNJets){
      TJet_Px[k]           = Jet[k].p.Px();
      TJet_Py[k]           = Jet[k].p.Py();
      TJet_Pz[k]           = Jet[k].p.Pz();
      TJet_E[k]            = Jet[k].p.E();
      TJet_isBJet[k]       = Jet[k].isbtag;
    }
    else{
      TJet_Px[k]           = 0;
      TJet_Py[k]           = 0;
      TJet_Pz[k]           = 0;
      TJet_E[k]            = 0;
			TJet_isBJet[k]       = 0;
		}
	}

	// JES, JER
	ResetOriginalObjects();
	gSysSource = JESUp;
	SmearJetPts(1);
	SetEventObjects();
	if( (TIsElMu && IsElMuEvent()) || (TIsDoubleElec &&  IsElElEvent()) || (TIsDoubleMuon && IsMuMuEvent()) ){
		TNJetsJESUp = getNJets();
		TMETJESUp = Get<Float_t>("met_jecUp_pt"); //met
		for(int k = 0; k<40; k++){
			if(k<TNJetsJESUp) TJetJESUp_Pt[k] = JetPt.at(k); 
			else              TJetJESUp_Pt[k] = 0;
		}
  THTJESUp = getHT();
	}

	ResetOriginalObjects();
	gSysSource = JESDown;
	SmearJetPts(2);
	SetEventObjects();
	if( (TIsElMu && IsElMuEvent()) || (TIsDoubleElec &&  IsElElEvent()) || (TIsDoubleMuon && IsMuMuEvent()) ){
		TMETJESDown = Get<Float_t>("met_jecDown_pt");
		TNJetsJESDown = getNJets();
		for(int k = 0; k<40; k++){
			if(k<TNJetsJESDown) TJetJESDown_Pt[k] = JetPt.at(k); 
			else                TJetJESDown_Pt[k] = 0;
		}
  THTJESDown = getHT();
	}

	ResetOriginalObjects();
	SmearJetPts(3);
	gSysSource = JER;
	SetEventObjects();
	if( (TIsElMu && IsElMuEvent()) || (TIsDoubleElec &&  IsElElEvent()) || (TIsDoubleMuon && IsMuMuEvent()) ){
		TNJetsJER = getNJets();
		for(int k = 0; k<40; k++){
			if(k<TNJetsJER) TJetJER_Pt[k] = JetPt.at(k); 
			else            TJetJER_Pt[k] = 0;
		}
	}

	// Btag
	fChargeSwitch = false;

	ResetOriginalObjects();
	gSysSource = BtagUp;
	SetEventObjects();
	if( (TIsElMu && IsElMuEvent()) || (TIsDoubleElec &&  IsElElEvent()) || (TIsDoubleMuon && IsMuMuEvent()) ) TNJetsBtagUp  = getNBTags();

	ResetOriginalObjects();
	gSysSource = BtagDown;
	SetEventObjects();
  if( (TIsElMu && IsElMuEvent()) || (TIsDoubleElec &&  IsElElEvent()) || (TIsDoubleMuon && IsMuMuEvent()) ) TNJetsBtagDown  = getNBTags();

  ResetOriginalObjects();
  gSysSource = MisTagUp;
  SetEventObjects();
  if( (TIsElMu && IsElMuEvent()) || (TIsDoubleElec &&  IsElElEvent()) || (TIsDoubleMuon && IsMuMuEvent()) ) TNJetsBtagMisTagUp  = getNBTags();

  ResetOriginalObjects();
  gSysSource = MisTagDown;
  SetEventObjects();
  TNJetsBtagMisTagDown  = getNBTags();

  ResetOriginalObjects();
  gSysSource = Norm;
  SetEventObjects();
}


//-----------------------------------------------------------------------
// InsideLoop
//-----------------------------------------------------------------------
void StopMiniTrees::InsideLoop() {
  if(gIsT2tt)  if(fabs(fabs((Get<Int_t>("GenSusyMStop"))- gStopMass)) > 1 || fabs((fabs(Get<Int_t>("GenSusyMNeutralino")) - gLspMass)) > 1) return;

	fHDummy->Fill(0.5);
	if (!METFilter()) return;

	// Calculate PU Weight
	if (!gIsData) PUSF = fPUWeight->GetWeight(Get<Float_t>("nTrueInt")); //True       //nTruePU

	// Init data members ........................................................
	GetTreeVariables();
	SetOriginalObjects();
	SetEventObjects();
	if (!gIsData) SelectedGenLepton(); // Get number of generated leptons 

	FillYields(); // Fill yields
}

void StopMiniTrees::Summary(){}

//------------------------------------------------------------------------------
// TRIGGER INFORMATION
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopTrigger#Run2015C_D_25_ns_data_with_RunII
bool StopMiniTrees::PassTriggerMuMu() {
   if(!gIsData) return true;
   //return true;
   Bool_t pass         = false;
   Bool_t passDoubleMu = false;
   Bool_t passSingleMu = false;
   if (gIsData){
    passDoubleMu = (Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"  ) ||
                      Get<Int_t>("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") );
      if (gSampleName == "SingleMuon") passSingleMu = ( Get<Int_t>("HLT_BIT_HLT_IsoTkMu24_v") || 
                                                        Get<Int_t>("HLT_BIT_HLT_IsoMu24_v"  ) );
      if ( (gSampleName == "DoubleMuon" && passDoubleMu                 ) ||
           (gSampleName == "SingleMuon" && passSingleMu && !passDoubleMu) )
         pass = true;
   }else{
      pass = false;
   }
   return pass;
}

bool StopMiniTrees::PassTriggerEE(){
   if(!gIsData) return true;
   //return true;
   Bool_t pass         = false;
   Bool_t passDoubleEl = false;
   Bool_t passSingleEl = false;
   if (gIsData){
      passDoubleEl = Get<Int_t>("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      if (gSampleName == "SingleElectron") passSingleEl = Get<Int_t>("HLT_BIT_HLT_Ele27_WPTight_Gsf_v");
      if ( (gSampleName == "DoubleEG"  && passDoubleEl                 ) ||
           (gSampleName == "SingleElectron" && passSingleEl && !passDoubleEl) )
         pass = true;
   }else{
      pass = false;
   }
   return pass;
}

bool StopMiniTrees::PassTriggerEMu(){
   if(!gIsData) return true;
   //return true;
   Bool_t pass         = false;
   Bool_t passElMu     = false;
   Bool_t passSingleEl = false;
   Bool_t passSingleMu = false;
   if (gIsData){
      passElMu     = (Get<Int_t>("HLT_BIT_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") ||
            Get<Int_t>("HLT_BIT_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") );
      if (gSampleName == "SingleElectron" || gSampleName == "SingleMuon") passSingleEl =  Get<Int_t>("HLT_BIT_HLT_Ele27_WPTight_Gsf_v");
      if (gSampleName == "SingleElectron" || gSampleName == "SingleMuon") passSingleMu = (Get<Int_t>("HLT_BIT_HLT_IsoTkMu24_v") ||
                                                                                   Get<Int_t>("HLT_BIT_HLT_IsoMu24_v") );

      if ( (gSampleName == "MuonEG"         && passElMu                                  ) ||
           (gSampleName == "SingleMuon"     && passSingleMu && !passElMu                 ) ||
           (gSampleName == "SingleElectron" && passSingleEl && !passElMu && !passSingleMu) )
         pass = true;
   }else{
      pass = false;
   }
   return pass;
}

//------------------------------------------------------------------------------
// Get METHODS
//------------------------------------------------------------------------------
float StopMiniTrees::getHT(){
  float ht(0);
  for (unsigned int i=0; i<Jet.size(); i++) ht+=Jet[i].p.Pt();
  return ht;
}
float StopMiniTrees::getJetPtIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return Jet[ind].p.Pt();
}
float StopMiniTrees::getJetEtaIndex(unsigned int ind){
  if (Jet.size() <= ind) return -999.;
  return TMath::Abs(Jet[ind].p.Eta());
}
float StopMiniTrees::getMT(gChannel chan){
  float ptl1 = fHypLepton1.p.Pt();
	float ptl2 = fHypLepton2.p.Pt();
	float dphi = getDelPhill();
	return TMath::Sqrt(2*ptl1*ptl2*(1-TMath::Cos(dphi)));
}

float StopMiniTrees::getDPhiLepJet(){
	if (fHypLepton1.index == -1) return -999.; if (fHypLepton2.index == -1) return -999.;
	// Int_t ij = getLeadingJetbTag(); if (ij < 0) return -999.; 
	if(Jet.size()<1) return -999.;
	TLorentzVector jet = Jet[0].p;
	TLorentzVector plep = fHypLepton1.p;
	if (plep.Pt() < fHypLepton2.p.Pt()) plep = fHypLepton2.p;
	return TMath::Abs(plep.DeltaPhi(jet));
}

float StopMiniTrees::getMinDPhiMetJets(){
	if (fHypLepton1.index == -1) return -999.; if (fHypLepton2.index == -1) return -999.;
	if(Jet.size()<2) return -999.;
	TLorentzVector jet1 = Jet[0].p;
  TLorentzVector jet2 = Jet[1].p;
  TLorentzVector pmet;
  pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
  float MinDelta = TMath::Min(TMath::Abs(pmet.DeltaPhi(jet1)), TMath::Abs(pmet.DeltaPhi(jet2)) );
  return MinDelta;
}

float StopMiniTrees::getDelPhill(){ return fHypLepton1.p.DeltaPhi(fHypLepton2.p);}

float StopMiniTrees::getDPhiJetMet(){
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	return getDPhiClosestJet(pmet);
}

float StopMiniTrees::getDPhiLepMet(){
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	TLorentzVector plep = fHypLepton1.p;
	if (plep.Pt() < fHypLepton2.p.Pt()) plep = fHypLepton2.p;
	return TMath::Abs(plep.DeltaPhi(pmet));
}

float StopMiniTrees::getMT2(TLorentzVector plep1, TLorentzVector plep2, TLorentzVector pmet, float mass){
  double pa[3]; double pb[3]; double pmiss[3];
  pmiss[0] = 0.; // irrelevant
  pmiss[1] = pmet.Px(); pmiss[2] = pmet.Py();
  pa[0] = 0.; pa[1] = plep1.Px(); pa[2] = plep1.Py();
  pb[0] = 0.; pb[1] = plep2.Px(); pb[2] = plep2.Py();
  mt2 MT2bisect;
  MT2bisect.set_momenta(pa, pb, pmiss);
  MT2bisect.set_mn(mass); // testmass
  float MT2 = MT2bisect.get_mt2();
  return MT2;
}

float StopMiniTrees::getMT2ll(gChannel chan){
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
  return getMT2(plep1, plep2, pmet, 0.);
}

float StopMiniTrees::getMT2b(gChannel chan){
  if (getNJets() < 2) return -1;
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
  TLorentzVector pjet0 = Jet[0].p; TLorentzVector pjet1 = Jet[1].p;
  TLorentzVector lv;
  lv.SetPtEtaPhiM(TMath::Sqrt( pow(pmet.Px()+(plep1+plep2).Px(), 2) +  pow(pmet.Py()+(plep1+plep2).Py(), 2)), 0., TMath::ATan2(pmet.Py()+(plep1+plep2).Py(),pmet.Px()+(plep1+plep2).Px()), 0.);
  return getMT2(pjet0, pjet1, lv, 80.398);
}

float StopMiniTrees::getMT2lb(gChannel chan){
  if (getNJets() < 2) return -1;
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
  TLorentzVector pjet0 = Jet[0].p; TLorentzVector pjet1 = Jet[1].p;
  float MT2llbb00, MT2llbb01, MT2llbb;
  float METx = pmet.Px(); float METy = pmet.Py();
  float MET = pmet.Pt(); float MET_phi = pmet.Phi();
  TLorentzVector LepPlusBtagJet00 = plep1 + pjet0;
  TLorentzVector LepPlusBtagJet10 = plep2 + pjet0;
  TLorentzVector LepPlusBtagJet11 = plep2 + pjet1;
  TLorentzVector LepPlusBtagJet01 = plep1 + pjet1;
  if (LepPlusBtagJet11.M()<173 && LepPlusBtagJet00.M()<173 && (LepPlusBtagJet10.M()>173 || LepPlusBtagJet01.M()>173))
    MT2llbb=getMT2(LepPlusBtagJet00, LepPlusBtagJet11, pmet, 0.);
  else if ((LepPlusBtagJet11.M()>173 || LepPlusBtagJet00.M()>173) && LepPlusBtagJet10.M()<173 && LepPlusBtagJet01.M()<173)
    MT2llbb=getMT2(LepPlusBtagJet01, LepPlusBtagJet10, pmet, 0.);
  else if (LepPlusBtagJet11.M()<173 && LepPlusBtagJet00.M()<173 && LepPlusBtagJet10.M()<173 && LepPlusBtagJet01.M()<173) {
    if ( fabs(LepPlusBtagJet11.M()-LepPlusBtagJet00.M()) < fabs(LepPlusBtagJet10.M()-LepPlusBtagJet01.M()) )
      MT2llbb=getMT2(LepPlusBtagJet00, LepPlusBtagJet11, pmet, 0.);
    else
      MT2llbb=getMT2(LepPlusBtagJet01, LepPlusBtagJet10, pmet, 0.);
  }
  else
    MT2llbb=0;
  return MT2llbb;
}

float StopMiniTrees::getMeff(){
	if(Jet.size()<2) return -999.;
	TLorentzVector pmet;
	pmet.SetPtEtaPhiM(getMET(), 0, getMETPhi(), 0);
	return fHypLepton1.p.Pt() + fHypLepton2.p.Pt() + Jet[0].p.Pt() + Jet[1].p.Pt() + pmet.Pt();
}

TLorentzVector StopMiniTrees::getPtllb(){
	TLorentzVector pmet; pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	TLorentzVector pl1 = fHypLepton1.p; TLorentzVector pl2 = fHypLepton2.p;
	return pl1 + pl2 + pmet;
}

float StopMiniTrees::getDPhibMet(){
	TLorentzVector pmet; pmet.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	TLorentzVector Ptllb = getPtllb();
	return pmet.DeltaPhi(Ptllb); 
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#JER_Scaling_factors_and_Uncertai
float StopMiniTrees::getJERScale(int jet){
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
float StopMiniTrees::getJERScaleUp(int jet){
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
float StopMiniTrees::getJERScaleDown(int jet){
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

float StopMiniTrees::getErrPt(float Pt, float Eta) {
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

float StopMiniTrees::getLeptonError(gChannel chan){
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

float StopMiniTrees::getTriggerError(gChannel chan){
	float trig(0.);
	if (chan==Muon) trig = fLeptonSF->GetDoubleMuSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	if (chan==Elec) trig = fLeptonSF->GetDoubleElSF_err(fHypLepton1.p.Eta(),fHypLepton2.p.Eta());
	if (chan==ElMu){
		float leadingPt  = fHypLepton1.p.Pt() > fHypLepton2.p.Pt()? fHypLepton1.p.Pt() : fHypLepton2.p.Pt();
		float leadingEta = fHypLepton1.p.Pt() > fHypLepton2.p.Pt()? fHypLepton1.p.Eta() : fHypLepton2.p.Eta();
		trig = fLeptonSF->GetMuEGSF_err    (leadingPt,leadingEta);
	}
	return trig;
}

float StopMiniTrees::getFSError(gChannel chan){
	if(!gIsT2tt) return 0;
  float eFS1 = 1; float eFS2 = 1;
	if (chan == Muon){
    eFS1  = fLeptonSF->GetFastSimMuonSF_err(fHypLepton1.p.Pt(), TMath::Abs(fHypLepton1.p.Eta()));
    eFS2  = fLeptonSF->GetFastSimMuonSF_err(fHypLepton2.p.Pt(), TMath::Abs(fHypLepton2.p.Eta()));
  } 
  else if (chan == Elec){
    eFS1  = fLeptonSF->GetFastSimElectronSF_err(fHypLepton1.p.Pt(), TMath::Abs(fHypLepton1.p.Eta()));
    eFS2  = fLeptonSF->GetFastSimElectronSF_err(fHypLepton2.p.Pt(), TMath::Abs(fHypLepton2.p.Eta()));
  }
  else if (chan == ElMu){
    eFS1  = fLeptonSF->GetFastSimMuonSF_err(fHypLepton1.p.Pt(),TMath::Abs(fHypLepton1.p.Eta()));
    eFS2  = fLeptonSF->GetFastSimElectronSF_err(fHypLepton2.p.Pt(), TMath::Abs(fHypLepton2.p.Eta()));
	}
  return TMath::Sqrt(eFS1*eFS1 + eFS2*eFS2); 
}

float StopMiniTrees::getSF(gChannel chan) {
	if (gIsData)              return 1.; //Don't scale data
	float id1(1.),id2(1.), trig(1.);
	float err1(0.), err2(0.), err_trg(0.);
  float SF = 1; float FSSF = 1;
  float FS1 = 1; float FS2 = 1;
	if (chan == Muon){
		id1  = fLeptonSF->GetTightMuonSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta());
		id2  = fLeptonSF->GetTightMuonSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
		trig = fLeptonSF->GetDoubleMuSF (fHypLepton1.p.Pt(),fHypLepton1.p.Eta());
    FS1  = fLeptonSF->GetFastSimMuonSF(fHypLepton1.p.Pt(), TMath::Abs(fHypLepton1.p.Eta()));
    FS2  = fLeptonSF->GetFastSimMuonSF(fHypLepton2.p.Pt(), TMath::Abs(fHypLepton2.p.Eta()));
  } 
  else if (chan == Elec){
    id1  = fLeptonSF->GetTightElectronSF(fHypLepton1.p.Pt(), fHypLepton1.p.Eta()); 
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta()); 
    trig = fLeptonSF->GetDoubleElSF     (fHypLepton1.p.Pt(),fHypLepton1.p.Eta()); 
    FS1  = fLeptonSF->GetFastSimElectronSF(fHypLepton1.p.Pt(), TMath::Abs(fHypLepton1.p.Eta()));
    FS2  = fLeptonSF->GetFastSimElectronSF(fHypLepton2.p.Pt(), TMath::Abs(fHypLepton2.p.Eta()));
  }
  else if (chan == ElMu){
    float leadingPt  = fHypLepton1.p.Pt() > fHypLepton2.p.Pt()? fHypLepton1.p.Pt() : fHypLepton2.p.Pt();
    float leadingEta = fHypLepton1.p.Pt() > fHypLepton2.p.Pt()? fHypLepton1.p.Eta() : fHypLepton2.p.Eta();
    id1  = fLeptonSF->GetTightMuonSF    (fHypLepton1.p.Pt(), fHypLepton1.p.Eta()); 
    id2  = fLeptonSF->GetTightElectronSF(fHypLepton2.p.Pt(), fHypLepton2.p.Eta());
    trig = fLeptonSF->GetMuEGSF         (leadingPt,leadingEta);
    FS1  = fLeptonSF->GetFastSimMuonSF(fHypLepton1.p.Pt(),TMath::Abs(fHypLepton1.p.Eta()));
    FS2  = fLeptonSF->GetFastSimElectronSF(fHypLepton2.p.Pt(), TMath::Abs(fHypLepton2.p.Eta()));
	}
  SF = PUSF*id1*id2*trig;
  FSSF = FS1*FS2;
	if(gIsT2tt) SF*= FSSF;
  //cout << "FS1  FS2 = " << FS1 << "  " << FS2 << endl;
  return (SF);
}

float StopMiniTrees::getISRJetsWeight(Int_t nISRJet){
  Float_t SF = 0;
  if     (nISRJet == 0) SF = 1;
  else if(nISRJet == 1) SF = 0.882;
  else if(nISRJet == 2) SF = 0.792;
  else if(nISRJet == 3) SF = 0.702;
  else if(nISRJet == 4) SF = 0.648;
  else if(nISRJet == 5) SF = 0.601;
  else                  SF = 0.515;
  return SF;
}

//--------------------------------------------------------------------------
// Fill histograms      
//------------------------------------------------------------------------
void StopMiniTrees::FillYields(gSystFlag sys){
	ResetHypLeptons();  

	if (PassTriggerEMu()  && IsElMuEvent()){
		// Define Hypothesis Leptons...
		EventWeight = gWeight * getSF(ElMu);

		EventWeight_LepEffUp   = gWeight * (getSF(ElMu) + getLeptonError(ElMu));
		EventWeight_LepEffDown = gWeight * (getSF(ElMu) - getLeptonError(ElMu));
		EventWeight_FSUp       = gWeight * (getSF(ElMu) + getFSError(ElMu));
		EventWeight_FSDown     = gWeight * (getSF(ElMu) - getFSError(ElMu));
		EventWeight_TrigUp     = gWeight * (getSF(ElMu) + getTriggerError(ElMu));
		EventWeight_TrigDown   = gWeight * (getSF(ElMu) - getTriggerError(ElMu));
		hWeight -> Fill(EventWeight,1.);

		if (!gIsData){
			PUSF = fPUWeightUp->GetWeight(Get<Float_t>("nTrueInt"));
			EventWeight_PUUp       = gWeight * getSF(ElMu);
			PUSF = fPUWeightDown->GetWeight(Get<Float_t>("nTrueInt"));
			EventWeight_PUDown     = gWeight * getSF(ElMu);
			PUSF = fPUWeight->GetWeight(Get<Float_t>("nTrueInt")); 
		}

		if(gIsMCatNLO){
			EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); // *(1.-2.*0.115)); 

			EventWeight_LepEffUp   *= genWeight; 
			EventWeight_LepEffDown *= genWeight; 
			EventWeight_FSUp       *= genWeight; 
			EventWeight_FSDown     *= genWeight; 
			EventWeight_TrigUp     *= genWeight; 
			EventWeight_TrigDown   *= genWeight; 
			EventWeight_PUUp       *= genWeight; 
			EventWeight_PUDown     *= genWeight; 
		}
		if((gCreateTree) && (sys==Norm) && !(fChargeSwitch) && PassesMllVeto()){
			SetTreeVariables(ElMu);
			if( (TNJets >= 2) || (TNJetsJESUp >= 2) || (TNJetsJESDown >= 2) || (TNJetsJER >= 2))  fTree->Fill();
		}

	}
	ResetHypLeptons(); 

	if (PassTriggerMuMu() && IsMuMuEvent()){
		EventWeight = gWeight * getSF(Muon); 
		//EventWeight = 1.;
		EventWeight_LepEffUp   = gWeight * (getSF(Muon) + getLeptonError(Muon));
		EventWeight_LepEffDown = gWeight * (getSF(Muon) - getLeptonError(Muon));
		EventWeight_FSUp       = gWeight * (getSF(Muon) + getFSError(Muon));
		EventWeight_FSDown     = gWeight * (getSF(Muon) - getFSError(Muon));
		EventWeight_TrigUp     = gWeight * (getSF(Muon) + getTriggerError(Muon));
		EventWeight_TrigDown   = gWeight * (getSF(Muon) - getTriggerError(Muon));

		if (!gIsData){
			PUSF = fPUWeightUp->GetWeight(Get<Float_t>("nTrueInt"));
			EventWeight_PUUp       = gWeight * getSF(Muon);
			PUSF = fPUWeightDown->GetWeight(Get<Float_t>("nTrueInt"));
			EventWeight_PUDown     = gWeight * getSF(Muon);
			PUSF = fPUWeight->GetWeight(Get<Float_t>("nTrueInt")); 
		}

		if(gIsMCatNLO){
			EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115)); 

			EventWeight_LepEffUp   *= genWeight; 
			EventWeight_LepEffDown *= genWeight; 
			EventWeight_FSUp       *= genWeight; 
			EventWeight_FSDown     *= genWeight; 
			EventWeight_TrigUp     *= genWeight; 
			EventWeight_TrigDown   *= genWeight; 
			EventWeight_PUUp       *= genWeight; 
			EventWeight_PUDown     *= genWeight; 
		}
		if( (gCreateTree) &&	(sys==Norm) && !(fChargeSwitch) && PassesZVeto()      && PassesMllVeto()){
		SetTreeVariables(Muon);
			if((TNJets >= 2) || (TNJetsJESUp >= 2) || (TNJetsJESDown >= 2) || (TNJetsJER >= 2) ) 	fTree->Fill();
		}
  }

  ResetHypLeptons(); 
  if (PassTriggerEE() && IsElElEvent()){
		EventWeight = gWeight * getSF(Elec);
		EventWeight_LepEffUp   = gWeight * (getSF(Elec) + getLeptonError(Elec));
		EventWeight_LepEffDown = gWeight * (getSF(Elec) - getLeptonError(Elec));
		EventWeight_FSUp       = gWeight * (getSF(Elec) + getFSError(Elec));
		EventWeight_FSDown     = gWeight * (getSF(Elec) - getFSError(Elec));
		EventWeight_TrigUp     = gWeight * (getSF(Elec) + getTriggerError(Elec));
		EventWeight_TrigDown   = gWeight * (getSF(Elec) - getTriggerError(Elec));

		if (!gIsData){
			PUSF = fPUWeightUp->GetWeight(Get<Float_t>("nTrueInt"));
			EventWeight_PUUp       = gWeight * getSF(Elec);
			PUSF = fPUWeightDown->GetWeight(Get<Float_t>("nTrueInt"));
			EventWeight_PUDown     = gWeight * getSF(Elec);
			PUSF = fPUWeight->GetWeight(Get<Float_t>("nTrueInt")); 
		}

		if(gIsMCatNLO){
			EventWeight = EventWeight * genWeight;// /(TMath::Abs(T_Event_weight)); //*(1.-2.*0.115)); 
			EventWeight_LepEffUp   *= genWeight; 
			EventWeight_LepEffDown *= genWeight; 
			EventWeight_FSUp       *= genWeight; 
			EventWeight_FSDown     *= genWeight; 
			EventWeight_PUUp       *= genWeight; 
			EventWeight_PUDown     *= genWeight; 
			EventWeight_TrigUp     *= genWeight; 
			EventWeight_TrigDown   *= genWeight; 
		}
		if( (gCreateTree) &&	(sys==Norm)        && !(fChargeSwitch) && PassesZVeto()      && PassesMllVeto()){
			SetTreeVariables(Elec);
			if((TNJets >= 2) || (TNJetsJESUp >= 2) || (TNJetsJESDown >= 2) || (TNJetsJER >= 2))  fTree->Fill();
		}

  }
  ResetHypLeptons();
}

//----------------------------------------------------------------------
// Passes
//----------------------------------------------------------------------
bool StopMiniTrees::PassesMllVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
	if (InvMass < 20.)            return false; 
	return true;
}

bool StopMiniTrees::PassesZVeto(){
	// Check consistency.
	if (fHypLepton1.index == -1) return false;
	if (fHypLepton2.index == -1) return false;
	float InvMass = (fHypLepton1.p+fHypLepton2.p).M();
  if (InvMass > 76. && InvMass < 106.) return false;
  return true;
}

bool StopMiniTrees::PassesNJetsCut(){
	if (getNJets() <= 1) return false;
	return true;
}

bool StopMiniTrees::PassesMETCut(){
	if (getMET() < 50.) return false;
	return true;
}

bool StopMiniTrees::PassesNBtagCut(){
	if (getNBTags() < 1) return false;
	return true;
}

bool StopMiniTrees::PassesDYVetoCut(){
  if(getMinDPhiMetJets() < 0.25) return false;
  if (getHT() == 0) return false;
  if(getMET()/TMath::Sqrt(getHT()) < 5.0) return false;
  return true;
}

bool StopMiniTrees::IsElMuEvent(){
	if (fChargeSwitch){      return (IsDileptonEvent()  == 3);   }
	return (IsDileptonEvent() == -3);
}

bool StopMiniTrees::IsMuMuEvent(){
	if (fChargeSwitch){  return (IsDileptonEvent()  == 1); }
	return (IsDileptonEvent() == -1);
}

bool StopMiniTrees::IsElElEvent(){
  if (fChargeSwitch){    return (IsDileptonEvent()  == 2); }
  return (IsDileptonEvent() == -2);
}

int StopMiniTrees::IsDileptonEvent(){
	if(Lepton.size() != 2) return 0; // third 3rd lepton veto
  if(Lepton[0].p.Pt()<25) return 0;
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
	return result;
}

//------------------------------------------------------------------------------
// LEPTON SELECTORS
//------------------------------------------------------------------------------
bool momentumComparator(lepton i, lepton j){ return (i.p.Pt()>j.p.Pt()); }

vector<lepton> StopMiniTrees::SortLeptonsByPt(vector<lepton>& leptons){
  vector<lepton> theLep = leptons;
  sort (theLep.begin(), theLep.end(), momentumComparator);
  return theLep;
}

int StopMiniTrees::getSelectedLeptons(){
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
		if(IsTightMuon(i)){ thetype = 0; nMuon ++; }
    else if(IsTightElectron(i)){ thetype = 1; nElec++; }
    else  continue;
 
    lep.SetPxPyPzE(LepGood_px[i], LepGood_py[i], LepGood_pz[i], LepGood_energy[i]); 
    lepton tmpLepton(lep, LepGood_charge[i], thetype, i); 
    tmp_lepton.push_back(tmpLepton);
  }
  Lepton = SortLeptonsByPt(tmp_lepton);
  return Lepton.size();
}

bool StopMiniTrees::METFilter(){
  if(gIsT2tt) return true;
  return true;
  if (Get<Int_t>("Flag_HBHENoiseFilter") && 
      Get<Int_t>("Flag_HBHENoiseIsoFilter") && 
      Get<Int_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") && 
      Get<Int_t>("Flag_goodVertices") && 
			Get<Int_t>("Flag_eeBadScFilter") &&
			Get<Int_t>("Flag_badMuonFilter") && 
			Get<Int_t>("Flag_badChargedHadronFilter") &&  
			Get<Int_t>("Flag_globalTightHalo2016Filter")
		 ){
    return true;
  }
  else return false;
}

//------------------------------------------------------------------------------
// Muon Selectors  
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/view/CMS/TopMUO
// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
bool StopMiniTrees::IsTightMuon(unsigned int iMuon,float ptcut){
   if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return false;
   if (LepGood_pt[iMuon]              < ptcut) return false;
   if (TMath::Abs(LepGood_eta[iMuon]) > 2.4)   return false;
   //if (LepGood_relIso04[iMuon]        > 0.15)  return false; 
	 if (!getMultiIso(iMuon)) return false;
   if (TMath::Abs(LepGood_dxy[iMuon]) >= 0.05)  return false; 
   if (TMath::Abs(LepGood_dz[iMuon])  >= 0.1)  return false;
   if (Get<Float_t>("LepGood_sip3d", iMuon) > 4.0) return false; 
   //if (Get<Int_t>("LepGood_tightId",iMuon)           < 1   )  return false;
   if (Get<Int_t>("LepGood_mediumMuonId",iMuon)    < 1   )  return false; // medium ID?
   return true;
}

float StopMiniTrees::getMuonIso(int iMuon){
	if (iMuon < 0) return 9999.;
	if ((TMath::Abs(LepGood_pdgId[iMuon])) != 13) return 9999.;
	return LepGood_relIso04[iMuon];
}

//------------------------------------------------------------------------------
// Electron Selectors
//------------------------------------------------------------------------------
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
bool StopMiniTrees::IsTightElectron(unsigned int iElec, float ptcut){
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return false; 
	if (LepGood_pt[iElec] < ptcut) return false;
	if (TMath::Abs(LepGood_eta[iElec]) > 2.4) return false;
  if (!getMultiIso(iElec)) return false;
  //if (getElecIso(iElec) > 0.0766) return false;
	if (TMath::Abs(LepGood_dxy[iElec]) >= 0.05) return false; 
	if (TMath::Abs(LepGood_dz[ iElec]) >= 0.1 ) return false;
	if (TMath::Abs(LepGood_etaSc[iElec]) > 1.4442 && 
			TMath::Abs(LepGood_etaSc[iElec]) < 1.566) return false;
	if (Get<Float_t>("LepGood_sip3d", iElec) > 4.0) return false; 
  // Tight Electron Id:
  if (Get<Int_t>("LepGood_tightId",iElec) < 3) return false;
  /* We use the string tightId because it contains the right definition of LepGood_dEtaScTrkIn...
  if(TMath::Abs(Get<Float_t>("LepGood_etaSc", iElec)) < 1.479){ // central
    if(TMath::Abs(Get<Float_t>("LepGood_sigmaIEtaIEta", iElec)) > 0.00998) return false;
    if(TMath::Abs(Get<Float_t>("LepGood_dEtaScTrkIn", iElec)) > 0.00308  ) return false; 
    if(TMath::Abs(Get<Float_t>("LepGood_dPhiScTrkIn", iElec)) > 0.0816) return false;
    if(Get<Float_t>("LepGood_hadronicOverEm", iElec) > 0.0414) return false;
    if(TMath::Abs(Get<Float_t>("LepGood_eInvMinusPInv", iElec)) > 0.0129) return false;
	}
  else{ // forward
    if(TMath::Abs(Get<Float_t>("LepGood_sigmaIEtaIEta", iElec)) > 0.0292) return false;
    if(TMath::Abs(Get<Float_t>("LepGood_dEtaScTrkIn", iElec)) > 0.00605) return false; 
    if(TMath::Abs(Get<Float_t>("LepGood_dPhiScTrkIn", iElec)) > 0.0394) return false;
    if(Get<Float_t>("LepGood_hadronicOverEm", iElec) > 0.0641) return false;
    if(TMath::Abs(Get<Float_t>("LepGood_eInvMinusPInv", iElec)) >  	0.0129) return false;
  }*/
	if (Get<Int_t>("LepGood_lostHits", iElec)         > 1      ) return false; 
	if (Get<Int_t>("LepGood_convVeto", iElec)        == 0      ) return false; 

	return true;
}

float StopMiniTrees::getElecIso(int iElec){
	if (iElec < 0) return 9999.;
	if ((TMath::Abs(LepGood_pdgId[iElec])) != 11) return 9999.;
	return LepGood_relIso03[iElec];
}

bool StopMiniTrees::getMultiIso(unsigned int index){
  //https://www.dropbox.com/s/fsfw0gummwsc61v/lepawareJECv2_bkg_wp_300915.pdf?dl=0
  float max_mRelIso = 0.09;  
  float min_ptRatio = 0.84; // Very tight
  float min_ptRel   = 7.2 ; // Tight
	//if      ((TMath::Abs(LepGood_pdgId[index])) == 11) max_mRelIso = 0.13;  // Tight for electrons
	//else if ((TMath::Abs(LepGood_pdgId[index])) == 13) max_mRelIso = 0.2;   // Medium for muons
  float mRelIso = Get<Float_t>("LepGood_miniRelIso", index);
  float ptRatio = Get<Float_t>("LepGood_jetPtRatiov2", index);
  float ptRel   = Get<Float_t>("LepGood_jetPtRelv2", index);
  // min_ptRatio = 0; min_ptRel = 0; // No MultiIso
  return (mRelIso < max_mRelIso && (ptRatio > min_ptRatio || ptRel > min_ptRel));
}

void StopMiniTrees::setMET(float newmet){ MET = newmet;}
float StopMiniTrees::getMET(){ return MET; }
float StopMiniTrees::getMETPhi(){ return MET_Phi;}
int StopMiniTrees::getNJets(){ return nJets;}

float StopMiniTrees::getDPhiClosestJet(TLorentzVector lep){
	float minDphi = 9999.;
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (minDphi > TMath::Abs(lep.DeltaPhi(Jet[i].p))) minDphi = TMath::Abs(lep.DeltaPhi(Jet[i].p));
	}
	return minDphi;
}

int StopMiniTrees::getLeadingJetbTag(){
	for (unsigned int i=0; i<Jet.size(); i++) {
		if (Jet[i].isbtag) return i;
	}
	return  -1;
}

int StopMiniTrees::getNBTags(){
	int ntags(0);
	for(UInt_t i = 0; i <Jet.size(); i++){
		if (Jet[i].isbtag) ntags++;
	}
	return ntags;
}

float StopMiniTrees::getDeltaPhillJet(){
	if (fHypLepton1.index == -1) return -999.;
	if (fHypLepton2.index == -1) return -999.;
	Int_t ij = getLeadingJetbTag();
	if (ij < 0) return -999.; 
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	TLorentzVector jet = Jet[ij].p; 
	return TMath::Abs(dilep.DeltaPhi(jet));
}

float StopMiniTrees::getTopD(){
	if (fHypLepton1.index == -1) return -999;
	if (fHypLepton2.index == -1) return -999;
	// Make Dilepton system
	TLorentzVector dilep = fHypLepton1.p+fHypLepton2.p;
	Float_t DeltaPhi(0.),TopD(0.);
	TLorentzVector jet;
	if (nJets == 0) return -999.;
	TopD     = 1 - (DeltaPhi/TMath::Pi()) * (1 - Jet_btagCSV[(Jet[0].index)]);
	return TopD;
}

int StopMiniTrees::getSelectedJets(){
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
    if      (gSysSource == BtagUp)      isbtag = fBTagSFbUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai); 
    else if (gSysSource == BtagDown)    isbtag = fBTagSFbDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
    else if (gSysSource == MisTagUp)    isbtag = fBTagSFlUp->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
    else if (gSysSource == MisTagDown)  isbtag = fBTagSFlDo->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
    else                                isbtag = fBTagSFnom->IsTagged(jetbtagi, jetmcflavouri, JetPt.at(i), jetetai);
    }
    jet tmpjet(jt, isbtag, i);
    Jet.push_back(tmpjet);
    nj++;
  }
  return nj;
}

bool StopMiniTrees::IsGoodJet(unsigned int ijet, float ptcut){
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
  for(unsigned int i = 0; i < Lepton.size(); i++){
    if(jet.DeltaR(Lepton[i].p) < minDR) return false;
  }
  return true;
}

//------------------------------------------------------------------------------
// SelectedGenLepton
//------------------------------------------------------------------------------
void StopMiniTrees::SelectedGenLepton() {
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
	}
}

void StopMiniTrees::propagateMET(TLorentzVector nVec, TLorentzVector oVec){
	TLorentzVector met;
	met.SetPtEtaPhiM(getMET(), 0., getMETPhi(), 0.);
	// set the pfMET to the old MET minus original vector plus new vector 
	setMET( (met+oVec-nVec).Pt() );
}

std::vector<int> StopMiniTrees::CleanedJetIndices(float pt){
	std::vector<int> cleanJetsInd;
	for(Int_t i = 0; i <nJet; i++){
		if (IsGoodJet(i,pt)) cleanJetsInd.push_back(i);
	}
	return cleanJetsInd;
}

void StopMiniTrees::SmearJetPts(int flag){
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
			if(flag == 1) JetPt.at(*it) = Get<Float_t>("Jet_rawPt", *it)*Get<Float_t>("Jet_corr_JECUp",*it);   // vary up   for flag 1 
			if(flag == 2) JetPt.at(*it) = Get<Float_t>("Jet_rawPt", *it)*Get<Float_t>("Jet_corr_JECDown",*it); // vary down for flag 2;
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

void StopMiniTrees::ScaleLeptons(int flag){
	// Shift the lepton pts for systematics studies
	if(gIsData) return; // don't smear data
	if(flag == 0) return;
	//float scale = 0.003; // 0.3% for muons
	float scale = 0.005; 
	TLorentzVector oleps, leps, tmp;
	for(Int_t k = 0; k < nMuon; ++k){ //xxx Muon
		if(TMath::Abs(LepGood_pdgId[k]) != 13) continue; 
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

void StopMiniTrees::ScaleMET(int flag){
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

