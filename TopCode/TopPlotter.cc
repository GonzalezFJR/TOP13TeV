#include "TopPlotter.h"

#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLatex.h"

using namespace std;

bool gUseTTMadSpin = true;
const float fLumiNorm = 2173*1.023; //2222.979

TopPlotter::TopPlotter(){}

void TopPlotter::Init(TString pathtofiles){
  cout << "ResetDataMembers" << endl;
  ResetDataMembers();
  //  ttbar_TLWG = 252.89;
  //ttbar_TLWG = 245.8;
  ttbar_TLWG = 831.8;
  
  LoadSamples(pathtofiles);

  DoDF = false;  if (pathtofiles.Contains("DF")) DoDF = true;
  DoSF = false;  if (pathtofiles.Contains("SF")) DoSF = true; 
}
void TopPlotter::Loop(){
  //  cout << "Printing Yields with MC estimations..." << endl;
  //  tA->PrintYieldsWithMC();

  cout << "Calculate Data-Driven backgrounds" << endl;
  CalculateNonWZLeptonsBkg();
  CalculateDYBkg();

  //  cout << "Draw Kinematic Plots with MC/DD estimations..." << endl;
  // tA->DrawKinematicPlotsWithMC(-1, NBTagsNJets, -1);
  //  DrawKinematicPlots(true, ElMu, NBtagJets, i2jets);
  //DrawKinematicPlots(true, ElMu);
  //DrawKinematicPlots(false, ElMu);
  
  cout << "Draw Plots for Likelihood..." << endl;
  //DrawNbjetsNjets(false);
//  //DrawNbjetsNjets(true);
//  
  cout << "Saving Plots for Likelihood..." << endl;
  //SaveHistosForLH(false);
//  //SaveHistosForLH(true);
//  
  cout << "Calculating cross section ... " << endl;
  CalculateCrossSection(false);
  CalculateCrossSection(true );
  
}
void TopPlotter::LoadSamples(TString pathtofiles){
  TFile *_file ;
  TH1F *hOSyields;
  TH1F *hOSyields_sys;
  TH1F *hSSyields;
  
  cout << "Loading Samples.... " << endl;
  Float_t Weight = fLumiNorm; 
  for (size_t sample=0; sample<gNSAMPLES; sample++){
    TString samplename = pathtofiles + "/Tree_13TeV_EA_"+SampleName[sample]+".root";
    cout << "Loading " + samplename +" ..." << endl;
    _file = new TFile(samplename);
    //Bool_t IsData = (sample == DoubleEG || sample == DoubleMuon || sample == MuEG || sample ==SingleMuon || sample ==SingleElectron );
    Bool_t IsData = (sample == DoubleEG || sample == DoubleMuon || sample == MuEG );
    
    S[sample].name = SampleName[sample];
    
    // CALCULATE WEIGHT
    if (sample==TTJets_MadSpin    || 
	//sample==TTJets_matchingup || sample==TTJets_matchingdown ||
	sample==TTJets_scaleup    || sample==TTJets_scaledown        
        )   
      Weight = fLumiNorm ;//* (0.108*9.)*(0.108*9.);   
    else 
      Weight = fLumiNorm; 
    
    if(sample==TTJets_MadSpin)    toppt_weight = ((TH1F*) _file->Get("H_TopPtWeight"))->GetMean();
    //if(sample==TTJets_MadSpinPDF) {
    //  S[sample].pdfWeights    = (TH1F*)_file->Get("H_pdfWeight");
    //  S[sample].pdfWeightsSum = (TH1F*)_file->Get("H_pdfWeightSum");
    //}
    // Load numbers... 
    for (size_t ch=0; ch<gNCHANNELS; ch++){
      hOSyields = (TH1F*) _file->Get("H_Yields_"  +gChanLabel[ch]);
      hSSyields = (TH1F*) _file->Get("H_SSYields_"+gChanLabel[ch]);
      
      for (size_t cut=0; cut<iNCUTS; cut++){
	if (IsData){
	  S[sample].Yields	 [ch][cut] = hOSyields->GetBinContent(cut+1);
	  S[sample].Yields_stat  [ch][cut] = hOSyields->GetBinError(cut+1);
	  S[sample].SSYields	 [ch][cut] = hSSyields->GetBinContent(cut+1);
	  S[sample].SSYields_stat[ch][cut] = hSSyields->GetBinError(cut+1);
	}
	else {  
	  S[sample].Yields       [ch][cut] = hOSyields->GetBinContent(cut+1) * Weight;
	  S[sample].Yields_stat  [ch][cut] = hOSyields->GetBinError(cut+1)   * Weight;
	  S[sample].SSYields     [ch][cut] = hSSyields->GetBinContent(cut+1) * Weight;
	  S[sample].SSYields_stat[ch][cut] = hSSyields->GetBinError(cut+1)   * Weight;
	}
      }
    }      
    
    //if (sample==TTJets_MadSpinPDF   || sample==TTbar_Powheg            ||
	//sample==TTbar_Powheg_Herwig || sample==TTJetsFullLeptMGTuneP11 ||
	//sample==TTJetsFullLeptMGTuneP11noCR) continue;

    // Load Systematics (ONLY FOR MC...) 
    for (size_t ch=0; ch<gNCHANNELS; ch++){                                                            
     if (!IsData){
	for (size_t sys=1; sys<gNSYST; sys++){
	  hOSyields_sys = (TH1F*) _file->Get("H_Yields_"  +gChanLabel[ch]+"_"+SystName[sys]);
	  for (size_t cut=0; cut<iNCUTS; cut++)
	    S[sample].Yields_syst[ch][cut][sys] = hOSyields_sys->GetBinContent(cut+1) * Weight;
	}
	
	
	S[sample].SystError[0][SFIDISO] = 0.0318;                      // mu-mu   //((TH1F*)_file->Get("H_LepSys_"+ gChanLabel[ch]+"_dilepton"))->GetMean();
	S[sample].SystError[1][SFIDISO] = 0.0462;                      // e-e	      //((TH1F*)_file->Get("H_LepSys_"+ gChanLabel[ch]+"_dilepton"))->GetMean();
	S[sample].SystError[2][SFIDISO] = sqrt(0.0227*0.0227 + 0.0155*0.0155); // e-mu        //((TH1F*)_file->Get("H_LepSys_"+ gChanLabel[ch]+"_dilepton"))->GetMean();
	
	
	
	S[sample].SystError[ch][SFTrig] = 0.044;  //((TH1F*)_file->Get("H_TrigSys_"+ gChanLabel[ch]+"_dilepton"))->GetMean();
	S[sample].SystError[0][SFTrig]  = 0.0108; // mu-mu      
	S[sample].SystError[1][SFTrig]  = 0.0125; // e-e	     
	S[sample].SystError[2][SFTrig]  = 0.0123; // e-mu       
  
      }
    }
     
    // Load kinematic histograms of the samples. 
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t cut=0; cut<iNCUTS; cut++){
	if (SampleName[sample] == "DoubleMuonSum"         ||       
	    SampleName[sample] == "DoubleEGsum"           || 
	    SampleName[sample] == "MuonEGsum"             ||       
	    SampleName[sample] == "DYJetsToLL_M10to50_aMCatNLO_ext" ||
	    SampleName[sample] == "DYJetsToLL_M50_aMCatNLO")  {
	    //SampleName[sample] == "ZJets_Madgraph")  {
	  
	  /// LOAD DY DD estimation:
	  TString histoname = "H_DY_InvMass_" + gChanLabel[chan] +  "_" + sCut[cut];
 	  //S[sample].MllHistos[chan][cut] = GetHisto1D(_file, histoname);
	  S[sample].MllHistos[chan][cut] = GetHisto1D(_file, histoname);
 	  if (!IsData)  S[sample].MllHistos[chan][cut]->Scale(Weight);
	}

	TString histoname = "";
	for (size_t var=0; var<gNVARS; var++){
	  histoname = "H_" + KinVarName[var] + "_" + gChanLabel[chan] + "_" + sCut[cut];
	  S[sample].KinHistos[chan][cut][var] = GetHisto1D(_file, histoname);
	  if (!IsData)                S[sample].KinHistos[chan][cut][var]->Scale(Weight);
	}

	// SYSTEMATIC UNCERTAINTIES...
	for (size_t sys=0; sys<gNSYST; sys++){
	  if (sys==0) histoname = "H_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "H_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  // cout << "Reading... " << histoname << " " << endl;
	  S[sample].NBtagsNJets[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].NBtagsNJets[chan][cut][sys]->Scale(Weight);
	  
	  if (sys==0) histoname = "HSS_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "HSS_NBtagsNJets_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  
	  S[sample].SSNBtagsNJets[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].SSNBtagsNJets[chan][cut][sys]->Scale(Weight);
	  
	  
	  if (sys==0) histoname = "H_InvMass_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "H_InvMass_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  //	  cout << "Reading... " << histoname << endl;
	  S[sample].InvMass[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].InvMass[chan][cut][sys]->Scale(Weight);
	  
	  if (sys==0) histoname = "HSS_InvMass_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "HSS_InvMass_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  
	  S[sample].SSInvMass[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].SSInvMass[chan][cut][sys]->Scale(Weight);
	  
	  
	  if (sys==0) histoname = "H_AbsDelPhiLeps_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "H_AbsDelPhiLeps_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  //	  cout << "Reading... " << histoname << endl;
	  S[sample].AbsDelPhiLeps[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].AbsDelPhiLeps[chan][cut][sys]->Scale(Weight);
	  
	  if (sys==0) histoname = "HSS_AbsDelPhiLeps_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "HSS_AbsDelPhiLeps_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  
	  S[sample].SSAbsDelPhiLeps[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].SSAbsDelPhiLeps[chan][cut][sys]->Scale(Weight);
	  
	  
	  if (sys==0) histoname = "H_delPhi2LeadJets_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "H_delPhi2LeadJets_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  //	  cout << "Reading... " << histoname << endl;
	  S[sample].delPhi2LeadJets[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].delPhi2LeadJets[chan][cut][sys]->Scale(Weight);
	  
	  if (sys==0) histoname = "HSS_delPhi2LeadJets_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "HSS_delPhi2LeadJets_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  
	  S[sample].SSdelPhi2LeadJets[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].SSdelPhi2LeadJets[chan][cut][sys]->Scale(Weight);
	  
	  
	  if (sys==0) histoname = "H_minDelRJetsLeps_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "H_minDelRJetsLeps_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  //	  cout << "Reading... " << histoname << endl;
	  S[sample].minDelRJetsLeps[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].minDelRJetsLeps[chan][cut][sys]->Scale(Weight);
	  
	  if (sys==0) histoname = "HSS_minDelRJetsLeps_" + gChanLabel[chan] + "_" + sCut[cut];
	  else        histoname = "HSS_minDelRJetsLeps_" + gChanLabel[chan] + "_" + sCut[cut] + "_" + SystName[sys];
	  
	  S[sample].SSminDelRJetsLeps[chan][cut][sys] = GetHisto1D(_file, histoname);
	  if (!IsData) S[sample].SSminDelRJetsLeps[chan][cut][sys]->Scale(Weight);
	}
      }
    }
    _file->Close();
  }
  LoadCategories();
}
/*
void TopPlotter::LoadCategory(Categories &Cat, TString sample){
  
  // CREATE LIST WITH SELECTED SAMPLES...
  std::vector<Samples> SampleList;
  if (sample=="ttbar"){
    SampleList.push_back(TTJets_MadSpin);
  }
  else if (sample=="stop"){  
    SampleList.push_back(TbarWDilep);
    SampleList.push_back(TWDilep);
  }
  else if (sample=="dy"){
    SampleList.push_back(ZJets_Madgraph);
    SampleList.push_back(DYJets_Madgraph);
  }
  else if (sample=="vv"){
    SampleList.push_back(WZ);
    SampleList.push_back(ZZ);
    SampleList.push_back(WWTo2L2Nu_Madgraph);
  }
  else if (sample=="rare"){
    //    SampleList.push_back(TTWJets );
    //    SampleList.push_back(TTZJets );
    //    SampleList.push_back(TTGJets );
    SampleList.push_back(TTWWJets);
    SampleList.push_back(WWWJets );
    SampleList.push_back(WWZJets );
    SampleList.push_back(WZZJets );
    SampleList.push_back(ZZZJets );
  }
  else if (sample=="fake"){
    SampleList.push_back(TTJetsSemiLeptMGtauola);
    SampleList.push_back(Wbb_Madgraph          );
    SampleList.push_back(WgammaToLNuG          );
  }
  else if (sample=="Data"){
    SampleList.push_back(DoubleMu);
    SampleList.push_back(DoubleElectron);
    SampleList.push_back(MuEG);
  }
  else if (sample=="total"){
    SampleList.push_back(TbarWDilep);
    SampleList.push_back(TWDilep);
    SampleList.push_back(ZJets_Madgraph);
    SampleList.push_back(DYJets_Madgraph);
    SampleList.push_back(WZ);
    SampleList.push_back(ZZ);
    SampleList.push_back(WWTo2L2Nu_Madgraph);
    //    SampleList.push_back(TTWJets );
    //    SampleList.push_back(TTZJets );
    //    SampleList.push_back(TTGJets );
    SampleList.push_back(TTWWJets);
    SampleList.push_back(WWWJets );
    SampleList.push_back(WWZJets );
    SampleList.push_back(WZZJets );
    SampleList.push_back(ZZZJets );
    SampleList.push_back(TTJetsSemiLeptMGtauola);
    SampleList.push_back(Wbb_Madgraph          );
    SampleList.push_back(WgammaToLNuG          );
  }
  else { 
    cout << "[ERROR]: you have not selected a valid category..." << endl;
    break;
  }
  
  // NOW LOAD YIELDS FOR THE SELECTED CATEGORY:
  for (size_t cut=0; cut<iNCUTS; cut++){
    for (size_t sl=0; sl<SampleList.size(); sl++){
      for (size_t chan=0; chan<gNCHANNELS; chan++){
	Cat.name = sample;
	
	// YIELDS
	if (sl==0) Cat.Yields[chan][cut]  = S[SampleList[sl]].Yields[chan][cut];
	else       Cat.Yields[chan][cut] += S[SampleList[sl]].Yields[chan][cut];
	
	if (sl==0) Cat.Yields_stat[chan][cut]  = S[SampleList[sl]].Yields_stat[chan][cut];
	else       Cat.Yields_stat[chan][cut] += S[SampleList[sl]].Yields_stat[chan][cut];

	for (size_t sys=0; sys<gNSYST; sys++){
	  if (sl==0) Cat.Yields_syst[chan][cut]  = S[SampleList[sl]].Yields_syst[chan][cut];
	  else       Cat.Yields_syst[chan][cut] += S[SampleList[sl]].Yields_syst[chan][cut];
	}

	// SSYIELDS
	if (sl==0) Cat.SSYields[chan][cut]  = S[SampleList[sl]].SSYields[chan][cut];
	else       Cat.SSYields[chan][cut] += S[SampleList[sl]].SSYields[chan][cut];

	if (sl==0) Cat.SSYields_stat[chan][cut]  = S[SampleList[sl]].SSYields_stat[chan][cut];
	else       Cat.SSYields_stat[chan][cut] += S[SampleList[sl]].SSYields_stat[chan][cut];
      }
    }
  }
  // NOW LOAD KINEMATIC HISTOGRAMS
  
}
*/
void TopPlotter::LoadCategories(){
  
  // Now Load everything (Yields, SSYields and KinHistos) to the different bkg. categories...
  for (size_t cut=0; cut<iNCUTS; cut++){
   Data.name = "Data";

    Data .Yields[Muon][cut] = S[DoubleMuon].Yields[Muon][cut];
    Data .Yields[Elec][cut] = S[DoubleEG]  .Yields[Elec][cut];
    Data .Yields[ElMu][cut] = S[MuEG]      .Yields[ElMu][cut];
    
    Data .Yields_stat[Muon][cut] = S[DoubleMuon].Yields_stat[Muon][cut];
    Data .Yields_stat[Elec][cut] = S[DoubleEG]  .Yields_stat[Elec][cut];
    Data .Yields_stat[ElMu][cut] = S[MuEG]      .Yields_stat[ElMu][cut];
    
    Data .SSYields[Muon][cut] = S[DoubleMuon].SSYields[Muon][cut];
    Data .SSYields[Elec][cut] = S[DoubleEG]  .SSYields[Elec][cut];
    Data .SSYields[ElMu][cut] = S[MuEG]      .SSYields[ElMu][cut];
   
    Data .SSYields_stat[Muon][cut] = S[DoubleMuon]      .SSYields_stat[Muon][cut];
    Data .SSYields_stat[Elec][cut] = S[DoubleEG].SSYields_stat[Elec][cut];
    Data .SSYields_stat[ElMu][cut] = S[MuEG]          .SSYields_stat[ElMu][cut];
 
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      TTbar.name = "ttbar";
      if (gUseTTMadSpin) TTbar.Yields[chan][cut] = S[TTJets_MadSpin]        .Yields[chan][cut]; 
      //else               TTbar.Yields[chan][cut] = S[TTJetsFullLeptMGtauola].Yields[chan][cut]; 
      
      SUSYstop.name = "SUSYstop";
      //SUSYstop.Yields[chan][cut]  = S[T2tt_150to250LSP1to100_LeptonFilter].Yields[chan][cut];
      
      STop.name = "stop";
      STop .Yields[chan][cut]  = S[TbarWDilep].Yields[chan][cut];
      STop .Yields[chan][cut] += S[TWDilep]   .Yields[chan][cut];
      
      DY.name = "dy";
      DY   .Yields[chan][cut]  = S[ZJets_Madgraph].Yields[chan][cut];
      DY   .Yields[chan][cut] += S[DYJets_Madgraph].Yields[chan][cut];
      
      VV.name = "vv";
      VV   .Yields[chan][cut]  = S[WZ]                .Yields[chan][cut];
      VV   .Yields[chan][cut] += S[ZZ]                .Yields[chan][cut];
      VV   .Yields[chan][cut] += S[WWTo2L2Nu_Madgraph].Yields[chan][cut];
      
      Rare.name = "rare";
      //Rare .Yields[chan][cut]  = S[TTWJets] .Yields[chan][cut];
      //Rare .Yields[chan][cut] += S[TTZJets] .Yields[chan][cut];
      //      Rare .Yields[chan][cut] += S[TTGJets] .Yields[chan][cut];
      //Rare .Yields[chan][cut] += S[TTWWJets].Yields[chan][cut];
      //Rare .Yields[chan][cut] += S[WWWJets] .Yields[chan][cut];
      //Rare .Yields[chan][cut] += S[WWZJets] .Yields[chan][cut];
      //Rare .Yields[chan][cut] += S[WZZJets] .Yields[chan][cut];
      //Rare .Yields[chan][cut] += S[ZZZJets] .Yields[chan][cut];

      Fake.name = "fake";
      Fake .Yields[chan][cut]  = S[TTJetsSemiLeptMGtauola].Yields[chan][cut];
      Fake .Yields[chan][cut]  += S[Wbb_Madgraph]          .Yields[chan][cut];
      //Fake .Yields[chan][cut] += S[WgammaToLNuG]          .Yields[chan][cut];
//if(cut==3 && chan==2){
//cout << chan << " " << cut << " tt: " << S[TTJetsSemiLeptMGtauola].Yields[chan][cut] << endl;
//cout << chan << " " << cut << " wj: " << S[Wbb_Madgraph].Yields[chan][cut] << endl;
//cout << chan << " " << cut << " fake" << Fake .Yields[chan][cut] << endl;
//}

      Total.name = "total";
      Total.Yields[chan][cut]  = STop.Yields[chan][cut];
      Total.Yields[chan][cut] += DY  .Yields[chan][cut];
      Total.Yields[chan][cut] += VV  .Yields[chan][cut];
      Total.Yields[chan][cut] += Rare.Yields[chan][cut];
      Total.Yields[chan][cut] += Fake.Yields[chan][cut];

      // SYSTEMATIC ERROR
      for (size_t sys=0; sys<gNSYST; sys++){
	if (gUseTTMadSpin) {
	  TTbar.Yields_syst[chan][cut][sys]=S[TTJets_MadSpin].        Yields_syst[chan][cut][sys];
	}
	//else               {
	//  TTbar.Yields_syst[chan][cut][sys]=S[TTJetsFullLeptMGtauola].Yields_syst[chan][cut][sys];
	//}
      
	//SUSYstop.Yields_syst[chan][cut][sys]  = S[T2tt_150to250LSP1to100_LeptonFilter].Yields_syst[chan][cut][sys];
      //cout << chan << " " << cut << " "<< TTbar.Yields[chan][cut] << endl;
	
	STop .Yields_syst[chan][cut][sys]  = S[TbarWDilep].Yields_syst[chan][cut][sys];
	STop .Yields_syst[chan][cut][sys] += S[TWDilep]   .Yields_syst[chan][cut][sys];
	
	DY   .Yields_syst[chan][cut][sys]  = S[ZJets_Madgraph] .Yields_syst[chan][cut][sys];
	DY   .Yields_syst[chan][cut][sys] += S[DYJets_Madgraph].Yields_syst[chan][cut][sys];
	VV   .Yields_syst[chan][cut][sys]  = S[WZ]                .Yields_syst[chan][cut][sys];
	VV   .Yields_syst[chan][cut][sys] += S[ZZ]                .Yields_syst[chan][cut][sys];
	VV   .Yields_syst[chan][cut][sys] += S[WWTo2L2Nu_Madgraph].Yields_syst[chan][cut][sys];

	Rare .Yields_syst[chan][cut][sys]  = 0.;//S[TTWJets] .Yields_syst[chan][cut][sys];
	//Rare .Yields_syst[chan][cut][sys] += S[TTZJets] .Yields_syst[chan][cut][sys];
	//Rare .Yields_syst[chan][cut][sys] += S[TTWWJets].Yields_syst[chan][cut][sys];
	//Rare .Yields_syst[chan][cut][sys] += S[WWWJets] .Yields_syst[chan][cut][sys];
	//Rare .Yields_syst[chan][cut][sys] += S[WWZJets] .Yields_syst[chan][cut][sys];
	//Rare .Yields_syst[chan][cut][sys] += S[WZZJets] .Yields_syst[chan ][cut][sys];
	//Rare .Yields_syst[chan][cut][sys] += S[ZZZJets] .Yields_syst[chan][cut][sys];

	Fake .Yields_syst[chan][cut][sys]  = S[TTJetsSemiLeptMGtauola].Yields_syst[chan][cut][sys];
	Fake .Yields_syst[chan][cut][sys]  += S[Wbb_Madgraph]         .Yields_syst[chan][cut][sys];
	//Fake .Yields_syst[chan][cut][sys] += S[WgammaToLNuG]          .Yields_syst[chan][cut][sys];
      }
      
      // STATISTICAL ERROR
      if (gUseTTMadSpin) TTbar.Yields_stat[chan][cut]  = S[TTJets_MadSpin].        Yields_stat[chan][cut] * S[TTJets_MadSpin].        Yields_stat[chan][cut];
      //else               TTbar.Yields_stat[chan][cut]  = S[TTJetsFullLeptMGtauola].Yields_stat[chan][cut] * S[TTJetsFullLeptMGtauola].Yields_stat[chan][cut];
      
      //SUSYstop .Yields_stat[chan][cut]  = S[T2tt_150to250LSP1to100_LeptonFilter].Yields_stat[chan][cut] * S[T2tt_150to250LSP1to100_LeptonFilter].Yields_stat[chan][cut];
      
      STop .Yields_stat[chan][cut]  = S[TbarWDilep].Yields_stat[chan][cut] * S[TbarWDilep].Yields_stat[chan][cut];
      STop .Yields_stat[chan][cut] += S[TWDilep]   .Yields_stat[chan][cut] * S[TWDilep]   .Yields_stat[chan][cut];
      
      DY   .Yields_stat[chan][cut]  = S[ZJets_Madgraph] .Yields_stat[chan][cut] * S[ZJets_Madgraph] .Yields_stat[chan][cut];
      DY   .Yields_stat[chan][cut] += S[DYJets_Madgraph].Yields_stat[chan][cut] * S[DYJets_Madgraph].Yields_stat[chan][cut];
      
      VV   .Yields_stat[chan][cut]  = S[WZ]                .Yields_stat[chan][cut] * S[WZ]                .Yields_stat[chan][cut];
      VV   .Yields_stat[chan][cut] += S[ZZ]                .Yields_stat[chan][cut] * S[ZZ]                .Yields_stat[chan][cut];
      VV   .Yields_stat[chan][cut] += S[WWTo2L2Nu_Madgraph].Yields_stat[chan][cut] * S[WWTo2L2Nu_Madgraph].Yields_stat[chan][cut];
	
      //Rare .Yields_stat[chan][cut]  = S[TTWJets] .Yields_stat[chan][cut] * S[TTWJets] .Yields_stat[chan][cut];
      //Rare .Yields_stat[chan][cut] += S[TTZJets] .Yields_stat[chan][cut] * S[TTZJets] .Yields_stat[chan][cut];
      //Rare .Yields_stat[chan][cut] += S[TTWWJets].Yields_stat[chan][cut] * S[TTWWJets].Yields_stat[chan][cut];
      //Rare .Yields_stat[chan][cut] += S[WWWJets] .Yields_stat[chan][cut] * S[WWWJets] .Yields_stat[chan][cut];
      //Rare .Yields_stat[chan][cut] += S[WWZJets] .Yields_stat[chan][cut] * S[WWZJets] .Yields_stat[chan][cut];
      //Rare .Yields_stat[chan][cut] += S[WZZJets] .Yields_stat[chan][cut] * S[WZZJets] .Yields_stat[chan][cut];
      //Rare .Yields_stat[chan][cut] += S[ZZZJets] .Yields_stat[chan][cut] * S[ZZZJets] .Yields_stat[chan][cut];
	
      Fake .Yields_stat[chan][cut]   = S[TTJetsSemiLeptMGtauola].Yields_stat[chan][cut] * S[TTJetsSemiLeptMGtauola].Yields_stat[chan][cut];
      Fake .Yields_stat[chan][cut]  += S[Wbb_Madgraph]          .Yields_stat[chan][cut] * S[Wbb_Madgraph]          .Yields_stat[chan][cut];
      //Fake .Yields_stat[chan][cut] += S[WgammaToLNuG]          .Yields_stat[chan][cut] * S[WgammaToLNuG]          .Yields_stat[chan][cut];
       
      SUSYstop  .Yields_stat[chan][cut] = TMath::Sqrt(SUSYstop  .Yields_stat[chan][cut]);
      TTbar.Yields_stat[chan][cut] = TMath::Sqrt(TTbar.Yields_stat[chan][cut]);
      STop .Yields_stat[chan][cut] = TMath::Sqrt(STop .Yields_stat[chan][cut]);
      DY   .Yields_stat[chan][cut] = TMath::Sqrt(DY   .Yields_stat[chan][cut]);
      VV   .Yields_stat[chan][cut] = TMath::Sqrt(VV   .Yields_stat[chan][cut]);
      Rare .Yields_stat[chan][cut] = 0.;//TMath::Sqrt(Rare .Yields_stat[chan][cut]);
      Fake .Yields_stat[chan][cut] = TMath::Sqrt(Fake .Yields_stat[chan][cut]);
              
      Total.Yields_stat[chan][cut]  = Fake .Yields_stat[chan][cut] * Fake .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += Rare .Yields_stat[chan][cut] * Rare .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += VV   .Yields_stat[chan][cut] * VV   .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += DY   .Yields_stat[chan][cut] * DY   .Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut] += STop .Yields_stat[chan][cut] * STop .Yields_stat[chan][cut];
      //Total.Yields_stat[chan][cut] += TTbar.Yields_stat[chan][cut] * TTbar.Yields_stat[chan][cut];
      Total.Yields_stat[chan][cut]  = TMath::Sqrt(Total.Yields_stat[chan][cut]);

      if(cut==4 && chan==2){
         cout << "Fake = " << Fake .Yields_stat[chan][cut]<< ", " 
              << "Rare = " << Rare .Yields_stat[chan][cut]<< ", "
              << "VV   = " << VV   .Yields_stat[chan][cut]<< ", "
              << "DY   = " << DY   .Yields_stat[chan][cut]<< ", "
              //<< "DY DD= " << DD_DY.Yields_stat[chan][cut]<< ", "
              << "STop = " << STop .Yields_stat[chan][cut]<< ", " << endl;
         cout << "Tot  = " << Total.Yields_stat[chan][cut]<< endl;       
      }
      
      // SS Yields
      if (gUseTTMadSpin) TTbar.SSYields[chan][cut]  = S[TTJets_MadSpin].        SSYields[chan][cut];
      //else               TTbar.SSYields[chan][cut]  = S[TTJetsFullLeptMGtauola].SSYields[chan][cut];
      
      STop .SSYields[chan][cut]  = S[TbarWDilep].SSYields[chan][cut];
      STop .SSYields[chan][cut] += S[TWDilep]   .SSYields[chan][cut];
            
      DY   .SSYields[chan][cut]  = S[ZJets_Madgraph].SSYields[chan][cut];
      DY   .SSYields[chan][cut] += S[DYJets_Madgraph].SSYields[chan][cut];
      
      VV   .SSYields[chan][cut]  = S[WZ]                .SSYields[chan][cut];
      VV   .SSYields[chan][cut] += S[ZZ]                .SSYields[chan][cut];
      VV   .SSYields[chan][cut] += S[WWTo2L2Nu_Madgraph].SSYields[chan][cut];
      
      //Rare .SSYields[chan][cut]  = S[TTWJets] .SSYields[chan][cut];
      //Rare .SSYields[chan][cut] += S[TTZJets] .SSYields[chan][cut];
      //Rare .SSYields[chan][cut] += S[TTWWJets].SSYields[chan][cut];
      //Rare .SSYields[chan][cut] += S[WWWJets] .SSYields[chan][cut];
      //Rare .SSYields[chan][cut] += S[WWZJets] .SSYields[chan][cut];
      //Rare .SSYields[chan][cut] += S[WZZJets] .SSYields[chan][cut];
      //Rare .SSYields[chan][cut] += S[ZZZJets] .SSYields[chan][cut];

      Fake .SSYields[chan][cut]   = S[TTJetsSemiLeptMGtauola].SSYields[chan][cut];
      Fake .SSYields[chan][cut]  += S[Wbb_Madgraph]          .SSYields[chan][cut];
      //Fake .SSYields[chan][cut] += S[WgammaToLNuG]          .SSYields[chan][cut];
      
//if(cut==3 && chan==2){
//cout << chan << " " << cut << " tt: " << S[TTJetsSemiLeptMGtauola].SSYields[chan][cut] << endl;
//cout << chan << " " << cut << " wj: " << S[Wbb_Madgraph].SSYields[chan][cut] << endl;
//cout << chan << " " << cut << " fake" << Fake .SSYields[chan][cut] << endl;
//}
      Total.SSYields[chan][cut]  = STop .SSYields[chan][cut];
      Total.SSYields[chan][cut] += DY   .SSYields[chan][cut];
      Total.SSYields[chan][cut] += VV   .SSYields[chan][cut];
      Total.SSYields[chan][cut] += Rare .SSYields[chan][cut];
      Total.SSYields[chan][cut] += TTbar.SSYields[chan][cut];

      if (gUseTTMadSpin) TTbar.SSYields_stat[chan][cut]  = S[TTJets_MadSpin].        SSYields_stat[chan][cut] * S[TTJets_MadSpin].        SSYields_stat[chan][cut];
      //else               TTbar.SSYields_stat[chan][cut]  = S[TTJetsFullLeptMGtauola].SSYields_stat[chan][cut] * S[TTJetsFullLeptMGtauola].SSYields_stat[chan][cut];
      
      STop .SSYields_stat[chan][cut]  = S[TbarWDilep].SSYields_stat[chan][cut] * S[TbarWDilep].SSYields_stat[chan][cut];
      STop .SSYields_stat[chan][cut] += S[TWDilep]   .SSYields_stat[chan][cut] * S[TWDilep]   .SSYields_stat[chan][cut];
      
      DY   .SSYields_stat[chan][cut]  = S[ZJets_Madgraph] .SSYields_stat[chan][cut] * S[ZJets_Madgraph] .SSYields_stat[chan][cut];
      DY   .SSYields_stat[chan][cut] += S[DYJets_Madgraph].SSYields_stat[chan][cut] * S[DYJets_Madgraph].SSYields_stat[chan][cut];
      
      VV   .SSYields_stat[chan][cut]  = S[WZ]                .SSYields_stat[chan][cut] * S[WZ]                .SSYields_stat[chan][cut];
      VV   .SSYields_stat[chan][cut] += S[ZZ]                .SSYields_stat[chan][cut] * S[ZZ]                .SSYields_stat[chan][cut];
      VV   .SSYields_stat[chan][cut] += S[WWTo2L2Nu_Madgraph].SSYields_stat[chan][cut] * S[WWTo2L2Nu_Madgraph].SSYields_stat[chan][cut];
	
      //Rare .SSYields_stat[chan][cut]  = S[TTWJets] .SSYields_stat[chan][cut] * S[TTWJets] .SSYields_stat[chan][cut];
      //Rare .SSYields_stat[chan][cut] += S[TTZJets] .SSYields_stat[chan][cut] * S[TTZJets] .SSYields_stat[chan][cut];
      //      Rare .SSYields_stat[chan][cut] += S[TTGJets] .SSYields_stat[chan][cut] * S[TTGJets] .SSYields_stat[chan][cut];
      //Rare .SSYields_stat[chan][cut] += S[TTWWJets].SSYields_stat[chan][cut] * S[TTWWJets].SSYields_stat[chan][cut];
      //Rare .SSYields_stat[chan][cut] += S[WWWJets] .SSYields_stat[chan][cut] * S[WWWJets] .SSYields_stat[chan][cut];
      //Rare .SSYields_stat[chan][cut] += S[WWZJets] .SSYields_stat[chan][cut] * S[WWZJets] .SSYields_stat[chan][cut];
      //Rare .SSYields_stat[chan][cut] += S[WZZJets] .SSYields_stat[chan][cut] * S[WZZJets] .SSYields_stat[chan][cut];
      //Rare .SSYields_stat[chan][cut] += S[ZZZJets] .SSYields_stat[chan][cut] * S[ZZZJets] .SSYields_stat[chan][cut];
	
      Fake .SSYields_stat[chan][cut]  = S[TTJetsSemiLeptMGtauola].SSYields_stat[chan][cut] * S[TTJetsSemiLeptMGtauola].SSYields_stat[chan][cut];
      Fake .SSYields_stat[chan][cut]  += S[Wbb_Madgraph]          .SSYields_stat[chan][cut] * S[Wbb_Madgraph]          .SSYields_stat[chan][cut];
      //Fake .SSYields_stat[chan][cut] += S[WgammaToLNuG]          .SSYields_stat[chan][cut] * S[WgammaToLNuG]          .SSYields_stat[chan][cut];
      
      TTbar.SSYields_stat[chan][cut] = TMath::Sqrt(TTbar.SSYields_stat[chan][cut]);
      STop .SSYields_stat[chan][cut] = TMath::Sqrt(STop .SSYields_stat[chan][cut]);
      DY   .SSYields_stat[chan][cut] = TMath::Sqrt(DY   .SSYields_stat[chan][cut]);
      VV   .SSYields_stat[chan][cut] = TMath::Sqrt(VV   .SSYields_stat[chan][cut]);
      Rare .SSYields_stat[chan][cut] = TMath::Sqrt(Rare .SSYields_stat[chan][cut]);
      Fake .SSYields_stat[chan][cut] = TMath::Sqrt(Fake .SSYields_stat[chan][cut]);

      Total.SSYields_stat[chan][cut]  = Rare .SSYields_stat[chan][cut] * Rare .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += VV   .SSYields_stat[chan][cut] * VV   .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += DY   .SSYields_stat[chan][cut] * DY   .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += STop .SSYields_stat[chan][cut] * STop .SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut] += TTbar.SSYields_stat[chan][cut] * TTbar.SSYields_stat[chan][cut];
      Total.SSYields_stat[chan][cut]  = TMath::Sqrt(Total.SSYields_stat[chan][cut]);
      
    }
  
  
 
    /// FOR DY ESTIMATION:
    Data.MllHistos[Muon][cut] = (TH1F*)S[DoubleMuon].MllHistos[Muon][cut]->Clone();
    Data.MllHistos[Elec][cut] = (TH1F*)S[DoubleEG]  .MllHistos[Elec][cut]->Clone();
    Data.MllHistos[ElMu][cut] = (TH1F*)S[MuEG]      .MllHistos[ElMu][cut]->Clone();
  
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      DY.MllHistos[chan][cut] = (TH1F*)S[ZJets_Madgraph].MllHistos[chan][cut]->Clone();
      DY.MllHistos[chan][cut] ->Add(   S[DYJets_Madgraph] .MllHistos[chan][cut]);	
    }
  }
  cout << " ---> Kinematic Histograms..." << endl;
  for (size_t chan=0; chan<gNCHANNELS+1; chan++){
    for (size_t cut=0; cut<iNCUTS; cut++){
     for (size_t var=0; var<gNVARS; var++){
	if (chan == gNCHANNELS){
	  Data .KinHistos[chan][cut][var] = (TH1F*)S[DoubleMuon]              .KinHistos[Muon][cut][var]->Clone();
	  Data .KinHistos[chan][cut][var] ->Add(   S[DoubleEG]        .KinHistos[Elec][cut][var]);
	  if (gUseTTMadSpin) { 
	    TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJets_MadSpin]        .KinHistos[Muon][cut][var]->Clone();
	    TTbar.KinHistos[chan][cut][var] ->Add(   S[TTJets_MadSpin]        .KinHistos[Elec][cut][var]);
	  }
	  //else {
	  //  TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJetsFullLeptMGtauola].KinHistos[Muon][cut][var]->Clone();
	  //  TTbar.KinHistos[chan][cut][var] ->Add(   S[TTJetsFullLeptMGtauola].KinHistos[Elec][cut][var]);
	  //}
	  
	  STop .KinHistos[chan][cut][var] = (TH1F*)S[TbarWDilep]            .KinHistos[Muon][cut][var]->Clone();
	  STop .KinHistos[chan][cut][var] ->Add(   S[TbarWDilep]            .KinHistos[Elec][cut][var]);
	  STop .KinHistos[chan][cut][var] ->Add(   S[TWDilep]               .KinHistos[Muon][cut][var]);
	  STop .KinHistos[chan][cut][var] ->Add(   S[TWDilep]               .KinHistos[Elec][cut][var]);
	  DY   .KinHistos[chan][cut][var] = (TH1F*)S[ZJets_Madgraph]        .KinHistos[Muon][cut][var]->Clone();
	  DY   .KinHistos[chan][cut][var] ->Add(   S[DYJets_Madgraph]       .KinHistos[Muon][cut][var]);
	  DY   .KinHistos[chan][cut][var] ->Add(   S[ZJets_Madgraph]        .KinHistos[Elec][cut][var]);
	  DY   .KinHistos[chan][cut][var] ->Add(   S[DYJets_Madgraph]       .KinHistos[Elec][cut][var]);
	  VV   .KinHistos[chan][cut][var] = (TH1F*)S[WZ]    .KinHistos[Muon][cut][var]->Clone();
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WWTo2L2Nu_Madgraph]    .KinHistos[Muon][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[ZZ]                    .KinHistos[Muon][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WWTo2L2Nu_Madgraph]    .KinHistos[Elec][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WZ]                    .KinHistos[Elec][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[ZZ]                    .KinHistos[Elec][cut][var]);
	  Fake .KinHistos[chan][cut][var] = (TH1F*)S[TTJetsSemiLeptMGtauola].KinHistos[Muon][cut][var]->Clone();  
	  Fake .KinHistos[chan][cut][var] ->Add(   S[Wbb_Madgraph]          .KinHistos[Muon][cut][var]);
	  //Fake .KinHistos[chan][cut][var] ->Add(   S[WgammaToLNuG]          .KinHistos[Muon][cut][var]);
	  Fake .KinHistos[chan][cut][var] = (TH1F*)S[TTJetsSemiLeptMGtauola].KinHistos[Elec][cut][var]->Clone();
	  Fake .KinHistos[chan][cut][var] ->Add(   S[Wbb_Madgraph]          .KinHistos[Elec][cut][var]);
	  //Fake .KinHistos[chan][cut][var] ->Add(   S[WgammaToLNuG]          .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] = (TH1F*)S[TTWJets]               .KinHistos[Muon][cut][var]->Clone();  
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTWWJets]              .KinHistos[Muon][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTZJets]               .KinHistos[Muon][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WWWJets]               .KinHistos[Muon][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WWZJets]               .KinHistos[Muon][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WZZJets]               .KinHistos[Muon][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[ZZZJets]               .KinHistos[Muon][cut][var]);
	  //	  Rare .KinHistos[chan][cut][var] ->Add(   S[TTGJets]               .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTWJets]               .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTWWJets]              .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTZJets]               .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WWWJets]               .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WWZJets]               .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WZZJets]               .KinHistos[Elec][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[ZZZJets]               .KinHistos[Elec][cut][var]);
	}
	else {
	  if (chan == Muon){
	    Data .KinHistos[chan][cut][var] = (TH1F*)S[DoubleMuon]            .KinHistos[chan][cut][var]->Clone();
	  }
	  else if (chan == Elec){
	    Data .KinHistos[chan][cut][var] = (TH1F*)S[DoubleEG]      .KinHistos[chan][cut][var]->Clone();
	  }
	  else if (chan == ElMu){
	    Data .KinHistos[chan][cut][var] = (TH1F*)S[MuEG]                .KinHistos[chan][cut][var]->Clone();
	  }
	  if (gUseTTMadSpin) 
	    TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJets_MadSpin]        .KinHistos[chan][cut][var]->Clone();
	  //else 
	  //  TTbar.KinHistos[chan][cut][var] = (TH1F*)S[TTJetsFullLeptMGtauola].KinHistos[chan][cut][var]->Clone();
	  
	  STop .KinHistos[chan][cut][var] = (TH1F*)S[TbarWDilep]            .KinHistos[chan][cut][var]->Clone();
	  STop .KinHistos[chan][cut][var] ->Add(   S[TWDilep]               .KinHistos[chan][cut][var]);
	  DY   .KinHistos[chan][cut][var] = (TH1F*)S[ZJets_Madgraph]        .KinHistos[chan][cut][var]->Clone();
	  DY   .KinHistos[chan][cut][var] ->Add(   S[DYJets_Madgraph]       .KinHistos[chan][cut][var]);
	  VV   .KinHistos[chan][cut][var] = (TH1F*)S[WZ]    .KinHistos[chan][cut][var]->Clone();
	  VV   .KinHistos[chan][cut][var] ->Add(   S[WWTo2L2Nu_Madgraph]    .KinHistos[chan][cut][var]);
	  VV   .KinHistos[chan][cut][var] ->Add(   S[ZZ]                    .KinHistos[chan][cut][var]);
	  Fake .KinHistos[chan][cut][var] = (TH1F*)S[Wbb_Madgraph]          .KinHistos[chan][cut][var]->Clone();  
	  Fake .KinHistos[chan][cut][var] ->Add(   S[TTJetsSemiLeptMGtauola].KinHistos[chan][cut][var]);
	  //Fake .KinHistos[chan][cut][var] ->Add(   S[WgammaToLNuG]          .KinHistos[chan][cut][var]);
	  //Rare .KinHistos[chan][cut][var] = (TH1F*)S[TTWJets]               .KinHistos[chan][cut][var]->Clone();  
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTWWJets]              .KinHistos[chan][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[TTZJets]               .KinHistos[chan][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WWWJets]               .KinHistos[chan][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WWZJets]               .KinHistos[chan][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[WZZJets]               .KinHistos[chan][cut][var]);
	  //Rare .KinHistos[chan][cut][var] ->Add(   S[ZZZJets]               .KinHistos[chan][cut][var]);

	  
	}
	SetupDraw(Data .KinHistos[chan][cut][var], kBlack, var); 
	SetupDraw(TTbar.KinHistos[chan][cut][var], 633   , var);  
	SetupDraw(DY   .KinHistos[chan][cut][var], 852   , var);  
	SetupDraw(VV   .KinHistos[chan][cut][var], 390   , var);  
	SetupDraw(STop .KinHistos[chan][cut][var], 616   , var);  
	SetupDraw(Fake .KinHistos[chan][cut][var], 413   , var);  
	//SetupDraw(Rare .KinHistos[chan][cut][var], kYellow  , var);  
     }
    }
    TTbar.SystError[chan][SFIDISO]  = TMath::Sqrt(S[TTJets_MadSpin].SystError[chan][SFIDISO]*S[TTJets_MadSpin].SystError[chan][SFIDISO]);

    //SUSYstop .SystError[chan][SFIDISO]  = TMath::Sqrt(S[T2tt_150to250LSP1to100_LeptonFilter].SystError[chan][SFIDISO]*S[T2tt_150to250LSP1to100_LeptonFilter].SystError[chan][SFIDISO]);

    STop .SystError[chan][SFIDISO]  = S[TbarWDilep].SystError[chan][SFIDISO]*S[TbarWDilep].SystError[chan][SFIDISO];
    STop .SystError[chan][SFIDISO] += S[TWDilep]   .SystError[chan][SFIDISO]*S[TWDilep]   .SystError[chan][SFIDISO];
    STop .SystError[chan][SFIDISO]  = TMath::Sqrt(STop .SystError[chan][SFIDISO]);
    
    VV   .SystError[chan][SFIDISO]  = S[WWTo2L2Nu_Madgraph].SystError[chan][SFIDISO]*S[WWTo2L2Nu_Madgraph].SystError[chan][SFIDISO];
    VV   .SystError[chan][SFIDISO] += S[WZ].SystError[chan][SFIDISO]*S[WZ].SystError[chan][SFIDISO];
    VV   .SystError[chan][SFIDISO] += S[ZZ].SystError[chan][SFIDISO]*S[ZZ].SystError[chan][SFIDISO];
    VV   .SystError[chan][SFIDISO]  = TMath::Sqrt(VV.SystError[chan][SFIDISO]);

    Fake .SystError[chan][SFIDISO]  = S[Wbb_Madgraph]          .SystError[chan][SFIDISO]*S[Wbb_Madgraph]          .SystError[chan][SFIDISO];
    Fake .SystError[chan][SFIDISO] += S[TTJetsSemiLeptMGtauola].SystError[chan][SFIDISO]*S[TTJetsSemiLeptMGtauola].SystError[chan][SFIDISO];
    //Fake .SystError[chan][SFIDISO] += S[WgammaToLNuG]          .SystError[chan][SFIDISO]*S[WgammaToLNuG]          .SystError[chan][SFIDISO];
    Fake .SystError[chan][SFIDISO]  = TMath::Sqrt(Fake .SystError[chan][SFIDISO]);

    //  Rare .SystError[chan][SFIDISO]  = S[TTGJets]               .SystError[chan][SFIDISO]*S[TTGJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO]  = S[TTWJets]               .SystError[chan][SFIDISO]*S[TTWJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO] += S[TTWWJets]              .SystError[chan][SFIDISO]*S[TTWWJets]              .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO] += S[TTZJets]               .SystError[chan][SFIDISO]*S[TTZJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO] += S[WWWJets]               .SystError[chan][SFIDISO]*S[WWWJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO] += S[WWZJets]               .SystError[chan][SFIDISO]*S[WWZJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO] += S[WZZJets]               .SystError[chan][SFIDISO]*S[WZZJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO] += S[ZZZJets]               .SystError[chan][SFIDISO]*S[ZZZJets]               .SystError[chan][SFIDISO];
    //Rare .SystError[chan][SFIDISO]  = TMath::Sqrt(Rare .SystError[chan][SFIDISO]);

    //SUSYstop .SystError[chan][SFTrig]  = TMath::Sqrt(S[T2tt_150to250LSP1to100_LeptonFilter].SystError[chan][SFTrig]*S[T2tt_150to250LSP1to100_LeptonFilter].SystError[chan][SFTrig]);

    TTbar.SystError[chan][SFTrig]  = TMath::Sqrt(S[TTJets_MadSpin].SystError[chan][SFTrig]*S[TTJets_MadSpin].SystError[chan][SFTrig]);

    STop .SystError[chan][SFTrig]  = S[TbarWDilep].SystError[chan][SFTrig]*S[TbarWDilep].SystError[chan][SFTrig];
    STop .SystError[chan][SFTrig] += S[TWDilep]   .SystError[chan][SFTrig]*S[TWDilep]   .SystError[chan][SFTrig];
    STop .SystError[chan][SFTrig]  = TMath::Sqrt(STop .SystError[chan][SFTrig]);
    
    VV   .SystError[chan][SFTrig]  = S[WWTo2L2Nu_Madgraph].SystError[chan][SFTrig]*S[WWTo2L2Nu_Madgraph].SystError[chan][SFTrig];
    VV   .SystError[chan][SFTrig] += S[WZ].SystError[chan][SFTrig]*S[WZ].SystError[chan][SFTrig];
    VV   .SystError[chan][SFTrig] += S[ZZ].SystError[chan][SFTrig]*S[ZZ].SystError[chan][SFTrig];
    VV   .SystError[chan][SFTrig]  = TMath::Sqrt(VV.SystError[chan][SFTrig]);

    Fake .SystError[chan][SFTrig]  = S[Wbb_Madgraph]          .SystError[chan][SFTrig]*S[Wbb_Madgraph]          .SystError[chan][SFTrig];
    Fake .SystError[chan][SFTrig] += S[TTJetsSemiLeptMGtauola].SystError[chan][SFTrig]*S[TTJetsSemiLeptMGtauola].SystError[chan][SFTrig];
    //Fake .SystError[chan][SFTrig] += S[WgammaToLNuG]          .SystError[chan][SFTrig]*S[WgammaToLNuG]          .SystError[chan][SFTrig];
    Fake .SystError[chan][SFTrig]  = TMath::Sqrt(Fake .SystError[chan][SFTrig]);

   // Rare .SystError[chan][SFTrig]  = S[TTWJets]               .SystError[chan][SFTrig]*S[TTWJets]               .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig] += S[TTWWJets]              .SystError[chan][SFTrig]*S[TTWWJets]              .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig] += S[TTZJets]               .SystError[chan][SFTrig]*S[TTZJets]               .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig] += S[WWWJets]               .SystError[chan][SFTrig]*S[WWWJets]               .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig] += S[WWZJets]               .SystError[chan][SFTrig]*S[WWZJets]               .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig] += S[WZZJets]               .SystError[chan][SFTrig]*S[WZZJets]               .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig] += S[ZZZJets]               .SystError[chan][SFTrig]*S[ZZZJets]               .SystError[chan][SFTrig];
    //Rare .SystError[chan][SFTrig]  = TMath::Sqrt(Rare .SystError[chan][SFTrig]);
  }
  
  // SYSTEMATIC ERRORS HISTOS
  cout << " ---> Systematic Histograms" << endl;
  for (size_t cut=0; cut<iNCUTS; cut++){
    for (size_t chan=0; chan<gNCHANNELS; chan++){

      for (size_t sys=0; sys<gNSYST; sys++){
	if (chan == Muon){
	  Data .NBtagsNJets[chan][cut][sys]   = (TH1F*)S[DoubleMuon]        .NBtagsNJets[chan][cut][0]->Clone();
	  Data .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[DoubleMuon]      .SSNBtagsNJets[chan][cut][0]->Clone();
	}
	else if (chan == Elec){
	  Data .NBtagsNJets[chan][cut][sys]   = (TH1F*)S[DoubleEG]  .NBtagsNJets[chan][cut][0]->Clone();
	  Data .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[DoubleEG].SSNBtagsNJets[chan][cut][0]->Clone();
	}
	else if (chan == ElMu){
	  Data .NBtagsNJets[chan][cut][sys]   = (TH1F*)S[MuEG]            .NBtagsNJets[chan][cut][0]->Clone();
	  Data .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[MuEG]          .SSNBtagsNJets[chan][cut][0]->Clone();
	}
	
	if (gUseTTMadSpin) {
	  TTbar.NBtagsNJets[chan][cut][sys]   = (TH1F*)S[TTJets_MadSpin]        .NBtagsNJets[chan][cut][sys]->Clone();
	  TTbar.SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[TTJets_MadSpin]      .SSNBtagsNJets[chan][cut][sys]->Clone();
	}
	//else {
	//  TTbar.NBtagsNJets[chan][cut][sys]   = (TH1F*)S[TTJetsFullLeptMGtauola].NBtagsNJets[chan][cut][sys]->Clone();
	//  TTbar.SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[TTJetsFullLeptMGtauola].SSNBtagsNJets[chan][cut][sys]->Clone();
	//}
        //SUSYstop.NBtagsNJets[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].NBtagsNJets[chan][cut][sys]->Clone();
	STop .NBtagsNJets[chan][cut][sys] = (TH1F*)S[TbarWDilep]            .NBtagsNJets[chan][cut][sys]->Clone();
	STop .NBtagsNJets[chan][cut][sys] ->Add(   S[TWDilep]               .NBtagsNJets[chan][cut][sys]);
	DY   .NBtagsNJets[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]        .NBtagsNJets[chan][cut][sys]->Clone();
	DY   .NBtagsNJets[chan][cut][sys] ->Add(   S[DYJets_Madgraph]        .NBtagsNJets[chan][cut][sys]);
	VV   .NBtagsNJets[chan][cut][sys] = (TH1F*)S[WZ]    .NBtagsNJets[chan][cut][sys]->Clone();
	VV   .NBtagsNJets[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]    .NBtagsNJets[chan][cut][sys]);
	VV   .NBtagsNJets[chan][cut][sys] ->Add(   S[ZZ]                    .NBtagsNJets[chan][cut][sys]);
	Fake .NBtagsNJets[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .NBtagsNJets[chan][cut][sys]->Clone();  
	Fake .NBtagsNJets[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].NBtagsNJets[chan][cut][sys]);
	//Fake .NBtagsNJets[chan][cut][sys] ->Add(   S[WgammaToLNuG]          .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] = (TH1F*)S[TTWJets]               .NBtagsNJets[chan][cut][sys]->Clone();  
	//	Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[TTWJets]               .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[TTWWJets]              .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[TTZJets]               .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[WWWJets]               .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[WWZJets]               .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[WZZJets]               .NBtagsNJets[chan][cut][sys]);
	//Rare .NBtagsNJets[chan][cut][sys] ->Add(   S[ZZZJets]               .NBtagsNJets[chan][cut][sys]);
	
	SetupDraw(Data .NBtagsNJets[chan][cut][sys], kBlack, NBTagsNJets); 
	SetupDraw(TTbar.NBtagsNJets[chan][cut][sys], 633   , NBTagsNJets);  
	SetupDraw(DY   .NBtagsNJets[chan][cut][sys], 852   , NBTagsNJets);  
	SetupDraw(VV   .NBtagsNJets[chan][cut][sys], 390   , NBTagsNJets);  
	SetupDraw(STop .NBtagsNJets[chan][cut][sys], 616   , NBTagsNJets);  
	SetupDraw(Fake .NBtagsNJets[chan][cut][sys], 413   , NBTagsNJets);  
	//SetupDraw(Rare .NBtagsNJets[chan][cut][sys], kYellow  , NBTagsNJets);  

        //SUSYstop.SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].SSNBtagsNJets[chan][cut][sys]->Clone();
	STop .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[TbarWDilep]            .SSNBtagsNJets[chan][cut][sys]->Clone();
	STop .SSNBtagsNJets[chan][cut][sys] ->Add(   S[TWDilep]               .SSNBtagsNJets[chan][cut][sys]);
	DY   .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]       .SSNBtagsNJets[chan][cut][sys]->Clone();
	DY   .SSNBtagsNJets[chan][cut][sys] ->Add(   S[DYJets_Madgraph]        .SSNBtagsNJets[chan][cut][sys]);
	VV   .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[WZ]    .SSNBtagsNJets[chan][cut][sys]->Clone();
	VV   .SSNBtagsNJets[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]    .SSNBtagsNJets[chan][cut][sys]);
	VV   .SSNBtagsNJets[chan][cut][sys] ->Add(   S[ZZ]                    .SSNBtagsNJets[chan][cut][sys]);
	Fake .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .SSNBtagsNJets[chan][cut][sys]->Clone();  
	Fake .SSNBtagsNJets[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].SSNBtagsNJets[chan][cut][sys]);
	//Fake .SSNBtagsNJets[chan][cut][sys] ->Add(   S[WgammaToLNuG]          .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] = (TH1F*)S[TTWJets]               .SSNBtagsNJets[chan][cut][sys]->Clone();  
	//	Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[TTWJets]               .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[TTWWJets]              .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[TTZJets]               .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[WWWJets]               .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[WWZJets]               .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[WZZJets]               .SSNBtagsNJets[chan][cut][sys]);
	//Rare .SSNBtagsNJets[chan][cut][sys] ->Add(   S[ZZZJets]               .SSNBtagsNJets[chan][cut][sys]);
	
	SetupDraw(Data .SSNBtagsNJets[chan][cut][sys], kBlack, NBTagsNJets); 
	SetupDraw(TTbar.SSNBtagsNJets[chan][cut][sys], 633   , NBTagsNJets);  
	SetupDraw(DY   .SSNBtagsNJets[chan][cut][sys], 852   , NBTagsNJets);  
	SetupDraw(VV   .SSNBtagsNJets[chan][cut][sys], 390   , NBTagsNJets);  
	SetupDraw(STop .SSNBtagsNJets[chan][cut][sys], 616   , NBTagsNJets);  
	SetupDraw(Fake .SSNBtagsNJets[chan][cut][sys], 413   , NBTagsNJets);  
	//SetupDraw(Rare .SSNBtagsNJets[chan][cut][sys], kYellow  , NBTagsNJets);  
      }
      
      if (gUseTTMadSpin) {
	TTbar.NBtagsNJets[chan][cut][Q2ScaleUp   ] = (TH1F*)S[TTJets_scaleup]     .NBtagsNJets[chan][cut][0]->Clone();
	TTbar.NBtagsNJets[chan][cut][Q2ScaleDown ] = (TH1F*)S[TTJets_scaledown]   .NBtagsNJets[chan][cut][0]->Clone();
	//TTbar.NBtagsNJets[chan][cut][MatchingUp  ] = (TH1F*)S[TTJets_matchingup]  .NBtagsNJets[chan][cut][0]->Clone();
	//TTbar.NBtagsNJets[chan][cut][MatchingDown] = (TH1F*)S[TTJets_matchingdown].NBtagsNJets[chan][cut][0]->Clone();
	
	SetupDraw(TTbar.NBtagsNJets[chan][cut][Q2ScaleUp   ], 633   , NBTagsNJets);  
	SetupDraw(TTbar.NBtagsNJets[chan][cut][Q2ScaleDown ], 633   , NBTagsNJets);  
	//SetupDraw(TTbar.NBtagsNJets[chan][cut][MatchingUp  ], 633   , NBTagsNJets);  
	//SetupDraw(TTbar.NBtagsNJets[chan][cut][MatchingDown], 633   , NBTagsNJets);  	

	//// SS 
	TTbar.SSNBtagsNJets[chan][cut][Q2ScaleUp   ]=(TH1F*)S[TTJets_scaleup]     .SSNBtagsNJets[chan][cut][0]->Clone();
	TTbar.SSNBtagsNJets[chan][cut][Q2ScaleDown ]=(TH1F*)S[TTJets_scaledown]   .SSNBtagsNJets[chan][cut][0]->Clone();
	//TTbar.SSNBtagsNJets[chan][cut][MatchingUp  ]=(TH1F*)S[TTJets_matchingup]  .SSNBtagsNJets[chan][cut][0]->Clone();
	//TTbar.SSNBtagsNJets[chan][cut][MatchingDown]=(TH1F*)S[TTJets_matchingdown].SSNBtagsNJets[chan][cut][0]->Clone();
	
	SetupDraw(TTbar.SSNBtagsNJets[chan][cut][Q2ScaleUp   ], 633   , NBTagsNJets);  
	SetupDraw(TTbar.SSNBtagsNJets[chan][cut][Q2ScaleDown ], 633   , NBTagsNJets);  
	//SetupDraw(TTbar.SSNBtagsNJets[chan][cut][MatchingUp  ], 633   , NBTagsNJets);  
	//SetupDraw(TTbar.SSNBtagsNJets[chan][cut][MatchingDown], 633   , NBTagsNJets);  	
      }
	
	
      for (size_t sys=0; sys<gNSYST; sys++){
        if (chan == Muon){
          Data .InvMass[chan][cut][sys]   = (TH1F*)S[DoubleMuon]        .InvMass[chan][cut][0]->Clone();
          Data .SSInvMass[chan][cut][sys] = (TH1F*)S[DoubleMuon]      .SSInvMass[chan][cut][0]->Clone();
        }
        else if (chan == Elec){
          Data .InvMass[chan][cut][sys]   = (TH1F*)S[DoubleEG]  .InvMass[chan][cut][0]->Clone();
          Data .SSInvMass[chan][cut][sys] = (TH1F*)S[DoubleEG].SSInvMass[chan][cut][0]->Clone();
        }
        else if (chan == ElMu){
          Data .InvMass[chan][cut][sys]   = (TH1F*)S[MuEG]	      .InvMass[chan][cut][0]->Clone();
          Data .SSInvMass[chan][cut][sys] = (TH1F*)S[MuEG]	    .SSInvMass[chan][cut][0]->Clone();
        }
        
        if (gUseTTMadSpin) {
          TTbar.InvMass[chan][cut][sys]   = (TH1F*)S[TTJets_MadSpin]	    .InvMass[chan][cut][sys]->Clone();
          TTbar.SSInvMass[chan][cut][sys] = (TH1F*)S[TTJets_MadSpin]	  .SSInvMass[chan][cut][sys]->Clone();
        }
        //else {
        //  TTbar.InvMass[chan][cut][sys]   = (TH1F*)S[TTJetsFullLeptMGtauola].InvMass[chan][cut][sys]->Clone();
        //  TTbar.SSInvMass[chan][cut][sys] = (TH1F*)S[TTJetsFullLeptMGtauola].SSInvMass[chan][cut][sys]->Clone();
        //}
        //SUSYstop .InvMass[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].InvMass[chan][cut][sys]->Clone();
        STop .InvMass[chan][cut][sys] = (TH1F*)S[TbarWDilep]		.InvMass[chan][cut][sys]->Clone();
        STop .InvMass[chan][cut][sys] ->Add(   S[TWDilep]		.InvMass[chan][cut][sys]);
        DY   .InvMass[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	.InvMass[chan][cut][sys]->Clone();
        DY   .InvMass[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	.InvMass[chan][cut][sys]);
        VV   .InvMass[chan][cut][sys] = (TH1F*)S[WZ]	.InvMass[chan][cut][sys]->Clone();
        VV   .InvMass[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]	.InvMass[chan][cut][sys]);
        VV   .InvMass[chan][cut][sys] ->Add(   S[ZZ]			.InvMass[chan][cut][sys]);
        Fake .InvMass[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .InvMass[chan][cut][sys]->Clone();  
        Fake .InvMass[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].InvMass[chan][cut][sys]);
        //Fake .InvMass[chan][cut][sys] ->Add(   S[WgammaToLNuG]  	.InvMass[chan][cut][sys]);
        //Rare .InvMass[chan][cut][sys] = (TH1F*)S[TTWJets]		.InvMass[chan][cut][sys]->Clone();  
        //	Rare .InvMass[chan][cut][sys] ->Add(   S[TTWJets]		.InvMass[chan][cut][sys]);
        //Rare .InvMass[chan][cut][sys] ->Add(   S[TTWWJets]		.InvMass[chan][cut][sys]);
        //Rare .InvMass[chan][cut][sys] ->Add(   S[TTZJets]		.InvMass[chan][cut][sys]);
        //Rare .InvMass[chan][cut][sys] ->Add(   S[WWWJets]		.InvMass[chan][cut][sys]);
        //Rare .InvMass[chan][cut][sys] ->Add(   S[WWZJets]		.InvMass[chan][cut][sys]);
       // Rare .InvMass[chan][cut][sys] ->Add(   S[WZZJets]		.InvMass[chan][cut][sys]);
        //Rare .InvMass[chan][cut][sys] ->Add(   S[ZZZJets]		.InvMass[chan][cut][sys]);
        
        SetupDraw(Data .InvMass[chan][cut][sys], kBlack, InvMass); 
        SetupDraw(TTbar.InvMass[chan][cut][sys], 633   , InvMass);  
        SetupDraw(DY   .InvMass[chan][cut][sys], 852   , InvMass);  
        SetupDraw(VV   .InvMass[chan][cut][sys], 390   , InvMass);  
        SetupDraw(STop .InvMass[chan][cut][sys], 616   , InvMass);  
        SetupDraw(Fake .InvMass[chan][cut][sys], 413   , InvMass);  
        //SetupDraw(Rare .InvMass[chan][cut][sys], kYellow  , InvMass);  
        
        //SUSYstop.SSInvMass[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].SSInvMass[chan][cut][sys]->Clone();
        STop .SSInvMass[chan][cut][sys] = (TH1F*)S[TbarWDilep]  	  .SSInvMass[chan][cut][sys]->Clone();
        STop .SSInvMass[chan][cut][sys] ->Add(   S[TWDilep]		  .SSInvMass[chan][cut][sys]);
        DY   .SSInvMass[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	  .SSInvMass[chan][cut][sys]->Clone();
        DY   .SSInvMass[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	  .SSInvMass[chan][cut][sys]);
        VV   .SSInvMass[chan][cut][sys] = (TH1F*)S[WZ]                    .SSInvMass[chan][cut][sys]->Clone();
        VV   .SSInvMass[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]    .SSInvMass[chan][cut][sys]);
        VV   .SSInvMass[chan][cut][sys] ->Add(   S[ZZ]  		  .SSInvMass[chan][cut][sys]);
        Fake .SSInvMass[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .SSInvMass[chan][cut][sys]->Clone();  
        Fake .SSInvMass[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].SSInvMass[chan][cut][sys]);
        //Fake .SSInvMass[chan][cut][sys] ->Add(   S[WgammaToLNuG]	  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] = (TH1F*)S[TTWJets]		  .SSInvMass[chan][cut][sys]->Clone();  
        //	Rare .SSInvMass[chan][cut][sys] ->Add(   S[TTWJets]		  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] ->Add(   S[TTWWJets]		  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] ->Add(   S[TTZJets]		  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] ->Add(   S[WWWJets]		  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] ->Add(   S[WWZJets]		  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] ->Add(   S[WZZJets]		  .SSInvMass[chan][cut][sys]);
        //Rare .SSInvMass[chan][cut][sys] ->Add(   S[ZZZJets]		  .SSInvMass[chan][cut][sys]);

        SetupDraw(Data .SSInvMass[chan][cut][sys], kBlack, InvMass); 
        SetupDraw(TTbar.SSInvMass[chan][cut][sys], 633   , InvMass);  
        SetupDraw(DY   .SSInvMass[chan][cut][sys], 852   , InvMass);  
        SetupDraw(VV   .SSInvMass[chan][cut][sys], 390   , InvMass);  
        SetupDraw(STop .SSInvMass[chan][cut][sys], 616   , InvMass);  
        SetupDraw(Fake .SSInvMass[chan][cut][sys], 413   , InvMass);  
        //SetupDraw(Rare .SSInvMass[chan][cut][sys], kYellow  , InvMass);  
      }
      if (gUseTTMadSpin) {
        TTbar.InvMass[chan][cut][Q2ScaleUp   ] = (TH1F*)S[TTJets_scaleup]     .InvMass[chan][cut][0]->Clone();
        TTbar.InvMass[chan][cut][Q2ScaleDown ] = (TH1F*)S[TTJets_scaledown]   .InvMass[chan][cut][0]->Clone();
        //TTbar.InvMass[chan][cut][MatchingUp  ] = (TH1F*)S[TTJets_matchingup]  .InvMass[chan][cut][0]->Clone();
        //TTbar.InvMass[chan][cut][MatchingDown] = (TH1F*)S[TTJets_matchingdown].InvMass[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.InvMass[chan][cut][Q2ScaleUp   ], 633   , InvMass);  
        SetupDraw(TTbar.InvMass[chan][cut][Q2ScaleDown ], 633   , InvMass);  
        //SetupDraw(TTbar.InvMass[chan][cut][MatchingUp  ], 633   , InvMass);  
        //SetupDraw(TTbar.InvMass[chan][cut][MatchingDown], 633   , InvMass);        
        
        //// SS 
        TTbar.SSInvMass[chan][cut][Q2ScaleUp   ]=(TH1F*)S[TTJets_scaleup]     .SSInvMass[chan][cut][0]->Clone();
        TTbar.SSInvMass[chan][cut][Q2ScaleDown ]=(TH1F*)S[TTJets_scaledown]   .SSInvMass[chan][cut][0]->Clone();
        //TTbar.SSInvMass[chan][cut][MatchingUp  ]=(TH1F*)S[TTJets_matchingup]  .SSInvMass[chan][cut][0]->Clone();
        //TTbar.SSInvMass[chan][cut][MatchingDown]=(TH1F*)S[TTJets_matchingdown].SSInvMass[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.SSInvMass[chan][cut][Q2ScaleUp   ], 633   , InvMass);  
        SetupDraw(TTbar.SSInvMass[chan][cut][Q2ScaleDown ], 633   , InvMass);  
        //SetupDraw(TTbar.SSInvMass[chan][cut][MatchingUp  ], 633   , InvMass);  
        //SetupDraw(TTbar.SSInvMass[chan][cut][MatchingDown], 633   , InvMass);      
      } 	
 	
	
      for (size_t sys=0; sys<gNSYST; sys++){
        if (chan == Muon){
          Data .AbsDelPhiLeps[chan][cut][sys]   = (TH1F*)S[DoubleMuon]        .AbsDelPhiLeps[chan][cut][0]->Clone();
          Data .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[DoubleMuon]      .SSAbsDelPhiLeps[chan][cut][0]->Clone();
        }
        else if (chan == Elec){
          Data .AbsDelPhiLeps[chan][cut][sys]   = (TH1F*)S[DoubleEG]  .AbsDelPhiLeps[chan][cut][0]->Clone();
          Data .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[DoubleEG].SSAbsDelPhiLeps[chan][cut][0]->Clone();
        }
        else if (chan == ElMu){
          Data .AbsDelPhiLeps[chan][cut][sys]   = (TH1F*)S[MuEG]	      .AbsDelPhiLeps[chan][cut][0]->Clone();
          Data .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[MuEG]	    .SSAbsDelPhiLeps[chan][cut][0]->Clone();
        }
        
        if (gUseTTMadSpin) {
          TTbar.AbsDelPhiLeps[chan][cut][sys]   = (TH1F*)S[TTJets_MadSpin]	    .AbsDelPhiLeps[chan][cut][sys]->Clone();
          TTbar.SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[TTJets_MadSpin]	  .SSAbsDelPhiLeps[chan][cut][sys]->Clone();
        }
        //else {
        //  TTbar.AbsDelPhiLeps[chan][cut][sys]   = (TH1F*)S[TTJetsFullLeptMGtauola].AbsDelPhiLeps[chan][cut][sys]->Clone();
        //  TTbar.SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[TTJetsFullLeptMGtauola].SSAbsDelPhiLeps[chan][cut][sys]->Clone();
        //}
        //SUSYstop.AbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].AbsDelPhiLeps[chan][cut][sys]->Clone();
        STop .AbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[TbarWDilep]		.AbsDelPhiLeps[chan][cut][sys]->Clone();
        STop .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[TWDilep]		.AbsDelPhiLeps[chan][cut][sys]);
        DY   .AbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	.AbsDelPhiLeps[chan][cut][sys]->Clone();
        DY   .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	.AbsDelPhiLeps[chan][cut][sys]);
        VV   .AbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[WZ]	.AbsDelPhiLeps[chan][cut][sys]->Clone();
        VV   .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]    .AbsDelPhiLeps[chan][cut][sys]);
        VV   .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[ZZ]		      .AbsDelPhiLeps[chan][cut][sys]);
        Fake .AbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .AbsDelPhiLeps[chan][cut][sys]->Clone();  
        Fake .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].AbsDelPhiLeps[chan][cut][sys]);
        //Fake .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[WgammaToLNuG]  	.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[TTWJets]		.AbsDelPhiLeps[chan][cut][sys]->Clone();  
        //	Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTWJets]		.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTWWJets]		.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTZJets]		.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[WWWJets]		.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[WWZJets]		.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[WZZJets]		.AbsDelPhiLeps[chan][cut][sys]);
        //Rare .AbsDelPhiLeps[chan][cut][sys] ->Add(   S[ZZZJets]		.AbsDelPhiLeps[chan][cut][sys]);
        
        SetupDraw(Data .AbsDelPhiLeps[chan][cut][sys], kBlack, AbsDelPhiLeps); 
        SetupDraw(TTbar.AbsDelPhiLeps[chan][cut][sys], 633   , AbsDelPhiLeps);  
        SetupDraw(DY   .AbsDelPhiLeps[chan][cut][sys], 852   , AbsDelPhiLeps);  
        SetupDraw(VV   .AbsDelPhiLeps[chan][cut][sys], 390   , AbsDelPhiLeps);  
        SetupDraw(STop .AbsDelPhiLeps[chan][cut][sys], 616   , AbsDelPhiLeps);  
        SetupDraw(Fake .AbsDelPhiLeps[chan][cut][sys], 413   , AbsDelPhiLeps);  
        //SetupDraw(Rare .AbsDelPhiLeps[chan][cut][sys], kYellow  , AbsDelPhiLeps);  
        
        //SUSYstop.SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].SSAbsDelPhiLeps[chan][cut][sys]->Clone();
        STop .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[TbarWDilep]  	  .SSAbsDelPhiLeps[chan][cut][sys]->Clone();
        STop .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[TWDilep]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        DY   .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	  .SSAbsDelPhiLeps[chan][cut][sys]->Clone();
        DY   .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	  .SSAbsDelPhiLeps[chan][cut][sys]);
        VV   .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[WZ]                      .SSAbsDelPhiLeps[chan][cut][sys]->Clone();
        VV   .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]  	  .SSAbsDelPhiLeps[chan][cut][sys]);
        VV   .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[ZZ]  		          .SSAbsDelPhiLeps[chan][cut][sys]);
        Fake .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]            .SSAbsDelPhiLeps[chan][cut][sys]->Clone();  
        Fake .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola]  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Fake .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[WgammaToLNuG]	  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] = (TH1F*)S[TTWJets]		  .SSAbsDelPhiLeps[chan][cut][sys]->Clone();  
        //	Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTWJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTWWJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[TTZJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[WWWJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[WWZJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[WZZJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);
        //Rare .SSAbsDelPhiLeps[chan][cut][sys] ->Add(   S[ZZZJets]		  .SSAbsDelPhiLeps[chan][cut][sys]);

        SetupDraw(Data .SSAbsDelPhiLeps[chan][cut][sys], kBlack, AbsDelPhiLeps); 
        SetupDraw(TTbar.SSAbsDelPhiLeps[chan][cut][sys], 633   , AbsDelPhiLeps);  
        SetupDraw(DY   .SSAbsDelPhiLeps[chan][cut][sys], 852   , AbsDelPhiLeps);  
        SetupDraw(VV   .SSAbsDelPhiLeps[chan][cut][sys], 390   , AbsDelPhiLeps);  
        SetupDraw(STop .SSAbsDelPhiLeps[chan][cut][sys], 616   , AbsDelPhiLeps);  
        SetupDraw(Fake .SSAbsDelPhiLeps[chan][cut][sys], 413   , AbsDelPhiLeps);  
        //SetupDraw(Rare .SSAbsDelPhiLeps[chan][cut][sys], kYellow  , AbsDelPhiLeps);  
      }
      if (gUseTTMadSpin) {
        TTbar.AbsDelPhiLeps[chan][cut][Q2ScaleUp   ] = (TH1F*)S[TTJets_scaleup]     .AbsDelPhiLeps[chan][cut][0]->Clone();
        TTbar.AbsDelPhiLeps[chan][cut][Q2ScaleDown ] = (TH1F*)S[TTJets_scaledown]   .AbsDelPhiLeps[chan][cut][0]->Clone();
        //TTbar.AbsDelPhiLeps[chan][cut][MatchingUp  ] = (TH1F*)S[TTJets_matchingup]  .AbsDelPhiLeps[chan][cut][0]->Clone();
        //TTbar.AbsDelPhiLeps[chan][cut][MatchingDown] = (TH1F*)S[TTJets_matchingdown].AbsDelPhiLeps[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.AbsDelPhiLeps[chan][cut][Q2ScaleUp   ], 633   , AbsDelPhiLeps);  
        SetupDraw(TTbar.AbsDelPhiLeps[chan][cut][Q2ScaleDown ], 633   , AbsDelPhiLeps);  
        //SetupDraw(TTbar.AbsDelPhiLeps[chan][cut][MatchingUp  ], 633   , AbsDelPhiLeps);  
        //SetupDraw(TTbar.AbsDelPhiLeps[chan][cut][MatchingDown], 633   , AbsDelPhiLeps);        
        
        //// SS 
        TTbar.SSAbsDelPhiLeps[chan][cut][Q2ScaleUp   ]=(TH1F*)S[TTJets_scaleup]     .SSAbsDelPhiLeps[chan][cut][0]->Clone();
        TTbar.SSAbsDelPhiLeps[chan][cut][Q2ScaleDown ]=(TH1F*)S[TTJets_scaledown]   .SSAbsDelPhiLeps[chan][cut][0]->Clone();
        //TTbar.SSAbsDelPhiLeps[chan][cut][MatchingUp  ]=(TH1F*)S[TTJets_matchingup]  .SSAbsDelPhiLeps[chan][cut][0]->Clone();
        //TTbar.SSAbsDelPhiLeps[chan][cut][MatchingDown]=(TH1F*)S[TTJets_matchingdown].SSAbsDelPhiLeps[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.SSAbsDelPhiLeps[chan][cut][Q2ScaleUp   ], 633   , AbsDelPhiLeps);  
        SetupDraw(TTbar.SSAbsDelPhiLeps[chan][cut][Q2ScaleDown ], 633   , AbsDelPhiLeps);  
        //SetupDraw(TTbar.SSAbsDelPhiLeps[chan][cut][MatchingUp  ], 633   , AbsDelPhiLeps);  
        //SetupDraw(TTbar.SSAbsDelPhiLeps[chan][cut][MatchingDown], 633   , AbsDelPhiLeps);      
      } 	
 	
	
      for (size_t sys=0; sys<gNSYST; sys++){
        if (chan == Muon){
          Data .delPhi2LeadJets[chan][cut][sys]   = (TH1F*)S[DoubleMuon]        .delPhi2LeadJets[chan][cut][0]->Clone();
          Data .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[DoubleMuon]      .SSdelPhi2LeadJets[chan][cut][0]->Clone();
        }
        else if (chan == Elec){
          Data .delPhi2LeadJets[chan][cut][sys]   = (TH1F*)S[DoubleEG]  .delPhi2LeadJets[chan][cut][0]->Clone();
          Data .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[DoubleEG].SSdelPhi2LeadJets[chan][cut][0]->Clone();
        }
        else if (chan == ElMu){
          Data .delPhi2LeadJets[chan][cut][sys]   = (TH1F*)S[MuEG]	      .delPhi2LeadJets[chan][cut][0]->Clone();
          Data .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[MuEG]	    .SSdelPhi2LeadJets[chan][cut][0]->Clone();
        }
        
        if (gUseTTMadSpin) {
          TTbar.delPhi2LeadJets[chan][cut][sys]   = (TH1F*)S[TTJets_MadSpin]	    .delPhi2LeadJets[chan][cut][sys]->Clone();
          TTbar.SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[TTJets_MadSpin]	  .SSdelPhi2LeadJets[chan][cut][sys]->Clone();
        }
        //else {
        //  TTbar.delPhi2LeadJets[chan][cut][sys]   = (TH1F*)S[TTJetsFullLeptMGtauola].delPhi2LeadJets[chan][cut][sys]->Clone();
        //  TTbar.SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[TTJetsFullLeptMGtauola].SSdelPhi2LeadJets[chan][cut][sys]->Clone();
        //}
        //SUSYstop.delPhi2LeadJets[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].delPhi2LeadJets[chan][cut][sys]->Clone();
        STop .delPhi2LeadJets[chan][cut][sys] = (TH1F*)S[TbarWDilep]		.delPhi2LeadJets[chan][cut][sys]->Clone();
        STop .delPhi2LeadJets[chan][cut][sys] ->Add(   S[TWDilep]		.delPhi2LeadJets[chan][cut][sys]);
        DY   .delPhi2LeadJets[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	.delPhi2LeadJets[chan][cut][sys]->Clone();
        DY   .delPhi2LeadJets[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	.delPhi2LeadJets[chan][cut][sys]);
        VV   .delPhi2LeadJets[chan][cut][sys] = (TH1F*)S[WZ]	                .delPhi2LeadJets[chan][cut][sys]->Clone();
        VV   .delPhi2LeadJets[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]	.delPhi2LeadJets[chan][cut][sys]);
        VV   .delPhi2LeadJets[chan][cut][sys] ->Add(   S[ZZ]			.delPhi2LeadJets[chan][cut][sys]);
        Fake .delPhi2LeadJets[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .delPhi2LeadJets[chan][cut][sys]->Clone();  
        Fake .delPhi2LeadJets[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].delPhi2LeadJets[chan][cut][sys]);
        //Fake .delPhi2LeadJets[chan][cut][sys] ->Add(   S[WgammaToLNuG]  	.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] = (TH1F*)S[TTWJets]		.delPhi2LeadJets[chan][cut][sys]->Clone();  
        //	Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[TTWJets]		.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[TTWWJets]		.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[TTZJets]		.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[WWWJets]		.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[WWZJets]		.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[WZZJets]		.delPhi2LeadJets[chan][cut][sys]);
        //Rare .delPhi2LeadJets[chan][cut][sys] ->Add(   S[ZZZJets]		.delPhi2LeadJets[chan][cut][sys]);
        
        SetupDraw(Data .delPhi2LeadJets[chan][cut][sys], kBlack, delPhi2LeadJets); 
        SetupDraw(TTbar.delPhi2LeadJets[chan][cut][sys], 633   , delPhi2LeadJets);  
        SetupDraw(DY   .delPhi2LeadJets[chan][cut][sys], 852   , delPhi2LeadJets);  
        SetupDraw(VV   .delPhi2LeadJets[chan][cut][sys], 390   , delPhi2LeadJets);  
        SetupDraw(STop .delPhi2LeadJets[chan][cut][sys], 616   , delPhi2LeadJets);  
        SetupDraw(Fake .delPhi2LeadJets[chan][cut][sys], 413   , delPhi2LeadJets);  
        //SetupDraw(Rare .delPhi2LeadJets[chan][cut][sys], kYellow  , delPhi2LeadJets);  
        
        //SUSYstop.SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].SSdelPhi2LeadJets[chan][cut][sys]->Clone();
        STop .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[TbarWDilep]  	  .SSdelPhi2LeadJets[chan][cut][sys]->Clone();
        STop .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[TWDilep]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        DY   .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	  .SSdelPhi2LeadJets[chan][cut][sys]->Clone();
        DY   .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	  .SSdelPhi2LeadJets[chan][cut][sys]);
        VV   .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[WZ]                    .SSdelPhi2LeadJets[chan][cut][sys]->Clone();
        VV   .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]    .SSdelPhi2LeadJets[chan][cut][sys]);
        VV   .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[ZZ]  		  .SSdelPhi2LeadJets[chan][cut][sys]);
        Fake .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .SSdelPhi2LeadJets[chan][cut][sys]->Clone();  
        Fake .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].SSdelPhi2LeadJets[chan][cut][sys]);
        //Fake .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[WgammaToLNuG]	  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] = (TH1F*)S[TTWJets]		  .SSdelPhi2LeadJets[chan][cut][sys]->Clone();  
        //	Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[TTWJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[TTWWJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[TTZJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[WWWJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[WWZJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[WZZJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);
        //Rare .SSdelPhi2LeadJets[chan][cut][sys] ->Add(   S[ZZZJets]		  .SSdelPhi2LeadJets[chan][cut][sys]);

        SetupDraw(Data .SSdelPhi2LeadJets[chan][cut][sys], kBlack, delPhi2LeadJets); 
        SetupDraw(TTbar.SSdelPhi2LeadJets[chan][cut][sys], 633   , delPhi2LeadJets);  
        SetupDraw(DY   .SSdelPhi2LeadJets[chan][cut][sys], 852   , delPhi2LeadJets);  
        SetupDraw(VV   .SSdelPhi2LeadJets[chan][cut][sys], 390   , delPhi2LeadJets);  
        SetupDraw(STop .SSdelPhi2LeadJets[chan][cut][sys], 616   , delPhi2LeadJets);  
        SetupDraw(Fake .SSdelPhi2LeadJets[chan][cut][sys], 413   , delPhi2LeadJets);  
        //SetupDraw(Rare .SSdelPhi2LeadJets[chan][cut][sys], kYellow  , delPhi2LeadJets);  
      }
      if (gUseTTMadSpin) {
        TTbar.delPhi2LeadJets[chan][cut][Q2ScaleUp   ] = (TH1F*)S[TTJets_scaleup]     .delPhi2LeadJets[chan][cut][0]->Clone();
        TTbar.delPhi2LeadJets[chan][cut][Q2ScaleDown ] = (TH1F*)S[TTJets_scaledown]   .delPhi2LeadJets[chan][cut][0]->Clone();
        //TTbar.delPhi2LeadJets[chan][cut][MatchingUp  ] = (TH1F*)S[TTJets_matchingup]  .delPhi2LeadJets[chan][cut][0]->Clone();
        //TTbar.delPhi2LeadJets[chan][cut][MatchingDown] = (TH1F*)S[TTJets_matchingdown].delPhi2LeadJets[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.delPhi2LeadJets[chan][cut][Q2ScaleUp   ], 633   , delPhi2LeadJets);  
        SetupDraw(TTbar.delPhi2LeadJets[chan][cut][Q2ScaleDown ], 633   , delPhi2LeadJets);  
        //SetupDraw(TTbar.delPhi2LeadJets[chan][cut][MatchingUp  ], 633   , delPhi2LeadJets);  
        //SetupDraw(TTbar.delPhi2LeadJets[chan][cut][MatchingDown], 633   , delPhi2LeadJets);        
        
        //// SS 
        TTbar.SSdelPhi2LeadJets[chan][cut][Q2ScaleUp   ]=(TH1F*)S[TTJets_scaleup]     .SSdelPhi2LeadJets[chan][cut][0]->Clone();
        TTbar.SSdelPhi2LeadJets[chan][cut][Q2ScaleDown ]=(TH1F*)S[TTJets_scaledown]   .SSdelPhi2LeadJets[chan][cut][0]->Clone();
        //TTbar.SSdelPhi2LeadJets[chan][cut][MatchingUp  ]=(TH1F*)S[TTJets_matchingup]  .SSdelPhi2LeadJets[chan][cut][0]->Clone();
        //TTbar.SSdelPhi2LeadJets[chan][cut][MatchingDown]=(TH1F*)S[TTJets_matchingdown].SSdelPhi2LeadJets[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.SSdelPhi2LeadJets[chan][cut][Q2ScaleUp   ], 633   , delPhi2LeadJets);  
        SetupDraw(TTbar.SSdelPhi2LeadJets[chan][cut][Q2ScaleDown ], 633   , delPhi2LeadJets);  
        //SetupDraw(TTbar.SSdelPhi2LeadJets[chan][cut][MatchingUp  ], 633   , delPhi2LeadJets);  
        //SetupDraw(TTbar.SSdelPhi2LeadJets[chan][cut][MatchingDown], 633   , delPhi2LeadJets);      
      } 	
 	
	
      for (size_t sys=0; sys<gNSYST; sys++){
        if (chan == Muon){
          Data .minDelRJetsLeps[chan][cut][sys]   = (TH1F*)S[DoubleMuon]        .minDelRJetsLeps[chan][cut][0]->Clone();
          Data .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[DoubleMuon]      .SSminDelRJetsLeps[chan][cut][0]->Clone();
        }
        else if (chan == Elec){
          Data .minDelRJetsLeps[chan][cut][sys]   = (TH1F*)S[DoubleEG]  .minDelRJetsLeps[chan][cut][0]->Clone();
          Data .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[DoubleEG].SSminDelRJetsLeps[chan][cut][0]->Clone();
        }
        else if (chan == ElMu){
          Data .minDelRJetsLeps[chan][cut][sys]   = (TH1F*)S[MuEG]	      .minDelRJetsLeps[chan][cut][0]->Clone();
          Data .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[MuEG]	    .SSminDelRJetsLeps[chan][cut][0]->Clone();
        }
        
        if (gUseTTMadSpin) {
          TTbar.minDelRJetsLeps[chan][cut][sys]   = (TH1F*)S[TTJets_MadSpin]	    .minDelRJetsLeps[chan][cut][sys]->Clone();
          TTbar.SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[TTJets_MadSpin]	  .SSminDelRJetsLeps[chan][cut][sys]->Clone();
        }
        //else {
        //  TTbar.minDelRJetsLeps[chan][cut][sys]   = (TH1F*)S[TTJetsFullLeptMGtauola].minDelRJetsLeps[chan][cut][sys]->Clone();
        //  TTbar.SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[TTJetsFullLeptMGtauola].SSminDelRJetsLeps[chan][cut][sys]->Clone();
        //}
        //SUSYstop.minDelRJetsLeps[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].minDelRJetsLeps[chan][cut][sys]->Clone();
        STop .minDelRJetsLeps[chan][cut][sys] = (TH1F*)S[TbarWDilep]		.minDelRJetsLeps[chan][cut][sys]->Clone();
        STop .minDelRJetsLeps[chan][cut][sys] ->Add(   S[TWDilep]		.minDelRJetsLeps[chan][cut][sys]);
        DY   .minDelRJetsLeps[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	.minDelRJetsLeps[chan][cut][sys]->Clone();
        DY   .minDelRJetsLeps[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	.minDelRJetsLeps[chan][cut][sys]);
        VV   .minDelRJetsLeps[chan][cut][sys] = (TH1F*)S[WZ]	                .minDelRJetsLeps[chan][cut][sys]->Clone();
        VV   .minDelRJetsLeps[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]	.minDelRJetsLeps[chan][cut][sys]);
        VV   .minDelRJetsLeps[chan][cut][sys] ->Add(   S[ZZ]			.minDelRJetsLeps[chan][cut][sys]);
        Fake .minDelRJetsLeps[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .minDelRJetsLeps[chan][cut][sys]->Clone();  
        Fake .minDelRJetsLeps[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].minDelRJetsLeps[chan][cut][sys]);
        //Fake .minDelRJetsLeps[chan][cut][sys] ->Add(   S[WgammaToLNuG]  	.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] = (TH1F*)S[TTWJets]		.minDelRJetsLeps[chan][cut][sys]->Clone();  
        //	Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[TTWJets]		.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[TTWWJets]		.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[TTZJets]		.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[WWWJets]		.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[WWZJets]		.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[WZZJets]		.minDelRJetsLeps[chan][cut][sys]);
        //Rare .minDelRJetsLeps[chan][cut][sys] ->Add(   S[ZZZJets]		.minDelRJetsLeps[chan][cut][sys]);
        
        SetupDraw(Data .minDelRJetsLeps[chan][cut][sys], kBlack, minDelRJetsLeps); 
        SetupDraw(TTbar.minDelRJetsLeps[chan][cut][sys], 633   , minDelRJetsLeps);  
        SetupDraw(DY   .minDelRJetsLeps[chan][cut][sys], 852   , minDelRJetsLeps);  
        SetupDraw(VV   .minDelRJetsLeps[chan][cut][sys], 390   , minDelRJetsLeps);  
        SetupDraw(STop .minDelRJetsLeps[chan][cut][sys], 616   , minDelRJetsLeps);  
        SetupDraw(Fake .minDelRJetsLeps[chan][cut][sys], 413   , minDelRJetsLeps);  
        //SetupDraw(Rare .minDelRJetsLeps[chan][cut][sys], kYellow  , minDelRJetsLeps);  
        
        //SUSYstop.SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[T2tt_150to250LSP1to100_LeptonFilter].SSminDelRJetsLeps[chan][cut][sys]->Clone();
        STop .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[TbarWDilep]  	  .SSminDelRJetsLeps[chan][cut][sys]->Clone();
        STop .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[TWDilep]		  .SSminDelRJetsLeps[chan][cut][sys]);
        DY   .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[ZJets_Madgraph]	  .SSminDelRJetsLeps[chan][cut][sys]->Clone();
        DY   .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[DYJets_Madgraph]	  .SSminDelRJetsLeps[chan][cut][sys]);
        VV   .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[WZ]                    .SSminDelRJetsLeps[chan][cut][sys]->Clone();
        VV   .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[WWTo2L2Nu_Madgraph]    .SSminDelRJetsLeps[chan][cut][sys]);
        VV   .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[ZZ]  		  .SSminDelRJetsLeps[chan][cut][sys]);
        Fake .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[Wbb_Madgraph]          .SSminDelRJetsLeps[chan][cut][sys]->Clone();  
        Fake .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[TTJetsSemiLeptMGtauola].SSminDelRJetsLeps[chan][cut][sys]);
        //Fake .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[WgammaToLNuG]	  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] = (TH1F*)S[TTWJets]		  .SSminDelRJetsLeps[chan][cut][sys]->Clone();  
        //	Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[TTWJets]		  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[TTWWJets]		  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[TTZJets]		  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[WWWJets]		  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[WWZJets]		  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[WZZJets]		  .SSminDelRJetsLeps[chan][cut][sys]);
        //Rare .SSminDelRJetsLeps[chan][cut][sys] ->Add(   S[ZZZJets]		  .SSminDelRJetsLeps[chan][cut][sys]);

        SetupDraw(Data .SSminDelRJetsLeps[chan][cut][sys], kBlack, minDelRJetsLeps); 
        SetupDraw(TTbar.SSminDelRJetsLeps[chan][cut][sys], 633   , minDelRJetsLeps);  
        SetupDraw(DY   .SSminDelRJetsLeps[chan][cut][sys], 852   , minDelRJetsLeps);  
        SetupDraw(VV   .SSminDelRJetsLeps[chan][cut][sys], 390   , minDelRJetsLeps);  
        SetupDraw(STop .SSminDelRJetsLeps[chan][cut][sys], 616   , minDelRJetsLeps);  
        SetupDraw(Fake .SSminDelRJetsLeps[chan][cut][sys], 413   , minDelRJetsLeps);  
        //SetupDraw(Rare .SSminDelRJetsLeps[chan][cut][sys], kYellow  , minDelRJetsLeps);  
      }

      if (gUseTTMadSpin) {
        TTbar.minDelRJetsLeps[chan][cut][Q2ScaleUp   ] = (TH1F*)S[TTJets_scaleup]     .minDelRJetsLeps[chan][cut][0]->Clone();
        TTbar.minDelRJetsLeps[chan][cut][Q2ScaleDown ] = (TH1F*)S[TTJets_scaledown]   .minDelRJetsLeps[chan][cut][0]->Clone();
        //TTbar.minDelRJetsLeps[chan][cut][MatchingUp  ] = (TH1F*)S[TTJets_matchingup]  .minDelRJetsLeps[chan][cut][0]->Clone();
        //TTbar.minDelRJetsLeps[chan][cut][MatchingDown] = (TH1F*)S[TTJets_matchingdown].minDelRJetsLeps[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.minDelRJetsLeps[chan][cut][Q2ScaleUp   ], 633   , minDelRJetsLeps);  
        SetupDraw(TTbar.minDelRJetsLeps[chan][cut][Q2ScaleDown ], 633   , minDelRJetsLeps);  
        //SetupDraw(TTbar.minDelRJetsLeps[chan][cut][MatchingUp  ], 633   , minDelRJetsLeps);  
        //SetupDraw(TTbar.minDelRJetsLeps[chan][cut][MatchingDown], 633   , minDelRJetsLeps);        
        
        //// SS 
        TTbar.SSminDelRJetsLeps[chan][cut][Q2ScaleUp   ]=(TH1F*)S[TTJets_scaleup]     .SSminDelRJetsLeps[chan][cut][0]->Clone();
        TTbar.SSminDelRJetsLeps[chan][cut][Q2ScaleDown ]=(TH1F*)S[TTJets_scaledown]   .SSminDelRJetsLeps[chan][cut][0]->Clone();
        //TTbar.SSminDelRJetsLeps[chan][cut][MatchingUp  ]=(TH1F*)S[TTJets_matchingup]  .SSminDelRJetsLeps[chan][cut][0]->Clone();
        //TTbar.SSminDelRJetsLeps[chan][cut][MatchingDown]=(TH1F*)S[TTJets_matchingdown].SSminDelRJetsLeps[chan][cut][0]->Clone();
        
        SetupDraw(TTbar.SSminDelRJetsLeps[chan][cut][Q2ScaleUp   ], 633   , minDelRJetsLeps);  
        SetupDraw(TTbar.SSminDelRJetsLeps[chan][cut][Q2ScaleDown ], 633   , minDelRJetsLeps);  
        //SetupDraw(TTbar.SSminDelRJetsLeps[chan][cut][MatchingUp  ], 633   , minDelRJetsLeps);  
        //SetupDraw(TTbar.SSminDelRJetsLeps[chan][cut][MatchingDown], 633   , minDelRJetsLeps);      
      } 	
    }
  }
  cout << "DONE! " << endl;
  
}
void TopPlotter::ResetDataMembers(){
  
  // Resetting samples holder...
  for (size_t sample=0; sample<gNSAMPLES; sample++){
    for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t cut=0; cut<iNCUTS; cut++){
	S[sample].Yields     [chan][cut] = 0.;
	S[sample].Yields_stat[chan][cut] = 0.;
	S[sample].SSYields   [chan][cut] = 0.;
	S[sample].SSYields_stat[chan][cut] = 0.;

	for (size_t sys=0; sys<gNSYST; sys++){
	  S[sample].Yields_syst[chan][cut][sys] = 0.;
	}
      }
      for (size_t sys=0; sys<gNSYSTERRTypes; sys++){
	S[sample].SystError[chan][sys] = 0.;
      }
    }
  }
  
  // Resetting background / signal categories...
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    for (size_t cut=0; cut<iNCUTS; cut++){
      Data .  Yields[chan][cut] = 0.;    Data .  Yields_stat[chan][cut] = 0.;
      TTbar.  Yields[chan][cut] = 0.;	 TTbar.  Yields_stat[chan][cut] = 0.;
      STop .  Yields[chan][cut] = 0.;	 STop .  Yields_stat[chan][cut] = 0.;
      DY   .  Yields[chan][cut] = 0.;	 DY   .  Yields_stat[chan][cut] = 0.;
      VV   .  Yields[chan][cut] = 0.;	 VV   .  Yields_stat[chan][cut] = 0.;
      Rare .  Yields[chan][cut] = 0.;	 Rare .  Yields_stat[chan][cut] = 0.;
      Fake .  Yields[chan][cut] = 0.;	 Fake .  Yields_stat[chan][cut] = 0.;
      Total.  Yields[chan][cut] = 0.;	 Total.  Yields_stat[chan][cut] = 0.;
      SUSYstop.Yields[chan][cut] = 0.;	 SUSYstop.  Yields_stat[chan][cut] = 0.;

      for (size_t sys=0; sys<gNSYST; sys++){
	Data .  Yields_syst[chan][cut][sys] = 0.;    
	TTbar.  Yields_syst[chan][cut][sys] = 0.;    
	STop .  Yields_syst[chan][cut][sys] = 0.;    
	DY   .  Yields_syst[chan][cut][sys] = 0.;    
	VV   .  Yields_syst[chan][cut][sys] = 0.;    
	Rare .  Yields_syst[chan][cut][sys] = 0.;    
	Fake .  Yields_syst[chan][cut][sys] = 0.;    
	Total.  Yields_syst[chan][cut][sys] = 0.;    
	SUSYstop.Yields_syst[chan][cut][sys] = 0.;    

	DD_DY  .Yields_syst[chan][cut][sys] = 0.;    
	DD_NonW.Yields_syst[chan][cut][sys] = 0.;    
      }

      Data .SSYields[chan][cut] = 0.; Data .SSYields_stat[chan][cut] = 0.;
      TTbar.SSYields[chan][cut] = 0.; TTbar.SSYields_stat[chan][cut] = 0.;
      STop .SSYields[chan][cut] = 0.; STop .SSYields_stat[chan][cut] = 0.;
      DY   .SSYields[chan][cut] = 0.; DY   .SSYields_stat[chan][cut] = 0.;
      VV   .SSYields[chan][cut] = 0.; VV   .SSYields_stat[chan][cut] = 0.;
      Rare .SSYields[chan][cut] = 0.; Rare .SSYields_stat[chan][cut] = 0.;
      Fake .SSYields[chan][cut] = 0.; Fake .SSYields_stat[chan][cut] = 0.;
      Total.SSYields[chan][cut] = 0.; Total.SSYields_stat[chan][cut] = 0.;

      // Data Driven
      DD_DY  .  Yields[chan][cut] = 0.;   DD_DY  .Yields_stat[chan][cut] = 0.;
      DD_NonW.  Yields[chan][cut] = 0.;	  DD_NonW.Yields_stat[chan][cut] = 0.;
      DD_DY  .SSYields[chan][cut] = 0.;
      DD_NonW.SSYields[chan][cut] = 0.;
    }
    for (size_t sys=0; sys<gNSYSTERRTypesALL; sys++){
      Data .SystError[chan][sys] = 0.;
      TTbar.SystError[chan][sys] = 0.;
      STop .SystError[chan][sys] = 0.;
      DY   .SystError[chan][sys] = 0.;
      VV   .SystError[chan][sys] = 0.;
      Rare .SystError[chan][sys] = 0.;
      Fake .SystError[chan][sys] = 0.;
      Total.SystError[chan][sys] = 0.;
      SUSYstop.SystError[chan][sys] = 0.;
    }
  }
  
  // Reset XSection...
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    ttbar.xsec     [ch] = 0.;
    ttbar.xsec_syst[ch] = 0.;
    ttbar.xsec_stat[ch] = 0.;
    ttbar.xsec_lumi[ch] = 0.;
    ttbar.acc      [ch] = 0.;
    ttbar.acc_stat [ch] = 0.;
    ttbar.acc_syst [ch] = 0.;

    ttbar.err_VV   [ch] = 0.;
    ttbar.err_DY   [ch] = 0.;
    ttbar.err_STop [ch] = 0.;
    ttbar.err_Fake [ch] = 0.;
    ttbar.err_Rare [ch] = 0.;
    ttbar.err_IDIso[ch] = 0.;
    ttbar.err_Trig [ch] = 0.;
    ttbar.err_LES  [ch] = 0.;
    ttbar.err_JES  [ch] = 0.;
    ttbar.err_JER  [ch] = 0.;
    ttbar.err_Btag [ch] = 0.;
    ttbar.err_PU   [ch] = 0.;
    ttbar.err_Q2   [ch] = 0.;
    ttbar.err_Match[ch] = 0.;
  }
}
void TopPlotter::ResetSystematicErrors(){
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    for (size_t sys=les; sys<gNSYSTERRTypesALL; sys++){
      //      Data .SystError[chan][sys] = 0.;
      TTbar.SystError[chan][sys] = 0.;
      STop .SystError[chan][sys] = 0.;
      DY   .SystError[chan][sys] = 0.;
      VV   .SystError[chan][sys] = 0.;
      Rare .SystError[chan][sys] = 0.;
      Fake .SystError[chan][sys] = 0.;
      Total.SystError[chan][sys] = 0.;
      SUSYstop.SystError[chan][sys] = 0.;
    }
    
    for (size_t cut=0; cut<iNCUTS; cut++){
      //      Data .  Yields_syst[chan][cut][0] = 0.;    
      TTbar.  Yields_syst[chan][cut][0] = 0.;    
      STop .  Yields_syst[chan][cut][0] = 0.;    
      DY   .  Yields_syst[chan][cut][0] = 0.;    
      VV   .  Yields_syst[chan][cut][0] = 0.;    
      Rare .  Yields_syst[chan][cut][0] = 0.;    
      Fake .  Yields_syst[chan][cut][0] = 0.;    
      Total.  Yields_syst[chan][cut][0] = 0.;    
      SUSYstop  .Yields_syst[chan][cut][0] = 0.;    
    }
  }
  // Reset XSection...
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    ttbar.xsec     [ch] = 0.;
    ttbar.xsec_syst[ch] = 0.;
    ttbar.xsec_stat[ch] = 0.;
    ttbar.xsec_lumi[ch] = 0.;
    ttbar.acc      [ch] = 0.;
    ttbar.acc_stat [ch] = 0.;
    ttbar.acc_syst [ch] = 0.;

    ttbar.err_VV   [ch] = 0.;
    ttbar.err_DY   [ch] = 0.;
    ttbar.err_STop [ch] = 0.;
    ttbar.err_Fake [ch] = 0.;
    ttbar.err_Rare [ch] = 0.;
    ttbar.err_IDIso[ch] = 0.;
    ttbar.err_Trig [ch] = 0.;
    ttbar.err_LES  [ch] = 0.;
    ttbar.err_JES  [ch] = 0.;
    ttbar.err_JER  [ch] = 0.;
    ttbar.err_Btag [ch] = 0.;
    ttbar.err_PU   [ch] = 0.;
    ttbar.err_Q2   [ch] = 0.;
    ttbar.err_Match[ch] = 0.;
  }

}
void TopPlotter::PrintSystematicErrors(){
  fOutputSubDir = "XSection/";
  
  TString filename = "";
  if (gUseTTMadSpin) filename = fOutputDir + fOutputSubDir + Form("SystErrors_MadSpin.txt");
  else               filename = fOutputDir + fOutputSubDir + Form("SystErrors_Madgraph.txt");
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  fOUTSTREAM2.open(filename, ios::trunc);

  fOUTSTREAM2 << "===========================================================================" << endl;
  fOUTSTREAM2 << "                      Systematic Uncertainties in pb (%)                   " << endl;
  fOUTSTREAM2 << " Source                |      El/El      |      Mu/Mu     |     El/Mu      " << endl;
  fOUTSTREAM2 << "---------------------------------------------------------------------------" << endl;
  fOUTSTREAM2 << Form(" VV                    | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_VV[Elec],  100 * ttbar.err_VV[Elec] / ttbar.xsec[Elec],
		      ttbar.err_VV[Muon],  100 * ttbar.err_VV[Muon] / ttbar.xsec[Muon],
		      ttbar.err_VV[ElMu],  100 * ttbar.err_VV[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" DY                    | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_DY[Elec],  100 * ttbar.err_DY[Elec] / ttbar.xsec[Elec],
		      ttbar.err_DY[Muon],  100 * ttbar.err_DY[Muon] / ttbar.xsec[Muon],
		      ttbar.err_DY[ElMu],  100 * ttbar.err_DY[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" STop                  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_STop[Elec],  100 * ttbar.err_STop[Elec] / ttbar.xsec[Elec],
		      ttbar.err_STop[Muon],  100 * ttbar.err_STop[Muon] / ttbar.xsec[Muon],
		      ttbar.err_STop[ElMu],  100 * ttbar.err_STop[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Fake                  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Fake[Elec],  100 * ttbar.err_Fake[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Fake[Muon],  100 * ttbar.err_Fake[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Fake[ElMu],  100 * ttbar.err_Fake[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Rare                  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Rare[Elec],  100 * ttbar.err_Rare[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Rare[Muon],  100 * ttbar.err_Rare[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Rare[ElMu],  100 * ttbar.err_Rare[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Lepton Efficiencies   | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_IDIso[Elec],  100 * ttbar.err_IDIso[Elec] / ttbar.xsec[Elec],
		      ttbar.err_IDIso[Muon],  100 * ttbar.err_IDIso[Muon] / ttbar.xsec[Muon],
		      ttbar.err_IDIso[ElMu],  100 * ttbar.err_IDIso[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Trigger Efficiencies  | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Trig[Elec],  100 * ttbar.err_Trig[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Trig[Muon],  100 * ttbar.err_Trig[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Trig[ElMu],  100 * ttbar.err_Trig[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Lepton Energy Scale   | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_LES[Elec],  100 * ttbar.err_LES[Elec] / ttbar.xsec[Elec],
		      ttbar.err_LES[Muon],  100 * ttbar.err_LES[Muon] / ttbar.xsec[Muon],
		      ttbar.err_LES[ElMu],  100 * ttbar.err_LES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Jet Energy Scale      | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_JES[Elec],  100 * ttbar.err_JES[Elec] / ttbar.xsec[Elec],
		      ttbar.err_JES[Muon],  100 * ttbar.err_JES[Muon] / ttbar.xsec[Muon],
		      ttbar.err_JES[ElMu],  100 * ttbar.err_JES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Jet Energy Resolution | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_JER[Elec],  100 * ttbar.err_JER[Elec] / ttbar.xsec[Elec],
		      ttbar.err_JER[Muon],  100 * ttbar.err_JER[Muon] / ttbar.xsec[Muon],
		      ttbar.err_JER[ElMu],  100 * ttbar.err_JER[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" b-tagging             | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Btag[Elec],  100 * ttbar.err_Btag[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Btag[Muon],  100 * ttbar.err_Btag[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Btag[ElMu],  100 * ttbar.err_Btag[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Pile Up               | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_PU[Elec],  100 * ttbar.err_PU[Elec] / ttbar.xsec[Elec],
		      ttbar.err_PU[Muon],  100 * ttbar.err_PU[Muon] / ttbar.xsec[Muon],
		      ttbar.err_PU[ElMu],  100 * ttbar.err_PU[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" QCD scale             | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Q2[Elec],  100 * ttbar.err_Q2[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Q2[Muon],  100 * ttbar.err_Q2[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Q2[ElMu],  100 * ttbar.err_Q2[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2 << Form(" Matching partons      | %4.2f (%3.2f) | %4.2f (%3.2f) | %4.2f (%3.2f) ", 
		      ttbar.err_Match[Elec],  100 * ttbar.err_Match[Elec] / ttbar.xsec[Elec],
		      ttbar.err_Match[Muon],  100 * ttbar.err_Match[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Match[ElMu],  100 * ttbar.err_Match[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM2.close();
  
  gSystem->Exec("cat "+filename);
}
void TopPlotter::PrintYieldsWithDD(){
  if (gUseTTMadSpin) fOutputSubDir = "Yields_DD_Madspin/";
  else               fOutputSubDir = "Yields_DD_Madgraph/"; 

  TString yieldsfilename = "";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  for (size_t cut=0; cut<iNCUTS; cut++){
    yieldsfilename = fOutputDir + fOutputSubDir + "Yields_MC_"+sCut[cut]+".txt";
    fOUTSTREAM.open(yieldsfilename, ios::trunc);

    fOUTSTREAM << "///////////////////////////////////////////////////////////////////////////\\\\" << endl;
    fOUTSTREAM << " Producing MC yields for cut " << sCut[cut] << endl;
    fOUTSTREAM << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
    
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << "                YIELDS  |    El/El   |    Mu/Mu   |    El/Mu   ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    
    float sum_mm(0.), sum_em(0.), sum_ee(0.);
    for (size_t sample=TTJets_MadSpin; sample<gNSAMPLES; sample++){
      //if ( gUseTTMadSpin && sample == TTJetsFullLeptMGtauola) continue;
      if (!gUseTTMadSpin && sample == TTJets_MadSpin        ) continue;
      
      float temp_mm  = S[sample].Yields[Muon][cut]; sum_mm  += temp_mm ;
      float temp_em  = S[sample].Yields[ElMu][cut]; sum_em  += temp_em ;
      float temp_ee  = S[sample].Yields[Elec][cut]; sum_ee  += temp_ee ;
      
      fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", SampleName[sample].Data(),
			 temp_ee , temp_mm , temp_em );
    }
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", "MC Sum", 
		       sum_ee , sum_mm , sum_em );
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %5.0f     |  %5.0f     |  %5.0f     ||\n", "Data", 
		       S[DoubleEG].Yields[Elec][cut], S[DoubleMuon].Yields[Muon][cut], S[MuEG].Yields[ElMu][cut]); 
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    // Now POST a Summary...
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "                            SUMMARY (stat error only):                       " << endl; 
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "          Source  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||", 
		       DD_DY.Yields[Elec][cut], DD_DY.Yields_stat[Elec][cut], 
		       DD_DY.Yields[Muon][cut], DD_DY.Yields_stat[Muon][cut], 
		       DD_DY.Yields[ElMu][cut], DD_DY.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       DD_NonW.Yields[Elec][cut], DD_NonW.Yields_stat[Elec][cut], 
		       DD_NonW.Yields[Muon][cut], DD_NonW.Yields_stat[Muon][cut], 
		       DD_NonW.Yields[ElMu][cut], DD_NonW.Yields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" Single top quark | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       STop.Yields[Elec][cut], STop.Yields_stat[Elec][cut],
		       STop.Yields[Muon][cut], STop.Yields_stat[Muon][cut], 
		       STop.Yields[ElMu][cut], STop.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Dibosons         | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       VV.Yields[Elec][cut], VV.Yields_stat[Elec][cut], 
		       VV.Yields[Muon][cut], VV.Yields_stat[Muon][cut], 
		       VV.Yields[ElMu][cut], VV.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Rare.Yields[Elec][cut], Rare.Yields_stat[Elec][cut], 
		       Rare.Yields[Muon][cut], Rare.Yields_stat[Muon][cut], 
		       Rare.Yields[ElMu][cut], Rare.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Total background | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Total.Yields[Elec][cut], Total.Yields_stat[Elec][cut], 
		       Total.Yields[Muon][cut], Total.Yields_stat[Muon][cut], 
		       Total.Yields[ElMu][cut], Total.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       TTbar.Yields[Elec][cut], TTbar.Yields_stat[Elec][cut], 
		       TTbar.Yields[Muon][cut], TTbar.Yields_stat[Muon][cut], 
		       TTbar.Yields[ElMu][cut], TTbar.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Data             | %5.0f            | %5.0f            | %5.0f            ||", 
		       Data.Yields[Elec][cut], 
		       Data.Yields[Muon][cut], 
		       Data.Yields[ElMu][cut] )<< endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    fOUTSTREAM.close();
  }
}
void TopPlotter::PrintYieldsWithMC(){
  if (gUseTTMadSpin) fOutputSubDir = "Yields_MC_Madspin/";
  else               fOutputSubDir = "Yields_MC_Madgraph/"; 

  TString yieldsfilename = "";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  for (size_t cut=0; cut<iNCUTS; cut++){
    yieldsfilename = fOutputDir + fOutputSubDir + "Yields_MC_"+sCut[cut]+".txt";
    fOUTSTREAM.open(yieldsfilename, ios::trunc);

    fOUTSTREAM << "///////////////////////////////////////////////////////////////////////////\\\\" << endl;
    fOUTSTREAM << " Producing MC yields for cut " << sCut[cut] << endl;
    fOUTSTREAM << "  scaling MC to " << fLumiNorm << " /pb" << endl << endl;
    
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << "                YIELDS  |    El/El   |    Mu/Mu   |    El/Mu   ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    
    float sum_mm(0.), sum_em(0.), sum_ee(0.);
    for (size_t sample=TTJets_MadSpin; sample<gNSAMPLES; sample++){
      //if ( gUseTTMadSpin && sample == TTJetsFullLeptMGtauola) continue;
      if (!gUseTTMadSpin && sample == TTJets_MadSpin        ) continue;
      
      float temp_mm  = S[sample].Yields[Muon][cut]; sum_mm  += temp_mm ;
      float temp_em  = S[sample].Yields[ElMu][cut]; sum_em  += temp_em ;
      float temp_ee  = S[sample].Yields[Elec][cut]; sum_ee  += temp_ee ;
      
      fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", SampleName[sample].Data(),
			 temp_ee , temp_mm , temp_em );
    }
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %8.2f  |  %8.2f  |  %8.2f  ||\n", "MC Sum", 
		       sum_ee , sum_mm , sum_em );
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << Form("%23s |  %5.0f     |  %5.0f     |  %5.0f     ||\n", "Data", 
		       S[DoubleEG].Yields[Elec][cut], S[DoubleMuon].Yields[Muon][cut], S[MuEG].Yields[ElMu][cut]); 
    fOUTSTREAM << "-----------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    // Now POST a Summary...
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "                            SUMMARY (stat error only):                       " << endl; 
    fOUTSTREAM << "=============================================================================" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "          Source  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Drell-Yan        | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||", 
		       DY.Yields[Elec][cut], DY.Yields_stat[Elec][cut], 
		       DY.Yields[Muon][cut], DY.Yields_stat[Muon][cut], 
		       DY.Yields[ElMu][cut], DY.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Fake.Yields[Elec][cut], Fake.Yields_stat[Elec][cut], 
		       Fake.Yields[Muon][cut], Fake.Yields_stat[Muon][cut], 
		       Fake.Yields[ElMu][cut], Fake.Yields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" Single top quark | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       STop.Yields[Elec][cut], STop.Yields_stat[Elec][cut],
		       STop.Yields[Muon][cut], STop.Yields_stat[Muon][cut], 
		       STop.Yields[ElMu][cut], STop.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Dibosons         | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       VV.Yields[Elec][cut], VV.Yields_stat[Elec][cut], 
		       VV.Yields[Muon][cut], VV.Yields_stat[Muon][cut], 
		       VV.Yields[ElMu][cut], VV.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Rare.Yields[Elec][cut], Rare.Yields_stat[Elec][cut], 
		       Rare.Yields[Muon][cut], Rare.Yields_stat[Muon][cut], 
		       Rare.Yields[ElMu][cut], Rare.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Total background | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       Total.Yields[Elec][cut], Total.Yields_stat[Elec][cut], 
		       Total.Yields[Muon][cut], Total.Yields_stat[Muon][cut], 
		       Total.Yields[ElMu][cut], Total.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %7.1f +/- %4.1f ||",  
		       TTbar.Yields[Elec][cut], TTbar.Yields_stat[Elec][cut], 
		       TTbar.Yields[Muon][cut], TTbar.Yields_stat[Muon][cut], 
		       TTbar.Yields[ElMu][cut], TTbar.Yields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Data             | %5.0f            | %5.0f            | %5.0f            ||", 
		       Data.Yields[Elec][cut], 
		       Data.Yields[Muon][cut], 
		       Data.Yields[ElMu][cut] )<< endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    fOUTSTREAM.close();
  }
}

void TopPlotter::DrawKinematicPlots(Bool_t DD, Int_t onechan, Int_t onevar, Int_t onecut){

  if (DD) {
    if (gUseTTMadSpin) fOutputSubDir = "KinematicHistos_DD_MadSpin/";
    else               fOutputSubDir = "KinematicHistos_DD_Madgraph/";
  }
  else {
    if (gUseTTMadSpin) fOutputSubDir = "KinematicHistos_MC_MadSpin/";
    else               fOutputSubDir = "KinematicHistos_MC_Madgraph/";
  }

  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);
  TString figname = "";

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  gROOT->LoadMacro("~folgueras/TOP/TopCode/tdrstyle.h"); 
  setTDRStyle();
  
  TLegend *leg = new TLegend(0.73,0.58,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(Data .KinHistos[0][0][0], "Data", "PL");
  //leg->AddEntry(Rare .KinHistos[0][0][0], "t#bar{t}V, VVV", "F");
  leg->AddEntry(VV   .KinHistos[0][0][0], "VV"         , "F");
  leg->AddEntry(Fake .KinHistos[0][0][0], "Non W/Z"    , "F");
  leg->AddEntry(DY   .KinHistos[0][0][0], "DY"         , "F");
  leg->AddEntry(STop .KinHistos[0][0][0], "tW/#bar{t}W", "F");
  leg->AddEntry(TTbar.KinHistos[0][0][0], "t#bar{t}"   , "F");
  
  TCanvas *c1 = new TCanvas("c1","Plot");
  c1->Divide(1,2);
   
  //Plot Pad
  TPad *plot = (TPad*)c1->GetPad(1); 
  plot->SetPad(0.01, 0.23, 0.99, 0.99);
  plot->SetTopMargin(0.1);
  plot->SetRightMargin(0.04);
  
  TPad *ratio = (TPad*)c1->GetPad(2);
  ratio->SetPad(0.01, 0.02, 0.99, 0.3);
  ratio->SetGridx();
  ratio->SetGridy();
  ratio->SetTopMargin(0.05);
  ratio->SetBottomMargin(0.4);
  ratio->SetRightMargin(0.04);
  
  THStack* MC[gNCHANNELS+1][iNCUTS][gNVARS];
  for (Int_t var=0; var<gNVARS; var++){
    if (onevar > -1 && onevar != var) continue;  //print only one 
    for (Int_t ch=0; ch<gNCHANNELS+1; ch++){
      if (onechan > -1 && onechan != ch) continue; //print only one
      for (Int_t cut=0; cut<iNCUTS; cut++){
	if (onecut > -1 && onecut != cut) continue; //print only one
	
	if (DD){
	  if (ch==gNCHANNELS) {
	    DY  .KinHistos[ch][cut][var]->Scale((DD_DY.Yields[Muon][cut]+DD_DY.Yields[Elec][cut])/(DY   .Yields[Muon][cut]+DY   .Yields[Elec][cut]));
	    Fake.KinHistos[ch][cut][var]->Scale((DD_NonW.Yields[Muon][cut]+DD_NonW.Yields[Elec][cut])/(Fake .Yields[Muon][cut]+Fake .Yields[Elec][cut]));
	  }
	  else {
	    DY   .KinHistos[ch][cut][var]->Scale(DD_DY.Yields[ch][cut]/DY.Yields[ch][cut]);
	    Fake .KinHistos[ch][cut][var]->Scale(DD_NonW.Yields[ch][cut]/Fake.Yields[ch][cut]);
	  }
	}
	
	MC[ch][cut][var] = new THStack("MC_"+gChanLabel[ch]+sCut[cut]+KinVarName[var],"");
	
	Total.KinHistos[ch][cut][var] = (TH1F*)TTbar.KinHistos[ch][cut][var]->Clone(); 
	Total.KinHistos[ch][cut][var]->Add(    STop .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    DY   .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    Fake .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    VV   .KinHistos[ch][cut][var]);
	Total.KinHistos[ch][cut][var]->Add(    Rare .KinHistos[ch][cut][var]);
		
	MC[ch][cut][var]->Add(TTbar.KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(STop .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(DY   .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(Fake .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(VV   .KinHistos[ch][cut][var]);  
	MC[ch][cut][var]->Add(Rare .KinHistos[ch][cut][var]);  
	
	plot->cd();
	
	MC[ch][cut][var]->Draw("hist");
	
	MC[ch][cut][var]->GetYaxis()->SetTitle("Events");
	MC[ch][cut][var]->GetYaxis()->SetTitleOffset(1.2);
	MC[ch][cut][var]->GetYaxis()->SetTitleSize(0.07);
	MC[ch][cut][var]->GetYaxis()->SetLabelSize(0.055);
	MC[ch][cut][var]->GetYaxis()->SetNdivisions(607);
	
	MC[ch][cut][var]->GetXaxis()->SetLabelSize(0.0);
	MC[ch][cut][var]->GetXaxis()->SetTitle("");
	
	//	Data.KinHistos[ch][cut][var]->Sumw2();
	Data.KinHistos[ch][cut][var]->SetMarkerStyle(20);
	Data.KinHistos[ch][cut][var]->Draw("SAME");
	
	leg->Draw("SAME");
	
	float maxh = Data.KinHistos[ch][cut][var]->GetMaximum();
	if (maxh < MC[ch][cut][var]->GetMaximum()) maxh = MC[ch][cut][var]->GetMaximum();
	MC[ch][cut][var]->SetMaximum(1.7*maxh);
	//	MC[ch][cut][var]->SetMinimum(50);
	
	// set logY
	//if (cut == iDilepton) plot->SetLogy();
	//if (cut == iZVeto   ) plot->SetLogy();
	
	ratio->cd();
	
	TH1F *H_Ratio;
	H_Ratio = (TH1F*)Data.KinHistos[ch][cut][var]->Clone();
	H_Ratio->Divide(Total.KinHistos[ch][cut][var]);
	
	H_Ratio->SetFillStyle(1001);
	H_Ratio->SetLineWidth(1);
	H_Ratio->SetFillColor(  kGray+1);
	H_Ratio->SetLineColor(  kGray+1);
	H_Ratio->SetMarkerColor(kGray+1);
	H_Ratio->SetMarkerStyle(20);
	//	  H_Ratio->SetMarkerColor(1);
	//    H_Ratio->SetLineColor(0);
	H_Ratio->SetMaximum(2);
	H_Ratio->SetMinimum(0);
	H_Ratio->SetTitle("");
	
	H_Ratio->GetYaxis()->SetTitle("Obs/Exp");
	H_Ratio->GetYaxis()->CenterTitle();
	H_Ratio->GetYaxis()->SetTitleOffset(0.45);
	H_Ratio->GetYaxis()->SetTitleSize(0.16);
	H_Ratio->GetYaxis()->SetLabelSize(0.15);
	H_Ratio->GetYaxis()->SetNdivisions(402);
	H_Ratio->GetXaxis()->SetTitleOffset(1.1);
	H_Ratio->GetXaxis()->SetLabelSize(0.15);
	H_Ratio->GetXaxis()->SetTitleSize(0.16);
	
	H_Ratio->SetMinimum(0.6);
	H_Ratio->SetMaximum(1.4);
	
	H_Ratio->DrawCopy("E2");
	TGraphErrors *thegraphRatio = new TGraphErrors(H_Ratio);
	thegraphRatio->SetFillStyle(3001);
	thegraphRatio->SetFillColor(kGray+1);
	//   
	thegraphRatio->Draw("E2 SAME");
	
	plot->cd();
	DrawTopLine(ch);
	
	TString channel = "";
	if (ch==Muon)      channel = "_MM";
	else if (ch==Elec) channel = "_EE";
	else if (ch==ElMu) channel = "_DF";
	else               channel = "_SF";
    
	figname = fOutputDir+fOutputSubDir + KinVarName[var] + channel + "_" + sCut[cut];
	c1->SaveAs(figname + ".pdf" );
	c1->SaveAs(figname + ".png" );
	c1->SaveAs(figname + ".root");
	
	plot->SetLogy(0);
      }
    }
  }
}
void TopPlotter::DrawNbjetsNjets(bool DD){
  if (DD) {
    if (gUseTTMadSpin) fOutputSubDir = "NbjetsNjets_DD_MadSpin/";
    else               fOutputSubDir = "NbjetsNjets_DD_Madgraph/";
  }
  else {
    if (gUseTTMadSpin) fOutputSubDir = "NbjetsNjets_MC_MadSpin/";
    else               fOutputSubDir = "NbjetsNjets_MC_Madgraph/";
  }
  
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);
  TString figname = "";

  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  
  gROOT->LoadMacro("~folgueras/TOP/TopCode/tdrstyle.h"); 
  setTDRStyle();
  
  TLegend *leg = new TLegend(0.73,0.58,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(Data .KinHistos[0][0][0], "Data",                 "PL");
  leg->AddEntry(Rare .KinHistos[0][0][0], "t#bar{t}V, VVV",        "F");
  leg->AddEntry(VV   .KinHistos[0][0][0], "VV",                    "F");
  leg->AddEntry(Fake .KinHistos[0][0][0], "Non W/Z",               "F");
  leg->AddEntry(STop .KinHistos[0][0][0], "Single Top",            "F");
  leg->AddEntry(DY   .KinHistos[0][0][0], "Z/\\gamma* l^{+}l^{-}", "F");
  leg->AddEntry(TTbar.KinHistos[0][0][0], "t#bar{t}",              "F");
  
  TCanvas *c1 = new TCanvas("c1","Plot");
  c1->Divide(1,2);
   
  //Plot Pad
  TPad *plot = (TPad*)c1->GetPad(1); 
  plot->SetPad(0.01, 0.23, 0.99, 0.99);
  plot->SetTopMargin(0.1);
  plot->SetRightMargin(0.04);
  
  TPad *ratio = (TPad*)c1->GetPad(2);
  ratio->SetPad(0.01, 0.02, 0.99, 0.3);
  ratio->SetGridx();
  ratio->SetGridy();
  ratio->SetTopMargin(0.05);
  ratio->SetBottomMargin(0.4);
  ratio->SetRightMargin(0.04);

  THStack* MC[gNCHANNELS][iNCUTS][gNALLSYST];
  for (Int_t ch=0; ch<gNCHANNELS; ch++){
    if (ch!=ElMu) continue;
    for (Int_t cut=0; cut<iNCUTS; cut++){
      if (cut != iDilepton) continue;
      for (size_t sys=0; sys<gNALLSYST; sys++){
	Int_t syst(sys);
	if (sys >= gNSYST) syst = 0;
	else               syst = sys;
	
	Total.NBtagsNJets[ch][cut][sys] = (TH1F*)TTbar.NBtagsNJets[ch][cut][sys ]->Clone();
	Total.NBtagsNJets[ch][cut][sys]->Add(    STop .NBtagsNJets[ch][cut][syst]);
	if (DD) {
	  DY   .NBtagsNJets[ch][cut][syst]->Scale(DD_DY.Yields[ch][cut]/DY.Yields[ch][cut]);  // syst???
	  Total.NBtagsNJets[ch][cut][sys]->Add(  DY   .NBtagsNJets[ch][cut][syst]);
	  Total.NBtagsNJets[ch][cut][sys]->Add(DD_NonW.NBtagsNJets[ch][cut][syst]);
 	}
	else    {
	  Total.NBtagsNJets[ch][cut][sys]->Add(  DY   .NBtagsNJets[ch][cut][syst]);
	  Total.NBtagsNJets[ch][cut][sys]->Add(  Fake .NBtagsNJets[ch][cut][syst]);
	}
	Total.NBtagsNJets[ch][cut][sys]->Add(    VV   .NBtagsNJets[ch][cut][syst]);
	Total.NBtagsNJets[ch][cut][sys]->Add(    Rare .NBtagsNJets[ch][cut][syst]);
	
	MC[ch][cut][sys] = new THStack("MC_"+gChanLabel[ch]+sCut[cut]+SystName[sys],"");
	MC[ch][cut][sys]->Add(TTbar.NBtagsNJets[ch][cut][sys ]);  
	MC[ch][cut][sys]->Add(STop .NBtagsNJets[ch][cut][syst]);  
	MC[ch][cut][sys]->Add(DY   .NBtagsNJets[ch][cut][syst]);
	if (DD) { 
	  MC[ch][cut][sys]->Add(DD_NonW.NBtagsNJets[ch][cut][syst]);
	}
	else    { 
	  MC[ch][cut][sys]->Add(Fake.NBtagsNJets[ch][cut][syst]);
	}
	MC[ch][cut][sys]->Add(VV   .NBtagsNJets[ch][cut][syst]);  
	MC[ch][cut][sys]->Add(Rare .NBtagsNJets[ch][cut][syst]);  

	plot->cd();

	MC[ch][cut][sys]->Draw("HIST");
	
	MC[ch][cut][sys]->GetYaxis()->SetTitle("Events");
	MC[ch][cut][sys]->GetYaxis()->SetTitleOffset(1.4);
	MC[ch][cut][sys]->GetYaxis()->SetTitleSize(0.05);
	MC[ch][cut][sys]->GetYaxis()->SetLabelSize(0.055);
	MC[ch][cut][sys]->GetYaxis()->SetNdivisions(607);
	
	MC[ch][cut][sys]->GetXaxis()->SetLabelSize(0.0);
	MC[ch][cut][sys]->GetXaxis()->SetTitle("");
	
	Data.NBtagsNJets[ch][cut][syst]->Sumw2();
	Data.NBtagsNJets[ch][cut][syst]->SetMarkerStyle(20);
	Data.NBtagsNJets[ch][cut][syst]->Draw("SAME");
	
	leg->Draw("SAME");
	
	float maxh = Data.NBtagsNJets[ch][cut][syst]->GetMaximum();
	if (maxh < MC[ch][cut][sys]->GetMaximum()) maxh = MC[ch][cut][sys]->GetMaximum();
	MC[ch][cut][sys]->SetMaximum(1.7*maxh);

	ratio->cd();
	
	TH1F *H_Ratio;
	H_Ratio = (TH1F*)Data.NBtagsNJets[ch][cut][syst]->Clone();
	H_Ratio->Divide(Total.NBtagsNJets[ch][cut][sys]);
	
	H_Ratio->SetFillStyle(1001);
	H_Ratio->SetLineWidth(1);
	H_Ratio->SetFillColor(  kGray+1);
	H_Ratio->SetLineColor(  kGray+1);
	H_Ratio->SetMarkerColor(kGray+1);
	H_Ratio->SetMarkerStyle(20);
	//	  H_Ratio->SetMarkerColor(1);
	//    H_Ratio->SetLineColor(0);
	H_Ratio->SetMaximum(2);
	H_Ratio->SetMinimum(0);
	H_Ratio->SetTitle("");
	
	H_Ratio->GetYaxis()->SetTitle("Obs/Exp");
	H_Ratio->GetYaxis()->CenterTitle();
	H_Ratio->GetYaxis()->SetTitleOffset(0.45);
	H_Ratio->GetYaxis()->SetTitleSize(0.16);
	H_Ratio->GetYaxis()->SetLabelSize(0.15);
	H_Ratio->GetYaxis()->SetNdivisions(402);
	H_Ratio->GetXaxis()->SetTitleOffset(1.1);
	H_Ratio->GetXaxis()->SetLabelSize(0.15);
	H_Ratio->GetXaxis()->SetTitleSize(0.16);
	
	H_Ratio->SetMinimum(0.6);
	H_Ratio->SetMaximum(1.4);
	
	H_Ratio->DrawCopy("E2");
	TGraphErrors *thegraphRatio = new TGraphErrors(H_Ratio);
	thegraphRatio->SetFillStyle(3001);
	thegraphRatio->SetFillColor(kGray+1);
	//   
	thegraphRatio->Draw("E2 SAME");
	
	plot->cd();
	DrawTopLine(ch);
   
	TString channel = "";
	if      (ch==Muon) channel = "_MM";
	else if (ch==Elec) channel = "_EE";
	else if (ch==ElMu) channel = "_DF";
	else               channel = "_SF";
	
	figname = fOutputDir+fOutputSubDir + "NbtagsNjets" + channel + "_" + sCut[cut] + "_" + SystName[sys];
	c1->SaveAs(figname + ".pdf" );
	c1->SaveAs(figname + ".png" );
	c1->SaveAs(figname + ".root");
      }
    }
  }
}
void TopPlotter::SaveHistosForLH(bool DD){
  TString filename = fOutputDir+"HistosForHL_MC.root";
  if (DD) filename = fOutputDir+"HistosForHL_DD.root";
  
  TFile *hfile = new TFile(filename,"RECREATE");
  
  TH1F *fHNBNJ_Data ;
  TH1F *fHNBNJ_TTbar[gNALLSYST];
  TH1F *fHNBNJ_Stop [gNALLSYST];
  TH1F *fHNBNJ_SUSYstop [gNALLSYST];
  TH1F *fHNBNJ_DY   [gNALLSYST];
  TH1F *fHNBNJ_VV   [gNALLSYST];
  TH1F *fHNBNJ_Rare [gNALLSYST];
  TH1F *fHNBNJ_Fake [gNALLSYST];
  
  TH1F *fHMll_Data ;
  TH1F *fHMll_TTbar[gNALLSYST];
  TH1F *fHMll_SUSYstop [gNALLSYST];
  TH1F *fHMll_Stop [gNALLSYST];
  TH1F *fHMll_DY   [gNALLSYST];
  TH1F *fHMll_VV   [gNALLSYST];
  TH1F *fHMll_Rare [gNALLSYST];
  TH1F *fHMll_Fake [gNALLSYST];
  
  TH1F *fHdll_Data ;
  TH1F *fHdll_TTbar[gNALLSYST];
  TH1F *fHdll_SUSYstop [gNALLSYST];
  TH1F *fHdll_Stop [gNALLSYST];
  TH1F *fHdll_DY   [gNALLSYST];
  TH1F *fHdll_VV   [gNALLSYST];
  TH1F *fHdll_Rare [gNALLSYST];
  TH1F *fHdll_Fake [gNALLSYST];
  
  TH1F *fHdjj_Data ;
  TH1F *fHdjj_TTbar[gNALLSYST];
  TH1F *fHdjj_SUSYstop [gNALLSYST];
  TH1F *fHdjj_Stop [gNALLSYST];
  TH1F *fHdjj_DY   [gNALLSYST];
  TH1F *fHdjj_VV   [gNALLSYST];
  TH1F *fHdjj_Rare [gNALLSYST];
  TH1F *fHdjj_Fake [gNALLSYST];
  
  TH1F *fHdjl_Data ;
  TH1F *fHdjl_TTbar[gNALLSYST];
  TH1F *fHdjl_SUSYstop [gNALLSYST];
  TH1F *fHdjl_Stop [gNALLSYST];
  TH1F *fHdjl_DY   [gNALLSYST];
  TH1F *fHdjl_VV   [gNALLSYST];
  TH1F *fHdjl_Rare [gNALLSYST];
  TH1F *fHdjl_Fake [gNALLSYST];
  
  fHNBNJ_Data = (TH1F*) Data.NBtagsNJets    [ElMu][i1btag][0]->Clone("NJetsNBjets__DATA");
  fHNBNJ_Data->Write();
  fHMll_Data  = (TH1F*) Data.InvMass        [ElMu][i1btag][0]->Clone("InvMass__DATA");
  fHMll_Data ->Write();
  fHdll_Data  = (TH1F*) Data.AbsDelPhiLeps  [ElMu][i1btag][0]->Clone("AbsDelPhiLeps__DATA");
  fHdll_Data ->Write();
  fHdjj_Data  = (TH1F*) Data.delPhi2LeadJets[ElMu][i1btag][0]->Clone("delPhi2LeadJets__DATA");
  fHdjj_Data ->Write();
  fHdjl_Data  = (TH1F*) Data.minDelRJetsLeps[ElMu][i1btag][0]->Clone("minDelRJetsLeps__DATA");
  fHdjl_Data ->Write();
  
  for (size_t sys=0; sys<gNALLSYST; sys++){
    Int_t syst = sys;
    if (sys < gNSYST) syst = sys;
    else              syst = 0; 
    
    fHNBNJ_TTbar   [sys ] = (TH1F*)TTbar   .NBtagsNJets[ElMu][i1btag][sys ]->Clone("NJetsNBjets__"+TTbar   .name + sysname[sys]);
    fHNBNJ_SUSYstop[syst] = (TH1F*)SUSYstop.NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+SUSYstop.name + sysname[sys]);
    fHNBNJ_Stop    [syst] = (TH1F*)STop    .NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+STop    .name + sysname[sys]);
    fHNBNJ_DY      [syst] = (TH1F*)DY      .NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+DY      .name + sysname[sys]);
    fHNBNJ_VV      [syst] = (TH1F*)VV      .NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+VV      .name + sysname[sys]);

//    fHNBNJ_WW   [syst] = (TH1F*)S[WWTo2L2Nu_Madgraph].NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__ww"+sysname[sys]);
//    fHNBNJ_WZ   [syst] = (TH1F*)S[WZ]                .NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__wz"+sysname[sys]);
//    fHNBNJ_ZZ   [syst] = (TH1F*)S[ZZ]                .NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__zz"+sysname[sys]);
    
    if (DD) {
      // here crahses SaveHistosForLH(true)
      fHNBNJ_Fake[syst] = (TH1F*)DD_NonW.NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+Fake.name+sysname[sys]);
    } else {
      fHNBNJ_Fake[syst] = (TH1F*)   Fake.NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+Fake.name+sysname[sys]);
    }
    fHNBNJ_Rare  [syst] = (TH1F*)   Rare.NBtagsNJets[ElMu][i1btag][syst]->Clone("NJetsNBjets__"+Rare.name+sysname[sys]);
    
    //fHNBNJ_TTbar   [sys ]->Sumw2();
    //fHNBNJ_SUSYstop[syst]->Sumw2();
    //fHNBNJ_Stop    [syst]->Sumw2();
    //fHNBNJ_DY      [syst]->Sumw2();
    //fHNBNJ_VV      [syst]->Sumw2();
    //fHNBNJ_Rare    [syst]->Sumw2();
    //fHNBNJ_Fake    [syst]->Sumw2();
  
    fHNBNJ_TTbar[sys ]->Write();
    fHNBNJ_SUSYstop[syst]->Write();
    fHNBNJ_Stop [syst]->Write();
    fHNBNJ_DY   [syst]->Write();
    fHNBNJ_VV   [syst]->Write();
//    fHNBNJ_WW   [syst]->Write();
//    fHNBNJ_WZ   [syst]->Write();
//    fHNBNJ_ZZ   [syst]->Write();
    fHNBNJ_Rare [syst]->Write();
    fHNBNJ_Fake [syst]->Write();
  
  
    fHMll_TTbar   [sys ] = (TH1F*)TTbar   .InvMass[ElMu][i1btag][sys ]->Clone("InvMass__"+TTbar   .name + sysname[sys]);
    fHMll_SUSYstop[syst] = (TH1F*)SUSYstop.InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+SUSYstop.name + sysname[sys]);
    fHMll_Stop    [syst] = (TH1F*)STop    .InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+STop    .name + sysname[sys]);
    fHMll_DY      [syst] = (TH1F*)DY      .InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+DY      .name + sysname[sys]);
    fHMll_VV      [syst] = (TH1F*)VV      .InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+VV      .name + sysname[sys]);

//    fHMll_WW   [syst] = (TH1F*)S[WWTo2L2Nu_Madgraph].InvMass[ElMu][i1btag][syst]->Clone("InvMass__ww"+sysname[sys]);
//    fHMll_WZ   [syst] = (TH1F*)S[WZ]                .InvMass[ElMu][i1btag][syst]->Clone("InvMass__wz"+sysname[sys]);
//    fHMll_ZZ   [syst] = (TH1F*)S[ZZ]                .InvMass[ElMu][i1btag][syst]->Clone("InvMass__zz"+sysname[sys]);
    
    if (DD) {
      fHMll_Fake[syst] = (TH1F*)DD_NonW.InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+Fake.name+sysname[sys]);
    } else {
      fHMll_Fake[syst] = (TH1F*)   Fake.InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+Fake.name+sysname[sys]);
    }
    fHMll_Rare [syst]  = (TH1F*)   Rare.InvMass[ElMu][i1btag][syst]->Clone("InvMass__"+Rare.name+sysname[sys]);
    
    fHMll_TTbar[sys ]->Write();
    fHMll_SUSYstop[syst]->Write();
    fHMll_Stop [syst]->Write();
    fHMll_DY   [syst]->Write();
    fHMll_VV   [syst]->Write();
//    fHMll_WW   [syst]->Write();
//    fHMll_WZ   [syst]->Write();
//    fHMll_ZZ   [syst]->Write();
    fHMll_Rare [syst]->Write();
    fHMll_Fake [syst]->Write();
  
  
    fHdll_TTbar   [sys ] = (TH1F*)TTbar   .AbsDelPhiLeps[ElMu][i1btag][sys ]->Clone("AbsDelPhiLeps__"+TTbar   .name + sysname[sys]);
    fHdll_SUSYstop[syst] = (TH1F*)SUSYstop.AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+SUSYstop.name + sysname[sys]);
    fHdll_Stop    [syst] = (TH1F*)STop    .AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+STop    .name + sysname[sys]);
    fHdll_DY      [syst] = (TH1F*)DY      .AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+DY	 .name + sysname[sys]);
    fHdll_VV      [syst] = (TH1F*)VV      .AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+VV	 .name + sysname[sys]);

//    fHMll_WW   [syst] = (TH1F*)S[WWTo2L2Nu_Madgraph].AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__ww"+sysname[sys]);
//    fHMll_WZ   [syst] = (TH1F*)S[WZ]                .AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__wz"+sysname[sys]);
//    fHMll_ZZ   [syst] = (TH1F*)S[ZZ]                .AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__zz"+sysname[sys]);
    
    if (DD) {
      fHdll_Fake[syst] = (TH1F*)DD_NonW.AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+Fake.name+sysname[sys]);
    } else {
      fHdll_Fake[syst] = (TH1F*)   Fake.AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+Fake.name+sysname[sys]);
    }
    fHdll_Rare  [syst] = (TH1F*)   Rare.AbsDelPhiLeps[ElMu][i1btag][syst]->Clone("AbsDelPhiLeps__"+Rare.name+sysname[sys]);
    
    fHdll_TTbar[sys ]->Write();
    fHdll_SUSYstop [syst]->Write();
    fHdll_Stop [syst]->Write();
    fHdll_DY   [syst]->Write();
    fHdll_VV   [syst]->Write();
//    fHdll_WW   [syst]->Write();
//    fHdll_WZ   [syst]->Write();
//    fHdll_ZZ   [syst]->Write();
    fHdll_Rare [syst]->Write();
    fHdll_Fake [syst]->Write();
  
  
    fHdjj_TTbar   [sys ] = (TH1F*)TTbar   .delPhi2LeadJets[ElMu][i1btag][sys ]->Clone("delPhi2LeadJets__"+TTbar   .name + sysname[sys]);
    fHdjj_SUSYstop[syst] = (TH1F*)SUSYstop.delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+SUSYstop.name + sysname[sys]);
    fHdjj_Stop    [syst] = (TH1F*)STop    .delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+STop    .name + sysname[sys]);
    fHdjj_DY      [syst] = (TH1F*)DY      .delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+DY      .name + sysname[sys]);
    fHdjj_VV      [syst] = (TH1F*)VV      .delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+VV      .name + sysname[sys]);

    if (DD) {
      fHdjj_Fake[syst] = (TH1F*)DD_NonW.delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+Fake.name+sysname[sys]);
    } else {
      fHdjj_Fake[syst] = (TH1F*)   Fake.delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+Fake.name+sysname[sys]);
    }
    fHdjj_Rare  [syst] = (TH1F*)   Rare.delPhi2LeadJets[ElMu][i1btag][syst]->Clone("delPhi2LeadJets__"+Rare.name+sysname[sys]);
    
    fHdjj_TTbar[sys ]->Write();
    fHdjj_SUSYstop [syst]->Write();
    fHdjj_Stop [syst]->Write();
    fHdjj_DY   [syst]->Write();
    fHdjj_VV   [syst]->Write();
    fHdjj_Rare [syst]->Write();
    fHdjj_Fake [syst]->Write();
  
  
    fHdjl_TTbar   [sys ] = (TH1F*)TTbar   .minDelRJetsLeps[ElMu][i1btag][sys ]->Clone("minDelRJetsLeps__"+TTbar   .name + sysname[sys]);
    fHdjl_SUSYstop[syst] = (TH1F*)SUSYstop.minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+SUSYstop.name + sysname[sys]);
    fHdjl_Stop    [syst] = (TH1F*)STop    .minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+STop    .name + sysname[sys]);
    fHdjl_DY      [syst] = (TH1F*)DY      .minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+DY      .name + sysname[sys]);
    fHdjl_VV      [syst] = (TH1F*)VV      .minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+VV      .name + sysname[sys]);
    
    if (DD) {
      fHdjl_Fake[syst] = (TH1F*)DD_NonW.minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+Fake.name+sysname[sys]);
    } else {
      fHdjl_Fake[syst] = (TH1F*)   Fake.minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+Fake.name+sysname[sys]);
    }
    fHdjl_Rare  [syst] = (TH1F*)   Rare.minDelRJetsLeps[ElMu][i1btag][syst]->Clone("minDelRJetsLeps__"+Rare.name+sysname[sys]);
    
    fHdjl_TTbar[sys ]->Write();
    fHdjl_SUSYstop [syst]->Write();
    fHdjl_Stop [syst]->Write();
    fHdjl_DY   [syst]->Write();
    fHdjl_VV   [syst]->Write();
    fHdjl_Rare [syst]->Write();
    fHdjl_Fake [syst]->Write();
  
  
  }
  hfile->Close();
}
void TopPlotter::CalculateSystematicErrors(Categories &C, Int_t cut){
  
  // SF ID and Iso systematics.. 
  cout << C.name << endl;
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    if (C.name == "ttbar"){
      C.SystError[ch][Q2]       = 0.;//015;//TMath::Max(TMath::Abs(S[TTJets_scaleup].Yields[ch][cut]   - 
					//		S[TTJets_MadSpin].Yields[ch][cut]),
					//     TMath::Abs(S[TTJets_scaledown].Yields[ch][cut] - 
					//		S[TTJets_MadSpin].Yields[ch][cut])
					//     ) / (S[TTJets_MadSpin].Yields[ch][cut]);
      
      C.SystError[ch][toppt]    = 0.; 
      if (toppt_weight != 0) {
	C.SystError[ch][toppt]    = 0.015;//TMath::Abs(C.Yields_syst[ch][cut][TopPt]*(1./toppt_weight) - C.Yields[ch][cut]
				       //	       ) / C.Yields[ch][cut];
      }
      //C.SystError[ch][Matching] = TMath::Max(TMath::Abs(S[TTJets_matchingup].Yields[ch][cut]   - 
	//						S[TTJets_MadSpin].Yields[ch][cut]),
	//				     TMath::Abs(S[TTJets_matchingdown].Yields[ch][cut] - 
	//						S[TTJets_MadSpin].Yields[ch][cut])
	//				     ) / (S[TTJets_MadSpin].Yields[ch][cut]);
      C.SystError[ch][cr]  = 0.0; // Scale ME unc. from Sergio; //TMath::Abs(S[TTJetsFullLeptMGTuneP11].Yields[ch][cut] - S[TTJetsFullLeptMGTuneP11noCR].Yields[ch][cut]) / S[TTJetsFullLeptMGTuneP11].Yields[ch][cut]; 
 
      C.SystError[ch][pdf] = 0.0;  //GetPDFUncertainty();
      C.SystError[ch][had] = 0.0;    // TMath::Abs(S[TTbar_Powheg].Yields[ch][cut] - S[TTbar_Powheg_Herwig].Yields[ch][cut]) / S[TTbar_Powheg].Yields[ch][cut]; 
   
      // Use here numbers from Sergio and Juan R.   
      if(cut==3){                          // 2 jets
         C.SystError[ch][cr]    = 0.0;  // Scale ME
	 C.SystError[ch][toppt] = 0.0;  // NLO generator
	 C.SystError[ch][pdf]   = 0.0;  
         C.SystError[ch][had]   = 0.0;     
      }else if(cut==4){                    // 1 b tag
         //Scale ME
         C.SystError[0][cr]    = 0.0040;  //
         C.SystError[1][cr]    = 0.0042;  //
         C.SystError[2][cr]    = 0.0010;  //

         // NLO generator
         C.SystError[0][toppt]    = 0.0193;  //
         C.SystError[1][toppt]    = 0.0316;  //
         C.SystError[2][toppt]    = 0.0212;  //

         // PDF
	 C.SystError[0][pdf]   = 0.0053;  // mu-mu
	 C.SystError[1][pdf]   = 0.0055;  // e-e
	 C.SystError[2][pdf]   = 0.0060;  // e-mu

         // JES
	 C.SystError[0][jes]   = 0.0328;  // mu-mu
	 C.SystError[1][jes]   = 0.0315;  // e-e
	 C.SystError[2][jes]   = 0.0215;  // e-mu

         // hadronization
         C.SystError[0][had]   = 0.015;    
         C.SystError[1][had]   = 0.028;    
         C.SystError[2][had]   = 0.0128;    
      }
    
    
    }
    C.SystError[ch][les]      = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][LESUp] - C.Yields[ch][cut]), 
					   TMath::Abs(C.Yields_syst[ch][cut][LESDown] - C.Yields[ch][cut])
					   ) / C.Yields[ch][cut]; 
    C.SystError[ch][jes]      = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][JESUp] - C.Yields[ch][cut]), 
					   TMath::Abs(C.Yields_syst[ch][cut][JESDown] - C.Yields[ch][cut])
					   ) / C.Yields[ch][cut];
    if(cut==4){                          
          // JES
	 C.SystError[0][jes] = 0.0258;  // mu-mu
	 C.SystError[1][jes] = 0.0254;  // e-e
	 C.SystError[2][jes] = 0.0159;  // e-mu
    }
					   
    C.SystError[ch][jer]      = (TMath::Abs(C.Yields_syst[ch][cut][JER] - C.Yields[ch][cut])) / C.Yields[ch][cut];
    
    C.SystError[ch][btag]     = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][BtagUp] - C.Yields[ch][cut]), 
					   TMath::Abs(C.Yields_syst[ch][cut][BtagDown] - C.Yields[ch][cut])
					   ) / C.Yields[ch][cut];
    C.SystError[ch][mistag]   = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][MisTagUp] - C.Yields[ch][cut]), 
					   TMath::Abs(C.Yields_syst[ch][cut][MisTagDown] - C.Yields[ch][cut])
					   ) / C.Yields[ch][cut];
    
    C.SystError[ch][pu]       = TMath::Max(TMath::Abs(C.Yields_syst[ch][cut][PUUp] - C.Yields[ch][cut]), 
					   TMath::Abs(C.Yields_syst[ch][cut][PUDown] - C.Yields[ch][cut])
					   ) / C.Yields[ch][cut];
    
    if (C.name == "vv"  ) C.SystError[ch][vv]   = 0.30;
    if (C.name == "rare") C.SystError[ch][rare] = 0.30;
    if (C.name == "stop") C.SystError[ch][stop] = 0.30;
    
    for (size_t sys=0; sys<gNSYSTERRTypesALL; sys++){ 
      C.Yields_syst[ch][cut][0] += C.SystError[ch][sys] *  C.SystError[ch][sys];
    }
    C.Yields_syst[ch][cut][0] = C.Yields[ch][cut] * TMath::Sqrt(C.Yields_syst[ch][cut][0]);
  }
  return;
}

float TopPlotter::GetPDFUncertainty(){
  /*
  Float_t pdfunc_max(0.),pdfunc(0.);
  Float_t Nevents = S[TTJets_MadSpinPDF].pdfWeightsSum->GetEntries() / S[TTJets_MadSpinPDF].pdfWeightsSum->GetNbinsX();
  Float_t pdfREF  = S[TTJets_MadSpinPDF].pdfWeights->GetBinContent(1)/S[TTJets_MadSpinPDF].pdfWeightsSum->GetBinContent(1)*Nevents;
  for (Int_t i=0; i<51; i++){
    Float_t pdfnew = S[TTJets_MadSpinPDF].pdfWeights->GetBinContent(i+1)/S[TTJets_MadSpinPDF].pdfWeightsSum->GetBinContent(i+1)*Nevents;
    pdfunc = TMath::Abs(pdfnew - pdfREF)/pdfREF;
    if (pdfunc >= pdfunc_max) pdfunc_max = pdfunc;
  }
  return pdfunc_max;
  */
  return 1.;
}
void TopPlotter::CalculateSystematicErrorsWithXSec(Categories &C, Int_t cut){
  /////////////////////////////////////////////////////////////////////////
  //   This method calculates the cross-section for each systematic and    
  //   compare it with the nominal one, the difference is assigned as a    
  //   systematic uncertainty (for systematics varied up/down)             
  /////////////////////////////////////////////////////////////////////////
  fOutputSubDir = "SystematicErrors/";
  
  XSection tmp_up, tmp_down; 
  Float_t obs(0.), bkg(0.), eff(0.);
  
  for (Int_t ch = 0; ch<gNCHANNELS; ch++){
    obs = Data .Yields[ch][cut];
    bkg = Total.Yields[ch][cut];
    eff = TTbar.Yields[ch][cut];
    
    tmp_up  .xsec[ch] = ttbar_TLWG * (obs - bkg) / eff;
    tmp_down.xsec[ch] = ttbar_TLWG * (obs - bkg) / eff;
    if (C.name == "ttbar") eff = 0.;
  }
  
}
void TopPlotter::CalculateDYBkg(){
  fOutputSubDir = "DataDriven/";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  TString yieldsfilename = "";
  Double_t R    [gNCHANNELS][iNCUTS];
  Double_t R_err[gNCHANNELS][iNCUTS];
  Double_t N_in [gNCHANNELS];
  Double_t N_out[gNCHANNELS];
  Double_t k_ll [gNCHANNELS];
  Double_t SF   [gNCHANNELS];

  Double_t N_in_err [gNCHANNELS];
  Double_t N_out_err[gNCHANNELS];
  Double_t k_ll_err [gNCHANNELS];
  Double_t SF_err   [gNCHANNELS];

  yieldsfilename = fOutputDir + fOutputSubDir + "DY.txt";
  fOUTSTREAM.open(yieldsfilename, ios::trunc);

  for (size_t cut=0; cut<iNCUTS; cut++){    
    if (cut == iZVeto) continue;
    // Calculate R:
    Int_t low_in = DY.MllHistos[Muon][cut]->FindBin(76.);
    Int_t  up_in = DY.MllHistos[Muon][cut]->FindBin(106.);
    
    Double_t nout_ee(0.),nin_ee(0.),nout_err_ee(0.),nin_err_ee(0.);
    Double_t nout_mm(0.),nin_mm(0.),nout_err_mm(0.),nin_err_mm(0.);
    nin_mm  = DY.MllHistos[Muon][cut]->Integral(low_in, up_in); nin_err_mm = TMath::Sqrt(nin_mm);
    nin_ee  = DY.MllHistos[Elec][cut]->Integral(low_in, up_in); nin_err_mm = TMath::Sqrt(nin_ee); 
    
    nout_mm = DY.MllHistos[Muon][cut]->Integral(0, 200)-nin_mm; nout_err_mm = TMath::Sqrt(nout_mm);
    nout_ee = DY.MllHistos[Elec][cut]->Integral(0, 200)-nin_ee; nout_err_ee = TMath::Sqrt(nout_ee);
    
    R    [Muon][cut] = nout_mm / nin_mm;
    R    [Elec][cut] = nout_ee / nin_ee;
    R_err[Muon][cut] = (nout_err_mm/nout_mm + nin_err_mm/nin_mm) * R[Muon][cut];
    R_err[Elec][cut] = (nout_err_ee/nout_ee + nin_err_ee/nin_ee) * R[Elec][cut];
    

    R_err[Muon][cut] = TMath::Sqrt(R_err[Muon][cut]*R_err[Muon][cut] + 0.5*R[Muon][cut]*0.5*R[Muon][cut]);    
    R_err[Elec][cut] = TMath::Sqrt(R_err[Elec][cut]*R_err[Elec][cut] + 0.5*R[Elec][cut]*0.5*R[Elec][cut]);

    N_in [Muon] = Data.MllHistos[Muon][cut]->IntegralAndError(low_in, up_in, N_in_err[Muon]); 
    N_in [Elec] = Data.MllHistos[Elec][cut]->IntegralAndError(low_in, up_in, N_in_err[Elec]);
    N_in [ElMu] = Data.MllHistos[ElMu][cut]->IntegralAndError(low_in, up_in, N_in_err[ElMu]);
  
    k_ll [Muon] = TMath::Sqrt(nin_mm / nin_ee  );  
    k_ll [Elec] = TMath::Sqrt(nin_ee / nin_mm  );
    
    k_ll_err[Muon] = 0.5*(nin_err_mm/nin_mm + nin_err_ee/nin_ee);
    k_ll_err[Elec] = 0.5*(nin_err_ee/nin_ee + nin_err_mm/nin_mm);
    
    // CALCULATE DD ESTIMATES... 
    N_out[Muon] = R[Muon][cut] * (N_in[Muon] - 0.5 * N_in[ElMu] * k_ll[Muon]);
    N_out[Elec] = R[Elec][cut] * (N_in[Elec] - 0.5 * N_in[ElMu] * k_ll[Elec]);
    
    Double_t tmperr1 =  N_in_err[Elec]/N_in[Elec];
    Double_t tmperr2 =  N_in_err[ElMu]/N_in[ElMu] + k_ll_err[Elec]/k_ll[Elec];
    N_out_err[Elec] = R_err[Elec][cut]/R[Elec][cut]*(tmperr1 + 0.5*tmperr2) * N_out[Elec];
    
    tmperr1 =  N_in_err[Muon]/N_in[Muon];
    tmperr2 =  N_in_err[ElMu]/N_in[ElMu] + k_ll_err[Muon]/k_ll[Muon];
    N_out_err[Muon] = R_err[Muon][cut]/R[Muon][cut]*(tmperr1 + 0.5*tmperr2) * N_out[Muon];
    
    SF   [Muon] = N_out[Muon] / nout_mm;  
    SF   [Elec] = N_out[Elec] / nout_ee;
    SF   [ElMu] = TMath::Sqrt(SF[Elec] * SF[Muon]);

    SF_err[Muon] = (N_out_err[Muon]/N_out[Muon] + nout_err_mm/nout_mm) * SF[Muon];  
    SF_err[Elec] = (N_out_err[Elec]/N_out[Elec] + nout_err_ee/nout_ee) * SF[Elec];  
    SF_err[ElMu] = 0.5*SF[ElMu]*(SF_err[Elec]/SF[Elec] + SF_err[Muon]/SF[Muon]);
    
    fOUTSTREAM << " Estimation for " + sCut[cut] << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "                  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" n_in (MC)        | %7.1f +/- %4.1f | %7.1f +/- %4.1f |                  ||",   		       
		       nin_ee, N_in_err[Elec], 
		       nin_mm, N_in_err[Muon]) << endl;
    fOUTSTREAM << Form(" n_out (MC)       | %7.1f +/- %4.1f | %7.1f +/- %4.1f |                  ||",   		       
		       nout_ee, N_out_err[Elec], 
		       nout_mm, N_out_err[Muon]) << endl;
    fOUTSTREAM << Form(" R(Nout/Nin)(MC)  | %6.3f +/- %4.3f | %6.3f +/- %4.3f |                  ||",
		       R[Elec][cut], R_err[Elec][cut], 
		       R[Muon][cut], R_err[Muon][cut]) << endl;
    fOUTSTREAM << Form(" k_ll             | %6.3f +/- %4.3f | %6.3f +/- %4.3f |                  ||",   		       
		       k_ll[Elec], k_ll_err[Elec], 
		       k_ll[Muon], k_ll_err[Muon]) << endl;
    fOUTSTREAM << Form(" N_in (D)         | %7.1f +/- %4.1f | %7.1f +/- %4.1f |  %7.1f +/- %4.1f                 ||",   		       
		       N_in[Elec], N_in_err[Elec], 
		       N_in[Muon], N_in_err[Muon],
		       N_in[ElMu], N_in_err[ElMu]) << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" <<  endl;  
    fOUTSTREAM << Form(" N_out            | %7.1f +/- %4.1f | %7.1f +/- %4.1f |                  ||",   		       
		       N_out[Elec], N_out_err[Elec], 
		       N_out[Muon], N_out_err[Muon]) << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" <<  endl;  
    fOUTSTREAM << Form(" SF (D/MC)        | %6.3f +/- %4.3f | %6.3f +/- %4.3f | %6.3f +/- %4.3f ||",
		       SF[Elec], SF_err[Elec], SF[Muon], SF_err[Muon], SF[ElMu], SF_err[ElMu]) << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" <<  endl;  
    fOUTSTREAM << Form(" Drell-Yan (MC)   | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       DY.Yields[Elec][cut], DY.Yields_stat[Elec][cut], 
		       DY.Yields[Muon][cut], DY.Yields_stat[Muon][cut], 
		       DY.Yields[ElMu][cut], DY.Yields_stat[ElMu][cut]) << endl;
    fOUTSTREAM << Form(" Drell-Yan (D)    | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       DY.Yields[Elec][cut] * SF[Elec], DY.Yields_stat[Elec][cut] * SF[Elec], 
		       DY.Yields[Muon][cut] * SF[Muon], DY.Yields_stat[Muon][cut] * SF[Muon], 
		       DY.Yields[ElMu][cut] * SF[ElMu], DY.Yields_stat[ElMu][cut] * SF[ElMu]) << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    fOUTSTREAM << endl;
    
    DY_SF[Elec][cut] = SF[Elec];
    DY_SF[Muon][cut] = SF[Muon];
    DY_SF[ElMu][cut] = SF[ElMu];
  }
  fOUTSTREAM.close();
  gSystem->Exec("cat "+yieldsfilename);

  /// DRAW R
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1000);
  gStyle->SetOptStat("emruo");
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  Double_t x    [iNCUTS] = {1, 2, 3, 4, 5};
  Double_t x_err[iNCUTS] = {0., 0., 0., 0., 0.};
  
  R[Elec][iZVeto] = R[Elec][iDilepton];
  R[Muon][iZVeto] = R[Muon][iDilepton];
  TGraphErrors *el = new TGraphErrors(iNCUTS, x, R[Elec], x_err, R_err[Elec]);
  TGraphErrors *mu = new TGraphErrors(iNCUTS, x, R[Muon], x_err, R_err[Muon]);
  
  el->SetMarkerStyle(20);
  el->SetMarkerColor(kRed);
  el->SetLineColor(kRed);
  el->SetMaximum(0.20);
  el->SetMinimum(0.00);
  mu->SetMarkerStyle(20);
  mu->SetMarkerColor(kBlue);
  mu->SetLineColor(kBlue);

  TLegend *leg = new TLegend(0.73,0.73,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(el, "R (ee) ", "PL");
  leg->AddEntry(mu, "R (#mu#mu)", "PL");
  
  TCanvas *c1 = new TCanvas("c1","", 800, 600);
  c1->cd();
  el->Draw("AP");
  mu->Draw("P SAME");
  leg->Draw("SAME");

  //  cout << "ERRORS: " << endl;
//  Float_t mumax =  TMath::MaxElement(mu->GetN(),mu->GetY());
//  Float_t mumin =  TMath::MinElement(mu->GetN(),mu->GetY());
//  Float_t elmax =  TMath::MaxElement(el->GetN(),el->GetY());
//  Float_t elmin =  TMath::MinElement(el->GetN(),el->GetY());

//  cout<< elmax <<" - "<< elmin <<" = "<< (elmax-elmin)/el->GetMean(2) << endl;
//  cout<< mumax <<" - "<< mumin <<" = "<< (mumax-mumin)/mu->GetMean(2) << endl;
//  cout<< el->GetMean(2) << "+/-" <<  el->GetRMS(2) << endl;
//  cout<< mu->GetMean(2) << "+/-" <<  mu->GetRMS(2) << endl;

  c1->SaveAs(fOutputDir + fOutputSubDir + "Routin.png");
  c1->SaveAs(fOutputDir + fOutputSubDir + "Routin.pdf");
  c1->SaveAs(fOutputDir + fOutputSubDir + "Routin.png");
    
  for (size_t chan=0; chan<gNCHANNELS; chan++){
    for (size_t ct=0; ct<iNCUTS; ct++){
      DD_DY.Yields     [chan][ct]    = DY.Yields     [chan][ct]*DY_SF[chan][ct];
      DD_DY.Yields_stat[chan][ct]    = DY.Yields_stat[chan][ct]*DY_SF[chan][ct];
      DD_DY.Yields_syst[chan][ct][0] = 0.15*DD_DY.Yields[chan][ct];
    }
  }
}


void TopPlotter::CalculateNonWZLeptonsBkg(){
  fOutputSubDir = "DataDriven/";
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);

  TString yieldsfilename = "";
  Double_t R    [gNCHANNELS][iNCUTS];
  Double_t R_err[gNCHANNELS][iNCUTS];
  for (size_t cut=0; cut<iNCUTS; cut++){
    yieldsfilename = fOutputDir + fOutputSubDir + "NonWZ_"+sCut[cut]+".txt";
    fOUTSTREAM.open(yieldsfilename, ios::trunc);
    
    // Calculate R:
    for (size_t ch=0; ch<gNCHANNELS; ch++){
      R    [ch][cut] = Fake.Yields[ch][cut] / Fake.SSYields[ch][cut];
      R_err[ch][cut] = Fake.Yields_stat[ch][cut]*Fake.Yields_stat[ch][cut] / (Fake.SSYields[ch][cut]*Fake.SSYields[ch][cut]) + Fake.SSYields_stat[ch][cut]*Fake.SSYields_stat[ch][cut]*(R[ch][cut] / Fake.SSYields[ch][cut])*(R[ch][cut] / Fake.SSYields[ch][cut]);
      R_err[ch][cut] = TMath::Sqrt(R_err[ch][cut]);

      DD_NonW.Yields[ch][cut]	      = R[ch][cut] * (Data.SSYields[ch][cut]-Total.SSYields[ch][cut]);
      //DD_NonW.Yields_stat[ch][cut]    = TMath::Sqrt(TMath::Power(R[ch][cut]*Total.SSYields_stat[ch][cut],2) + TMath::Power(R_err[ch][cut]*(Data.SSYields[ch][cut]-Total.SSYields[Elec][cut]),2));
      DD_NonW.Yields_stat[ch][cut]    = TMath::Sqrt(TMath::Power(R[ch][cut]*Total.SSYields_stat[ch][cut],2) + TMath::Power(R_err[ch][cut]*(Data.SSYields[ch][cut]-Total.SSYields[ch][cut]),2));
      DD_NonW.Yields_syst[ch][cut][0] = 0.3*DD_NonW.Yields[ch][cut];  // ASSING a 30% uncertainty.
    }
    
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << "          Source  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" TTbar dilepton   | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       TTbar.SSYields[Elec][cut], TTbar.SSYields_stat[Elec][cut], 
		       TTbar.SSYields[Muon][cut], TTbar.SSYields_stat[Muon][cut], 
		       TTbar.SSYields[ElMu][cut], TTbar.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Drell-Yan        | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||", 
		       DY.SSYields[Elec][cut], DY.SSYields_stat[Elec][cut], 
		       DY.SSYields[Muon][cut], DY.SSYields_stat[Muon][cut], 
		       DY.SSYields[ElMu][cut], DY.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Single top quark | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       STop.SSYields[Elec][cut], STop.SSYields_stat[Elec][cut],
		       STop.SSYields[Muon][cut], STop.SSYields_stat[Muon][cut], 
		       STop.SSYields[ElMu][cut], STop.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Dibosons         | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       VV.SSYields[Elec][cut], VV.SSYields_stat[Elec][cut], 
		       VV.SSYields[Muon][cut], VV.SSYields_stat[Muon][cut], 
		       VV.SSYields[ElMu][cut], VV.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Rare             | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       Rare.SSYields[Elec][cut], Rare.SSYields_stat[Elec][cut], 
		       Rare.SSYields[Muon][cut], Rare.SSYields_stat[Muon][cut], 
		       Rare.SSYields[ElMu][cut], Rare.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Total background | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       Total.SSYields[Elec][cut], Total.SSYields_stat[Elec][cut], 
		       Total.SSYields[Muon][cut], Total.SSYields_stat[Muon][cut], 
		       Total.SSYields[ElMu][cut], Total.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z lep (SS) | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       Fake.SSYields[Elec][cut], Fake.SSYields_stat[Elec][cut], 
		       Fake.SSYields[Muon][cut], Fake.SSYields_stat[Muon][cut], 
		       Fake.SSYields[ElMu][cut], Fake.SSYields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" Non W/Z lep (OS) | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       Fake.Yields[Elec][cut], Fake.Yields_stat[Elec][cut], 
		       Fake.Yields[Muon][cut], Fake.Yields_stat[Muon][cut], 
		       Fake.Yields[ElMu][cut], Fake.Yields_stat[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << Form(" R (OS/SS)        | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       R[Elec][cut], R_err[Elec][cut], 
		       R[Muon][cut], R_err[Muon][cut], 
		       R[ElMu][cut], R_err[ElMu][cut])
	       << endl;							  
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" Data             | %5.0f            | %5.0f            | %5.0f            ||", 
		       Data.SSYields[Elec][cut], 
		       Data.SSYields[Muon][cut], 
		       Data.SSYields[ElMu][cut] )<< endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << Form(" SS prediction    | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       Data.SSYields[Elec][cut]-Total.SSYields[Elec][cut], Total.SSYields_stat[Elec][cut],
		       Data.SSYields[Muon][cut]-Total.SSYields[Muon][cut], Total.SSYields_stat[Muon][cut],
		       Data.SSYields[ElMu][cut]-Total.SSYields[ElMu][cut], Total.SSYields_stat[ElMu][cut])
	       << endl;
    fOUTSTREAM << Form(" SS x R(OSSS)     | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||",  
		       DD_NonW.Yields[Elec][cut], DD_NonW.Yields_stat[Elec][cut],
		       DD_NonW.Yields[Muon][cut], DD_NonW.Yields_stat[Muon][cut],
		       DD_NonW.Yields[ElMu][cut], DD_NonW.Yields_stat[ElMu][cut])
      	       << endl;
    fOUTSTREAM << "-----------------------------------------------------------------------------" << endl;
    fOUTSTREAM << endl;
    
    fOUTSTREAM << endl;
    fOUTSTREAM.close();
    
    if (cut == i1btag) gSystem->Exec("cat "+yieldsfilename);
  }
  // DRAWING STUFF...
  Double_t x    [iNCUTS] = {1, 2, 3, 4, 5};
  Double_t x_err[iNCUTS] = {0., 0., 0., 0., 0.};
  
  TGraphErrors *elel = new TGraphErrors(iNCUTS, x, R[Elec], x_err, R_err[Elec]);
  TGraphErrors *mumu = new TGraphErrors(iNCUTS, x, R[Muon], x_err, R_err[Muon]);
  TGraphErrors *elmu = new TGraphErrors(iNCUTS, x, R[ElMu], x_err, R_err[ElMu]);
  
  elel->SetMarkerStyle(20);
  elel->SetMarkerColor(kRed);
  elel->SetLineColor(kRed);
  elel->SetMaximum(2.0);
  elel->SetMinimum(0.5);
  mumu->SetMarkerStyle(20);
  mumu->SetMarkerColor(kBlue);
  mumu->SetLineColor(kBlue);
  elmu->SetMarkerStyle(20);
  elmu->SetMarkerColor(kBlack);
  elmu->SetLineColor(kBlack);

  TLegend *leg = new TLegend(0.73,0.73,0.90,0.89);
  leg->SetFillColor(0);
  leg->SetLineColor(1);
  //leg->SetLineWidth(4);
  leg->SetTextFont(62); // Events in the leg!
  leg->SetTextSize(0.04);
  
  leg->AddEntry(elel, "R (ee) ", "PL");
  leg->AddEntry(mumu, "R (#mu#mu)", "PL");
  leg->AddEntry(elmu, "R (e#mu)", "PL");
  
  TCanvas *c1 = new TCanvas("c1","", 800, 600);
  c1->cd();
  elel->Draw("AP");
  mumu->Draw("P SAME");
  elmu->Draw("P SAME");
  leg->Draw("SAME");

  //  cout << "ERRORS: " << endl;
//  Float_t mumumax =  TMath::MaxElement(mumu->GetN(),mumu->GetY());
//  Float_t mumumin =  TMath::MinElement(mumu->GetN(),mumu->GetY());
//  Float_t elmumax =  TMath::MaxElement(elmu->GetN(),elmu->GetY());
//  Float_t elmumin =  TMath::MinElement(elmu->GetN(),elmu->GetY());
//  Float_t elelmax =  TMath::MaxElement(elel->GetN(),elel->GetY());
//  Float_t elelmin =  TMath::MinElement(elel->GetN(),elel->GetY());

//  cout<< elelmax <<" - "<< elelmin <<" = "<< (elelmax-elelmin)/elel->GetMean(2) << endl;
//  cout<< mumumax <<" - "<< mumumin <<" = "<< (mumumax-mumumin)/mumu->GetMean(2) << endl;
//  cout<< elmumax <<" - "<< elmumin <<" = "<< (elmumax-elmumin)/elmu->GetMean(2) << endl;

  c1->SaveAs(fOutputDir + fOutputSubDir + "R_NonW.png");
  c1->SaveAs(fOutputDir + fOutputSubDir + "R_NonW.pdf");
  c1->SaveAs(fOutputDir + fOutputSubDir + "R_NonW.png");
  
  
  /////////////////////////////////////////////////
  ////  TEMPLATES 
  /////////////////////////////////////////////////
  /*  for (size_t cut=0; cut<iNCUTS; cut++){
      for (size_t chan=0; chan<gNCHANNELS; chan++){
      for (size_t sys=0; sys<gNSYST; sys++){
      
      DD_NonW.NBtagsNJets[chan][cut][sys] = (TH1F*) Fake.SSNBtagsNJets[chan][cut][sys]->Clone();
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    Fake.SSNBtagsNJets[chan][cut][sys], -1);
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    Data.SSNBtagsNJets[chan][cut][sys],  1);
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    STop.SSNBtagsNJets[chan][cut][sys], -1);
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    VV  .SSNBtagsNJets[chan][cut][sys], -1);
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    DY  .SSNBtagsNJets[chan][cut][sys], -1);
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    Rare.SSNBtagsNJets[chan][cut][sys], -1);
      DD_NonW.NBtagsNJets[chan][cut][sys] ->Add(    STop.SSNBtagsNJets[chan][cut][sys], -1);
      
      
      DD_NonW.InvMass[chan][cut][sys] = (TH1F*) Fake.SSInvMass[chan][cut][sys]->Clone();
      DD_NonW.InvMass[chan][cut][sys] ->Add(    Fake.SSInvMass[chan][cut][sys], -1);
      DD_NonW.InvMass[chan][cut][sys] ->Add(    Data.SSInvMass[chan][cut][sys],  1);
      DD_NonW.InvMass[chan][cut][sys] ->Add(    STop.SSInvMass[chan][cut][sys], -1);
      DD_NonW.InvMass[chan][cut][sys] ->Add(    VV  .SSInvMass[chan][cut][sys], -1);
      DD_NonW.InvMass[chan][cut][sys] ->Add(    DY  .SSInvMass[chan][cut][sys], -1);
      DD_NonW.InvMass[chan][cut][sys] ->Add(    Rare.SSInvMass[chan][cut][sys], -1);
      DD_NonW.InvMass[chan][cut][sys] ->Add(    STop.SSInvMass[chan][cut][sys], -1);
      }
      }
      }
  */
}
//float TopPlotter::GetXSection(Categories &C, Int_t ch, Int_t cut, Int_t sys){
//  float xsec = 0.;
//  
//  xsec = ttbar_TLWG * (Data.Yields[ch][cut] - () 
//  
//}
void TopPlotter::CalculateCrossSection(Bool_t DD){
  // CALCULATE XSECTION without DD calculation (default);
  fOutputSubDir = "XSection/";
  TString filename = "";
  iCut cut = i1btag; // //iZVeto; //i2jets; //i1btag;
  //iCut cut = i2jets; 
  //iCut cut = iDilepton; 
  const char* scut = sCut[cut].Data();
  if (DD) {
    filename = fOutputDir+fOutputSubDir+Form("Cross_Section_%3.1f_%s_DD.txt",fLumiNorm/1000.,scut);
  }
  else {
    filename = fOutputDir+fOutputSubDir+Form("Cross_Section_%3.1f_%s_MC.txt",fLumiNorm/1000.,scut);
  }
  gSystem->mkdir(fOutputDir+fOutputSubDir, kTRUE);
  
  // FIRST CALCULATE REFERENCE...
  for (Int_t ch = 0; ch<gNCHANNELS; ch++){
    Total.Yields[ch][cut]  = STop.Yields[ch][cut];
    Total.Yields[ch][cut] += VV  .Yields[ch][cut];    
    Total.Yields[ch][cut] += Rare.Yields[ch][cut];
    if (DD){
      Total.Yields[ch][cut] += DD_DY  .Yields[ch][cut];
      Total.Yields[ch][cut] += DD_NonW.Yields[ch][cut];
    }
    else {
      Total.Yields[ch][cut] += DY  .Yields[ch][cut];
      Total.Yields[ch][cut] += Fake.Yields[ch][cut];
    }
    
    ttbar.xsec     [ch] = ttbar_TLWG * (Data.Yields[ch][cut] - Total.Yields[ch][cut]) / TTbar.Yields[ch][cut];// *604354./34542493.0; 
    
    // statistical error
    ttbar.xsec_stat[ch] = ttbar.xsec[ch] * TMath::Sqrt(Data.Yields[ch][cut]) / (Data.Yields[ch][cut] - Total.Yields[ch][cut]);
  }
  
  // Get Systematic Errors
  ResetSystematicErrors();
  
  
  cout << "Calculating Systematic Errors..." << endl;
  CalculateSystematicErrors(SUSYstop, cut);   
  CalculateSystematicErrors(TTbar, cut);
  CalculateSystematicErrors(STop , cut);   
  CalculateSystematicErrors(VV   , cut);
  CalculateSystematicErrors(Rare , cut);
  if (!DD) {
    CalculateSystematicErrors(Fake , cut);
    CalculateSystematicErrors(DY   , cut);
  }
  
  for (size_t ch=0; ch<gNCHANNELS; ch++){
    Total.Yields_syst[ch][cut][0]  = STop .Yields_syst[ch][cut][0] * STop .Yields_syst[ch][cut][0];
    Total.Yields_syst[ch][cut][0] += VV   .Yields_syst[ch][cut][0] * VV   .Yields_syst[ch][cut][0];
    //Total.Yields_syst[ch][cut][0] += Rare .Yields_syst[ch][cut][0] * Rare .Yields_syst[ch][cut][0];
    if (!DD){
      Total.Yields_syst[ch][cut][0] += DY  .Yields_syst[ch][cut][0] * DY  .Yields_syst[ch][cut][0];
      Total.Yields_syst[ch][cut][0] += Fake.Yields_syst[ch][cut][0] * Fake.Yields_syst[ch][cut][0];
    }
    else { 
      Total.Yields_syst[ch][cut][0]+=DD_NonW.Yields_syst[ch][cut][0]*DD_NonW.Yields_syst[ch][cut][0];
      Total.Yields_syst[ch][cut][0]+=DD_DY.Yields_syst[ch][cut][0]*DD_DY.Yields_syst[ch][cut][0];
    }
    Total.Yields_syst[ch][cut][0]  = TMath::Sqrt(Total.Yields_syst[ch][cut][0]);
  }
  
  // Get ttbar XSection
  for (Int_t ch = 0; ch<gNCHANNELS; ch++){
    ttbar.xsec     [ch] = ttbar_TLWG * (Data.Yields[ch][cut] - Total.Yields[ch][cut]) / TTbar.Yields[ch][cut];// *604354./34542493.0; 
    
    // statistical error
    ttbar.xsec_stat[ch] = ttbar.xsec[ch] * 
      TMath::Sqrt(Data.Yields[ch][cut]) / (Data.Yields[ch][cut] - Total.Yields[ch][cut]);
    
    // systematic error
    ttbar.xsec_syst[ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((Total.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		  (Total.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) + 
		  (TTbar.Yields_syst[ch][cut][0] / TTbar.Yields[ch][cut]) * 
		  (TTbar.Yields_syst[ch][cut][0] / TTbar.Yields[ch][cut]));
    ttbar.xsec_lumi[ch] = ttbar.xsec[ch]  * 0.027; //0.12; // Lumi uncertainty 
    
    // Acceptances...
    ttbar.acc      [ch] = TTbar.Yields     [ch][cut]    / (fLumiNorm * ttbar_TLWG);
    ttbar.acc_stat [ch] = TTbar.Yields_stat[ch][cut]    / (fLumiNorm * ttbar_TLWG);
    ttbar.acc_syst [ch] = TTbar.Yields_syst[ch][cut][0] / (fLumiNorm * ttbar_TLWG);
    
    // Loading systematics...
    ttbar.err_VV   [ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((VV   .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
     		  (VV   .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    ttbar.err_STop [ch] = ttbar.xsec[ch] * 
      TMath::Sqrt((STop .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
     		  (STop .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    if (DD) {
      ttbar.err_DY   [ch] = ttbar.xsec[ch] * 
	TMath::Sqrt((DD_DY  .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		    (DD_DY  .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
      ttbar.err_Fake [ch] = ttbar.xsec[ch] * 
	TMath::Sqrt((DD_NonW.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		    (DD_NonW.Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    }
    else {
      ttbar.err_DY   [ch] = ttbar.xsec[ch] * 
	TMath::Sqrt((DY   .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		    (DY   .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
      ttbar.err_Fake [ch] = ttbar.xsec[ch] * 
	TMath::Sqrt((Fake .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
		    (Fake .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    }
    ttbar.err_Rare [ch] = 0.; 
      //ttbar.xsec[ch] * 
      //TMath::Sqrt((Rare .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])) *
      //		  (Rare .Yields_syst[ch][cut][0] / (Data.Yields[ch][cut] - Total.Yields[ch][cut])));
    ttbar.err_IDIso[ch] = ttbar.xsec[ch] * TTbar.SystError[ch][SFIDISO];
    ttbar.err_Trig [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][SFTrig];
    ttbar.err_LES  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][les];
    ttbar.err_JES  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][jes];
        
    ttbar.err_JER  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][jer];
    ttbar.err_Btag [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][btag];
    ttbar.err_mtag [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][mistag];
    ttbar.err_PU   [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][pu];
    ttbar.err_TopPt[ch] = ttbar.xsec[ch] * TTbar.SystError[ch][toppt];
    ttbar.err_cr   [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][cr];
    ttbar.err_pdf  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][pdf];
    ttbar.err_had  [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][had];
    ttbar.err_Q2   [ch] = ttbar.xsec[ch] * TTbar.SystError[ch][Q2];
    //ttbar.err_Match[ch] = ttbar.xsec[ch] * TTbar.SystError[ch][Matching];
  }
  
  // print systematic errors
  //  PrintSystematicErrors();
  
  fOUTSTREAM.open(filename, ios::trunc);

  fOUTSTREAM << "\n \\hline\\hline" << endl;
  fOUTSTREAM << "           &                        Systematic Uncertainties (pb)       \\" << endl;
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << " Source                &     El/El     &     Mu/Mu     &    El/Mu      \\" << endl;
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << Form(" VV modelling          & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%)  \\\\", 
		     ttbar.err_VV[Elec],  100 * ttbar.err_VV[Elec] / ttbar.xsec[Elec],
		     ttbar.err_VV[Muon],  100 * ttbar.err_VV[Muon] / ttbar.xsec[Muon],
		     ttbar.err_VV[ElMu],  100 * ttbar.err_VV[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Single-Top modelling  & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_STop[Elec],  100 * ttbar.err_STop[Elec] / ttbar.xsec[Elec],
		     ttbar.err_STop[Muon],  100 * ttbar.err_STop[Muon] / ttbar.xsec[Muon],
		     ttbar.err_STop[ElMu],  100 * ttbar.err_STop[ElMu] / ttbar.xsec[ElMu])<< endl;
  //fOUTSTREAM << Form(" Rare-SM modelling     & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) ", 
//		     ttbar.err_Rare[Elec],  100 * ttbar.err_Rare[Elec] / ttbar.xsec[Elec],
//		     ttbar.err_Rare[Muon],  100 * ttbar.err_Rare[Muon] / ttbar.xsec[Muon],
//		     ttbar.err_Rare[ElMu],  100 * ttbar.err_Rare[ElMu] / ttbar.xsec[ElMu])<< endl;
		     //10.3,  1.7   )<< endl;
  fOUTSTREAM << Form(" Drell-Yan background  & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%)  \\\\", 
		     ttbar.err_DY[Elec],  100 * ttbar.err_DY[Elec] / ttbar.xsec[Elec],
		     ttbar.err_DY[Muon],  100 * ttbar.err_DY[Muon] / ttbar.xsec[Muon],
		     ttbar.err_DY[ElMu],  100 * ttbar.err_DY[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Non-Prompt background & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_Fake[Elec],  100 * ttbar.err_Fake[Elec] / ttbar.xsec[Elec],
		     ttbar.err_Fake[Muon],  100 * ttbar.err_Fake[Muon] / ttbar.xsec[Muon],
		     ttbar.err_Fake[ElMu],  100 * ttbar.err_Fake[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << Form(" Lepton Efficiencies   & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_IDIso[Elec],  100 * ttbar.err_IDIso[Elec] / ttbar.xsec[Elec],
		     ttbar.err_IDIso[Muon],  100 * ttbar.err_IDIso[Muon] / ttbar.xsec[Muon],
		     ttbar.err_IDIso[ElMu],  100 * ttbar.err_IDIso[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Trigger Efficiencies  & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_Trig[Elec],  100 * ttbar.err_Trig[Elec] / ttbar.xsec[Elec],
		     ttbar.err_Trig[Muon],  100 * ttbar.err_Trig[Muon] / ttbar.xsec[Muon],
		      ttbar.err_Trig[ElMu],  100 * ttbar.err_Trig[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Lepton Energy Scale   & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_LES[Elec],  100 * ttbar.err_LES[Elec] / ttbar.xsec[Elec],
		     ttbar.err_LES[Muon],  100 * ttbar.err_LES[Muon] / ttbar.xsec[Muon],
		     ttbar.err_LES[ElMu],  100 * ttbar.err_LES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Jet Energy Scale      & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_JES[Elec],  100 * ttbar.err_JES[Elec] / ttbar.xsec[Elec],
		     ttbar.err_JES[Muon],  100 * ttbar.err_JES[Muon] / ttbar.xsec[Muon],
		     ttbar.err_JES[ElMu],  100 * ttbar.err_JES[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Jet Energy Resolution & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_JER[Elec],  100 * ttbar.err_JER[Elec] / ttbar.xsec[Elec],
		     ttbar.err_JER[Muon],  100 * ttbar.err_JER[Muon] / ttbar.xsec[Muon],
		     ttbar.err_JER[ElMu],  100 * ttbar.err_JER[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" b-tagging             & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_Btag[Elec],  100 * ttbar.err_Btag[Elec] / ttbar.xsec[Elec],
		     ttbar.err_Btag[Muon],  100 * ttbar.err_Btag[Muon] / ttbar.xsec[Muon],
		     ttbar.err_Btag[ElMu],  100 * ttbar.err_Btag[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" mis-tagging           & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_mtag[Elec],  100 * ttbar.err_mtag[Elec] / ttbar.xsec[Elec],
		     ttbar.err_mtag[Muon],  100 * ttbar.err_mtag[Muon] / ttbar.xsec[Muon],
		     ttbar.err_mtag[ElMu],  100 * ttbar.err_mtag[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Pile-Up               & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_PU[Elec],  100 * ttbar.err_PU[Elec] / ttbar.xsec[Elec],
		     ttbar.err_PU[Muon],  100 * ttbar.err_PU[Muon] / ttbar.xsec[Muon],
		     ttbar.err_PU[ElMu],  100 * ttbar.err_PU[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << "\\hline" << endl;
//  fOUTSTREAM << Form(" Top Pt                & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) ", 
//		     ttbar.err_TopPt[Elec],  100 * ttbar.err_TopPt[Elec] / ttbar.xsec[Elec],
//		     ttbar.err_TopPt[Muon],  100 * ttbar.err_TopPt[Muon] / ttbar.xsec[Muon],
//		     ttbar.err_TopPt[ElMu],  100 * ttbar.err_TopPt[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" NLO Generator         & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_TopPt[Elec],  100 * ttbar.err_TopPt[Elec] / ttbar.xsec[Elec],
		     ttbar.err_TopPt[Muon],  100 * ttbar.err_TopPt[Muon] / ttbar.xsec[Muon],
		     ttbar.err_TopPt[ElMu],  100 * ttbar.err_TopPt[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" PDF Uncertainty       & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_pdf[Elec],  100 * ttbar.err_pdf[Elec] / ttbar.xsec[Elec],
		     ttbar.err_pdf[Muon],  100 * ttbar.err_pdf[Muon] / ttbar.xsec[Muon],
		     ttbar.err_pdf[ElMu],  100 * ttbar.err_pdf[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Hadronization         & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_had[Elec],  100 * ttbar.err_had[Elec] / ttbar.xsec[Elec],
		     ttbar.err_had[Muon],  100 * ttbar.err_had[Muon] / ttbar.xsec[Muon],
		     ttbar.err_had[ElMu],  100 * ttbar.err_had[ElMu] / ttbar.xsec[ElMu])<< endl;
  //fOUTSTREAM << Form(" Color Reconnection    & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) ", 
  fOUTSTREAM << Form(" Scale ME              & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_cr[Elec],  100 * ttbar.err_cr[Elec] / ttbar.xsec[Elec],
		     ttbar.err_cr[Muon],  100 * ttbar.err_cr[Muon] / ttbar.xsec[Muon],
		     ttbar.err_cr[ElMu],  100 * ttbar.err_cr[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Scale PS              & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) \\\\", 
		     ttbar.err_Q2[Elec],  100 * ttbar.err_Q2[Elec] / ttbar.xsec[Elec],
		     ttbar.err_Q2[Muon],  100 * ttbar.err_Q2[Muon] / ttbar.xsec[Muon],
		     ttbar.err_Q2[ElMu],  100 * ttbar.err_Q2[ElMu] / ttbar.xsec[ElMu])<< endl;
//  fOUTSTREAM << Form(" Matching partons      & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) & %4.2f (%4.2f \\%) ", 
//		     ttbar.err_Match[Elec],  100 * ttbar.err_Match[Elec] / ttbar.xsec[Elec],
//		     ttbar.err_Match[Muon],  100 * ttbar.err_Match[Muon] / ttbar.xsec[Muon],
//		     ttbar.err_Match[ElMu],  100 * ttbar.err_Match[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << Form(" Total systematic      & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) \\\\", 
		     ttbar.xsec_syst[Elec],  100 * ttbar.xsec_syst[Elec] / ttbar.xsec[Elec],
		     ttbar.xsec_syst[Muon],  100 * ttbar.xsec_syst[Muon] / ttbar.xsec[Muon],
		     ttbar.xsec_syst[ElMu],  100 * ttbar.xsec_syst[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Luminosity            & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) \\\\", 
		     ttbar.xsec_lumi[Elec],  100 * ttbar.xsec_lumi[Elec] / ttbar.xsec[Elec],
		     ttbar.xsec_lumi[Muon],  100 * ttbar.xsec_lumi[Muon] / ttbar.xsec[Muon],
		     ttbar.xsec_lumi[ElMu],  100 * ttbar.xsec_lumi[ElMu] / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << Form(" Statistical           & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) \\\\", 
		     ttbar.xsec_stat[Elec],  100 * ttbar.xsec_stat[Elec] / ttbar.xsec[Elec],
		     ttbar.xsec_stat[Muon],  100 * ttbar.xsec_stat[Muon] / ttbar.xsec[Muon],
		     ttbar.xsec_stat[ElMu],  100 * ttbar.xsec_stat[ElMu] / ttbar.xsec[ElMu])<< endl;
  float total_ee = TMath::Sqrt(ttbar.xsec_syst[Elec]*ttbar.xsec_syst[Elec] + ttbar.xsec_lumi[Elec]*ttbar.xsec_lumi[Elec] + ttbar.xsec_stat[Elec]*ttbar.xsec_stat[Elec]);
  float total_mm = TMath::Sqrt(ttbar.xsec_syst[Muon]*ttbar.xsec_syst[Muon] + ttbar.xsec_lumi[Muon]*ttbar.xsec_lumi[Muon] + ttbar.xsec_stat[Muon]*ttbar.xsec_stat[Muon]);
  float total_em = TMath::Sqrt(ttbar.xsec_syst[ElMu]*ttbar.xsec_syst[ElMu] + ttbar.xsec_lumi[ElMu]*ttbar.xsec_lumi[ElMu] + ttbar.xsec_stat[ElMu]*ttbar.xsec_stat[ElMu]);
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << Form(" TOTAL                 & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) & %4.2f (%3.1f \\%) \\\\", 
		     total_ee,  100 * total_ee / ttbar.xsec[Elec],
		     total_mm,  100 * total_mm / ttbar.xsec[Muon],
		     total_em,  100 * total_em / ttbar.xsec[ElMu])<< endl;
  fOUTSTREAM << "\\hline\\hline" << endl;



  
  fOUTSTREAM << endl;
  fOUTSTREAM << endl;

  fOUTSTREAM << "\\hline\\hline" << endl;
  fOUTSTREAM << "          Source  &            El/El            &            Mu/Mu            &            El/Mu            \\\\" << endl;
  fOUTSTREAM << "\\hline" << endl;
  if (DD) {
    fOUTSTREAM << Form(" Drell-Yan        & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\", 
		       DD_DY.Yields[Elec][cut], DD_DY.Yields_stat[Elec][cut], DD_DY.Yields_syst[Elec][cut][0], 
		       DD_DY.Yields[Muon][cut], DD_DY.Yields_stat[Muon][cut], DD_DY.Yields_syst[Muon][cut][0], 
		       DD_DY.Yields[ElMu][cut], DD_DY.Yields_stat[ElMu][cut], DD_DY.Yields_syst[ElMu][cut][0])
	       << endl;

    fOUTSTREAM << Form(" Non W/Z leptons  & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\",  
		      DD_NonW.Yields[Elec][cut], DD_NonW.Yields_stat[Elec][cut], DD_NonW.Yields_syst[Elec][cut][0], 
		      DD_NonW.Yields[Muon][cut], DD_NonW.Yields_stat[Muon][cut], DD_NonW.Yields_syst[Muon][cut][0], 
		      DD_NonW.Yields[ElMu][cut], DD_NonW.Yields_stat[ElMu][cut], DD_NonW.Yields_syst[ElMu][cut][0])
	       << endl;	
  }
  else {
    fOUTSTREAM << Form(" Drell-Yan        & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\", 
		       DY.Yields[Elec][cut], DY.Yields_stat[Elec][cut], DY.Yields_syst[Elec][cut][0], 
		       DY.Yields[Muon][cut], DY.Yields_stat[Muon][cut], DY.Yields_syst[Muon][cut][0], 
		     DY.Yields[ElMu][cut], DY.Yields_stat[ElMu][cut], DY.Yields_syst[ElMu][cut][0])
	       << endl;
    fOUTSTREAM << Form(" Non W/Z leptons  & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f ",  
		       Fake.Yields[Elec][cut], Fake.Yields_stat[Elec][cut], Fake.Yields_syst[Elec][cut][0], 
		       Fake.Yields[Muon][cut], Fake.Yields_stat[Muon][cut], Fake.Yields_syst[Muon][cut][0], 
		       Fake.Yields[ElMu][cut], Fake.Yields_stat[ElMu][cut], Fake.Yields_syst[ElMu][cut][0])
	       << endl;	
  }						  
  fOUTSTREAM << Form(" Single top quark & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\",  
		     STop.Yields[Elec][cut], STop.Yields_stat[Elec][cut], STop.Yields_syst[Elec][cut][0],
		     STop.Yields[Muon][cut], STop.Yields_stat[Muon][cut], STop.Yields_syst[Muon][cut][0], 
		     STop.Yields[ElMu][cut], STop.Yields_stat[ElMu][cut], STop.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << Form(" Dibosons         & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\",  
		     VV.Yields[Elec][cut], VV.Yields_stat[Elec][cut], VV.Yields_syst[Elec][cut][0], 
		     VV.Yields[Muon][cut], VV.Yields_stat[Muon][cut], VV.Yields_syst[Muon][cut][0], 
		     VV.Yields[ElMu][cut], VV.Yields_stat[ElMu][cut], VV.Yields_syst[ElMu][cut][0])
	     << endl;
  //fOUTSTREAM << Form(" Rare             & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f ",  
//		     Rare.Yields[Elec][cut], Rare.Yields_stat[Elec][cut], Rare.Yields_syst[Elec][cut][0], 
//		     Rare.Yields[Muon][cut], Rare.Yields_stat[Muon][cut], Rare.Yields_syst[Muon][cut][0], 
//		     Rare.Yields[ElMu][cut], Rare.Yields_stat[ElMu][cut], Rare.Yields_syst[ElMu][cut][0])
//	     << endl;
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << Form(" Total background & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\",  
		     Total.Yields[Elec][cut], Total.Yields_stat[Elec][cut], Total.Yields_syst[Elec][cut][0], 
		     Total.Yields[Muon][cut], Total.Yields_stat[Muon][cut], Total.Yields_syst[Muon][cut][0], 
		     Total.Yields[ElMu][cut], Total.Yields_stat[ElMu][cut], Total.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << Form(" TTbar dilepton   & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f & %7.1f $\\pm$ %4.1f $\\pm$ %6.1f \\\\",  
		     TTbar.Yields[Elec][cut], TTbar.Yields_stat[Elec][cut], TTbar.Yields_syst[Elec][cut][0], 
		     TTbar.Yields[Muon][cut], TTbar.Yields_stat[Muon][cut], TTbar.Yields_syst[Muon][cut][0], 
		     TTbar.Yields[ElMu][cut], TTbar.Yields_stat[ElMu][cut], TTbar.Yields_syst[ElMu][cut][0])
	     << endl;
  fOUTSTREAM << "\\hline" << endl;
  fOUTSTREAM << Form(" Data             & %5.0f                       & %5.0f                       & %5.0f                      \\\\", 
		     Data.Yields[Elec][cut], 
		     Data.Yields[Muon][cut], 
		     Data.Yields[ElMu][cut] )<< endl;
  fOUTSTREAM << "\\hline\\hline" << endl;
  fOUTSTREAM << endl;
  fOUTSTREAM << endl;
  fOUTSTREAM << "===================================================================================" << endl;
  fOUTSTREAM << "                           ttbar acceptance x eff. x BR (%)                        " << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << "            El/El         &           Mu/Mu          &            El/Mu            " << endl;
  fOUTSTREAM << "-----------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form("    %7.4f $\\pm$ %5.4f    &    %7.3f $\\pm$ %5.4f    &    %7.3f $\\pm$ %5.4f ",
		     100*ttbar.acc[Elec], 100*TMath::Sqrt(ttbar.acc_stat[Elec]*ttbar.acc_stat[Elec]+ttbar.acc_syst[Elec]*ttbar.acc_syst[Elec]),
		     100*ttbar.acc[Muon], 100*TMath::Sqrt(ttbar.acc_stat[Muon]*ttbar.acc_stat[Muon]+ttbar.acc_syst[Muon]*ttbar.acc_syst[Muon]),
		     100*ttbar.acc[ElMu], 100*TMath::Sqrt(ttbar.acc_stat[ElMu]*ttbar.acc_stat[ElMu]+ttbar.acc_syst[ElMu]*ttbar.acc_syst[ElMu])) 
	     << endl;
  fOUTSTREAM << "===================================================================================" << endl;
  
  fOUTSTREAM << endl;
  fOUTSTREAM << "==============================================================================================================" << endl;
  fOUTSTREAM << "                                         Cross-section measurement (pb)                                       " << endl;
  fOUTSTREAM << "--------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << "                El/El               &               Mu/Mu                &              El/Mu                 " << endl;
  fOUTSTREAM << "--------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" %4.1f $\\pm$ %3.1f $\\pm$ %4.1f $\\pm$ %3.1f  & %4.1f $\\pm$ %3.1f $\\pm$ %4.1f $\\pm$ %3.1f  & %4.1f $\\pm$ %3.1f $\\pm$ %4.1f $\\pm$ %3.1f ",
		     ttbar.xsec[Elec], ttbar.xsec_stat[Elec], ttbar.xsec_syst[Elec], ttbar.xsec_lumi[Elec],
		     ttbar.xsec[Muon], ttbar.xsec_stat[Muon], ttbar.xsec_syst[Muon], ttbar.xsec_lumi[Muon],
		     ttbar.xsec[ElMu], ttbar.xsec_stat[ElMu], ttbar.xsec_syst[ElMu], ttbar.xsec_lumi[ElMu]) 
	     << endl;
  fOUTSTREAM << "==============================================================================================================" << endl;
  fOUTSTREAM << endl;


/*
  fOUTSTREAM << "=========================================================================================================================================================================================================" << endl;
  fOUTSTREAM << "          El/Mu            " << endl;
  fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << "          Source  |  Yields +/- stat +/- syst   | vv    | stop  | rare  | DY    | fakes | IdIso | Trig  | les   | jes   | JER   | btag  | mistag| PU    | TopPt | PDF   | Had   | CR    | Q2    | Match | " << endl;
  fOUTSTREAM << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << endl;
  fOUTSTREAM << Form(" Drell-Yan	  | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     DD_DY.Yields[ElMu][cut]       , DD_DY.Yields_stat[ElMu][cut]   , DD_DY.Yields_syst[ElMu][cut][0], 
		     DD_DY.SystError[ElMu][vv]     , DD_DY.SystError[ElMu][stop]    , DD_DY.SystError[ElMu][rare]    , DD_DY.SystError[ElMu][dy], DD_DY.SystError[ElMu][fake],
		     DD_DY.SystError[ElMu][SFIDISO], DD_DY.SystError[ElMu][SFTrig]  , DD_DY.SystError[ElMu][les]     , DD_DY.SystError[ElMu][jes], 
		     DD_DY.SystError[ElMu][jer]    , DD_DY.SystError[ElMu][btag]    , DD_DY.SystError[ElMu][mistag]  , DD_DY.SystError[ElMu][pu],
		     DD_DY.SystError[ElMu][TopPt]  , DD_DY.SystError[ElMu][pdf]     , DD_DY.SystError[ElMu][had]     , DD_DY.SystError[ElMu][cr]  , 
		     DD_DY.SystError[ElMu][Q2]     , DD_DY.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << Form(" Non W/Z leptons  | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     DD_NonW.Yields[ElMu][cut]       , DD_NonW.Yields_stat[ElMu][cut]   , DD_NonW.Yields_syst[ElMu][cut][0], 
		     DD_NonW.SystError[ElMu][vv]     , DD_NonW.SystError[ElMu][stop]    , DD_NonW.SystError[ElMu][rare]    , DD_NonW.SystError[ElMu][dy], DD_NonW.SystError[ElMu][fake],
		     DD_NonW.SystError[ElMu][SFIDISO], DD_NonW.SystError[ElMu][SFTrig]  , DD_NonW.SystError[ElMu][les]     , DD_NonW.SystError[ElMu][jes], 
		     DD_NonW.SystError[ElMu][jer]    , DD_NonW.SystError[ElMu][btag]    , DD_NonW.SystError[ElMu][mistag]  , DD_NonW.SystError[ElMu][pu],
		     DD_NonW.SystError[ElMu][TopPt]  , DD_NonW.SystError[ElMu][pdf]     , DD_NonW.SystError[ElMu][had]     , DD_NonW.SystError[ElMu][cr]  , 
		     DD_NonW.SystError[ElMu][Q2]     , DD_NonW.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << Form(" Single Top       | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     STop.Yields[ElMu][cut]       , STop.Yields_stat[ElMu][cut]   , STop.Yields_syst[ElMu][cut][0], 
		     STop.SystError[ElMu][vv]     , STop.SystError[ElMu][stop]    , STop.SystError[ElMu][rare]    , STop.SystError[ElMu][dy], STop.SystError[ElMu][fake],
		     STop.SystError[ElMu][SFIDISO], STop.SystError[ElMu][SFTrig]  , STop.SystError[ElMu][les]     , STop.SystError[ElMu][jes], 
		     STop.SystError[ElMu][jer]    , STop.SystError[ElMu][btag]    , STop.SystError[ElMu][mistag]  , STop.SystError[ElMu][pu],
		     STop.SystError[ElMu][TopPt]  , STop.SystError[ElMu][pdf]     , STop.SystError[ElMu][had]     , STop.SystError[ElMu][cr]  , 
		     STop.SystError[ElMu][Q2]     , STop.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << Form(" Diboson          | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     VV.Yields[ElMu][cut]       , VV.Yields_stat[ElMu][cut]   , VV.Yields_syst[ElMu][cut][0], 
		     VV.SystError[ElMu][vv]     , VV.SystError[ElMu][stop]    , VV.SystError[ElMu][rare]    , VV.SystError[ElMu][dy], VV.SystError[ElMu][fake],
		     VV.SystError[ElMu][SFIDISO], VV.SystError[ElMu][SFTrig]  , VV.SystError[ElMu][les]     , VV.SystError[ElMu][jes], 
		     VV.SystError[ElMu][jer]    , VV.SystError[ElMu][btag]    , VV.SystError[ElMu][mistag]  , VV.SystError[ElMu][pu],
		     VV.SystError[ElMu][TopPt]  , VV.SystError[ElMu][pdf]     , VV.SystError[ElMu][had]     , VV.SystError[ElMu][cr]  , 
		     VV.SystError[ElMu][Q2]     , VV.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << Form(" Rare             | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     Rare.Yields[ElMu][cut]       , Rare.Yields_stat[ElMu][cut]   , Rare.Yields_syst[ElMu][cut][0], 
		     Rare.SystError[ElMu][vv]     , Rare.SystError[ElMu][stop]    , Rare.SystError[ElMu][rare]    , Rare.SystError[ElMu][dy], Rare.SystError[ElMu][fake],
		     Rare.SystError[ElMu][SFIDISO], Rare.SystError[ElMu][SFTrig]  , Rare.SystError[ElMu][les]     , Rare.SystError[ElMu][jes], 
		     Rare.SystError[ElMu][jer]    , Rare.SystError[ElMu][btag]    , Rare.SystError[ElMu][mistag]  , Rare.SystError[ElMu][pu],
		     Rare.SystError[ElMu][TopPt]  , Rare.SystError[ElMu][pdf]     , Rare.SystError[ElMu][had]     , Rare.SystError[ElMu][cr]  , 
		     Rare.SystError[ElMu][Q2]     , Rare.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << Form(" TTbar dilepton   | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     TTbar.Yields[ElMu][cut]       , TTbar.Yields_stat[ElMu][cut]   , TTbar.Yields_syst[ElMu][cut][0], 
		     TTbar.SystError[ElMu][vv]     , TTbar.SystError[ElMu][stop]    , TTbar.SystError[ElMu][rare]    , TTbar.SystError[ElMu][dy], TTbar.SystError[ElMu][fake],
		     TTbar.SystError[ElMu][SFIDISO], TTbar.SystError[ElMu][SFTrig]  , TTbar.SystError[ElMu][les]     , TTbar.SystError[ElMu][jes], 
		     TTbar.SystError[ElMu][jer]    , TTbar.SystError[ElMu][btag]    , TTbar.SystError[ElMu][mistag]  , TTbar.SystError[ElMu][pu],
		     TTbar.SystError[ElMu][TopPt]  , TTbar.SystError[ElMu][pdf]     , TTbar.SystError[ElMu][had]     , TTbar.SystError[ElMu][cr]  , 
		     TTbar.SystError[ElMu][Q2]     , TTbar.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << Form(" SUSY stop        | %7.1f +/- %4.1f +/- %6.1f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f | %1.3f |",   
        	     SUSYstop.Yields[ElMu][cut], SUSYstop.Yields_stat[ElMu][cut], SUSYstop.Yields_syst[ElMu][cut][0], 
		     SUSYstop.SystError[ElMu][vv]     , SUSYstop.SystError[ElMu][stop]    , SUSYstop.SystError[ElMu][rare]    , SUSYstop.SystError[ElMu][dy], SUSYstop.SystError[ElMu][fake],
		     SUSYstop.SystError[ElMu][SFIDISO], SUSYstop.SystError[ElMu][SFTrig]  , SUSYstop.SystError[ElMu][les]     , SUSYstop.SystError[ElMu][jes], 
		     SUSYstop.SystError[ElMu][jer]    , SUSYstop.SystError[ElMu][btag]    , SUSYstop.SystError[ElMu][mistag]  , SUSYstop.SystError[ElMu][pu],
		     SUSYstop.SystError[ElMu][TopPt]  , SUSYstop.SystError[ElMu][pdf]     , SUSYstop.SystError[ElMu][had]     , SUSYstop.SystError[ElMu][cr]  , 
		     SUSYstop.SystError[ElMu][Q2]     , SUSYstop.SystError[ElMu][Matching]    )
             << endl;
  fOUTSTREAM << "=========================================================================================================================================================================================================" << endl;

  
  fOUTSTREAM << Form(" tt syst = %3.3f ",   
                     TTbar.Yields[ElMu][cut] * sqrt(
		     TTbar.SystError[ElMu][vv]      * TTbar.SystError[ElMu][vv]      + 
		     TTbar.SystError[ElMu][stop]    * TTbar.SystError[ElMu][stop]    + 
		     TTbar.SystError[ElMu][rare]    * TTbar.SystError[ElMu][rare]    + 
		     TTbar.SystError[ElMu][dy]      * TTbar.SystError[ElMu][dy]      + 
		     TTbar.SystError[ElMu][fake]    * TTbar.SystError[ElMu][fake]    +
		     TTbar.SystError[ElMu][SFIDISO] * TTbar.SystError[ElMu][SFIDISO] + 
		     TTbar.SystError[ElMu][SFTrig]  * TTbar.SystError[ElMu][SFTrig]  +
		     TTbar.SystError[ElMu][les]     * TTbar.SystError[ElMu][les]     + 
		     TTbar.SystError[ElMu][jes]     * TTbar.SystError[ElMu][jes]     + 
		     TTbar.SystError[ElMu][jer]     * TTbar.SystError[ElMu][jer]     + 
		     TTbar.SystError[ElMu][btag]    * TTbar.SystError[ElMu][btag]    + 
		     TTbar.SystError[ElMu][mistag]  * TTbar.SystError[ElMu][mistag]  + 
		     TTbar.SystError[ElMu][pu]      * TTbar.SystError[ElMu][pu]      +
		     TTbar.SystError[ElMu][TopPt]   * TTbar.SystError[ElMu][TopPt]   + 
		     TTbar.SystError[ElMu][pdf]     * TTbar.SystError[ElMu][pdf]     + 
		     TTbar.SystError[ElMu][had]     * TTbar.SystError[ElMu][had]     + 
		     TTbar.SystError[ElMu][cr]      * TTbar.SystError[ElMu][cr]      + 
		     TTbar.SystError[ElMu][Q2]      * TTbar.SystError[ElMu][Q2]      + 
		     TTbar.SystError[ElMu][Matching]* TTbar.SystError[ElMu][Matching]
		     ) )
             << endl;
  
*/  
  fOUTSTREAM.close();
  gSystem->Exec("cat "+filename);
}
void TopPlotter::DrawTopLine(Int_t chan, Float_t y){

  TString htitleCMSChannel;
  if(chan==Muon)      htitleCMSChannel="#mu#mu";
  else if(chan==Elec) htitleCMSChannel="ee";
  else if(chan==ElMu) htitleCMSChannel="e#mu";
  else                htitleCMSChannel="ee#mu#mu";
  
  //TString titlelabel = Form("CMS Preliminary, #sqrt{s}=13 TeV, %4.1f pb^{-1}",fLumiNorm/*/1000.*/);
  TString titlelabel = Form("CMS, %4.2f pb^{-1} (13 TeV) ",fLumiNorm/*/1000.*/);
  TLatex *title  = new TLatex(-20.,50.,titlelabel);
  title->SetNDC();
  title->SetTextAlign(12);
  title->SetX(0.66);
  title->SetY(y);
  title->SetTextFont(42);
  title->SetTextSize(0.04);
  title->SetTextSizePixels(22);
  title->Draw("SAME");
  
  TLatex *chtitle  = new TLatex(-20.,50.,htitleCMSChannel);
  chtitle->SetNDC();
  chtitle->SetTextAlign(12);
  chtitle->SetX(0.85);
  chtitle->SetY(y);
  chtitle->SetTextFont(42);
  chtitle->SetTextSize(0.04);
  chtitle->SetTextSizePixels(22);
  //chtitle->Draw("SAME");
}
TH1F* TopPlotter::GetHisto1D(TFile *file, TString histoname) {
  
  if (!file) {
    std::cerr << "ERROR: Could not load file" << std::endl;
    return 0;
  }
  TH1F* h = (TH1F*) file->Get(histoname)->Clone(histoname);
  if (!h) {
    std::cerr << "ERROR[TopSelector]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
  }
  h->SetDirectory(0);
  //  file->Close();
  
  return h;
};
TH2F* TopPlotter::GetHisto2D(TFile *file, TString histoname) {
  
  if (!file) {
    std::cerr << "ERROR: Could not load file" << std::endl;
    return 0;
  }
  TH2F* h = (TH2F*) file->Get(histoname)->Clone(histoname);
  if (!h) {
    std::cerr << "ERROR[TopSelector]: Could not find histogram " 
	      << histoname << std::endl;
    return 0;
  }
  h->SetDirectory(0);
  //  file->Close();
  
  return h;
};
void TopPlotter::SetupDraw(TH1F* h, Int_t color, Int_t var){
  h->Rebin(GetRebin(h,var));
  
  TString xaxis = KinAxisLabel[var];
  
  h->SetTitle(h->GetTitle());
  h->GetXaxis()->SetTitle(xaxis);
  h->SetLineColor(1);
  h->SetFillColor(color);
  h->SetFillStyle(1001);
  
  h->SetBinContent(h->GetNbinsX(),(h->GetBinContent(h->GetNbinsX()+1)+h->GetBinContent(h->GetNbinsX())));
  h->SetBinContent(h->GetNbinsX()+1,0);
  
  if (var == NBTagsNJets) {  //change bin labels
    h->GetXaxis()->SetBinLabel( 1, "0");
    h->GetXaxis()->SetBinLabel( 2, "0");
    h->GetXaxis()->SetBinLabel( 3, "1");
    h->GetXaxis()->SetBinLabel( 4, "0");
    h->GetXaxis()->SetBinLabel( 5, "1");
    h->GetXaxis()->SetBinLabel( 6, "2");
    h->GetXaxis()->SetBinLabel( 7, "0");
    h->GetXaxis()->SetBinLabel( 8, "1");
    h->GetXaxis()->SetBinLabel( 9, "2");
    h->GetXaxis()->SetBinLabel(10, "3");
    h->GetXaxis()->SetBinLabel(11, "0");
    h->GetXaxis()->SetBinLabel(12, "1");
    h->GetXaxis()->SetBinLabel(13, "2");
    h->GetXaxis()->SetBinLabel(14, "3");
    h->GetXaxis()->SetBinLabel(15, "4");
    
    float bincontent(0.);
    for (int bin=15; bin<h->GetNbinsX()+1; bin++) bincontent += h->GetBinContent(bin);
    h->SetBinContent(15, bincontent);
    
    h->GetXaxis()->SetRangeUser(-0.5,14.5);
  }
  return;
}
Int_t TopPlotter::GetRebin(TH1F* h, Int_t var){
  TString histo = h->GetTitle();
  
  Int_t nbins   = h->GetNbinsX();
  float xmax  = h->GetXaxis()->GetBinUpEdge(nbins);
  float xmin  = h->GetXaxis()->GetBinLowEdge(1);

  Int_t rebin = nbins/(Int_t)(xmax-xmin); // normalize to 1GeV?
  if (var==MET        ) return 20 * rebin;
  if (var==InvMass    ) return 20 * rebin;
  if (var==Lep0Pt     ) return 10 * rebin;
  if (var==Lep1Pt     ) return 10 * rebin;
  if (var==DelLepPhi  ) return 5;
  if (var==AbsDelPhiLeps  ) return 2;
  if (var==delPhi2LeadJets) return 4;
  if (var==minDelRJetsLeps) return 20;
  if (var==NJets      ) return rebin;
  if (var==NBtagJets  ) return rebin;
  if (var==Jet0Pt     ) return 10 * rebin;
  if (var==Jet1Pt     ) return 10 * rebin;
  if (var==Vtx        ) return 1;
  if (var==DelPhillJet) return 100;
  if (var==DiLepPt)     return 100;
  		      
  //  if (var==CSVTag     ) return 25;
  //  if (var==TopD       ) return 25;
  if (var==DelPhillJet) return 10;
  
  return rebin;
}
