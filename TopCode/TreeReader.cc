/*
  TreeReader.cpp

  Compilation: gcc -I`root-config --incdir` -o TreeReader.run TreeReader.cc `root-config --libs` 
             
  Arguments: Name of the cut
  Options:   -o Name of the output root file
             -i Name of the input root file

*/

#ifndef __CINT__

#include<string>
#include<iostream>
#include<iomanip>      
#include<sstream>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<set>
#include<vector>

#include "TDirectory.h"
#include "TROOT.h"
#endif

//#include "TKey.h"
#include "TObject.h"
//#include "TPrint.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TSystem.h"
#include "TTree.h"
//#include <exception>
#include <sys/stat.h>
#include "SFIDISOTrigger.h"
#include "TOutputListSelectorDataMap.h"
#include "TStatus.h"
#ifndef __CINT__

using namespace std;

void display_usage()
{
  cout << "\033[1;37musage:\033[1;m skimfile cutindex [options]" << endl;
  cout << "" << endl;
  cout << "Options:" << endl;
  cout << "    -i inputfile.root  Input file" << endl;
  cout << "    -o name in the output file \"h_\"" << endl;
  cout << "    -d Input file directory. Default directory: InputTrees" << endl;
  cout << "    -h                 displays this help message and exits " << endl;
  cout << "" << endl;
}

TString namech[3]      = {"mumu","ee","mue"};
TString titlenamech[3] = {"#mu#mu","ee","e#mu"};
TString namecut[4]     = {"dilepton","2Jets","MET","1btag"};
const float G_Total_Lumi = 5288.;

void printProgress(int entry, const int nentries, TString header){
  int step = nentries/20;
  if( step < 200 )   step = 200;
  if( step > 10000 ) step = 10000;
  
  if(entry%step != 0 && (entry+1 != nentries) ) return;

  float progress_f = (float)(entry+1)/(float)(nentries)*100.;
  char progress[10];
  sprintf(progress, "%5.1f", progress_f);
  cout << " Processing " << setw(50) << left << header << setw(6) << right << progress << " %      \r" << flush;
  if(entry+1 == nentries) cout << endl;
}

int main(int argc, const char* argv[]){
  //  gSystem->Load("libTree");
  
  const char * _output   = 0;
  const char * _input    = 0;
  const char * _dir      = "ReReco_CSVL_Dec13/";
  const char * _prefix   = "Tree_Legacy_";
  
  // Arguments used
  //set<int> usedargs;
  //Parsing input options

  if(argc == 1){
    display_usage();
    return -1;
  }
  
  else{
    //Argumet 1 must be a valid input fileName
    for (int i = 1; i < argc; i++){
      if( strcmp(argv[i],"-i") == 0 ){
	_input = argv[i+1];
	i++;
      }
      if( strcmp(argv[i],"-d") == 0 ){
	_dir = argv[i+1];
	i++;
      }
      if( strcmp(argv[i],"-o") == 0 ){
	_output= argv[i+1];
	i++;
      }
      if( strcmp(argv[i],"-h") == 0 ||
	  strcmp(argv[i],"--help") == 0 ){
	display_usage();
	return 0;
      }
    }
  }//else
  if( _input ==0 ){
    cerr << "\033[1;31mskimfile ERROR:\033[1;m The '-i' option is mandatory!"
	      << endl;
    display_usage();
    return -1;
  }
  
  // reassigning
  TString fname(_input);
  TString hname(_output);
  TString fdir(_dir);
  TString prefix(_prefix);
  
  TFile* file_ = TFile::Open(fdir+prefix+fname+".root");
  TTree* theTree = (TTree*) file_->Get("AnalysisTree"); //Chain theTree("AnalysisTree"); 
  theTree->ResetBranchAddresses();

  cout << "Signal: ";
  cout << fdir+prefix+fname+".root"<< endl;
  
  int Event,Type,PV;
  float PUWeight;
  float MET,MET_Phi,METSig;
  float Lep0Px,Lep0Py,Lep0Pz,Lep0E,Lep1Px,Lep1Py,Lep1Pz,Lep1E;
  float NJets,HT,NBtagJets,Btag_j0,Btag_j1;
  float Jet0Px,Jet0Py,Jet0Pz,Jet0E,Jet0Et,Jet1Px,Jet1Py,Jet1Pz,Jet1E,Jet1Et;
  float BtagJet0Px,BtagJet0Py,BtagJet0Pz,BtagJet0E,BtagJet0Et,BtagJet1Px,BtagJet1Py,BtagJet1Pz,BtagJet1E,BtagJet1Et;

  float tPx,tPy,tPz,tE;
  float tbarPx,tbarPy,tbarPz,tbarE;

  float SysVarJet0,SysVarJet1,SysVarBtagJet0,SysVarBtagJet1;
  float UncJet0,UncJet1,UncBtagJet0,UncBtagJet1;
  
  theTree->SetBranchAddress( "TEvent",   &Event );
  theTree->SetBranchAddress( "TWeight",  &PUWeight );
  theTree->SetBranchAddress( "TChannel", &Type );
  theTree->SetBranchAddress( "TNPV",     &PV );

  theTree->SetBranchAddress( "TMET",     &MET );
  theTree->SetBranchAddress( "TMET_Phi", &MET_Phi );
  theTree->SetBranchAddress( "TMETSig",  &METSig );

  theTree->SetBranchAddress( "TLep0Px", &Lep0Px );
  theTree->SetBranchAddress( "TLep0Py", &Lep0Py );
  theTree->SetBranchAddress( "TLep0Pz", &Lep0Pz );
  theTree->SetBranchAddress( "TLep0E",  &Lep0E );

  theTree->SetBranchAddress( "TLep1Px", &Lep1Px );
  theTree->SetBranchAddress( "TLep1Py", &Lep1Py );
  theTree->SetBranchAddress( "TLep1Pz", &Lep1Pz );
  theTree->SetBranchAddress( "TLep1E",  &Lep1E );

  theTree->SetBranchAddress( "TNJets",    &NJets );
  theTree->SetBranchAddress( "THT",       &HT );
  theTree->SetBranchAddress( "TNJetsBtag",&NBtagJets );
  theTree->SetBranchAddress( "TBtagJet0", &Btag_j0 );
  theTree->SetBranchAddress( "TBtagJet1", &Btag_j1 );

  theTree->SetBranchAddress( "TJet0Px", &Jet0Px );
  theTree->SetBranchAddress( "TJet0Py", &Jet0Py );
  theTree->SetBranchAddress( "TJet0Pz", &Jet0Pz );
  theTree->SetBranchAddress( "TJet0E", &Jet0E );
  theTree->SetBranchAddress( "TJet0Et", &Jet0Et );

  theTree->SetBranchAddress( "TJet1Px", &Jet1Px );
  theTree->SetBranchAddress( "TJet1Py", &Jet1Py );
  theTree->SetBranchAddress( "TJet1Pz", &Jet1Pz );
  theTree->SetBranchAddress( "TJet1E", &Jet1E );  
  theTree->SetBranchAddress( "TJet1Et", &Jet1Et );

  theTree->SetBranchAddress( "TBtagJet0Px", &BtagJet0Px );
  theTree->SetBranchAddress( "TBtagJet0Py", &BtagJet0Py );
  theTree->SetBranchAddress( "TBtagJet0Pz", &BtagJet0Pz );
  theTree->SetBranchAddress( "TBtagJet0E", &BtagJet0E );
  theTree->SetBranchAddress( "TBtagJet0Et", &BtagJet0Et );

  theTree->SetBranchAddress( "TBtagJet1Px", &BtagJet1Px );
  theTree->SetBranchAddress( "TBtagJet1Py", &BtagJet1Py );
  theTree->SetBranchAddress( "TBtagJet1Pz", &BtagJet1Pz );
  theTree->SetBranchAddress( "TBtagJet1E", &BtagJet1E );
  theTree->SetBranchAddress( "TBtagJet1Et", &BtagJet1Et );

  if(fname.Contains("TT")){
    cout << "Reading pT reweight branches....." << endl;
    theTree->SetBranchAddress( "TtPx", &tPx );
    theTree->SetBranchAddress( "TtPy", &tPy );
    theTree->SetBranchAddress( "TtPz", &tPz );
    theTree->SetBranchAddress( "TtE", &tE );

    theTree->SetBranchAddress( "TtbarPx", &tbarPx );
    theTree->SetBranchAddress( "TtbarPy", &tbarPy );
    theTree->SetBranchAddress( "TtbarPz", &tbarPz );
    theTree->SetBranchAddress( "TtbarE",  &tbarE );
  }
  else{
    tPx=0.0;
    tPy=0.0;
    tPz=0.0;
    tE =0.0;

    tbarPx=0.0;
    tbarPy=0.0;
    tbarPz=0.0;
    tbarE =0.0;
  }

  if(fname.Contains("JES") || fname.Contains("JER")){
    theTree->SetBranchAddress( "TSysVarJet0", &SysVarJet0 );
    theTree->SetBranchAddress( "TSysVarJet1", &SysVarJet1 );
    theTree->SetBranchAddress( "TSysVarBtagJet0", &SysVarBtagJet0 );
    theTree->SetBranchAddress( "TSysVarBtagJet1", &SysVarBtagJet1 );

    theTree->SetBranchAddress( "TUncJet0", &UncJet0 );
    theTree->SetBranchAddress( "TUncJet1", &UncJet1 );
    theTree->SetBranchAddress( "TUncBtagJet0", &UncBtagJet0 );
    theTree->SetBranchAddress( "TUncBtagJet1", &UncBtagJet1 );
  }

  else {
    SysVarJet0=1.0;
    SysVarJet1=1.0;
    SysVarBtagJet0=1.0;
    SysVarBtagJet1=1.0;

    UncJet0=1.0;
    UncJet1=1.0;
    UncBtagJet0=1.0;
    UncBtagJet1=1.0;
  }  
  cout << "--- Processing: " << theTree->GetEntries() << " events" << endl;
  
  /*********************************
             Histograms
  **********************************/
  TH1F* hPV[4][3];
  TH1F* hMET[4][3];
  TH1F* hMET_Phi[4][3];
  TH1F* hMETSig[4][3];
  TH1F* hInvMass[4][3];
  TH1F* hDiLepPt[4][3]; 
  TH1F* hLep0Pt[4][3];
  TH1F* hLep1Pt[4][3];
  TH1F* hLep0Eta[4][3];
  TH1F* hLep1Eta[4][3];
  TH1F* hDelLepPhi[4][3];
  TH1F* hLep0Phi[4][3];
  TH1F* hLep1Phi[4][3];
  TH1F* hNJets[4][3];
  TH1F* hHT[4][3];
  TH1F* hNBtagJets[4][3];
  TH1F* hNBtagJets5[4][3];
  TH1F* hJet0Pt[4][3];
  TH1F* hJet1Pt[4][3];
  TH1F* hBtagJet0Pt[4][3];
  
  TH1F* hDYIn[4][3];
  TH1F* hDYOut[4][3];
  TH2F* hDY_NPV[4][3];
  TH2F* hDY_NJets[4][3];
  TH1F* hDY_TF_InvMass[4][3];
  
  TH1F* hSFIDISO[4][3];    TH1F* hSFIDISOError[4][3];
  TH1F* hSFIDISO_mu[4];    TH1F* hSFIDISOError_mu[4];
  TH1F* hSFIDISO_e[4];     TH1F* hSFIDISOError_e[4];
  TH1F* hSFTrigger[4][3];  TH1F* hSFTriggerError[4][3];
  
  TH1F* hSysVarJet0[4][3]        ;TH1F* hUncJet0[4][3];
  TH1F* hSysVarJet1[4][3]        ;TH1F* hUncJet1[4][3];
  TH1F* hSysVarBtagJet0[4][3]    ;TH1F* hUncBtagJet0[4][3];
  
  TH2F *hSysVarJet0Pt[4][3]      ;TH2F* hUncJet0Pt[4][3];
  TH2F *hSysVarJet1Pt[4][3]      ;TH2F* hUncJet1Pt[4][3];
  TH2F *hSysVarBtagJet0Pt[4][3]  ;TH2F* hUncBtagJet0Pt[4][3];
  TH2F *h2DSFIDISO_mu[4]  ;       TH2F* h2DSFIDISO_e[4];
  
  for(int j=0; j<4; j++){
    for(int i=0; i<3; i++){
      hPV[j][i]         = new TH1F("hPV_"+namech[i]+"_"+namecut[j],"PV Distribution  " + titlenamech[i],30,0,30);
      
      //hMET[j][i]        = new TH1F("hMET_"+namech[i]+"_"+namecut[j],"#slash{E}_{T} " + titlenamech[i],80,0,400);
      hMET[j][i]        = new TH1F("hMET_"+namech[i]+"_"+namecut[j],"#slash{E}_{T} " + titlenamech[i],30,0,150);
      hMET_Phi[j][i]    = new TH1F("hMET_Phi_"+namech[i]+"_"+namecut[j],"#Phi_{#slash{E}_{T}} " + titlenamech[i],160,-4,4);
      hMETSig[j][i]     = new TH1F("hMETSig_"+namech[i]+"_"+namecut[j],"#slash{E}_{T} Significance " + titlenamech[i],50,0,150);
      
      if(j==0) hInvMass[j][i]    = new TH1F("hInvMass_"+namech[i]+"_"+namecut[j],"Invariant Mass " + titlenamech[i],60,0.0,300.0);
      else hInvMass[j][i]    = new TH1F("hInvMass_"+namech[i]+"_"+namecut[j],"Invariant Mass " + titlenamech[i],60,6.0,306.0);
      hDiLepPt[j][i]    = new TH1F("hDiLepPt_"+namech[i]+"_"+namecut[j],"Dilepton p_{T} " + titlenamech[i],50,0,250);
      hLep0Pt[j][i]     = new TH1F("hLep0Pt_"+namech[i]+"_"+namecut[j],"p_{T} leading lepton " + titlenamech[i],50,0.0,250.0);
      hLep1Pt[j][i]     = new TH1F("hLep1Pt_"+namech[i]+"_"+namecut[j],"p_{T} second leading lepton " + titlenamech[i],50,0.0,250.0);
      hLep0Eta[j][i]    = new TH1F("hLep0Eta_"+namech[i]+"_"+namecut[j],"#eta_{Lep_{0}} " + titlenamech[i],160,-8,8);
      hLep1Eta[j][i]    = new TH1F("hLep1Eta_"+namech[i]+"_"+namecut[j],"#eta_{Lep_{1}} " + titlenamech[i],160,-8,8);
      hDelLepPhi[j][i]  = new TH1F("hDelLepPhi_"+namech[i]+"_"+namecut[j],"#Delta#phi_{lep} " + titlenamech[i],100,-5,5);
      hLep0Phi[j][i]    = new TH1F("hLep0Phi_"+namech[i]+"_"+namecut[j],"#phi_{Lep_{0}} " + titlenamech[i],100,-5,5);
      hLep1Phi[j][i]    = new TH1F("hLep1Phi_"+namech[i]+"_"+namecut[j],"#phi_{Lep_{1}} " + titlenamech[i],100,-5,5);
      
      hNJets[j][i]      = new TH1F("hNJets_"+namech[i]+"_"+namecut[j],"Jet multiplicity " + titlenamech[i],5,-0.5,4.5);
      hNJets[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNJets[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNJets[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNJets[j][i]->GetXaxis()->SetBinLabel(4,"3");
      hNJets[j][i]->GetXaxis()->SetBinLabel(5,"#geq 4");
 
      hNBtagJets5[j][i]  = new TH1F("hNBtagJets5_"+namech[i]+"_"+namecut[j],"b-tag jet multiplicity " + titlenamech[i],5,-0.5,4.5);
      hNBtagJets5[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNBtagJets5[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNBtagJets5[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNBtagJets5[j][i]->GetXaxis()->SetBinLabel(4,"3");
      hNBtagJets5[j][i]->GetXaxis()->SetBinLabel(5,"#geq 4");
 
      hNBtagJets[j][i]  = new TH1F("hNBtagJets_"+namech[i]+"_"+namecut[j],"b-tag jet multiplicity " + titlenamech[i],4,-0.5,3.5);
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(1,"0");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(2,"1");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(3,"2");
      hNBtagJets[j][i]->GetXaxis()->SetBinLabel(4,"#geq 3");
 
 
      hHT[j][i]         = new TH1F("hHT_"+namech[i]+"_"+namecut[j],"H_{T} " + titlenamech[i],300,0,600);
      hJet0Pt[j][i]     = new TH1F("hJet0Pt_"+namech[i]+"_"+namecut[j],"p_{T} leading jet " + titlenamech[i],50,0,250);
      hJet1Pt[j][i]     = new TH1F("hJet1Pt_"+namech[i]+"_"+namecut[j],"p_{T} Second leading jet " + titlenamech[i],50,0,250);
      hBtagJet0Pt[j][i] = new TH1F("hBtagJet0Pt_"+namech[i]+"_"+namecut[j],"p_{T} leading b-tagged jet " + titlenamech[i],50,0,250);
      
      
      /***************************
                 DY-DD
      ***************************/
      hDYIn[j][i]           = new TH1F("hDYIn_"+namech[i]+"_"+namecut[j],"DY In " + titlenamech[i],100,0,1000);    
      hDYOut[j][i]          = new TH1F("hDYOut_"+namech[i]+"_"+namecut[j],"DY Out " + titlenamech[i],100,0,1000);    
 
      hDY_NPV[j][i]         = new TH2F("hDY_vs_NPV_"+namech[i]+"_"+namecut[j],"DY vs NPV" + titlenamech[i],750,0.,750.,31,-0.5,30.5);    
      hDY_NJets[j][i]       = new TH2F("hDY_vs_NJets_"+namech[i]+"_"+namecut[j],"DY vs NJets" + titlenamech[i],750,0.,750.,5,-0.5,4.5);
 
      hDY_TF_InvMass[j][i]  = new TH1F("hDY_TF_InvMass_"+namech[i]+"_"+namecut[j],"(DY-TF) Invariant Mass " + titlenamech[i],60,6.0,306.0);
 
      /***************************
          SF(ID,ISO & Trigger)
      ***************************/
      hSFIDISO[j][i]           = new TH1F("hSFIDISO_"+namech[i]+"_"+namecut[j],"SF_{ID,ISO} " + titlenamech[i],400,0.8,1.2);    
      hSFIDISOError[j][i]      = new TH1F("hSFIDISOError_"+namech[i]+"_"+namecut[j],"#Delta SF_{ID,ISO} " + titlenamech[i],400,0,0.05); 
      hSFTrigger[j][i]         = new TH1F("hSFTrigger_"+namech[i]+"_"+namecut[j],"SF^{Trigger} " + titlenamech[i],400,0.8,1.2);    
      hSFTriggerError[j][i]    = new TH1F("hSFTriggerError_"+namech[i]+"_"+namecut[j],"#Delta SF^{Trigger} " + titlenamech[i],400,0,0.05);
 
 
      /***************************
              Systematics
      ***************************/
      hSysVarJet0[j][i]       = new TH1F("hSysVarJet0_"+namech[i]+"_"+namecut[j],"Systematic Variation: Jet0 " + titlenamech[i],100,0.9,1.1);    
      hSysVarJet1[j][i]       = new TH1F("hSysVarJet1_"+namech[i]+"_"+namecut[j],"Systematic Variation: Jet1 " + titlenamech[i],100,0.9,1.1);    
      hSysVarBtagJet0[j][i]   = new TH1F("hSysVarBtagJet0_"+namech[i]+"_"+namecut[j],"Systematic Variation: BtagJet0 " + titlenamech[i],100,0.9,1.1);    
      
      hUncJet0[j][i]          = new TH1F("hUncJet0_"+namech[i]+"_"+namecut[j],"Jet0 Uncertainty " + titlenamech[i],100,0.0,0.1);    
      hUncJet1[j][i]          = new TH1F("hUncJet1_"+namech[i]+"_"+namecut[j],"Jet1 Uncertainty " + titlenamech[i],100,0.0,0.1);    
      hUncBtagJet0[j][i]      = new TH1F("hUncBtagJet0_"+namech[i]+"_"+namecut[j],"BtagJet0 Uncertainty " + titlenamech[i],100,0.0,0.1);
 
      hSysVarJet0Pt[j][i]     = new TH2F("hSysVarJet0Pt_"+namech[i]+"_"+namecut[j],"Systematic Variation Vs P_{t}^{Jet0} " + titlenamech[i],120,0.0,120.0,100,0.9,1.1);    
      hSysVarJet1Pt[j][i]     = new TH2F("hSysVarJet1Pt_"+namech[i]+"_"+namecut[j],"Systematic Variation Vs P_{t}^{Jet1} " + titlenamech[i],120,0.0,120.0,100,0.9,1.1);    
      hSysVarBtagJet0Pt[j][i] = new TH2F("hSysVarBtagJet0Pt_"+namech[i]+"_"+namecut[j],"Systematic Variation Vs P_{t}^{BtagJet0} " + titlenamech[i],120,0.0,120,100,0.9,1.1);    
 
      hUncJet0Pt[j][i]        = new TH2F("hUncJet0Pt_"+namech[i]+"_"+namecut[j],"Jet Uncertainty Vs P_{t}^{Jet0} " + titlenamech[i],120,0.0,120.0,100,0.0,0.1);    
      hUncJet1Pt[j][i]        = new TH2F("hUncJet1Pt_"+namech[i]+"_"+namecut[j],"Jet Uncertainty Vs P_{t}^{Jet1} " + titlenamech[i],120,0.0,120.0,100,0.0,0.1);    
      hUncBtagJet0Pt[j][i]    = new TH2F("hUncBtagJet0Pt_"+namech[i]+"_"+namecut[j],"Jet Uncertainty Vs P_{t}^{BtagJet0} " + titlenamech[i],120,0.0,120.0,100,0.0,0.1);    
 
    }//for(i)
    
 
    hSFIDISO_mu[j]        = new TH1F("hSFIDISO_mu_"+namecut[j],"SF_{ID,ISO}^{#mu} ",400,0.9,1.1);    
    hSFIDISOError_mu[j]   = new TH1F("hSFIDISOError_mu_"+namecut[j],"#Delta SF_{ID,ISO}^{mu} ",300,0,0.03); 
 
    hSFIDISO_e[j]         = new TH1F("hSFIDISO_e_"+namecut[j],"SF_{ID,ISO}^{e} ",400,0.9,1.1);    
    hSFIDISOError_e[j]    = new TH1F("hSFIDISOError_e_"+namecut[j],"#Delta SF_{ID,ISO}^{e} ",300,0,0.03); 
 
 
    float hbinmuEta[7]={-2.4 , -1.2 , -0.9 , 0.0 , 0.9 , 1.2 , 2.4}; // muon Eta values
    float hbinmuPt[5] ={20 , 30 , 40 , 60 , 200}; // muon Pt values      
    
    h2DSFIDISO_mu[j]    = new TH2F("h2DSFIDISO_muon_"+namecut[j],"SF_{ID,ISO}^{#mu}(#eta,p_{t}) ",6,hbinmuEta,4,hbinmuPt);
    
    float hbineEta[11]={-2.5 , -2.0 , -1.556 , -1.442 , -0.8 , 0.0 , 0.8 , 1.442 , 1.556 , 2.0 , 2.5}; // electron Eta values
    float hbinePt[5] ={20 , 30 , 40 , 50 , 200}; // electron Pt values
    
    h2DSFIDISO_e[j]    = new TH2F("h2DSFIDISO_e_"+namecut[j],"SF_{ID,ISO}^{e}(#eta,p_{t}) ",10,hbineEta,4,hbinePt);
    
  }//for(j)

  TLorentzVector Lep0,Lep1;
  TLorentzVector Jet0,Jet1,BtagJet0,BtagJet1;

  TStopwatch sw;
  sw.Start(kTRUE);

  ///////////////////////////////////////
  // Please, IGNORE. Temporal solution //
  ///////////////////////////////////////
  //  TCanvas *mydummycanvas=new TCanvas();// 
  ///////////////////////////////////////
  
  
  /*******************
   SF Parametrization
  ******************/
  TH2F *hmuIDSF, *hmumuTriggerSF;
  TH2F *heIDSF, *heeTriggerSF;
  TH2F *hmueTriggerSF;  
 
  // Lepton SFs: ID and Iso w/ stat. + syst. Errors
  TFile *MuonSF = TFile::Open("/mnt_pool/fanae105/user/folgueras/TOP/LeptonSF/ReReco_SF_Mu.root"); // ReReco
  TFile *ElecSF = TFile::Open("/mnt_pool/fanae105/user/folgueras/TOP/LeptonSF/ReReco_SF_EG.root"); // ReReco
  TFile *ElMuSF = TFile::Open("/mnt_pool/fanae105/user/folgueras/TOP/LeptonSF/ReReco_trigger_SF_emu.root"); // ReReco  
  if(!MuonSF){ cerr << "ERROR [MuonSF]: Could not open file " << MuonSF << "!"  << endl;  }
  if(!ElecSF){ cerr << "ERROR [ElecSF]: Could not open file " << ElecSF << "!"  << endl;  }
  if(!ElMuSF){ cerr << "ERROR [MuEGSF]: Could not open file " << ElMuSF << "!"  << endl;  }
  
  hmuIDSF = (TH2F*) MuonSF->Get("GlobalSF")->Clone("muIDSF");
  hmumuTriggerSF = (TH2F*) MuonSF->Get("scalefactor eta2d with syst")->Clone("mumuTriggerSF"); // ReReco
  if(!hmuIDSF || !hmumuTriggerSF){
    cerr << "ERROR [MuonSF]: Could not find histogram for SF reweighting" << endl;
  }
  heIDSF = (TH2F*) ElecSF->Get("GlobalSF")->Clone("eIDSF");
  //heeTriggerSF = (TH2F*) ElecSF->Get("scalefactor_eta2d_with_syst")->Clone("eeTriggerSF"); // PromptReco
  heeTriggerSF = (TH2F*) ElecSF->Get("scalefactor eta2d with syst")->Clone("eeTriggerSF"); // ReReco
  if(!heIDSF || !heeTriggerSF){
    cerr << "ERROR [ElectronSF]: Could not find histogram for SF reweighting" << endl;
  }
  //hmueTriggerSF = (TH2F*) ElMuSF->Get("scalefactor_eta2d_with_syst")->Clone("mueTriggerSF"); // PromptReco
  hmueTriggerSF = (TH2F*) ElMuSF->Get("scalefactor eta2d with syst")->Clone("mueTriggerSF"); // ReReco
  if(!hmueTriggerSF){
    cerr << "ERROR [MuonElectronSF]: Could not find histogram for SF reweighting" << endl;
  }
  
  /********************************************
   Non W/Z Scale Factors/bin(btag) (MET)
  *********************************************/
  //             [bin][cut][channel]
  float NonWZ_DD_bin[4][4][3];
  
  for(int jbin=0; jbin<5; jbin++){
    for(int jcut=0; jcut<5; jcut++){
      for(int jchannel=0; jchannel<4; jchannel++){
	NonWZ_DD_bin[jbin][jcut][jchannel]=1.0; // NO SF
      }
    }
  }

  // [bin][cut][channel]
  // v8 04/07/13 (mumu)
  // MET
  //NonWZ_DD_bin[0][2][0]=3.99;
  //NonWZ_DD_bin[1][2][0]=4.75;
  // v8 04/07/13 (ee)
  // MET
  //NonWZ_DD_bin[0][2][1]=0.73;
  //NonWZ_DD_bin[1][2][1]=0.25;
  // v8 04/07/13 (mue)
  // MET
  //NonWZ_DD_bin[0][2][2]=2.26;
  //NonWZ_DD_bin[1][2][2]=3.35;



  /********************************************
   DY-DD Scale Factors/bin(btag) and Cut level
  *********************************************/
  //             [bin][cut][channel]
  float DY_DD_bin[4][4][3];
  
  for(int jbin=0; jbin<5; jbin++){
    for(int jcut=0; jcut<5; jcut++){
      for(int jchannel=0; jchannel<4; jchannel++){
	if(jcut==0) DY_DD_bin[jbin][jcut][jchannel]=1.0; // NO DY-DD at Dilepton level
	else DY_DD_bin[jbin][jcut][jchannel]=0.0;
      }
    }
  }
  // [bin][cut][channel]
  // v7 03/06/13 CF Old (mumu  1.312 & 1.572 & 1.572)
  // 2Jets
  DY_DD_bin[0][1][0]=1.017;
  DY_DD_bin[1][1][0]=1.312;
  // MET
  DY_DD_bin[0][2][0]=1.320;
  DY_DD_bin[1][2][0]=1.572;
  // 1btag
  DY_DD_bin[0][3][0]=1.320;
  DY_DD_bin[1][3][0]=1.572;
  
  // v7 03/06/13 CF Old (ee   1.355 & 1.664 & 1.664)
  // 2Jets
  DY_DD_bin[0][1][1]=1.070;
  DY_DD_bin[1][1][1]=1.355;
  // MET
  DY_DD_bin[0][2][1]=1.440;
  DY_DD_bin[1][2][1]=1.664;
  // 1btag
  DY_DD_bin[0][3][1]=1.440;
  DY_DD_bin[1][3][1]=1.664;

  // // [bin][cut][channel]
  // // NewJEC 01/05/13 CF Old (mumu)
  // // 2Jets
  // DY_DD_bin[0][1][0]=0.990;
  // DY_DD_bin[1][1][0]=1.225;
  // // MET
  // DY_DD_bin[0][2][0]=1.353;
  // DY_DD_bin[1][2][0]=1.570;
  // // 1btag
  // DY_DD_bin[0][3][0]=1.353;
  // DY_DD_bin[1][3][0]=1.570;
  
  // // NewJEC 01/05/13 CF Old (ee)
  // // 2Jets
  // DY_DD_bin[0][1][1]=1.094;
  // DY_DD_bin[1][1][1]=1.313;
  // // MET
  // DY_DD_bin[0][2][1]=1.546;
  // DY_DD_bin[1][2][1]=1.740;
  // // 1btag
  // DY_DD_bin[0][3][1]=1.546;
  // DY_DD_bin[1][3][1]=1.740;
  
  // // Pure MC (mue)
  // // 2Jets
  // DY_DD_bin[0][1][2]=1.0;
  // DY_DD_bin[1][1][2]=1.0;
  // // MET
  // DY_DD_bin[0][2][2]=1.0;
  // DY_DD_bin[1][2][2]=1.0;
  // // // 1btag
  // DY_DD_bin[0][3][2]=1.0;
  // DY_DD_bin[1][3][2]=1.0;

  // // NewJEC 01/05/13 (mue)
  // // 2Jets
  // DY_DD_bin[0][1][2]=1.01178; // + 0.30;
  // DY_DD_bin[1][1][2]=2.99075; // + 0.90;
  // // MET
  // DY_DD_bin[0][2][2]=1.01178; // + 0.30;
  // DY_DD_bin[1][2][2]=2.99075; // + 0.90;
  // // // 1btag
  // DY_DD_bin[0][3][2]=1.01178; // + 0.30;
  // DY_DD_bin[1][3][2]=2.99075; // + 0.90;

  // // Andreas' TEST 21/05/13 (mue)
  // // 2Jets
  // DY_DD_bin[0][1][2]=1.01178; // + 0.30;
  // DY_DD_bin[1][1][2]=1.8; // + 0.90;
  // // MET
  // DY_DD_bin[0][2][2]=1.01178; // + 0.30;
  // DY_DD_bin[1][2][2]=1.8; // + 0.90;
  // // // 1btag
  // DY_DD_bin[0][3][2]=1.01178; // + 0.30;
  // DY_DD_bin[1][3][2]=1.8; // + 0.90;

  ////////////////////////////////////////////////

  // // From R_{out/in} 03/06/13 (mue)
  // // 2Jets
  // DY_DD_bin[0][1][2]=1.04; // + 0.0;
  // DY_DD_bin[1][1][2]=1.27; // + 0.0;
  // // MET
  // DY_DD_bin[0][2][2]=1.04; // + 0.0;
  // DY_DD_bin[1][2][2]=1.27; // + 0.0;
  // // // 1btag
  // DY_DD_bin[0][3][2]=1.04; // + 0.0;
  // DY_DD_bin[1][3][2]=1.27; // + 0.0;

  // From R_{out/in}/ReReco-ReRecoSF Partial Exercise 30/10/13 (mue)
  // 2Jets
  DY_DD_bin[0][1][2]=1.04; // + 0.0;
  DY_DD_bin[1][1][2]=1.30; // + 0.0;
  // MET
  DY_DD_bin[0][2][2]=1.04; // + 0.0;
  DY_DD_bin[1][2][2]=1.30; // + 0.0;
  // // 1btag
  DY_DD_bin[0][3][2]=1.04; // + 0.0;
  DY_DD_bin[1][3][2]=1.30; // + 0.0;

  // TTJets 15/05/13 (mue)
  // 2Jets
  // DY_DD_bin[0][1][2]=0.904529; // + 0.30;
  // DY_DD_bin[1][1][2]=1.25063; // + 0.90;
  // // MET
  // DY_DD_bin[0][2][2]=0.904529; // + 0.30;
  // DY_DD_bin[1][2][2]=1.25063; // + 0.90;
  // // // 1btag
  // DY_DD_bin[0][3][2]=0.904529; // + 0.30;
  // DY_DD_bin[1][3][2]=1.25063; // + 0.90;

  // // Global 01/05/13 (mue)
  // // 2Jets
  // DY_DD_bin[0][1][2]=1.70;
  // DY_DD_bin[1][1][2]=1.70;
  // // MET
  // DY_DD_bin[0][2][2]=1.70;
  // DY_DD_bin[1][2][2]=1.70;
  // // 1btag
  // DY_DD_bin[0][3][2]=1.70;
  // DY_DD_bin[1][3][2]=1.70;
  
  // //////////////////////////////
  // // Without 2 Jet Cut
  // // NewJEC 01/05/13 CF Old (mumu)
  // // 2Jets
  // DY_DD_bin[0][1][0]=1.075;
  // DY_DD_bin[1][1][0]=1.175;
  // // MET
  // DY_DD_bin[0][2][0]=1.580;
  // DY_DD_bin[1][2][0]=1.595;
  // // 1btag
  // DY_DD_bin[0][3][0]=1.580;
  // DY_DD_bin[1][3][0]=1.595;
  
  // // NewJEC 01/05/13 CF Old (ee)
  // // 2Jets
  // DY_DD_bin[0][1][1]=1.154;
  // DY_DD_bin[1][1][1]=1.273;
  // // MET
  // DY_DD_bin[0][2][1]=1.695;
  // DY_DD_bin[1][2][1]=1.820;
  // // 1btag
  // DY_DD_bin[0][3][1]=1.695;
  // DY_DD_bin[1][3][1]=1.820;
  
  // // NewJEC 01/05/13 (mue)
  // // 2Jets
  // DY_DD_bin[0][1][2]=1.048;
  // DY_DD_bin[1][1][2]=1.820;
  // // MET
  // DY_DD_bin[0][2][2]=1.048;
  // DY_DD_bin[1][2][2]=1.820;
  // // 1btag
  // DY_DD_bin[0][3][2]=1.048;
  // DY_DD_bin[1][3][2]=1.820;
  // //////////////////////////////
  
    
  // Number de events for acceptance
  //          [Cut][Channel]
  int AccEvent[4][3]={0,0,0,0,
		      0,0,0,0,
		      0,0,0,0};
  // Number de events for acceptance
  //          [Cut][Channel]
  float EffEvent[4][3]={0,0,0,0,
			0,0,0,0,
			0,0,0,0};
  
  Long64_t nentries = theTree->GetEntriesFast();
  for (Long64_t ievt=0; ievt<theTree->GetEntriesFast();ievt++) {
    printProgress(ievt,nentries, fname);
    //    if (ievt%10000 == 0) cout << "--- ... Processing event: " << ievt << endl;
    theTree->GetEntry(ievt);
    
    /*************************
      Acceptance: Generation
    **************************/
    if(Type==3 || Type==4 || Type==5 || Type==-5) continue; // RECO Process
    
    // if(Type==0 || Type==1 || Type==2 || Type==-2) continue; // GEN Process
    // else if(Type== 3) Type=0;
    // else if(Type== 4) Type=1;
    // else if(Type== 5) Type=2;
    // else if(Type==-5) Type=-2;
    
    //emu=2;mue=-2
    int channel=Type; // SF/bin
    if(Type==-2) Type=2;// Fill Histograms MuonElectron=ElectronMuon    
    
    Lep0.SetPxPyPzE(Lep0Px,Lep0Py,Lep0Pz,Lep0E);
    Lep1.SetPxPyPzE(Lep1Px,Lep1Py,Lep1Pz,Lep1E);
    
    float InvMass=(Lep0+Lep1).M();       
 
    Jet0.SetPxPyPzE(Jet0Px*SysVarJet0,Jet0Py*SysVarJet0,Jet0Pz,Jet0E);
    Jet1.SetPxPyPzE(Jet1Px*SysVarJet1,Jet1Py*SysVarJet1,Jet1Pz,Jet1E);
    
    BtagJet0.SetPxPyPzE(BtagJet0Px*SysVarBtagJet0,BtagJet0Py*SysVarBtagJet0,BtagJet0Pz,BtagJet0E);
    BtagJet1.SetPxPyPzE(BtagJet1Px*SysVarBtagJet1,BtagJet1Py*SysVarBtagJet1,BtagJet1Pz,BtagJet1E);
    
    bool ZVeto=false;
    bool ZMass=false;
    if(InvMass<76. || InvMass>106.) ZVeto=true;
    if(InvMass>20.)                 ZMass=true;
      
    
    /*******************************************
   Trigger,ID & ISO Scale Factors/bin(Pt,Eta)
    *******************************************/
    
    vector<float> SF_ID_ISO_Tr;
  
    if (fname.Contains("DoubleMu") || fname.Contains("DoubleElectron") || fname.Contains("MuEG")){
      SF_ID_ISO_Tr.push_back(1.0);//SF_ID_ISO_Tr    [0] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_ID_ISO       [1] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_ID_ISO_Error [2] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Tr           [3] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Tr_Error     [4] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Lep0         [5] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Lep0_Error   [6] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Lep1         [7] 
      SF_ID_ISO_Tr.push_back(1.0);//SF_Lep1_Error   [8] 
    
      PUWeight = 1.;
    }
    else { // MC 
      SFIDISOTrigger(SF_ID_ISO_Tr,
 		     Lep0, Lep1,channel,
 		     hmuIDSF, hmumuTriggerSF,
 		     heIDSF, heeTriggerSF,
 		     hmueTriggerSF);
    
 
      PUWeight=PUWeight*SF_ID_ISO_Tr[0]*G_Total_Lumi; 
    
      // Systematic Uncertainties muon-electron: ELECTRON
      // if(channel==2)  PUWeight=PUWeight*( (SF_ID_ISO_Tr[5]-SF_ID_ISO_Tr[6])*(SF_ID_ISO_Tr[7]) ); //emu
      // if(channel==-2) PUWeight=PUWeight*( (SF_ID_ISO_Tr[5])*(SF_ID_ISO_Tr[7]-SF_ID_ISO_Tr[8]) ); //emu
 
      // Systematic Uncertainties muon-electron: MUON
      // if(channel==2)  PUWeight=PUWeight*( (SF_ID_ISO_Tr[5])*(SF_ID_ISO_Tr[7]-SF_ID_ISO_Tr[8]) ); //emu
      // if(channel==-2) PUWeight=PUWeight*( (SF_ID_ISO_Tr[5]-SF_ID_ISO_Tr[6])*(SF_ID_ISO_Tr[7]) ); //emu
 
    }// else(Contain("Data"))
  
 
    /*******************************************
      DY Data Driven: Rout/in and Template Fit
    *******************************************/
    int cut=-1;
 
    //(DY_DD)/bin
    float btagDYUp=100.0;
    float btagDYDown=-1.0;
    if(hname.Contains("Bin0"))   btagDYUp=1.0;
    if(hname.Contains("Binge1")) btagDYDown=0.0;
 
    if(ZMass && (NBtagJets>btagDYDown && NBtagJets<btagDYUp))                                                    cut=0; // Dileptons
    if(ZMass && NJets>1 && (NBtagJets>btagDYDown && NBtagJets<btagDYUp))                                         cut=1; // + 2 Jets
    if(ZMass && NJets>1 && (Type==2 || MET>40.0) && (NBtagJets>btagDYDown && NBtagJets<btagDYUp))                cut=2; // + MET
    if(ZMass && NJets>1 && (Type==2 || MET>40.0) && NBtagJets>0 && (NBtagJets>btagDYDown && NBtagJets<btagDYUp)) cut=3; // + 1 Btag
    
    
    for(int icut=0; icut<cut+1; icut++){
      // CF mumu and ee
      if(icut == 0)            hDY_NPV[icut][Type]  ->Fill(InvMass,PV,PUWeight);
      if(icut == 0 && MET<10 ) hDY_NJets[icut][Type]->Fill(InvMass,NJets,PUWeight);
      
      // Histograms DY-DD (Rout/in mumu and ee)
      if(ZVeto)  hDYOut[icut][Type]->Fill(InvMass,PUWeight);
      else       hDYIn[icut][Type] ->Fill(InvMass,PUWeight);
      
      // Histograms DY-TF (All Channels)
      if(icut==0 || Type==2 || ZVeto) hDY_TF_InvMass[icut][Type]->Fill(InvMass,PUWeight);
    }
    
 
    /***************************
            Selection
    ***************************/
    cut=-1;    
    if(ZMass)                                                               cut=0; // Dilepton
    if(ZMass && NJets>1 && (Type==2 || (ZVeto)))                            cut=1; // + 2Jets and ZVeto 
    if(ZMass && NJets>1 && (Type==2 || (ZVeto && MET>40.0)))                cut=2; // + MET
    if(ZMass && NJets>1 && (Type==2 || (ZVeto && MET>40.0)) && NBtagJets>0) cut=3; // + 1 Btag
 
 
    float DiLepPt=(Lep0+Lep1).Pt();       
    
    float DeltaPhiLep=Lep0.DeltaPhi(Lep1);
    float DeltaEtaLep=Lep0.Eta()-Lep1.Eta();
 
 
    /*******************
      Fill Histograms
    ******************/
     
    float PUWeightNoDY=PUWeight;
    
    for(int icut=0; icut<cut+1; icut++){
      
      // PUWeight reset for each cut 
      PUWeight=PUWeightNoDY;
 
      /************************************************************
       DY-DD Scale Factors/bin(btag) and Cut level and Pt reweight
      *************************************************************/
      if(fname.Contains("ZDY")){
 	if(NBtagJets==0) PUWeight=PUWeight*DY_DD_bin[0][icut][Type];
 	if(NBtagJets>0)  PUWeight=PUWeight*DY_DD_bin[1][icut][Type];
      }
 
      if(fname.Contains("WJets") || fname.Contains("Bkg")){
 	//if(NBtagJets==0) PUWeight=PUWeight*NonWZ_DD_bin[0][icut][Type];
 	//if(NBtagJets>0)  PUWeight=PUWeight*NonWZ_DD_bin[1][icut][Type];
      }
 
      if(fname.Contains("TT")){
 	TLorentzVector t,tbar;
 	t.SetPxPyPzE(tPx,tPy,tPz,tE);
 	tbar.SetPxPyPzE(tbarPx,tbarPy,tbarPz,tbarE);
 
 	// Definition in Gorner's slides
 	// https://indico.cern.ch/materialDisplay.py?contribId=1&materialId=slides&confId=254297 
 	float SF_tPt=sqrt( exp(0.156-0.00137*t.Pt()) * exp(0.156-0.00137*tbar.Pt()) );
 
 	//PUWeight=PUWeight*SF_tPt;   
 
      }

      /*************************************************************/
 
      hSFIDISO[icut][Type]->Fill(SF_ID_ISO_Tr[1],PUWeight);
      hSFIDISOError[icut][Type]->Fill(SF_ID_ISO_Tr[2],PUWeight);
      hSFTrigger[icut][Type]->Fill(SF_ID_ISO_Tr[3],PUWeight);
      hSFTriggerError[icut][Type]->Fill(SF_ID_ISO_Tr[4],PUWeight);
 
      if(channel==0){
 	hSFIDISO_mu[icut]->Fill(SF_ID_ISO_Tr[5],PUWeight);
 	hSFIDISOError_mu[icut]->Fill(SF_ID_ISO_Tr[6],PUWeight);
 	hSFIDISO_mu[icut]->Fill(SF_ID_ISO_Tr[7],PUWeight);
 	hSFIDISOError_mu[icut]->Fill(SF_ID_ISO_Tr[8],PUWeight);
 
 	h2DSFIDISO_mu[icut]->Fill(Lep0.Eta(),Lep0.Pt(),PUWeight);
 	h2DSFIDISO_mu[icut]->Fill(Lep1.Eta(),Lep1.Pt(),PUWeight);
      }
      if(channel==1){
 	hSFIDISO_e[icut]->Fill(SF_ID_ISO_Tr[5],PUWeight);
 	hSFIDISOError_e[icut]->Fill(SF_ID_ISO_Tr[6],PUWeight);
 	hSFIDISO_e[icut]->Fill(SF_ID_ISO_Tr[7],PUWeight);
 	hSFIDISOError_e[icut]->Fill(SF_ID_ISO_Tr[8],PUWeight);
 
 	h2DSFIDISO_e[icut]->Fill(Lep0.Eta(),Lep0.Pt(),PUWeight);
 	h2DSFIDISO_e[icut]->Fill(Lep1.Eta(),Lep1.Pt(),PUWeight);
      }
      if(channel==2){
 	hSFIDISO_e[icut]->Fill(SF_ID_ISO_Tr[5],PUWeight);
 	hSFIDISOError_e[icut]->Fill(SF_ID_ISO_Tr[6],PUWeight);
 	hSFIDISO_mu[icut]->Fill(SF_ID_ISO_Tr[7],PUWeight);
 	hSFIDISOError_mu[icut]->Fill(SF_ID_ISO_Tr[8],PUWeight);
 
 	h2DSFIDISO_e[icut] ->Fill(Lep0.Eta(),Lep0.Pt(),PUWeight);
 	h2DSFIDISO_mu[icut]->Fill(Lep1.Eta(),Lep1.Pt(),PUWeight);
      }
      if(channel==-2){
 	hSFIDISO_mu[icut]->Fill(SF_ID_ISO_Tr[5],PUWeight);
 	hSFIDISOError_mu[icut]->Fill(SF_ID_ISO_Tr[6],PUWeight);
 	hSFIDISO_e[icut]->Fill(SF_ID_ISO_Tr[7],PUWeight);
 	hSFIDISOError_e[icut]->Fill(SF_ID_ISO_Tr[8],PUWeight);
 
 	h2DSFIDISO_mu[icut]->Fill(Lep0.Eta(),Lep0.Pt(),PUWeight);
 	h2DSFIDISO_e[icut] ->Fill(Lep1.Eta(),Lep1.Pt(),PUWeight);
      }
      
      /******************
          Acceptace
      ******************/
      AccEvent[icut][Type]++;

      //TTbar Madgraph
      //EffEvent[icut][Type]+=PUWeight*(6923750/(234*5300))*(0.108*9)*(0.108*9);//No Luminosity
      //EffEvent[icut][Type]+=PUWeight*(0.108*9)*(0.108*9);//N_events
      //TTbar MC@NLO
      //EffEvent[icut][Type]+=PUWeight*(32842589/(234*5300));//No Luminosity
      //TTbar MC@NLO NoSC
      //EffEvent[icut][Type]+=PUWeight*(21916326/(234*5300));//No Luminosity
      //TTbar Powheg
      //EffEvent[icut][Type]+=PUWeight*(21675970/(234*5300));//No Luminosity
      //TTJetsMG
      //EffEvent[icut][Type]+=PUWeight*((12119013/0.11)/(234*5300));//No Luminosity
      //EffEvent[icut][Type]+=PUWeight*(2067614.0/(24.6*5311.0));//(245.0/234.0)//N_Events
      //EffEvent[icut][Type]+=PUWeight*(12011428.0/(24.6*5311.0));//(245.0/234.0);//N_Events
//SANTI     if(fname.Contains("TT")) EffEvent[icut][Type]+=PUWeight*(245.0/234.0)*(0.108*9)*(0.108*9);
//SANTI     else 

      EffEvent[icut][Type]+=PUWeight;
 
      hPV[icut][Type]->Fill(PV,PUWeight);
      
      hMET[icut][Type]->Fill(MET,PUWeight);
      hMET_Phi[icut][Type]->Fill(MET_Phi,PUWeight);
      hMETSig[icut][Type]->Fill(METSig,PUWeight);
      
      hInvMass[icut][Type]->Fill(InvMass,PUWeight);
      hDiLepPt[icut][Type]->Fill(DiLepPt,PUWeight);
      hLep0Pt[icut][Type]->Fill(Lep0.Pt(),PUWeight);
      hLep1Pt[icut][Type]->Fill(Lep1.Pt(),PUWeight);
      hLep0Eta[icut][Type]->Fill(Lep0.Eta(),PUWeight);
      hLep1Eta[icut][Type]->Fill(Lep1.Eta(),PUWeight);
      hDelLepPhi[icut][Type]->Fill(DeltaPhiLep,PUWeight);
      hLep0Phi[icut][Type]->Fill(Lep0.Phi(),PUWeight);
      hLep1Phi[icut][Type]->Fill(Lep1.Phi(),PUWeight);
      
      if(NJets>4) hNJets[icut][Type]->Fill(4,PUWeight);
      else hNJets[icut][Type]->Fill(NJets,PUWeight); 
      
      if(NBtagJets>4) hNBtagJets5[icut][Type]->Fill(4,PUWeight);
      else hNBtagJets5[icut][Type]->Fill(NBtagJets,PUWeight);
 
      if(NBtagJets>3) hNBtagJets[icut][Type]->Fill(3,PUWeight);
      else hNBtagJets[icut][Type]->Fill(NBtagJets,PUWeight);
 
      hHT[icut][Type]->Fill(HT,PUWeight);
      hJet0Pt[icut][Type]->Fill(Jet0.Pt(),PUWeight);
      hJet1Pt[icut][Type]->Fill(Jet1.Pt(),PUWeight);
      hBtagJet0Pt[icut][Type]->Fill(BtagJet0.Pt(),PUWeight);
 
       /***************************
          Systematic Histograms
       ***************************/
      hUncJet0[icut][Type]->Fill(UncJet0,PUWeight);
      hUncJet1[icut][Type]->Fill(UncJet1,PUWeight);
      hUncBtagJet0[icut][Type]->Fill(UncBtagJet0,PUWeight);
  
      hSysVarJet0[icut][Type]->Fill(SysVarJet0,PUWeight);
      hSysVarJet1[icut][Type]->Fill(SysVarJet1,PUWeight);
      hSysVarBtagJet0[icut][Type]->Fill(SysVarBtagJet0,PUWeight);
  
      hUncJet0Pt[icut][Type]->Fill(Jet0.Pt(),UncJet0);
      hUncJet1Pt[icut][Type]->Fill(Jet1.Pt(),UncJet1);
      //      hUncBtagJet0Pt[icut][Type] ->Fill(BtagJet0.Pt(),UncBtagJet0);
      
      hSysVarJet0Pt[icut][Type]->Fill(Jet0.Pt(),SysVarJet0);
      hSysVarJet1Pt[icut][Type]->Fill(Jet1.Pt(),SysVarJet1);
      hSysVarBtagJet0Pt[icut][Type]->Fill(BtagJet0.Pt(),SysVarBtagJet0);
  
    }//for(icuts)     
  }//for(events)
 
  // Get elapsed time
  sw.Stop();
  cout << "--- End of event loop: "; sw.Print();
  
  //Acceptance-Efficiency
  cout << "--------  Acceptace  --------" << endl;
  cout << "Number of RAW-mumu events:" << endl;
  cout << "Dilepton: " << AccEvent[0][0] << endl;
  cout << "2 Jets: "   << AccEvent[1][0] << endl;
  cout << "MET: "      << AccEvent[2][0] << endl;
  cout << "1 btag: "   << AccEvent[3][0] << endl;

  cout << "--------  Efficiency  --------" << endl;
  cout << "Number of Weigthed-mumu events:" << endl;
  cout << "Dilepton: " << EffEvent[0][0] << endl;
  cout << "2 Jets: "   << EffEvent[1][0] << endl;
  cout << "MET: "      << EffEvent[2][0] << endl;
  cout << "1 btag: "   << EffEvent[3][0] << endl;
  
  cout << "--------  Acceptace  --------" << endl;
  cout << "Number of RAW-ee events:" << endl;
  cout << "Dilepton: " << AccEvent[0][1] << endl;
  cout << "2 Jets: "   << AccEvent[1][1] << endl;
  cout << "MET: "      << AccEvent[2][1] << endl;
  cout << "1 btag: "   << AccEvent[3][1] << endl;

  cout << "--------  Efficiency  --------" << endl;
  cout << "Number of Weigthed-ee events: " << endl;
  cout << "Dilepton: " << EffEvent[0][1] << endl;
  cout << "2 Jets: "   << EffEvent[1][1] << endl;
  cout << "MET: "      << EffEvent[2][1] << endl;
  cout << "1 btag: "   << EffEvent[3][1] << endl;
  
  cout << "--------  Acceptace  --------" << endl;
  cout << "Number of RAW-mue events:" << endl;
  cout << "Dilepton: " << AccEvent[0][2] << endl;
  cout << "2 Jets: "   << AccEvent[1][2] << endl;
  cout << "MET: "      << AccEvent[2][2] << endl;
  cout << "1 btag: "   << AccEvent[3][2] << endl;

  cout << "--------  Efficiency  --------" << endl;
  cout << "Number of Weighted-mue events:" << endl;
  cout << "Dilepton: " << EffEvent[0][2] << endl;
  cout << "2 Jets: "   << EffEvent[1][2] << endl;
  cout << "MET: "      << EffEvent[2][2] << endl;
  cout << "1 btag: "   << EffEvent[3][2] << endl;


  // --- Write histograms
  TString dirname="TopResults_5fb_CSVL";   
  // make a dir if it does not exist!!
  system("mkdir -p " + dirname);
  
  TString outfname=dirname + "/" + hname + "_" + fname + ".root";
  TFile *target  = new TFile(outfname,"RECREATE" );  
  
  for(int j=0; j<4; j++){
    for(int i=0; i<3; i++){
      hPV[j][i]       ->Write();
      hMET[j][i]      ->Write();
      hMET_Phi[j][i]  ->Write();
      hMETSig[j][i]   ->Write();

      hInvMass[j][i]  ->Write();
      hDiLepPt[j][i]  ->Write();
      hLep0Pt[j][i]   ->Write();
      hLep1Pt[j][i]   ->Write();
      hLep0Eta[j][i]  ->Write();
      hLep1Eta[j][i]  ->Write();
      hDelLepPhi[j][i]->Write();
      hLep0Phi[j][i]  ->Write();
      hLep1Phi[j][i]  ->Write();
  
       hNJets[j][i]     ->Write();
       hHT[j][i]        ->Write();
       hNBtagJets5[j][i]->Write();
       hNBtagJets[j][i] ->Write();
       hJet0Pt[j][i]    ->Write();    
       hJet1Pt[j][i]    ->Write();    
       hBtagJet0Pt[j][i]->Write();    
             
       hDYIn[j][i]->Write();
       hDYOut[j][i]->Write();
       hDY_NPV[j][i]->Write();
       hDY_NJets[j][i]->Write();
       hDY_TF_InvMass[j][i]->Write();
  
       hSFIDISO[j][i]->Write();
       hSFIDISOError[j][i]->Write();
       hSFTrigger[j][i]->Write();
       hSFTriggerError[j][i]->Write();
  
       hUncJet0[j][i]->Write();
       hUncJet1[j][i]->Write();
       hUncBtagJet0[j][i]->Write();
  
       hSysVarJet0[j][i]->Write();
       hSysVarJet1[j][i]->Write();
       hSysVarBtagJet0[j][i]->Write();
  
       hUncJet0Pt[j][i]->Write();
       hUncJet1Pt[j][i]->Write();
       //       hUncBtagJet0Pt[j][i]->Write();
  
       hSysVarJet0Pt[j][i]->Write();
       hSysVarJet1Pt[j][i]->Write();
       hSysVarBtagJet0Pt[j][i]->Write();
    }//for(i)
 
    hSFIDISO_mu[j]->Write();
    hSFIDISO_e[j]->Write();
    hSFIDISOError_mu[j]->Write();
    hSFIDISOError_e[j]->Write();
    
    h2DSFIDISO_mu[j]->Write();
    h2DSFIDISO_e[j]->Write();
  }//for(j)
 
  target->Write();
  target->Close();
  cout << "File saved as " << outfname << endl;
}
#endif

