//////////////////////////////////////////////////////////////////////////
// Use the definitions in SetPlotter.C and the functions in DrawPlots.C
// Reads the rootfiles and save the histograms in rootfiles prepared
// to run the code of the combine tool
// Write datacard for shape asymptotic limits
// The histograms are saved in LimitHistosFolder and the datacards are
// written in DatacardsFolder
//
// Usage:
// MT2ShapeHistos(TString chan, float somelumi);
// MT2ShapeHistos("ElMu", 1264);
//
////////////////////////////////////////////////////////////////////////////

#include "DrawPlots.C"

void MT2ShapeHistos(TString chan, float somelumi, int mStop, int mLsp);
float getYield(TString rootfile, TString process);
void WriteDatacard(TString rootfile, float systfac);
void WriteCombinationDatacard(TString rootfile, float systfac);
//void CreateDatacard(Int_t mStop, Int_t mLsp, TString chan, float somelumi, bool systfac, float cut);
//void CreateDatacardCombination(Int_t mStop, Int_t mLsp, float somelumi, bool systfac, float cut);

void MT2YieldHistos(TString chan, float somelumi, float cut);
float Get8TeVStopScale(int mStop);
float getxsec8TeV(int StopMass);
float getxsec(int StopMass);

void FillEmptySRvector(TString rootfile);
void WriteDatacardSR(TString rootfile, float systfac);
bool IsEmptySR(TString SR, TString ch, TString rootfile);
void IntegrateEmptySR(TString rootfile);
void merge(TString MergeSF, TString rootfile, TString ch);
void CreateDatacard(Int_t mStop, Int_t mLsp, float somelumi, bool systfac, float cut);


// Constants
const TString LimitHistosFolder = "/mnt_pool/fanae105/user/juanr/CMSSW_7_4_7_patch1/StopLimits2016/rootfiles/";
const TString DatacardFolder = "/mnt_pool/fanae105/user/juanr/CMSSW_7_4_7_patch1/StopLimits2016/datacards/";
const TString mt2 = "MT2";
const TString ulevel = "DYVeto";

const int nSignals = 2;
const TString signals[nSignals] = {"T2tt_500_325", "T2tt_850_100"};

//######################################################################################
//###### Create rootfiles with MT2 Shape Histos ########################################
//######################################################################################

void MT2ShapeHistos(TString chan, float somelumi, int mStop, int mLsp){
  bool docomb = false;
  if(chan.BeginsWith("comb_")){
    docomb = true;
    chan.ReplaceAll("comb_", "");
  }
 TH1F* htt; TH1F* htW; TH1F* hDY; TH1F* hVV; TH1F* httV; TH1F* hdata; TH1F* hWJets;
  //TH1F* hS500_325; TH1F* hS850_100;
  TH1F* hStop;
  TString stopsample = Form("T2tt_mStop%i_mLsp%i", mStop, mLsp);
  htt    = H_ttbar(mt2, chan, ulevel, "0"); htt   ->Scale(somelumi/Lumi); htt   ->SetName("h_tt"      ); StackOverflow(htt);
  htW    = H_tW(   mt2, chan, ulevel, "0"); htW   ->Scale(somelumi/Lumi); htW   ->SetName("h_tW"      ); StackOverflow(htW);
  httV   = H_ttV(  mt2, chan, ulevel, "0"); httV  ->Scale(somelumi/Lumi); httV  ->SetName("h_ttV"     ); StackOverflow(httV);
  hDY    = H_DY(   mt2, chan, ulevel, "0"); hDY   ->Scale(somelumi/Lumi); hDY   ->SetName("h_DY"      ); StackOverflow(hDY);
  hVV    = H_VV(   mt2, chan, ulevel, "0"); hVV   ->Scale(somelumi/Lumi); hVV   ->SetName("h_VV"      ); StackOverflow(hVV);
  hWJets = H_fake( mt2, chan, ulevel, "0"); hWJets->Scale(somelumi/Lumi); hWJets->SetName("h_WJets"   ); StackOverflow(hWJets);
  hdata  = H_data( mt2, chan, ulevel); hdata ->Scale(somelumi/Lumi); hdata ->SetName("h_data_obs"); StackOverflow(hdata);
  hStop = loadHistogram(stopsample, mt2, chan, ulevel); hStop->Scale(somelumi/Lumi); hStop->SetName("h_Stop"); StackOverflow(hStop);
  if(do8TeV)  hStop->Scale(Get8TeVStopScale(mStop));
  TH1F* hlumi = new TH1F("thelumi", "thelumi", 1, 0, 1);
  hlumi->SetBinContent(1, somelumi);

  TFile* t;
  if(docomb) t = new TFile(LimitHistosFolder + "MT2shape" + Form("_%i_%i_", mStop, mLsp) + chan + ".root" , "recreate");
  else  t = new TFile(LimitHistosFolder + Form("MT2shape_%i_%i.root", mStop, mLsp), "recreate");
  htt->Write(); htW->Write(); hDY->Write(); hVV->Write(); httV->Write(); hdata->Write(); hWJets->Write();
  hStop->Write(); hStop->Write();
  hlumi->Write(); delete hlumi;
  delete htt; delete htW; delete hDY; delete hVV; delete httV; delete hdata; delete hWJets;
  delete hStop;
  t->Close(); delete t;
}

// With Signal Regions:
/*
enum SR{ //ADG, 
 //ADG,
 ADH, ADI, AEG, AEH, AEI, AFG, AFH, AFI,
 BDG, BDH, BDI, BEG, BFG, BFH, BFI, BEH, BEI, 
 CDG, CDH, CDI, CEG, CEH, CEI, CFG, CFH, CFI,
 nSR
};
const TString SRlabel[nSR] = { //"ADG", 
 //"ADG", 
 "ADH", "ADI", "AEG", "AEH", "AEI", "AFG", "AFH", "AFI", 
 "BDG", "BDH", "BDI", "BEG", "BFG", "BFH", "BFI", "BEH", "BEI", 
 "CDG", "CDH", "CDI", "CEG", "CEH", "CEI", "CFG", "CFH", "CFI"
};*/

int ESR[200]; 
const int ESRpos[3] = {0, 50, 100}; 

void MT2ShapeSR(TString chan, float somelumi, int mStop, int mLsp, TString SR){
 bool docomb = true;
 TH1F* htt; TH1F* htW; TH1F* hDY; TH1F* hVV; TH1F* httV; TH1F* hdata; TH1F* hWJets;
  //TH1F* hS500_325; TH1F* hS850_100;
  TH1F* hStop;
  TString mt2 = "MT2_" + SR;
  TString stopsample = Form("T2tt_mStop%i_mLsp%i", mStop, mLsp);
  htt    = H_ttbar(mt2, chan, ulevel, "0"); htt   ->Scale(somelumi/Lumi); htt   ->SetName("h_tt"      ); StackOverflow(htt);
  htW    = H_tW(   mt2, chan, ulevel, "0"); htW   ->Scale(somelumi/Lumi); htW   ->SetName("h_tW"      ); StackOverflow(htW);
  httV   = H_ttV(  mt2, chan, ulevel, "0"); httV  ->Scale(somelumi/Lumi); httV  ->SetName("h_ttV"     ); StackOverflow(httV);
  hDY    = H_DY(   mt2, chan, ulevel, "0"); hDY   ->Scale(somelumi/Lumi); hDY   ->SetName("h_DY"      ); StackOverflow(hDY);
  hVV    = H_VV(   mt2, chan, ulevel, "0"); hVV   ->Scale(somelumi/Lumi); hVV   ->SetName("h_VV"      ); StackOverflow(hVV);
  hWJets = H_fake( mt2, chan, ulevel, "0"); hWJets->Scale(somelumi/Lumi); hWJets->SetName("h_WJets"   ); StackOverflow(hWJets);
  hdata  = H_data( mt2, chan, ulevel); hdata ->Scale(somelumi/Lumi); hdata ->SetName("h_data_obs"); StackOverflow(hdata);

  hStop = loadHistogram(stopsample, mt2, chan, ulevel); hStop->Scale(somelumi/Lumi); hStop->SetName("h_Stop"); StackOverflow(hStop);
  if(do8TeV)  hStop->Scale(Get8TeVStopScale(mStop));
  TH1F* hlumi = new TH1F("thelumi", "thelumi", 1, 0, 1);
  hlumi->SetBinContent(1, somelumi);

  TFile* t;
  if(docomb) t = new TFile(LimitHistosFolder + "MT2" + SR + Form("_%i_%i_", mStop, mLsp) + chan + ".root" , "recreate");
  else       t = new TFile(LimitHistosFolder + "MT2" + SR + Form("_%i_%i.root", mStop, mLsp), "recreate");
  htt->Write(); htW->Write(); hDY->Write(); hVV->Write(); httV->Write(); hdata->Write(); hWJets->Write(); hStop->Write(); 
  hlumi->Write(); delete hlumi;
  delete htt; delete htW; delete hDY; delete hVV; delete httV; delete hdata; delete hWJets;
  delete hStop;
  t->Close(); delete t;
}



//######################################################################################
//###### Create Datacards                       ########################################
//######################################################################################
float getYield(TString rootfile, TString process){
  TFile *inputfile = TFile::Open(LimitHistosFolder+rootfile);
  TH1F* hist; inputfile->GetObject("h_" + process, hist);
  float yield = hist->Integral();
  delete inputfile; 
  return yield;
}
/*
void WriteDatacard(TString rootfile, float systfac){
  TString(rootfile).ReplaceAll(".root","");
  rootfile += ".root";
  TFile *inputfile = TFile::Open(LimitHistosFolder+rootfile);
  
  TH1F* lumihist; inputfile->GetObject("thelumi", lumihist);
  float thelumi = lumihist->GetBinContent(1);
  inputfile->Close(); delete inputfile;

  TString filename = DatacardFolder + "/datacard_" + TString(rootfile).ReplaceAll(".root","") + ".txt";
  ofstream outputfile;
  outputfile.open(filename);

  int nchannels = 1;
  int nbkgs = 5;
  int nsyst = 7;
  TString pro[]   = {"Stop", "tt", "tW", "ttV", "DY", "VV"};
  float systunc[] = {20.   , 25. , 25. ,  25. , 25. , 25. };

  outputfile << Form("imax %i\n", nchannels);
  outputfile << Form("jmax %i\n", nbkgs);
  outputfile << Form("kmax %i\n", nsyst);
  outputfile << "-----------\n";
  outputfile << "shapes * * " + LimitHistosFolder + rootfile + " h_$PROCESS $PROCESS_$SYSTEMATIC\n";
  outputfile << "-----------\n"; 
  outputfile << "bin 1\n";
  outputfile << Form("observation %2.3f\n", getYield(rootfile, "data_obs"));
  outputfile << "-----------\n"; 

  TString bin =      "bin             ";
  TString process1 = "process         "; 
  TString process2 = "process         "; 
  TString rate =     "rate            ";
  TString lumisyst = "lumi lnN        ";
  float y;
  for(int i = 0; i<nbkgs+1; i++){
    bin      += "1     ";
    process1 += (pro[i] + "   ");
    process2 += Form("%i     ", i); 
    rate     += Form("%2.3f     ", getYield(rootfile, pro[i]));
    lumisyst += Form("%2.3f     ", 1.0);
  }
  outputfile << bin + "\n";
  outputfile << process1 + "\n";
  outputfile << process2 + "\n";
  outputfile << rate + "\n";
  outputfile << "-----------\n"; 

  outputfile << lumisyst + "\n";
  TString prosyst = "";
  for(int i = 0; i<nbkgs+1; i++){
    prosyst = "";
    for(int j = 0; j<nbkgs+1; j++){
      if(i == j) prosyst += TString(Form("%f",1+systunc[i]/100*systfac)) + "    ";
      else       prosyst += "-    ";
    }
    outputfile << pro[i] + "  lnN     " + prosyst + "\n";
  }
  outputfile.close();
}

void CreateDatacardMT2(Int_t mStop, Int_t mLsp, TString chan, float somelumi, bool systfac, float cut){
  TString rootfile = Form("MT2shape_%i_%i", mStop, mLsp);
  if(chan == "comb"){
     CreateDatacardCombination(mStop, mLsp, somelumi, systfac, cut); 
     return;
  }
  else{
    //if     (rootfile == "MT2shape") MT2ShapeHistos(chan, somelumi, mStop, mLsp);
    //else if(rootfile == "MT2yields" ) MT2YieldHistos(chan, somelumi, cut, mStop, mLsp);
    //else   {cout << "Wrong rootfile name!" << endl; return;} 
    MT2ShapeHistos(chan, somelumi, mStop, mLsp);
    float fact;
    if(systfac) fact = 1/TMath::Sqrt(somelumi/Lumi);
    else        fact = 1;
    WriteDatacard(rootfile, fact);
    return;
  }
}

void WriteCombinationDatacard(TString rootfile, float systfac){
  TString(rootfile).ReplaceAll(".root","");
  //rootfile += "_ElMu.root";
  TFile *inputfile = TFile::Open(LimitHistosFolder+rootfile + "_ElMu.root");
  
  TH1F* lumihist; inputfile->GetObject("thelumi", lumihist);
  float thelumi = lumihist->GetBinContent(1);
  inputfile->Close(); delete inputfile;

  TString filename = DatacardFolder + "/datacard_" + TString(rootfile).ReplaceAll(".root","") + "_comb.txt";
  ofstream outputfile;
  outputfile.open(filename);

  int nchannels = 3;
  int nbkgs = 5;
  int nsyst = 7;
  TString pro[]   = {"Stop", "tt", "tW", "ttV", "DY", "VV"};
  float systunc[] = {20.   , 25. , 25. ,  25. , 25. , 25. };

  outputfile << Form("imax %i\n", nchannels);
  outputfile << Form("jmax %i\n", nbkgs);
  outputfile << Form("kmax %i\n", nsyst);
  outputfile << "-----------\n";
  for(int k = 0; k<3; k++){
    TString ch = chan[k];
    outputfile << "shapes * " + ch + " " + LimitHistosFolder + rootfile + "_" + ch + ".root" + " h_$PROCESS $PROCESS_$SYSTEMATIC\n";
  }
  outputfile << "-----------\n"; 
  TString bin = "bin ";
  TString obs = "observation ";
  for(int k = 0; k<3; k++){
    bin += chan[k] + "  "; 
    obs += Form("   %2.3f  ", getYield(rootfile + "_" + chan[k] + ".root", "data_obs"));
  }
  outputfile << bin + " \n";
  outputfile << obs + " \n";
  outputfile << "-----------\n"; 

  TString bin2 =     "bin             ";
  TString process1 = "process         "; 
  TString process2 = "process         "; 
  TString rate =     "rate            ";
  TString lumisyst = "lumi lnN        ";
  float y;
  for(int k = 0; k<3; k++){
    TString ch = chan[k];
    for(int i = 0; i<nbkgs+1; i++){
      bin2      += ch + "    ";
      process1 += (pro[i] + "   ");
      process2 += Form("%i     ", i); 
      rate     += Form("%2.3f     ", getYield(rootfile + "_" + ch + ".root", pro[i]));
      lumisyst += Form("%2.3f     ", 1.0);
    }
  }
  outputfile << bin2 + "\n";
  outputfile << process1 + "\n";
  outputfile << process2 + "\n";
  outputfile << rate + "\n";
  outputfile << "-----------\n"; 

  outputfile << lumisyst + "\n";
  TString prosyst = "";
  for(int i = 0; i<nbkgs+1; i++){
    prosyst = "";
    for(int j = 0; j<nbkgs+1; j++){
      if(i == j) prosyst += TString(Form("%f",1+systunc[i]/100*systfac)) + "    ";
      else       prosyst += "-    ";
    }
    outputfile << pro[i] + "  lnN     " + prosyst + prosyst + prosyst + "\n";
  }
  outputfile.close();
}

void CreateDatacardCombination(Int_t mStop, Int_t mLsp, float somelumi, bool systfac, float cut){
  TString rootfile = Form("MT2shape_%i_%i", mStop, mLsp);
  for(int k = 0; k < 3; k++) MT2ShapeHistos("comb_" + chan[k], somelumi, mStop, mLsp);
  float fact;
  if(systfac) fact = 1/TMath::Sqrt(somelumi/Lumi);
	else        fact = 1;
	WriteCombinationDatacard(rootfile, fact);
	return;
}
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FillEmptySRvector(TString rootfile){
	for (int i = 0; i < 200; i++) ESR[i] = 0;
	for(int k=0; k<3; k++){
		TString ch = chan[k];
		for(int s = 0; s<nSR; s++){
			TString sr = SRlabel[s];
      if(sr == "AA") continue;
			if(IsEmptySR(sr, ch, rootfile)) ESR[ESRpos[k] + s] = 1;
		}
	}
}

bool IsEmptySR(TString SR, TString ch, TString rootfile){
 // if (SR != "ADG" && SR != "BDG" && SR != "CDG") return 1; else return 0;
	TString(rootfile).ReplaceAll(".root","");
	TString rootname = TString(rootfile);
	rootname.ReplaceAll("MT2", "MT2"+SR);
	int nbkgs = 5;
	TString pro[]   = {"Stop", "tt", "tW", "ttV", "DY", "VV"};
	float y_Stop = 0;
	float y_bkg = 0;
	y_Stop = getYield(rootname + "_" + ch + ".root", "Stop");
	for(int i = 1; i<nbkgs+1; i++){
		y_bkg += getYield(rootname + "_" + ch + ".root", pro[i]);
	}
	if( (y_Stop < 1e-3) || (y_bkg < 5e-2) ){
		//cout << "  Stop yield = " << y_Stop << ";    "; cout << "Bkg yield = " << y_bkg << "; ";
		return 1; //merge(MergeSF, rootname, ch);
	}
	else return 0;
}

void IntegrateEmptySR(TString rootfile){
	TString rootname = TString(rootfile);
	for(int k=0; k<3; k++){
		TString MergeSF = "AB";//"ADH"; // ADH
		TString ch = chan[k];
		for(int s = 0; s<nSR; s++){
      TString sr = SRlabel[s];
      if(sr == "AA") continue;
      if(sr == MergeSF) continue;
      rootname = TString(rootfile);
			rootname.ReplaceAll("MT2", "MT2"+sr);
			if(IsEmptySR(sr, ch, rootfile)){ 
        cout << "   Merging SR " << sr << " into SR " << MergeSF << " for channel " << ch << "... ";
				merge(MergeSF, rootname, ch);
        ESR[ESRpos[k] + s] = 1;
			}
      else MergeSF = sr; 
		}
	}
}

void merge(TString MergeSF, TString rootfile, TString ch){
	TH1F* htt0; TH1F* htW0; TH1F* hDY0; TH1F* hVV0; TH1F* httV0; TH1F* hdata0; TH1F* hWJets0; TH1F* hStop0;
	TH1F* htt;  TH1F* htW;  TH1F* hDY;  TH1F* hVV;  TH1F* httV;  TH1F* hdata;  TH1F* hWJets;  TH1F* hStop;
  TFile* t0 = new TFile(LimitHistosFolder + rootfile + "_" + ch + ".root");
  t0->GetObject("h_tt", htt0);         htt0   ->SetDirectory(0);
  t0->GetObject("h_tW", htW0);         htW0   ->SetDirectory(0);
  t0->GetObject("h_DY", hDY0);         hDY0   ->SetDirectory(0);
  t0->GetObject("h_VV", hVV0);         hVV0   ->SetDirectory(0);
  t0->GetObject("h_ttV", httV0);       httV0  ->SetDirectory(0);
  t0->GetObject("h_WJets", hWJets0);   hWJets0->SetDirectory(0);
  t0->GetObject("h_data_obs", hdata0); hdata0 ->SetDirectory(0);
  t0->GetObject("h_Stop", hStop0);     hStop0 ->SetDirectory(0);
  float y0_bkg  = htt0->Integral() + htW0->Integral() + hDY0->Integral() + hVV0->Integral() + httV0->Integral();
  float y0_stop = hStop0->Integral();
  cout << "   Stop yield = " << y0_stop << ";    "; cout << "Bkg yield = " << y0_bkg << "; ";
  
  if( (y0_bkg == 0) && (y0_stop == 0) ){
    cout << "  Totaly empty SR! Skipping..." << endl;
    return;
  }
  else if( (y0_bkg != 0) && (y0_stop == 0) ){
   // cout << "   No signal in this SR... the bkg will not be merged!" << endl;
     cout << "   Merging SR with no signal... " << endl;
   // return;
  }
  
  TString rootfile2 = TString(rootfile);
  rootfile2.ReplaceAll("MT2" + TString(rootfile[3]) + TString(rootfile[4]), "MT2"+MergeSF);
  TFile* t = new TFile(LimitHistosFolder + rootfile2 + "_" + ch + ".root");
  t->GetObject("h_tt", htt);    
  htt->Add(htt0); htt   ->SetDirectory(0);
  t->GetObject("h_tW", htW);            htW->Add(htW0); htW   ->SetDirectory(0);
  t->GetObject("h_DY", hDY);            hDY->Add(hDY0); hDY   ->SetDirectory(0);
  t->GetObject("h_VV", hVV);            hVV->Add(hVV0); hVV   ->SetDirectory(0);
  t->GetObject("h_ttV", httV);          httV->Add(httV0); httV   ->SetDirectory(0);
  t->GetObject("h_WJets", hWJets);      hWJets->Add(hWJets0); hWJets   ->SetDirectory(0);
  t->GetObject("h_data_obs", hdata);    hdata->Add(hdata0); hdata   ->SetDirectory(0);
  t->GetObject("h_Stop", hStop);        hStop->Add(hStop0); hStop   ->SetDirectory(0);
  TH1F* hlumi; t->GetObject("thelumi", hlumi); hlumi->SetDirectory(0);
  float y_bkg  = htt->Integral() + htW->Integral() + hDY->Integral() + hVV->Integral() + httV->Integral();
  float y_stop = hStop->Integral();
  cout << "   Merged stop yield = " << y_stop << ";    "; cout << "Merged bkg yield = " << y_bkg << endl;


  
  TFile* tf = new TFile(LimitHistosFolder + rootfile2 + "_" + ch + ".root", "recreate");
  htt->Write(); htW->Write(); hDY->Write(); hVV->Write(); httV->Write(); hdata->Write(); hWJets->Write(); hStop->Write(); 
  hlumi->Write();
  delete htt0; delete htW0; delete hDY0; delete hVV0; delete httV0; delete hdata0; delete hWJets0; delete hStop0;
  delete htt;  delete htW;  delete hDY;  delete hVV;  delete httV;  delete hdata;  delete hWJets;  delete hStop;
  t->Close(); delete t; t0->Close(); tf->Close(); delete t0;
}









void WriteDatacardSR(TString rootfile, float systfac){
	TString(rootfile).ReplaceAll(".root","");
	TString rootname = TString(rootfile); 

  FillEmptySRvector(rootfile);
	rootname.ReplaceAll("MT2", "MT2AB");
	TFile *inputfile = TFile::Open(LimitHistosFolder+rootname + "_ElMu.root");

	TH1F* lumihist; inputfile->GetObject("thelumi", lumihist);
	float thelumi = lumihist->GetBinContent(1);
	inputfile->Close(); delete inputfile;

  TString filename = DatacardFolder + "/datacard_SR_" + TString(rootfile).ReplaceAll(".root","") + ".txt";
  ofstream outputfile;
  outputfile.open(filename);

  int nsyst = 7;
  int nbkgs = 5;
  TString pro[]   = {"Stop", "tt", "tW", "ttV", "DY", "VV"};
  float systunc[] = {20.   , 25. , 25. ,  25. , 25. , 25. };

	int nchannels = 0;
	for(int s=0; s<nSR; s++){
    if(SRlabel[s] == "AA") continue;
		for(int k=0; k<3; k++){
			//if(!IsEmptySR(SRlabel[s],chan[k], rootfile)) nchannels++;
			if(!ESR[ESRpos[k]+s]) nchannels++; 
		}
	}
  cout << "   Number of SRs after merge = " << nchannels << endl; 

  outputfile << Form("imax %i\n", nchannels);
	outputfile << Form("jmax %i\n", nbkgs);
	outputfile << Form("kmax %i\n", nsyst);
	outputfile << "-----------\n";
	for(int s = 0; s<nSR; s++){
    if(SRlabel[s] == "AA") continue;
		rootname = TString(rootfile);
		rootname.ReplaceAll("MT2", "MT2"+SRlabel[s]);
		for(int k = 0; k<3; k++){
			TString ch = chan[k];
			//if(IsEmptySR(SRlabel[s],ch, rootfile)) continue;
			if(ESR[ESRpos[k]+s]) continue; 
			outputfile << "shapes * " + ch + SRlabel[s] + " " + LimitHistosFolder + rootname + "_" + ch + ".root" + " h_$PROCESS $PROCESS_$SYSTEMATIC\n";
		}
	}
	outputfile << "-----------\n"; 
	TString bin = "bin ";
  TString obs = "observation ";
	for(int s = 0; s<nSR; s++){
    if(SRlabel[s] == "AA") continue;
		rootname = TString(rootfile);
		rootname.ReplaceAll("MT2", "MT2"+SRlabel[s]);
		for(int k = 0; k<3; k++){
			//if(IsEmptySR(SRlabel[s],chan[k], rootfile)) continue;
			if(ESR[ESRpos[k]+s]) continue; 
			bin += chan[k] + SRlabel[s] + "  "; 
			obs += Form("   %2.3f  ", getYield(rootname + "_" + chan[k] + ".root", "data_obs"));
		}
	}
  outputfile << bin + " \n";
  outputfile << obs + " \n";
  outputfile << "-----------\n"; 

  TString bin2 =     "bin             ";
  TString process1 = "process         "; 
  TString process2 = "process         "; 
  TString rate =     "rate            ";
	TString lumisyst = "lumi lnN        ";
	float y;
	for(int s = 0; s<nSR; s++){
    if(SRlabel[s] == "AA") continue;
		rootname = TString(rootfile);
		rootname.ReplaceAll("MT2", "MT2"+SRlabel[s]);
		for(int k = 0; k<3; k++){
			//if(IsEmptySR(SRlabel[s],chan[k], rootfile)) continue;
			if(ESR[ESRpos[k]+s]) continue; 
			TString ch = chan[k];
			for(int i = 0; i<nbkgs+1; i++){
				bin2      += ch + SRlabel[s] + "    ";
				process1 += (pro[i] + "   ");
				process2 += Form("%i     ", i); 
				rate     += Form("%2.3f     ", getYield(rootname + "_" + ch + ".root", pro[i]));
				lumisyst += Form("%2.3f     ", 1.0);
			}
		}
	}
	outputfile << bin2 + "\n";
	outputfile << process1 + "\n";
	outputfile << process2 + "\n";
  outputfile << rate + "\n";
  outputfile << "-----------\n"; 

  outputfile << lumisyst + "\n";
  TString prosyst = "";
  for(int i = 0; i<nbkgs+1; i++){
    prosyst = "";
    for(int j = 0; j<nbkgs+1; j++){
      if(i == j) prosyst += TString(Form("%f",1+systunc[i]/100*systfac)) + "    ";
      else       prosyst += "-    ";
    }
    TString LongString = "";
    for(int k = 0; k < 3; k ++) for(int s = 0; s<nSR; s++){ 
			if(SRlabel[s] == "AA") continue;
			//if(IsEmptySR(SRlabel[s],chan[k], rootfile)) continue;
			if(ESR[ESRpos[k]+s]) continue; 
			LongString += prosyst;
		}
		outputfile << pro[i] + "  lnN     " + LongString  + "\n";
  }
  outputfile.close();
}

void CreateDatacard(Int_t mStop, Int_t mLsp, float somelumi, bool systfac, float cut){
	TString rootfile = Form("MT2_%i_%i", mStop, mLsp);
  int i = 0;
  cout << "\n  ############################" << endl;
  cout << "  ## (1) Creating rootfiles ##   " << endl;
  cout << "  ############################\n" << endl;
	for(int s = 0; s < nSR; s++){
    if(SRlabel[s] == "AA") continue;
		for(int k = 0; k < 3; k++){
    i++;
    cout << Form("   Rootfile created %i/%i...", i, nSR*3) << endl;
			MT2ShapeSR(chan[k], somelumi, mStop, mLsp, SRlabel[s]);
		}
	}
  cout << "\n  ############################" << endl;
  cout << "  ## (2) Merging empty SRs  ##" << endl;
  cout << "  ############################\n" << endl;
  IntegrateEmptySR(rootfile);
	float fact;
  if(systfac) fact = 1/TMath::Sqrt(somelumi/Lumi);
  else        fact = 1;
  cout << "\n  ############################" << endl;
  cout << "  ## (3)  Writing Datacard  ##" << endl;
  cout << "  ############################\n" << endl;
  WriteDatacardSR(rootfile, fact);
  cout << "\n\n  xxxxxxxxx DONE xxxxxxxxx \n\n" << endl;
  return;
}


























////////////////////////////////////////////////////////////////////////////////////

float Get8TeVStopScale(int mStop){
  float fact = 1;
  fact = (float) getxsec8TeV(mStop)/getxsec(mStop);
  return fact;
}

float getxsec8TeV(int StopMass){
  if (StopMass == 125) return 197.122;
  else if (StopMass == 150) return 80.268;
  else if (StopMass == 175) return 36.7994;
  else if (StopMass == 200) return 18.5245;
  else if (StopMass == 225) return 9.90959;
  else if (StopMass == 250) return 5.57596;
  else if (StopMass == 275) return 3.2781;
  else if (StopMass == 300) return 1.99608;
  else if (StopMass == 325) return 1.25277;
  else if (StopMass == 350) return 0.807323;
  else if (StopMass == 375) return 0.531443;
  else if (StopMass == 400) return 0.35683;
  else if (StopMass == 425) return 0.243755;
  else if (StopMass == 450) return 0.169668;
  else if (StopMass == 475) return 0.119275;
  else if (StopMass == 500) return 0.0855847;
  else if (StopMass == 525) return 0.0618641;
  else if (StopMass == 550) return 0.0452067;
  else if (StopMass == 575) return 0.0333988;
  else if (StopMass == 600) return 0.0248009;
  else if (StopMass == 625) return 0.0185257;
  else if (StopMass == 650) return 0.0139566;
  else if (StopMass == 675) return 0.0106123;
  else if (StopMass == 700) return 0.0081141;
  else if (StopMass == 725) return 0.00623244;
  else if (StopMass == 750) return 0.00506044;
  else if (StopMass == 775) return 0.00372717;
  else if (StopMass == 800) return 0.00289588;
  else if (StopMass == 825) return 0.00226163;
  else if (StopMass == 850) return 0.00176742;
  else if (StopMass == 875) return 0.0013878;
  else if (StopMass == 900) return 0.00109501;
  else if (StopMass == 925) return 0.000866391;
  else if (StopMass == 950) return 0.000687022;
  else if (StopMass == 975) return 0.000546728;
  else cout << "No 8TeV Cross Section for that mass!!" << endl;
  return 1;
}

float getxsec(int StopMass){
  if (StopMass == 125) return 574.981;
  else if (StopMass == 150) return 249.409;
  else if (StopMass == 175) return 121.416;
  else if (StopMass == 200) return 64.5085;
  else if (StopMass == 225) return 36.3818;
  else if (StopMass == 250) return 21.5949;
  else if (StopMass == 275) return 13.3231;
  else if (StopMass == 300) return 8.51615;
  else if (StopMass == 325) return 5.60471;
  else if (StopMass == 350) return 3.78661;
  else if (StopMass == 375) return 2.61162;
  else if (StopMass == 400) return 1.83537;
  else if (StopMass == 425) return 1.31169;
  else if (StopMass == 450) return 0.948333;
  else if (StopMass == 475) return 0.697075;
  else if (StopMass == 500) return 0.51848;
  else if (StopMass == 525) return 0.390303;
  else if (StopMass == 550) return 0.296128;
  else if (StopMass == 575) return 0.226118;
  else if (StopMass == 600) return 0.174599;
  else if (StopMass == 625) return 0.136372;
  else if (StopMass == 650) return 0.107045;
  else if (StopMass == 675) return 0.0844877;
  else if (StopMass == 700) return 0.0670476;
  else if (StopMass == 725) return 0.0536438;
  else if (StopMass == 750) return 0.0431418;
  else if (StopMass == 775) return 0.0348796;
  else if (StopMass == 800) return 0.0283338;
  else if (StopMass == 825) return 0.0241099;
  else if (StopMass == 850) return 0.0189612;
  else if (StopMass == 875) return 0.015625;
  else if (StopMass == 900) return 0.0128895;
  else if (StopMass == 925) return 0.0106631;
  else if (StopMass == 950) return 0.00883465;
  else if (StopMass == 975) return 0.00735655;
  else cout << "No Cross Section for that mass!!" << endl;
  return 1;
}










//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



















void MT2YieldHistos(TString chan, float somelumi, float cut){
  bool docomb = false;
  if(chan.BeginsWith("comb_")){
    docomb = true;
    chan.ReplaceAll("comb_", "");
  }
  TH1F* htt; TH1F* htW; TH1F* hDY; TH1F* hVV; TH1F* httV; TH1F* hdata; TH1F* hWJets;
  TH1F* hS500_325; TH1F* hS850_100;
  htt    = H_ttbar(mt2, chan, ulevel, "0"); htt   ->Scale(somelumi/Lumi); htt   ->SetName("h_tt"      ); StackOverflow(htt   ); int lastbin = htt->GetNbinsX(); float yield;
  htW    = H_tW(   mt2, chan, ulevel, "0"); htW   ->Scale(somelumi/Lumi); htW   ->SetName("h_tW"      ); StackOverflow(htW   ); 
  httV   = H_ttV(  mt2, chan, ulevel, "0"); httV  ->Scale(somelumi/Lumi); httV  ->SetName("h_ttV"     ); StackOverflow(httV  ); 
  hDY    = H_DY(   mt2, chan, ulevel, "0"); hDY   ->Scale(somelumi/Lumi); hDY   ->SetName("h_DY"      ); StackOverflow(hDY   ); 
  hVV    = H_VV(   mt2, chan, ulevel, "0"); hVV   ->Scale(somelumi/Lumi); hVV   ->SetName("h_VV"      ); StackOverflow(hVV   ); 
  hWJets = H_fake( mt2, chan, ulevel, "0"); hWJets->Scale(somelumi/Lumi); hWJets->SetName("h_WJets"   ); StackOverflow(hWJets); 
  hdata  = H_data( mt2, chan, ulevel); hdata ->Scale(somelumi/Lumi); hdata ->SetName("h_data_obs"); StackOverflow(hdata ); 
  hS500_325 = loadHistogram("T2tt_500_325", mt2, chan, ulevel); hS500_325->Scale(somelumi/Lumi); hS500_325->SetName("h_S500_325"); StackOverflow(hS500_325); 
  hS850_100 = loadHistogram("T2tt_850_100", mt2, chan, ulevel); hS850_100->Scale(somelumi/Lumi); hS850_100->SetName("h_S850_100"); StackOverflow(hS850_100); 

  yield = htt   ->Integral(htt   ->FindBin(cut), lastbin); htt   ->Rebin(htt   ->GetNbinsX()); htt   ->SetBinContent(1, yield);
  yield = htW   ->Integral(htW   ->FindBin(cut), lastbin); htW   ->Rebin(htW   ->GetNbinsX()); htW   ->SetBinContent(1, yield);
  yield = httV  ->Integral(httV  ->FindBin(cut), lastbin); httV  ->Rebin(httV  ->GetNbinsX()); httV  ->SetBinContent(1, yield);
  yield = hDY   ->Integral(hDY   ->FindBin(cut), lastbin); hDY   ->Rebin(hDY   ->GetNbinsX()); hDY   ->SetBinContent(1, yield);
  yield = hVV   ->Integral(hVV   ->FindBin(cut), lastbin); hVV   ->Rebin(hVV   ->GetNbinsX()); hVV   ->SetBinContent(1, yield);
  yield = hWJets->Integral(hWJets->FindBin(cut), lastbin); hWJets->Rebin(hWJets->GetNbinsX()); hWJets->SetBinContent(1, yield);
  yield = hdata ->Integral(hdata ->FindBin(cut), lastbin); hdata ->Rebin(hdata ->GetNbinsX()); hdata ->SetBinContent(1, yield);

  yield = hS500_325->Integral(hS500_325->FindBin(cut), lastbin); hS500_325->Rebin(hS500_325->GetNbinsX()); hS500_325->SetBinContent(1, yield);
  yield = hS850_100->Integral(hS850_100->FindBin(cut), lastbin); hS850_100->Rebin(hS850_100->GetNbinsX()); hS850_100->SetBinContent(1, yield);

  TH1F* hlumi = new TH1F("thelumi", "thelumi", 1, 0, 1);
  hlumi->SetBinContent(1, somelumi);

  TFile* t;
  if(docomb) t = new TFile(LimitHistosFolder + "MT2yields_" + chan + ".root", "recreate");
  else  t = new TFile(LimitHistosFolder + "MT2yields.root", "recreate");
  htt->Write(); htW->Write(); hDY->Write(); hVV->Write(); httV->Write(); hdata->Write(); hWJets->Write();
  hS500_325->Write(); hS850_100->Write();
  hlumi->Write(); delete hlumi;
  delete htt; delete htW; delete hDY; delete hVV; delete httV; delete hdata; delete hWJets;
  delete hS850_100; delete hS500_325;
  t->Close(); delete t;
}

