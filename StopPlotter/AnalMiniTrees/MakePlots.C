R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(AnalMiniTree.C+)
R__LOAD_LIBRARY(Plot.C+)
#include "Histo.h"
#include "AnalMiniTree.h"
#include "Plot.h"
void PlotAndDatacard(TString var, TString chan = "ElMu", TString signal = "T2tt_mStop175_mLsp1", TString tag = "0", int setLog = 1);
void MakePlots(TString var = "MT2", TString chan = "ElMu", TString signal = "T2tt_mStop175_mLsp1", TString tag = "0");

TString signals[] = {"T2tt_mStop200_mLsp25","T2tt_mStop225_mLsp50","T2tt_mStop258_mLsp75","T2tt_mStop183_mLsp1","T2tt_mStop208_mLsp25","T2tt_mStop242_mLsp75","T2tt_mStop192_mLsp25","T2tt_mStop217_mLsp50","T2tt_mStop250_mLsp75"};

TString vars[] = {"MT2", "MET", "HT", "DeltaPhi", "DeltaEta", "SR", "SR12", "cutandcount", "ghent", "DeltaPhiLepMet"};

TString channels[] = {"ElMu", "Muon", "Elec", "sameF", "All"};

void MakePlots(TString var, TString chan, TString signal, TString tag){
//	for(int i = 0; i < 10; i++){
  //TString t = signals[i]; t.ReplaceAll("T2tt_mStop",""); t.ReplaceAll("mLsp", "");
//		for(int k = 0; k < 8; k++){
//			for(int j = 0; j < 5; j++){
//				PlotAndDatacard(vars[k], channels[j], signals[i], t);
				PlotAndDatacard("DeltaPhi", "ElMu", "T2tt_mStop175_mLsp1", "pru");
//  		}
//				PlotAndDatacard(vars[k], "ElMu", signals[i], t);
//		}
//	}
//	PlotAndDatacard(var, chan, signal, tag);
}


void PlotAndDatacard(TString var, TString chan, TString signal, TString tag, int setLog){
  Plot *p = new Plot(var, chan, 0, signal);
  p->doSetLogy = false;
  p->DrawStack(tag,1);
  p->MakeDatacard(tag);
  delete p;
}
