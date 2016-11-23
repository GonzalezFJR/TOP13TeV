R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(AnalMiniTree.C+)
R__LOAD_LIBRARY(Plot.C+)
#include "Histo.h"
#include "AnalMiniTree.h"
#include "Plot.h"
void PlotAndDatacard(TString var, TString chan = "ElMu", TString signal = "T2tt_mStop175_mLsp1", TString tag = "0", int setLog = 1);
void MakePlots(TString var = "MT2", TString chan = "ElMu", TString signal = "T2tt_mStop175_mLsp1", TString tag = "0");

void MakePlots(TString var, TString chan, TString signal, TString tag){
  TString vars[] = {"MT2", "MET", "HT", "DeltaPhi", "DeltaEta", "DeltaR", "SR"};
  for(int k = 0; k < 7; k++){
    PlotAndDatacard(vars[k], chan, signal, tag);
  }
 //   PlotAndDatacard(var, chan, signal, tag);
}

void PlotAndDatacard(TString var, TString chan, TString signal, TString tag, int setLog){
  Plot *p = new Plot(var, chan, setLog, signal);
  p->DrawStack(tag,1);
  p->MakeDatacard(tag);
  delete p;
}
