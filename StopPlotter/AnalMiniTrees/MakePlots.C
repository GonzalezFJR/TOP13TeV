#include "Histo.h"
#include "AnalMiniTree.h"
#include "Plot.h"

void MakePlots(){

	gROOT->LoadMacro("Histo.C"); gROOT->LoadMacro("AnalMiniTree.C"); gROOT->LoadMacro("Plot.C");
	Plot *t = new Plot("MET", "ElMu", 1, "T2tt_mStop175_mLsp1");

}
