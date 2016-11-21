#include "AnalMiniTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalMiniTree::Loop(TString plot){

	TString htitle; Int_t nbins; Float_t binmin = 0; Float_t binmax = 0;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();


	if     (plot == "MET"){
    htitle = "MET";
    nbins = 17;
    float xbins[] = {0, 50, 80, 110, 140, 170, 200, 230, 260, 290, 320, 350, 400, 450, 500, 550, 600}; 
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
    histo->SetTag(theSample, "MET", "MET [GeV]", "e#mu, M_{ll} #geq 20 GeV, #geq 2jets, #geq 1btag");
  }
	else if(plot == "MT2"){ 
    htitle = "MT2";
    nbins = 9;
    float xbins[] = {0, 20, 40, 50, 100, 140, 200, 300, 400}; 
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
  }
  else               { std::cout << " -------> Wrong plot name!! " << endl; return;}


	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		if(chan == "ElMu" && !TIsElMu) continue;
		if(chan == "Elec" && !TIsDoubleElec) continue;
		if(chan == "Muon" && !TIsDoubleMuon) continue;
		if( (chan == "SF" || chan == "sameF" || chan == "SameF") && (!TIsDoubleMuon || !TIsDoubleElec)) continue;

		if     (plot == "SR"){
      if(TMET < 100) histo->Fill(1, TWeight);
      if(TMET > 100) histo->Fill(2, TWeight);
		}
		else if(plot == "SR2"){ 
      if(TMET < 100) histo->Fill(1, TWeight);
      if(TMET > 100) histo->Fill(2, TWeight);
		}
    else if(plot == "MET"){
      if(TMET > 50 && TNJetsBtag > 0 && (TMET/TMath::Sqrt(THT)) > 5 && TMinDPhiMetJets > 0.25)
				histo->Fill(TMET, TWeight);
		}
    else if(plot == "MT2"){
      if(TMET > 50 && TNJetsBtag > 0 && ((chan == "ElMu") || ( (TMET/TMath::Sqrt(THT)) > 5 && TMinDPhiMetJets > 0.25 ) ))
				histo->Fill(TMT2ll, TWeight);
		}

	}
  histo->SetStyle();
}
