#include "AnalMiniTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalMiniTree::Loop(TString plot){

	TString htitle; Int_t nbins; Float_t binmin = 0; Float_t binmax = 0;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
  //############ SRs ############
  if(plot == "SRMT23D"){
    htitle = "SR_MT2_3D";
    nbins = 3*3*3;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 27, -0.5, 27-0.5);
    histo->SetTag(theSample, "SRs", "SR (MT2 3D)", "");
  }
  else if(plot == "SR"){
    htitle = "SR";
    nbins = 37;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 37, -0.5, 37-0.5);
    histo->SetTag(theSample, "SRs", "SR", "");
  }
  else if(plot == "SR12"){
    htitle = "SR12";
    nbins = 12;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 12, 0.5, 12+0.5);
    histo->SetTag(theSample, "SRs", "SR", "");
  }
  else if(plot == "cutandcount"){
    htitle = "CutAndCount";
    nbins = 1;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 1, 0.5, 1.5);
    histo->SetTag(theSample, "CutAndCount", "CutAndCount", "");
  }
  else if(plot == "ghent"){
    htitle = "Ghent_Analysis";
    nbins = 1;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 1, 0.5, 1.5);
    histo->SetTag(theSample, "Ghent", "Ghent", "");
  }
  //############ Variables ############
 else if (plot == "MET"){
    htitle = "MET";
    nbins = 17;
    float xbins[] = {0, 50, 80, 110, 140, 170, 200, 230, 260, 290, 320, 350, 400, 450, 500, 550, 600}; 
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
    histo->SetTag(theSample, "MET", "MET [GeV]", "e#mu, M_{ll} #geq 20 GeV, #geq 2jets, #geq 1btag");
  }
 else if (plot == "HT"){
    htitle = "HT";
    //nbins = 17;
    //float xbins[] = {0, 60, 100, , 140, 170, 200, 230, 260, 290, 320, 350, 400, 450, 500, 550, 600}; 
    //histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 23, 0, 1380);
    histo->SetTag(theSample, "HT", "HT [GeV]", "e#mu, M_{ll} #geq 20 GeV, #geq 2jets, #geq 1btag");
  }
	else if(plot == "MT2"){ 
    htitle = "MT2";
    nbins = 18;
    float xbins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 170, 200, 250, 300, 400}; 
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
  }
	else if(plot == "METSR"){ 
    htitle = "MET";
    nbins = 9;
    float xbins[] = {0, 20, 40, 50, 100, 140, 200, 300, 400}; 
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
  }
  else if(plot == "DeltaPhi"){
    htitle = "DeltaPhi";
    nbins = 20;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20, -3.2, 3.2);
  }
  else if(plot == "DeltaPhiLepMet"){
    htitle = "DeltaPhiLepMet";
    nbins = 20;
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20, -3.2, 3.2);
  }
  else if(plot == "DeltaEta"){
    htitle = "DeltaEta";
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 15, 0, 3.2);
  }
  else if(plot == "DeltaR"){
    htitle = "DeltaR";
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 15, 0, 8);
  }
  else               { std::cout << " -------> Wrong plot name!! " << endl; return;}
  TLorentzVector l0;
  TLorentzVector l1;
  TLorentzVector jet0;
  TLorentzVector jet1;
  TLorentzVector vmet;


	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		Float_t met = TMET; Float_t ht = THT; 
	Float_t weight = TWeight; 
    Int_t njets = TNJets; Int_t nbjets = TNJetsBtag; // + JetPt...

		if     (sys == "0" || sys == "nom") weight = TWeight;
    else if(sys == "LepEffUp")          weight = TWeight_LepEffUp;
    else if(sys == "LepÈffDown")         weight = TWeight_LepEffDown;
    else if(sys == "TrigUp")            weight = TWeight_TrigUp;
    else if(sys == "TrigDown")          weight = TWeight_TrigDown;
    else if(sys == "FSUp")              weight = TWeight_FSUp;
    else if(sys == "FSDown")            weight = TWeight_FSDown;
    else if(sys == "PUUp")              weight = TWeight_PUUp;
    else if(sys == "PUDown")            weight = TWeight_PUDown;

    else if(sys == "BtagUp")         nbjets = TNJetsBtagUp;
    else if(sys == "BtagDown")       nbjets = TNJetsBtagDown;
    else if(sys == "MisTagUp")       nbjets = TNJetsBtagMisTagUp;
    else if(sys == "MisTagDown")     nbjets = TNJetsBtagMisTagDown;

		else if(sys == "JESUp"){
			njets = TNJetsJESUp;
			met = TMETJESUp;
			ht = 0;
			for(int k = 0; k<njets; k++) ht += TJetJESUp_Pt[k];
		}
    else if(sys == "JESDown"){
			met = TMETJESDown;
      njets = TNJetsJESDown;
			ht = 0;
			for(int k = 0; k<njets; k++) ht += TJetJESUp_Pt[k];
    }
    else if(sys == "JER" || sys == "JERUp"){
      njets = TNJetsJER;
    }


    else{
      weight = TWeight;
      njets = TNJets;
      nbjets = TNJetsBtag;
    }

    if(njets < 2) continue;

    Float_t lb = TMT2lblb; Float_t bb = TMT2bb; Float_t ll = TMT2ll;
    l0.SetPxPyPzE(TLep1_Px, TLep1_Py, TLep1_Pz, TLep1_E);
    l1.SetPxPyPzE(TLep2_Px, TLep2_Py, TLep2_Pz, TLep2_E);
    jet0.SetPxPyPzE(TJet_Px[0], TJet_Py[0], TJet_Pz[0], TJet_E[0]);
    jet1.SetPxPyPzE(TJet_Px[1], TJet_Py[1], TJet_Pz[1], TJet_E[1]);

    vmet.SetPtEtaPhiM(met, 0., TMET_Phi, 0);
    Float_t dEta = TMath::Abs(l0.Eta() - l1.Eta());
    Float_t dPhi = TMath::Abs(l0.DeltaPhi(l1));

		if(chan == "ElMu" && !TIsElMu) continue;
		if(chan == "Elec" && !TIsDoubleElec) continue;
		if(chan == "Muon" && !TIsDoubleMuon) continue;
		if( (chan == "SF" || chan == "sameF" || chan == "SameF") && (!TIsDoubleMuon && !TIsDoubleElec)) continue;

		//############ SRs ############
		if     (plot == "SRMT23D"){
			if(ll < 100){
          if     (lb<100 && bb<170)               histo->Fill(1, weight);
          else if(lb<100 && bb>170 && bb<270)     histo->Fill(2, weight);
          else if(lb<100 && bb>270)               histo->Fill(3, weight);
          else if(lb<200 && lb>100 && bb<170)     histo->Fill(4, weight);
					else if(lb<200 && lb>100 && bb<270)     histo->Fill(5, weight);
					else if(lb<200 && lb>100 && bb>270)     histo->Fill(6, weight);
					else if(lb>200 && bb<170)               histo->Fill(7, weight);
					else if(lb>200 && bb>170 && bb<270)     histo->Fill(8, weight);
					else if(lb>200 && bb>270)               histo->Fill(9, weight);
			}

			else if(ll > 100 && ll < 200){ 
          if     (lb<100 && bb<170)               histo->Fill(10, weight);
          else if(lb<100 && bb>170 && bb<270)     histo->Fill(11, weight);
          else if(lb<100 && bb>270)               histo->Fill(12, weight);
          else if(lb<200 && lb>100 && bb<170)     histo->Fill(13, weight);
					else if(lb<200 && lb>100 && bb<270)     histo->Fill(14, weight);
					else if(lb<200 && lb>100 && bb>270)     histo->Fill(15, weight);
					else if(lb>200 && bb<170)               histo->Fill(16, weight);
					else if(lb>200 && bb>170 && bb<270)     histo->Fill(17, weight);
					else if(lb>200 && bb>270)               histo->Fill(18, weight);
			}
			else if(ll  > 200){ 
          if     (lb<100 && bb<170)               histo->Fill(19, weight);
          else if(lb<100 && bb>170 && bb<270)     histo->Fill(20, weight);
          else if(lb<100 && bb>270)               histo->Fill(21, weight);
          else if(lb<200 && lb>100 && bb<170)     histo->Fill(22, weight);
					else if(lb<200 && lb>100 && bb<270)     histo->Fill(23, weight);
					else if(lb<200 && lb>100 && bb>270)     histo->Fill(24, weight);
					else if(lb>200 && bb<170)               histo->Fill(25, weight);
					else if(lb>200 && bb>170 && bb<270)     histo->Fill(26, weight);
					else if(lb>200 && bb>270)               histo->Fill(27, weight);
			}
		}

		else if(plot == "SR"){ 
			if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) )){


				if     (dEta < 0.5 && dPhi < 1 && met< 100) histo->Fill(1,weight);
				if     (dEta < 0.5 && dPhi < 1 && met> 100 && met< 150) histo->Fill(2,weight); 
				if     (dEta < 0.5 && dPhi < 1 && met> 150 && met< 200) histo->Fill(3,weight); 
				if     (dEta < 0.5 && dPhi < 1 && met> 200 && met< 300) histo->Fill(4,weight); 
				if     (dEta < 0.5 && dPhi < 1 && met> 300) histo->Fill(5,weight); 

				else if(dEta < 0.5 && dPhi > 1 && dPhi < 2.5 && met< 100) histo->Fill(6,weight); 
				else if(dEta < 0.5 && dPhi > 1 && dPhi < 2.5 && met> 100 && met< 150) histo->Fill(7,weight); 
				else if(dEta < 0.5 && dPhi > 1 && dPhi < 2.5 && met> 150 && met< 200) histo->Fill(8,weight); 
				else if(dEta < 0.5 && dPhi > 1 && dPhi < 2.5 && met> 200 && met< 300) histo->Fill(9,weight); 
				else if(dEta < 0.5 && dPhi > 1 && dPhi < 2.5 && met> 300) histo->Fill(10,weight); 

				else if(dEta < 0.5 && dPhi > 2.5 && met< 100) histo->Fill(11,weight); 
				else if(dEta < 0.5 && dPhi > 2.5 && met> 100 && met< 150) histo->Fill(12,weight); 
				else if(dEta < 0.5 && dPhi > 2.5 && met> 150 && met< 200) histo->Fill(13,weight); 
				else if(dEta < 0.5 && dPhi > 2.5 && met> 200 && met< 300) histo->Fill(14,weight); 
				else if(dEta < 0.5 && dPhi > 2.5 && met> 300) histo->Fill(15,weight); 

				else if(dEta > 0.5 && dEta<1 && dPhi < 1 && met< 100) histo->Fill(16,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi < 1 && met> 100 && met< 150) histo->Fill(17,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi < 1 && met> 150 && met< 200) histo->Fill(18,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi < 1 && met> 200 && met< 300) histo->Fill(19,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi < 1 && met> 300) histo->Fill(20,weight); 

				else if(dEta > 0.5 && dEta<1 && dPhi > 1 && dPhi < 2.5 && met< 100) histo->Fill(21,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 1 && dPhi < 2.5 && met> 100 && met< 150) histo->Fill(22,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 1 && dPhi < 2.5 && met> 150 && met< 200) histo->Fill(23,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 1 && dPhi < 2.5 && met> 200 && met< 300) histo->Fill(24,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 1 && dPhi < 2.5 && met> 300) histo->Fill(25,weight); 

				else if(dEta > 0.5 && dEta<1 && dPhi > 2.5 && met< 100) histo->Fill(26,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 2.5 && met> 100 && met< 150) histo->Fill(27,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 2.5 && met> 150 && met< 200) histo->Fill(28,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 2.5 && met> 200 && met< 300) histo->Fill(29,weight); 
				else if(dEta > 0.5 && dEta<1 && dPhi > 2.5 && met> 300) histo->Fill(30,weight); 

        // MET > 100, MET < 100
				else if(dEta > 1 && dEta<2 && dPhi < 1 && met< 100)  histo->Fill(31,weight); 
				else if(dEta > 1 && dEta<2 && dPhi < 1 && met> 100)  histo->Fill(32,weight); 
				else if(dEta > 1 && dEta<2 && dPhi > 1 && dPhi < 2.5 && met< 100)  histo->Fill(33,weight); 
				else if(dEta > 1 && dEta<2 && dPhi > 1 && dPhi > 2.5 && met> 100)  histo->Fill(34,weight); 
				else if(dEta > 1 && dEta<2 && dPhi > 2.5 && met< 100) histo->Fill(35,weight); 
				else if(dEta > 1 && dEta<2 && dPhi > 2.5 && met> 100) histo->Fill(36,weight); 
				else if(dEta > 2) histo->Fill(37,weight); 
			}
    }
		else if(plot == "SR12"){ 
			if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) )){


				if     (dEta < 0.5 && dPhi < 2 && met< 150)               histo->Fill(1,weight); 
        else if(dEta < 0.5 && dPhi < 2 && met> 150 && met< 300) histo->Fill(2,weight); 
				else if(dEta < 0.5 && dPhi < 2 && met> 300)               histo->Fill(3,weight); 
				else if(dEta < 0.5 && dPhi > 2 && met< 150)               histo->Fill(4,weight); 
        else if(dEta < 0.5 && dPhi > 2 && met> 150 && met< 300) histo->Fill(5,weight); 
				else if(dEta < 0.5 && dPhi > 2 && met> 300)               histo->Fill(6,weight); 
				else if(dEta > 0.5 && dPhi < 2 && met< 150)               histo->Fill(7,weight); 
        else if(dEta > 0.5 && dPhi < 2 && met> 150 && met< 300) histo->Fill(8,weight); 
				else if(dEta > 0.5 && dPhi < 2 && met> 300)               histo->Fill(9,weight); 
				else if(dEta > 0.5 && dPhi > 2 && met< 150)               histo->Fill(10,weight); 
        else if(dEta > 0.5 && dPhi > 2 && met> 150 && met< 300) histo->Fill(11,weight); 
				else if(dEta > 0.5 && dPhi > 2 && met> 300)               histo->Fill(12,weight); 
			}
		}
		else if(plot == "cutandcount"){
			if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
				histo->Fill(1, weight);
		}
		else if(plot == "ghent"){
			if(met> 80 && nbjets> 0 && (met/TMath::Sqrt(ht)) > 5 ) 
				histo->Fill(1, weight);
		}


		//############ Variables ############
		else if(plot == "MET"){
			if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
				histo->Fill(met, weight);
		}
		else if(plot == "MT2"){
      if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
				histo->Fill(TMT2ll, weight);
		}
    else if(plot == "DeltaPhi"){
      if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
        histo->Fill(l0.DeltaPhi(l1), weight);
    }
    else if(plot == "DeltaPhiLepMet"){
      if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
        histo->Fill(l0.DeltaPhi(vmet), weight);
    }
    else if(plot == "DeltaEta"){
      if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
        histo->Fill(TMath::Abs(l0.Eta() - l1.Eta()), weight);
    }
    else if(plot == "DeltaR"){
      if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
        histo->Fill(l0.DeltaR(l1), weight);
    }
    else if(plot == "HT"){
      if(met> 50 && nbjets> 0 && ((chan == "ElMu") || ( (met/TMath::Sqrt(ht)) > 5 && TMinDPhiMetJets > 0.25 ) ))
        histo->Fill(ht, weight);
    }
	}
  histo->SetStyle();
}
