#include "AnalMiniTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void AnalMiniTree::Loop(TString plot){
	TString htitle; Int_t nbins; Float_t binmin = 0; Float_t binmax = 0;
	if (fChain == 0) return;
	Long64_t nentries = fChain->GetEntriesFast();
	if(plot == "cutandcount"){
		htitle = "CutAndCount";
		nbins = 1;
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 1, 0.5, 1.5);
		histo->SetTag(theSample, "CutAndCount", "CutAndCount", "");
	}
	else if (plot == "MET"){
		htitle = "MET";
		nbins = 17;
		float xbins[] = {0, 50, 80, 110, 140, 170, 200, 230, 260, 290, 320, 350, 400, 450, 500, 550, 600}; 
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins-1, xbins);
		histo->SetTag(theSample, "MET", "MET [GeV]", "e#mu, M_{ll} #geq 20 GeV, #geq 2jets, #geq 1btag");
	}
	else if (plot == "InvMass" || plot == "Mll"){
		htitle = "Mll";
		nbins = 40;
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, nbins, 20, 220);
		histo->SetTag(theSample, "Mll", "Mll [GeV]", "e#mu");
	}
	else if (plot == "HT"){
		htitle = "HT";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 23, 0, 1380);
		histo->SetTag(theSample, "HT", "HT [GeV]", "e#mu, M_{ll} #geq 20 GeV, #geq 2jets, #geq 1btag");
	}
	else if(plot == "MT2"){ 
		htitle = "MT2";
		nbins = 20;
		float xbins[] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 250, 300, 400}; 
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 30, 0, 300);
	}
	else if(plot == "DeltaPhi"){
		htitle = "DeltaPhi";
		nbins = 20;
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20, 0, 3.2);
	}
  else if(plot == "NBtagsNJets"){
    htitle = "NBtagNJets";
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 15 , -0.5, 14.5);
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
	else if(plot == "Lep0Pt"){
		htitle = "Lep0Pt";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20, 0, 200);
	}
	else if(plot == "Lep1Pt"){
		htitle = "Lep1Pt";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20, 0, 200);
	}
	else if(plot == "Jet0Pt"){
		htitle = "Jet0Pt";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 30, 0, 300);
	}
	else if(plot == "Jet1Pt"){
		htitle = "Jet1Pt";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 30, 0, 300);
	}
	else if(plot == "NJets"){
		htitle = "NJets";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 8, -0.5, 7.5);
	}
	else if(plot == "Lep0Eta"){
		htitle = "Lep0Eta";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20,-2.4, 2.4);
	}
	else if(plot == "Lep1Eta"){
		htitle = "Lep1Eta";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20,-2.4, 2.4);
	}
	else if(plot == "Jet0Eta"){
		htitle = "Jet0Eta";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 30,-2.4, 2.4);
	}
	else if(plot == "Jet1Eta"){
		htitle = "Jet1Eta";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 20, -2.4, 2.4);
	}
	else if(plot  =="NVert"){
		htitle = "NVert";
		histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 100, -0.5, 99.5);
	}
  else if(plot == "InvMassWW"){
    htitle = "InvMassWW";
    histo = new Histo(theSample + "_" + htitle, theSample + "_" + htitle, 30, 0., 150);
  }
	else               { std::cout << " -------> Wrong plot name!! " << endl; return;}


	TLorentzVector l0;   TLorentzVector l1;
	TLorentzVector jet0; TLorentzVector jet1;
	TLorentzVector vmet;

	Float_t weight;
	Float_t met; Float_t ht; Float_t mt2;
	Int_t njets; Int_t nbjets;
	Float_t jetpt0; Float_t jetpt1;

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Values for no variations
		met = TMET; ht = THT; mt2 = TMT2ll;
		weight = TWeight; 
		njets = TNJets; nbjets = TNJetsBtag; 
		jetpt0 = 0; jetpt1 = 0;

    //##### Weights...
    // Lepton efficiencies, trigger efficiencies, FullSim/FastSim, PU reweighting
		if     (sys == "0" || sys == "nom") weight = TWeight;
		else if(sys == "LepEffUp")          weight = TWeight_LepEffUp;
		else if(sys == "LepEffDown")       weight = TWeight_LepEffDown;
		else if(sys == "TrigUp")            weight = TWeight_TrigUp;
		else if(sys == "TrigDown")          weight = TWeight_TrigDown;
		else if(sys == "FSUp")              weight = TWeight_FSUp;
		else if(sys == "FSDown")            weight = TWeight_FSDown;
		else if(sys == "PUUp")              weight = TWeight_PUUp;
		else if(sys == "PUDown")            weight = TWeight_PUDown;

    // ######### B-tagging
    // Btag weights, Mistag
		else if(sys == "BtagUp")         nbjets = TNJetsBtagUp;
		else if(sys == "BtagDown")       nbjets = TNJetsBtagDown;
		else if(sys == "MistagUp")       nbjets = TNJetsBtagMisTagUp;
		else if(sys == "MistagDown")     nbjets = TNJetsBtagMisTagDown;

    //##########  Jet energy scale and resolution
    // Propagate to: nJets, JetPt, MET, MT2, HT
		else if(sys == "JESUp"){
			njets = TNJetsJESUp;
			met = TMETJESUp;
			mt2 = TMT2llJESUp;
			jetpt0 = TJetJESUp_Pt[0];
			jetpt1 = TJetJESUp_Pt[1];
			ht = THTJESUp;
		}
		else if(sys == "JESDown"){
			met = TMETJESDown;
			mt2 = TMT2llJESDown;
			njets = TNJetsJESDown;
			jetpt0 = TJetJESDown_Pt[0];
			jetpt1 = TJetJESDown_Pt[1];
			ht = THTJESDown;
		}
		else if(sys == "JER" || sys == "JERUp"){
			njets = TNJetsJER;
		}

		// ######## Special FullSim/FastSim MET systematic uncertainty


		else if(sys == "FSMETUp"  ){ if(theSample.BeginsWith("T2tt"))  met = TMET + TMath::Abs(TGenMET - TMET)/2;}
		else if(sys == "FSMETDown"){ if(theSample.BeginsWith("T2tt"))  met = TMET - TMath::Abs(TGenMET - TMET)/2;}
		else{    std::cout << " [AnalMiniTree::Error] No variation found: " << sys << endl;}

		// Set leptons, jets, met
		l0.SetPxPyPzE(TLep1_Px, TLep1_Py, TLep1_Pz, TLep1_E);
		l1.SetPxPyPzE(TLep2_Px, TLep2_Py, TLep2_Pz, TLep2_E);
		jet0.SetPxPyPzE(TJet_Px[0], TJet_Py[0], TJet_Pz[0], TJet_E[0]); //jet0.SetPtEtaPhiM(jetpt0, jet0.Eta(), jet0.Phi(), jet0.M());
		jet1.SetPxPyPzE(TJet_Px[1], TJet_Py[1], TJet_Pz[1], TJet_E[1]); //jet1.SetPtEtaPhiM(jetpt1, jet1.Eta(), jet1.Phi(), jet1.M());
		vmet.SetPtEtaPhiM(met, 0., TMET_Phi, 0);

		Float_t dEta = TMath::Abs(l0.Eta() - l1.Eta());
		Float_t dPhi = TMath::Abs(l0.DeltaPhi(l1));

    // Just to make sure, check charge

    //########### Select channel ################
		if(chan == "ElMu" && (TChannel != 1)) continue;
		if(chan == "Elec" && (TChannel != 3)) continue;
		if(chan == "Muon" && (TChannel != 2)) continue;
		if( (chan == "SF" || chan == "sameF" || chan == "SameF") && (TChannel == 1) ) continue;

		//############ ############
		else if(plot == "InvMassWW" && (met > 20 && njets == 0 && (l0+l1).Pt() > 30)) histo->Fill((l0+l1).M());
		else if(plot == "NBtagsNJets"){
      if     (njets == 0)  histo->Fill(nbjets     , weight);
      else if(njets == 1)  histo->Fill(nbjets + 1 , weight);
      else if(njets == 2)  histo->Fill(nbjets + 3 , weight);
      else if(njets == 3)  histo->Fill(nbjets + 6 , weight);
      else if(njets >= 4)  histo->Fill(nbjets + 10, weight);
    }

		//############ Stop selection ############
		if(njets >= 2 && nbjets >= 1 && met > 50){
			     if(plot == "cutandcount")                histo->Fill(1, weight);
			else if(plot == "MET")   		     	            histo->Fill(met, weight);
			else if(plot == "MT2")                        histo->Fill(mt2, weight);
			else if(plot == "DeltaPhi")                   histo->Fill(l0.DeltaPhi(l1), weight);
			else if(plot == "DeltaPhiLepMet") 	          histo->Fill(l0.DeltaPhi(vmet), weight);
			else if(plot == "DeltaEta") 		          		histo->Fill(TMath::Abs(l0.Eta() - l1.Eta()), weight);
			else if(plot == "DeltaR")    			           	histo->Fill(l0.DeltaR(l1), weight);
			else if(plot == "HT")     				           	histo->Fill(ht, weight);
			else if(plot == "InvMass" || plot == "Mll") 	histo->Fill(TMll, weight);
			else if(plot == "Lep0Pt")               			histo->Fill(l0.Pt(), weight);
			else if(plot == "Lep1Pt")			               	histo->Fill(l1.Pt(), weight);
			else if(plot == "Jet0Pt")			               	histo->Fill(jet0.Pt(), weight);
			else if(plot == "Jet1Pt")			              	histo->Fill(jet1.Pt(), weight);
			else if(plot == "NJets")		              		histo->Fill(njets, weight);
			else if(plot == "Lep0Eta")		            		histo->Fill(l0.Eta(), weight);
			else if(plot == "Lep1Eta")		            		histo->Fill(l1.Eta(), weight);
			else if(plot == "Jet0Eta")		            		histo->Fill(jet0.Eta(), weight);
			else if(plot == "Jet1Eta")			            	histo->Fill(jet1.Eta(), weight);
			else if(plot == "NVert")	  	            		histo->Fill(TNVert, weight);
		}
	}
	histo->SetStyle();
}
