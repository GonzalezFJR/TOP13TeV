#include "DrawPlots.C"

TH1F* loadDYHisto(TString lab = "MC", TString ch = "ElMu", TString lev = "2jets");

double GetDYDD(  TString ch = "ElMu", TString lev = "2jets", bool IsErr = false, bool docout = false);
double GetNonWDD(TString ch = "ElMu", TString lev = "2jets", bool IsErr = false, bool docout = false);

TH1F* loadDYHisto(TString lab, TString ch, TString lev){
	TFile* f1; TFile* f2;

	TH1F*  h1; TH1F*  h2;
	if(lab == "MC" || lab == "DY"){
		f1 = TFile::Open(path + "Tree_DYJetsToLL_M50_aMCatNLO.root");
		f2 = TFile::Open(path + "Tree_DYJetsToLL_M10to50_aMCatNLO.root"); 
	}
	else{
		f1 = TFile::Open(path + "Tree_Data_SingleElec.root");
		f2 = TFile::Open(path + "Tree_Data_SingleMu.root"); 
	}
	f1->GetObject("H_DY_InvMass_" + ch + "_" + lev, h1);
	f2->GetObject("H_DY_InvMass_" + ch + "_" + lev, h2);
	h1->Add(h2);
  if(lab == "MC" || lab == "DY") h1->Scale(Lumi);
  h1->SetDirectory(0);
  f1->Close(); f2->Close(); delete f1; delete f2;
  return h1;
}

double GetDYDD(TString ch, TString lev, bool IsErr, bool docout){
  double Re     = 0; double N_ine  = 0; double N_oute = 0; double k_lle = 0; double SFel = 0;
	double R_erre = 0; double N_in_erre = 0; double N_out_erre = 0; double k_ll_erre = 0; double SFel_err = 0;

  double Rm     = 0; double N_inm  = 0; double N_outm = 0; double k_llm = 0; double SFmu = 0;
	double R_errm = 0; double N_in_errm = 0; double N_out_errm = 0; double k_ll_errm = 0; double SFmu_err = 0;

  double N_inem = 0; double N_in_errem = 0; double SFem = 0; double SFem_err = 0;

	Double_t nout_ee(0.),nin_ee(0.),nout_err_ee(0.),nin_err_ee(0.);
	Double_t nout_mm(0.),nin_mm(0.),nout_err_mm(0.),nin_err_mm(0.);
  
  TH1F* DYmm   = loadDYHisto("DY","Muon", lev);
  TH1F* DYee   = loadDYHisto("DY","Elec", lev);
  TH1F* Datamm = loadDYHisto("Da","Muon", lev);
  TH1F* Dataee = loadDYHisto("Da","Elec", lev);
  TH1F* Dataem = loadDYHisto("Da","ElMu", lev);

  Int_t low_in = DYmm->FindBin(76. );
  Int_t  up_in = DYmm->FindBin(106.);

  nin_mm = DYmm->Integral(low_in, up_in); nin_err_mm = TMath::Sqrt(nin_mm);
  nin_ee = DYee->Integral(low_in, up_in); nin_err_ee = TMath::Sqrt(nin_ee);
  nout_mm = DYmm->Integral(0,200)-nin_mm; nout_err_mm = TMath::Sqrt(nout_err_mm);
  nout_ee = DYee->Integral(0,200)-nin_ee; nout_err_ee = TMath::Sqrt(nout_err_ee);

  Rm     = nout_mm/nin_mm;
  R_errm = (nout_err_mm/nout_mm + nin_err_mm/nin_mm)*Rm;
  Re     = nout_ee/nin_ee;
  R_erre = (nout_err_ee/nout_ee + nin_err_ee/nin_ee)*Re;

  N_inm  = Datamm->IntegralAndError(low_in, up_in, N_in_errm);
  N_ine  = Dataee->IntegralAndError(low_in, up_in, N_in_erre);
  N_inem = Dataem->IntegralAndError(low_in, up_in, N_in_errem);

  k_llm = TMath::Sqrt(nin_mm / nin_ee  );
  k_lle = TMath::Sqrt(nin_ee / nin_mm  );

	k_ll_errm = 0.5*(nin_err_mm/nin_mm + nin_err_ee/nin_ee);
	k_ll_erre = 0.5*(nin_err_ee/nin_ee + nin_err_mm/nin_mm);

	N_outm = Rm * (N_inm - 0.5 * N_inem * k_llm);
	N_oute = Re * (N_ine - 0.5 * N_inem * k_lle);

  double tmperr1; double tmperr2;
	tmperr1 =  N_in_erre/N_ine;
	tmperr2 =  N_in_errem/N_inem + k_ll_erre/k_lle;
	N_out_erre = R_erre/Re*(tmperr1 + 0.5*tmperr2) * N_oute;

	tmperr1 =  N_in_errm/N_inm;
	tmperr2 =  N_in_errem/N_inem + k_ll_errm/k_llm;
	N_out_errm = R_errm/Rm*(tmperr1 + 0.5*tmperr2) * N_outm;

  SFmu   = N_outm / nout_mm;
  SFel   = N_oute / nout_ee;
  SFem   = TMath::Sqrt(SFel * SFmu);

  SFmu_err = (N_out_errm/N_outm + nout_err_mm/nout_mm) * SFmu;
  SFel_err = (N_out_erre/N_oute + nout_err_ee/nout_ee) * SFel;
  SFem_err = 0.5*SFem*(SFel_err/SFel + SFmu_err/SFmu);

  if(docout){
   cout << " DY DD estimation for " << lev << endl;
   cout << "=============================================================================" << endl;
   cout << "                  |       El/El      |       Mu/Mu      |       El/Mu      ||" << endl;
	 cout << "-----------------------------------------------------------------------------" << endl;
	 cout<< Form(" n_in (MC)        | %7.1f +/- %4.1f | %7.1f +/- %4.1f |                  ||", nin_ee, N_in_erre, nin_mm, N_in_errm) << endl;
	 cout<< Form(" n_out (MC)       | %7.1f +/- %4.1f | %7.1f +/- %4.1f |                  ||", nout_ee, N_out_erre, nout_mm, N_out_errm) << endl;
	 cout<< Form(" R(Nout/Nin)(MC)  | %6.3f +/- %4.3f | %6.3f +/- %4.3f |                  ||", Re, R_erre, Rm, R_errm) << endl;
	 cout<< Form(" k_ll             | %6.3f +/- %4.3f | %6.3f +/- %4.3f |                  ||", k_lle, k_ll_erre, k_llm, k_ll_errm) << endl;
	 cout<< Form(" N_in (D)         | %7.1f +/- %4.1f | %7.1f +/- %4.1f | %5.1f +/- %4.1f   ||",	N_ine, N_in_erre, N_inm, N_in_errm, N_inem, N_in_errem) << endl;
	 cout<< "-----------------------------------------------------------------------------" <<  endl;
	 cout<< Form(" N_out            | %7.1f +/- %4.1f | %7.1f +/- %4.1f |                  ||",	 N_oute, N_out_erre, N_outm, N_out_errm) << endl;
	 cout<< "-----------------------------------------------------------------------------" <<  endl;
	 cout<< Form(" SF (D/MC)        | %6.3f +/- %4.3f | %6.3f +/- %4.3f | %6.3f +/- %4.3f ||", SFel, SFel_err, SFmu, SFmu_err, SFem, SFem_err) << endl;
	 cout<< "-----------------------------------------------------------------------------" <<  endl;
	 cout<< Form(" Drell-Yan (MC)   | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||", yield("DY", "Elec", lev), getStatError("DY", "Elec", lev), yield("DY", "Muon", lev), getStatError("DY", "Muon", lev), yield("DY", "ElMu", lev), getStatError("DY", "ElMu", lev)) << endl;
	 cout<< Form(" Drell-Yan (D)    | %7.2f +/- %4.2f | %7.2f +/- %4.2f | %7.2f +/- %4.2f ||", yield("DY", "Elec", lev)*SFel, getStatError("DY", "Elec", lev)*SFel, yield("DY", "Muon", lev)*SFmu, getStatError("DY", "Muon", lev)*SFmu, yield("DY", "ElMu", lev)*SFem, getStatError("DY", "ElMu", lev)*SFem) << endl;
   cout << "=============================================================================" << endl;
	 cout<< endl;
  }

  if(IsErr){
    if     (ch == "ElMu") return SFem_err;
    else if(ch == "Elec") return SFel_err;
    else if(ch == "Muon") return SFmu_err;
  }
  else{
    if     (ch == "ElMu") return SFem;
    else if(ch == "Elec") return SFel;
    else if(ch == "Muon") return SFmu;
  }
  cout << " DY SF -----> Wrong input values!!" << endl;
  return 1.;
}

double GetNonWDD(TString ch, TString lev, bool IsErr, bool docout){
	double fake   = yield("fake", ch, lev, "0", 0);
	double fakeSS = yield("fake", ch, lev, "0", 1);
	double dataSS = yield("data", ch, lev, "0", 1);
	double bkgSS  = yield("bkg",  ch, lev, "0", 1) - fakeSS;  // prompt MC in SS
  
	double efake   = getStatError("fake", ch, lev, "0", 0);
	double efakeSS = getStatError("fake", ch, lev, "0", 1);
	double edataSS = getStatError("data", ch, lev, "0", 1);
	double ebkgSS  = getStatError("bkg",  ch, lev, "0", 1) - efakeSS;
  
	double  R = fake/fakeSS; 
  double eR = R*(efake/fake + efakeSS/fakeSS);

  double  DDfake = R*(dataSS - bkgSS);
  double eDDfake = DDfake*(eR/R  + (edataSS/dataSS - ebkgSS/bkgSS));

  double  SF = DDfake/fake;
  double eSF = SF*(eDDfake/DDfake + efake/fake);

	if(docout){
   cout <<      " NonW leptons DD estimation for " << lev << endl;
   cout <<      "==================================================================" << endl;
	 cout << Form(" MC fake estimation (OS WJets + ttbar semilep)   = %2.2f ± %0.2f", fake, efake) << endl;
	 cout <<      "------------------------------------------------------------------" << endl;
	 cout << Form(" MC fake SS (WJets + ttbar semilep)   = %2.2f ± %0.2f", fakeSS, efakeSS)  << endl;
   cout << Form(" R = fakeOS/fakeSS                    = %2.2f ± %0.2f", R, eR)            << endl;
	 cout << Form(" MC prompt SS (other sources)         = %2.2f ± %0.2f", bkgSS , ebkgSS )  << endl;
	 cout << Form(" Data SS events                       = %2.2f ± %0.2f", dataSS , edataSS )  << endl;
	 cout <<      "------------------------------------------------------------------" << endl;
   cout << Form(" DD fake estimation                   = %2.2f ± %0.2f", DDfake , eDDfake )  << endl;
   cout << Form(" Scake Factor (DD/MC)                 = %2.2f ± %0.2f", SF , eSF )  << endl;
   cout <<      "==================================================================" << endl;
	}

	if(IsErr) return eSF;
	else      return SF;
}

void GetDataDriven(){
  cout << " ###### Drell-Yan DD (R out/in) #####" << endl;
	cout << endl;
	GetDYDD("ElMu", "dilepton", 0, 1);
	cout << endl;
	GetDYDD("ElMu", "2jets", 0, 1);
	cout << endl;
	cout << endl;
  cout << " ###### NonW leptons estimation #####" << endl;
	cout << endl;
  GetNonWDD("ElMu", "dilepton", 0, 1);
	cout << endl;
  GetNonWDD("ElMu", "2jets", 0, 1);
	cout << endl;
}



