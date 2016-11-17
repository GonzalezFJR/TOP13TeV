const TString path = "/nfs/fanae/user/juanr/stop80/genAnal/genVars.root";
const TString process[]     = {"ttbar", "stop_200_25", "stop_250_75"};
const int nVars = 19;
const TString vars[nVars] = {"NJets","NBJets","Lep0Pt","Lep1Pt","Jet0Pt","Jet1Pt","HT","MET","Meff", "Mll", "DilepPt", "DeltaPhill", "DeltaPhijj", "DeltaPhilj", "DeltaPhilm", "DeltaPhijm","MT2", "MT2lblb", "MT2bb"};



TH1F* httbar; TH1F* httbarMLM; TH1F* hstop1; TH1F* hstop2;
TFile* file = TFile::Open(path);

void plotvar(TString var){
  file->GetObject( var + "_" + "ttbar_powheg",  httbar);
  file->GetObject( var + "_" + "ttbar_MLM",     httbarMLM);
  file->GetObject( var + "_" + "stop_200_25", hstop1);
  file->GetObject( var + "_" + "stop_250_75", hstop2);

  httbar->SetLineColor(kRed);
  httbar->SetMarkerColor(kRed);
  httbarMLM->SetLineColor(kPink);
  httbarMLM->SetMarkerColor(kPink);
  hstop1->SetLineColor(kBlue);
  hstop1->SetMarkerColor(kBlue);
  hstop2->SetLineColor(kGreen);
  hstop2->SetMarkerColor(kGreen);


  httbar->SetLineWidth(2); httbarMLM->SetLineWidth(2); hstop1->SetLineWidth(2); hstop2->SetLineWidth(2);
  httbar->SetStats(0); httbarMLM->SetStats(0); hstop1->SetStats(0); hstop2->SetStats(0);
  httbar->Scale(1/httbar->Integral()); httbarMLM->Scale(1/httbar->Integral()); hstop1->Scale(1/hstop1->Integral()); hstop2->Scale(1/hstop2->Integral());
  httbar->SetMarkerStyle(21); hstop1->SetMarkerStyle(21); hstop2->SetMarkerStyle(21);
  httbar->SetMarkerSize(2.0); hstop1->SetMarkerSize(2.0); hstop2->SetMarkerSize(2.0);

  httbar->SetTitle(var);
  httbar->GetXaxis()->SetTitle(var);

  TLegend* leg;
  leg = new TLegend(0.80,0.79,0.93,0.93);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(httbar, "ttbar powheg", "lp");
  leg->AddEntry(httbarMLM, "ttbar MLM", "lp");
  leg->AddEntry(hstop1, "stop_200_25", "lp");
  leg->AddEntry(hstop2, "stop_250_75", "lp");

  TCanvas* c= new TCanvas("c","c",10,10,800,600);
  
  httbar->SetMaximum(httbar->GetMaximum()*1.4);
  httbar->Draw("");
  httbarMLM->Draw("same");
  hstop1->Draw("same");
  hstop2->Draw("same");
  leg->Draw("same");

  c->Print(var + ".png", "png");
  //file->Close();
  delete c; 
}


void genPlotter(){
	for(int i = 0; i < nVars; i++){
		plotvar(vars[i]);
	}
}
