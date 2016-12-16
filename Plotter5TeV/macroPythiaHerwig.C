{
//const TString path = "/mnt_pool/fanae105/user/juanr/TOP13TeV/temp/";
const TString path = "/mnt_pool/fanae105/user/juanr/TOP13TeV/Trees5TeV/may23/";

const TString p8   = "TTbar_PowhegFidu";
const TString hpp  = "TTbar_Powheg_HerwigFidu";

const TString ngenLeptons = "fHnGenLeptons";
const TString nmu         = "fHnGenMuo";
const TString nel         = "fHnGenEle";

	TFile* f1 = TFile::Open(path + "Tree_" + p8  + ".root");
	TFile* f2 = TFile::Open(path + "Tree_" + hpp + ".root");
	TH1F* h1; TH1F* h2;

	TH1F* hp81; TH1F* hp82;
	TH1F* hhp1; TH1F* hhp2;

	f1->GetObject(ngenLeptons, h1);
	f2->GetObject(ngenLeptons, h2);
	h1->SetStats(0); h2->SetStats(0);

	f1->GetObject(nmu, hp81);
	f1->GetObject(nel, hp82);
	hp81->Add(hp82);
	hp81->SetStats(0);


	f2->GetObject(nmu, hhp1);
	f2->GetObject(nel, hhp2);
	hhp1->Add(hhp2);
	hhp1->SetStats(0);


	TCanvas *c = new TCanvas("c", "c", 10, 10, 1400, 600);
	c->Divide(2,1);
	TPad *plot1 = (TPad*) c->GetPad(1);
  plot1->SetLogy();
	c->cd(1);
  hhp1->SetTitle("# Gen Leptons in the sample");
  hhp1->GetXaxis()->SetTitle("Number of gen leptons");
	hhp1->SetLineColor(4); hhp1->SetLineWidth(2); hhp1->SetMarkerColor(4); hhp1->SetMarkerStyle(8);
	hp81->SetLineColor(2); hp81->SetLineWidth(2); hp81->SetMarkerColor(2); hp81->SetMarkerStyle(4);
	hhp1->Draw("lp");
 
	hp81->Draw("lpsame");

  TLegend *le = new TLegend(0.55,0.78,0.85,0.85);
	le->SetTextSize(0.025);
  le->SetBorderSize(2);
  le->SetFillColor(10);
  le->AddEntry(hp81, Form("P8,  nEvents = %1.0f", hp82->GetEntries()), "l");
  le->AddEntry(hhp1, Form("H++, nEvents = %1.0f", hhp2->GetEntries()), "l");
  le->Draw();

	TPad *plot2 = (TPad*) c->GetPad(2);
	c->cd(2);
  plot2->SetLogy();
	h1->SetLineColor(2); h1->SetLineWidth(2); h1->SetMarkerColor(2); h1->SetMarkerStyle(4);
	h2->SetLineColor(4); h2->SetLineWidth(2); h2->SetMarkerColor(4); h2->SetMarkerStyle(8);
  h2->SetTitle("# Gen Leptons for the fiducial selection");
  h2->GetXaxis()->SetTitle("Number of gen leptons");
	h2->Draw("lp");
	h1->Draw("lpsame");
	h1->SetStats(0);

  TLegend *leg = new TLegend(0.55,0.78,0.85,0.85);
	leg->SetTextSize(0.025);
  leg->SetBorderSize(2);
  leg->SetFillColor(10);
  leg->AddEntry(h1, Form("P8,  nEvents = %1.0f", h1->GetEntries()), "l");
  leg->AddEntry(h2, Form("H++, nEvents = %1.0f", h2->GetEntries()), "l");
  leg->Draw();

  c->Print("nGenLeptons_p8_h++.png", "png");
}
