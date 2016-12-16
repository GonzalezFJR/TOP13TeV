{
TString path = "/mnt_pool/fanae105/user/juanr/TOP13TeV/temp/";
	TFile *f = TFile::Open(path + "Tree_" + "Data_SingleMu" + ".root");
	TH1F* h1; TH1F* h2; f->GetObject("fJetDataPt", h1); f->GetObject("fJetDataL2L3Pt", h2);
//	h1->SetDirectory(1); h2->SetDirectory(1); f->Close(); delete f;

	h1->SetLineColor(2); h1->SetLineWidth(2);
	h2->SetLineColor(4); h2->SetLineWidth(2);

	h1->Draw(); h2->Draw("same");

	TLegend *leg = new TLegend(0.60,0.75,0.93,0.93);
	leg->AddEntry(h1, "Jet_pt");
	leg->AddEntry(h2, "Jet_pt after L2L3 correction");
	leg->Draw("same");
}
