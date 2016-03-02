{//const TString pathv1 = "/nfs/fanae/user/palencia/testHeppy/TOP/TopTrees/temp/Tree_13TeV_EA_TestHeppy_v1.root";
//const TString path24 = "/nfs/fanae/user/palencia/testHeppy/TOP/TopTrees/feb24/";
//const TString pathv2 = "/nfs/fanae/user/palencia/testHeppy/TOP/TopTrees/temp/Tree_13TeV_EA_TestHeppy_v2.root";
//const TString pathv3 = "/nfs/fanae/user/palencia/testHeppy/TOP/TopTrees/temp/Tree_13TeV_EA_TestHeppy_v3.root";
const TString filename = "Tree_13TeV_EA_TTbar_Powheg.root";
const TString heppytree = "Tree_TTbar_Powheg_0.root";
//const TString pathv1 = "/pool/ciencias/users/user/palencia/tree_v1.root";
//const TString pathv2 = "/pool/ciencias/users/user/palencia/tree_v2.root";
const TString pathv1 = "/pool/ciencias/heppyTreesDR76X/v1/Tree_TTbar_Powheg_0.root";
const TString pathv2 = "/pool/ciencias/heppyTreesDR76X/v2/Tree_TTbar_Powheg_0.root";

TChain * fChain1 = new TChain("tree","");
TChain * fChain2 = new TChain("tree","");


	fChain1->AddFile(pathv1 );
  fChain2->AddFile(pathv2 );

	TCanvas *c1 = new TCanvas();

	fChain1->Draw("Jet_jecUp_pt * Jet_corr_JER", "Jet_pt < 200 && Jet_corr_JER > 0");
  TH1F h1 = (*htemp);
  h1.SetDirectory(0);
  delete htemp;
  h1.SetLineColor(kRed+1);
	fChain1->Draw("Jet_jecDown_pt * Jet_corr_JER", "Jet_pt < 200 && Jet_corr_JER > 0");
  TH1F h11 = (*htemp);
  h11.SetDirectory(0);
  delete htemp;
  h11.SetLineColor(kBlue);

  fChain2->Draw("Jet_jecUp_pt * Jet_corr_JER", "Jet_pt < 200 && Jet_corr_JER > 0");
  TH1F h2 = *htemp;
  h2.SetDirectory(0);
  delete htemp;
  h2.SetLineColor(kGreen+1);
  fChain2->Draw("Jet_jecDown_pt * Jet_corr_JER", "Jet_pt < 200 && Jet_corr_JER > 0");
  TH1F h22 = *htemp;
  h22.SetDirectory(0);
  delete htemp;
  h22.SetLineColor(1);

	//fChain1->Draw("Jet_corr_JER", "Jet_pt < 200 && Jet_corr_JER > 0.9 && Jet_corr_JER < 1.1");
	//fChain1->Draw("Jet_jecDown_pt", "Jet_pt < 200", "sames");
  c1->Update();
  c1->Clear();

  h1.Draw();
  h1.SetStats(0);
  h2.Draw("same");
  h2.SetStats(0);
  h11.Draw("same");
  h11.SetStats(0);
  h22.Draw("same");
  h22.SetStats(0);
  TLegend* leg = new TLegend(0.70,0.80,0.80,0.89);
  leg->SetTextSize(0.035);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->AddEntry(&h1,  "v1 Up", "l");
  leg->AddEntry(&h11, "v1 Down", "l");
  leg->AddEntry(&h2,  "v2 Up", "l");
  leg->AddEntry(&h22, "v2 Down", "l");
  leg->Draw("same");


  c1->Print("jetptvarcorrv1v2.pdf", "pdf");
  c1->Print("jetptvarcorrv1v2.png", "png");

//	fChain2->Draw("Jet_pt", "Jet_pt < 500", "sames");
}
