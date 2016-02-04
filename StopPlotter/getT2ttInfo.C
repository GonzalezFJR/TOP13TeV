#include <vector>

const TString inputpath = "/pool/ciencias/TreesDR74X/heppyTrees/v2/";
const TString outputdir = "/nfs/fanae/user/juanr/StopTOP/plots/";
const Int_t nSamples = 9;
const TString samplename[nSamples] = {
"T2tt_150to175LSP1to100_0", "T2tt_250LSP1to175_0", "T2tt_275LSP75to200_0",
"T2tt_200LSP1to125_0",  "T2tt_300to375LSP1to300_0",    "T2tt_400to475LSP1to400_0",
"T2tt_225LSP25to150_0",     "T2tt_500to550LSP1to475_0", "T2tt_600to950LSP1to450_0"};

std::vector<std::vector<float>> Points;
std::vector<float> NeutralinoMass;
std::vector<float> StopMass;
std::vector<int> Events;

TH3* Counts;

float mstop; float mchi; int nEvents;
void GetValues(TString sample){
	TFile* inputfile = TFile::Open(inputpath + "Tree_" + sample + ".root");
	inputfile->GetObject("CountSMS", Counts);
	float val = 0;
	for(int i = 0; i < Counts->GetNbinsX(); i++){
		for(int j = 0; j< Counts->GetNbinsY(); j++){
			val = Counts->GetBinContent(i,j,1);
			if(val!=0){
				mstop = Counts->GetXaxis()->GetBinCenter(i);
				mchi  = Counts->GetYaxis()->GetBinCenter(j); 
        nEvents = Counts->GetBinContent(Counts->FindBin(mstop, mchi, 0));
				if(mchi == 0) mchi = 1;
				NeutralinoMass.push_back(mchi);
				StopMass.push_back(mstop);
        Events.push_back(nEvents); 
			}
		}
	} 
	inputfile->Close(); delete inputfile; 
}

void GetAllMasses(){
  for(int k = 0; k < nSamples; k++) GetValues( samplename[k] );
}

void PlotPoints(){
	GetAllMasses();
  TCanvas *c1 = new TCanvas("c1","MassPlot",200,10,700,500);
  c1 -> DrawFrame(0, 0, 850, 600);
	TGraph *MassPlot = new TGraph(StopMass.size(), (float*) &StopMass[0], (float*) &NeutralinoMass[0]);
  MassPlot->SetMarkerStyle(20);
  MassPlot->SetMarkerSize(1.1);
  MassPlot->SetMarkerColor(2);
  MassPlot->Draw("P");
  TLatex Title = TLatex(); Title.DrawLatexNDC(.40, .94, "Mass Points");
  TLatex Xaxis = TLatex(); Xaxis.DrawLatexNDC(0.65, 0.03, "Stop Mass [GeV]");
  TLatex Yaxis = TLatex(); Yaxis.SetTextAngle(90); Yaxis.DrawLatexNDC(0.03, 0.31, "Neutralino Mass [GeV]");
  c1 -> Print(outputdir + "MassPoints.pdf", "pdf");
  c1 -> Print(outputdir + "MassPoints.png", "png");
}

void Plotcolz(){
  if(StopMass.size() == 0) GetAllMasses();
	TCanvas *c1 = new TCanvas("c1","c1",600,400);
	TH2F* hist = new TH2F("hist", "nEvents for the different mass points", 40, 12.5, 1000+12.5, 24, -12.5, 600-12.5);
	for(int k = 0; k < StopMass.size(); k++){
		hist->SetBinContent(hist->FindBin(StopMass[k], NeutralinoMass[k]), Events[k]);
	}
	hist -> SetStats(0);
  hist->GetXaxis()->SetTitle("Stop Mass [GeV]");
  hist->GetYaxis()->SetTitle("Neutralino Mass [GeV]");
  hist->GetXaxis()->SetNdivisions(40308);
  hist->GetYaxis()->SetNdivisions(40306);
  c1->SetGrid();
	hist->Draw("colz");
  c1 -> Print(outputdir + "MassPointsEvents.pdf", "pdf");
  c1 -> Print(outputdir + "MassPointsEvents.png", "png");

}


void PrintEvents(){
  if(StopMass.size() == 0) GetAllMasses();
	float c = 0;
	cout << "=====================================================" << endl;
	cout << " Stop Mass     ||  Neutralino Mass  ||    nEvents   " << endl;
	cout << "-----------------------------------------------------" << endl;
	for(int k = 0; k < StopMass.size(); k++){
		cout << Form("    %3.2f      ||       %3.2f       ||   %i      ", StopMass[k], NeutralinoMass[k], Events[k]) << endl;
	}
	cout << "=====================================================" << endl;
}

