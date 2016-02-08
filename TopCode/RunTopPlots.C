//#include "TopPlotter.h"

//void RunTopPlots(TString pathtofile="/nfs/fanae/user/palencia/april14Run2/TOP/TopTrees/temp/", Int_t verbose=0){
//void RunTopPlots(TString pathtofile="/nfs/fanae/user/palencia/april14Run2/TOP/TopTrees/may15/", Int_t verbose=2){
void RunTopPlots(TString pathtofile="../temp/", Int_t verbose=2){
  TString outputdir = pathtofile + "/TopPlots/";
  
  cout << "--------------" << endl;
  cout << "OutputDir is:      " << outputdir << endl;
  cout << "Verbose level is:  " << verbose << endl;
  cout << "--------------" << endl;
  
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
  
  cout << "Loading TopPlotter.cc..."<< endl;
  gROOT->LoadMacro("TopPlotter.cc+");
  
  cout << "Creating TopPlotter class..."<<endl;
  TopPlotter *tA = new TopPlotter();

  cout << "Init TopPlotter..."<<endl;
  tA->Init(pathtofile);

  cout << "Set OutputDir and Vervose level..."<<endl;
  tA->SetOutputDir(outputdir);
  tA->SetVerbose(verbose);
  
  tA->Loop();
  
  delete tA;
}
