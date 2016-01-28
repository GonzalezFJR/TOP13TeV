///////////////////////////////////////////////////////////////////////
//
//    FILE: UseExample.C
//  AUTHOR: I. Gonzalez Caballero
//    DATE: January, 2011
//
// CONTENT: Illustrative macro. Example on how to use the DatasetManager
//          utility class.
///////////////////////////////////////////////////////////////////////

void UseExample() {
  //Load DatasetManager
  cout << ">> Loading Dataset Manager..." << endl;
  gROOT->LoadMacro("DatasetManager.C+");




  //Create a new DatasetManager for U.O.
  cout << ">> Creating DatasetManager..." << endl;
  
  DatasetManager* dm = new DatasetManager("Legacy_Summer12_53X");
  //If there is no skim you may use the following format
  //DatasetManager* dm = new DatasetManager("TescoSpring11");

  //Use this if you know that the information on the google doc table has
  //changed and you need to update the information... otherwise the class
  //itself will redownload the information automatically if there has been
  //more than 24 hours since the last synchronization.
  //
  //dm->RedownloadFiles();

  //Select your dataset and load its information
  dm->LoadDataset("T2tt_150to350LSP0to250");


  //Now print some information
  cout << ">> Let's print some information..." << endl;
  //   + cross section
  cout << "   + X Section = " << dm->GetCrossSection() << endl;
  //   + Events in the sample
  cout << "   + N Events  = " << dm->GetEventsInTheSample() << endl;
  //   + Local Folder
  cout << "   + Local folder = " << dm->GetLocalFolder() << endl;
  //   + List of files
  cout << ">> Getting files..." << endl;
  std::vector<TString> files = dm->GetFiles();
  if (files.size() == 0) {
    cerr << "ERROR: Could not retrieve files!" << endl;
    return;
  }

  //   + Dump all information to stdout
  cout << ">> Dumping..." << endl;
  dm->Dump();

  //Check if the static method for real data works
  cout << ">> Now finding real data..." << endl;
  vector<TString> realdata= dm->GetRealDataFiles("MC_Summer12_53X/Legacy", 
						 "Tree_DoubleMuParkedB_4412");
  for (unsigned int i = 0; i < realdata.size(); i++)
    cout << "   [" << i << "] " << realdata[i] << endl; 

}
