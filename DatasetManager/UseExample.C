///////////////////////////////////////////////////////////////////////
//
//    FILE: UseExample.C
//  AUTHOR: I. Gonzalez Caballero
//    DATE: January, 2011
//
// CONTENT: Illustrative macro. Example on how to use the DatasetManager
//          utility class.
///////////////////////////////////////////////////////////////////////

// Load DatasetManager in ROOT 6
R__LOAD_LIBRARY(DatasetManager.C+)
#include "DatasetManager.h"


void UseExample() {
  //Create a new DatasetManager for U.O.
  cout << ">> Retrieving DatasetManager..." << endl;
  DatasetManager* dm = DatasetManager::GetInstance();
  // Set the tab in the google doc to be used
  dm->SetTab("DR74X25nsMiniAODv2");

  //Use this if you know that the information on the google doc table has
  //changed and you need to update the information... otherwise the class
  //itself will redownload the information automatically if there has been
  //more than 24 hours since the last synchronization.
  //
  //dm->RedownloadFiles();

  //Select your dataset and load its information
  dm->LoadDataset("WZTo3LNu");


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
    cerr << "ERROR: Could not find files!" << endl;
    return;
  }

  //   + Dump all information to stdout
  cout << ">> Dumping..." << endl;
  dm->Dump();

  //Check if the static method for real data works
  cout << ">> Now finding real data..." << endl;
  vector<TString> realdata= DatasetManager::GetRealDataFiles("DoubleMuon_Run2015D_v4");
  for (unsigned int i = 0; i < realdata.size(); i++)
    cout << "   [" << i << "] " << realdata[i] << endl; 

}
