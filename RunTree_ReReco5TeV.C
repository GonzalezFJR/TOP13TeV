// Load DatasetManager in ROOT 6
R__LOAD_LIBRARY(DatasetManager/DatasetManager.C+)
#include "DatasetManager/DatasetManager.h"

/********************************************************
 * Main function
 ********************************************************/
void RunTree_ReReco5TeV(TString  sampleName     = "TTbar_Madgraph",
			Int_t    nSlots         = 1,
      Int_t    Selection      = 0,
			Bool_t   DoSystStudies  = false,
			Long64_t nEvents        = 0){

  // Selection:
  // 0: dilepton
  // 1: ngenLep < 2
  // 2: emu
  // 3: ee + mumu
  // 4: emu + fiducial
  
  // VARIABLES TO BE USED AS PARAMETERS...
  Float_t G_Total_Lumi    = 19664.225;   
  Float_t G_Event_Weight  = 1.0;         
  Bool_t  G_IsData        = false;       
  Float_t G_LumiForPUData = 19468.3;     // luminosity in http://www.hep.uniovi.es/jfernan/PUhistos
  Bool_t  G_IsMCatNLO     = false;

  TString sample = sampleName;

  // PAF mode
  //----------------------------------------------------------------------------
  cout << endl; 
  PAFIExecutionEnvironment* pafmode = 0;
  if (nSlots <=1 ) {
    PAF_INFO("RunTree_ReReco5TeV", "Sequential mode chosen");
    pafmode = new PAFSequentialEnvironment();
  }
  else if (nSlots <=8) {
    PAF_INFO("RunTree_ReReco5TeV", "PROOF Lite mode chosen");
    pafmode = new PAFPROOFLiteEnvironment(nSlots);
  }
  else {
    PAF_INFO("RunTree_ReReco5TeV", "PoD mode chosen");
    pafmode = new PAFPoDEnvironment(nSlots);
  }

  // Create PAF Project whith that environment
  //----------------------------------------------------------------------------
  PAFProject* myProject = new PAFProject(pafmode);

  // Base path to input files
  //----------------------------------------------------------------------------
  TString dataPath = "/pool/ciencias/";

  // INPUT DATA SAMPLE
  //----------------------------------------------------------------------------
  TString userhome = "/mnt_pool/fanae105/user/$USER/";
  DatasetManager* dm = DatasetManager::GetInstance();
  dm->SetTab("5TeVDR76X25nsAOD");
  //dm->RedownloadFiles();

  // Deal with data samples
  if (sampleName.Contains("Data") || sampleName.Contains("data")){
    cout << "   + Data..." << endl;
		TString asample = Form("Tree_%s",sampleName.Data());
		cout << "   + Looking for " << asample << " trees..." << endl;
		myProject->AddDataFiles(dm->GetRealDataFiles(asample));
    G_Event_Weight = 1.;
    G_IsData = true;
  }
  else{ // Deal with MC samples
    cout << "MC sample... " << endl;
    G_IsData = false;
    dm->LoadDataset(sampleName);
   if(!sampleName.BeginsWith("5TeV")) myProject->AddDataFiles(dm->GetFiles());

    if((sampleName.Contains("aMCatNLO") || sampleName.Contains("amcatnlo")) && !sampleName.BeginsWith("5TeV") ){
      G_Event_Weight = dm->GetCrossSection() / dm->GetSumWeights();
			PAF_INFO("RunTree_ReReco5TeV", "This is a MC@NLO sample!");
			G_IsMCatNLO = true;

      cout << endl;
      cout << " weightSum(MC@NLO) = " << dm->GetSumWeights()     << endl;
    }

    else if(sampleName.BeginsWith("5TeV")){
      //TString localpath="/pool/ciencias/users/user/palencia/";
      TString localpath="/mnt_pool/fanae105/user/juanr/Trees_5TeV_april29/";

      sample.ReplaceAll("5TeV_", "");
      myProject->AddDataFile(localpath + "Tree_" + sample + ".root");
      G_Event_Weight = 1;
      G_IsData = false;
      if(sample == "TTbar_Powheg")             G_Event_Weight = 65.   /498000;
      if(sample == "TW")                       G_Event_Weight = 3.04  /429600;
      if(sample == "TbarW")                    G_Event_Weight = 3.04  /460900;
      if(sample == "WZTo3LNU")                 G_Event_Weight = 0.21  /100000;
      if(sample == "WWto2LNu")                 G_Event_Weight = 1.77  / 90000;
      if(sample == "DYJetsToLL_M50_aMCatNLO")  G_Event_Weight = 2055. /6.69891e+09;
      if(sample == "WJetsToLNU_aMCatNLO")      G_Event_Weight = 21159./1.29517e+10;

      if(sample == "TTbar_Powheg_ScaleUp")     G_Event_Weight = 65./462000;
      if(sample == "TTbar_Powheg_ScaleDown")   G_Event_Weight = 65./378500;

      if(sample == "Data_SingleMu"){  G_Event_Weight = 1; G_IsData = true;}
      
    }

    else  G_Event_Weight = dm->GetCrossSection() / dm->GetEventsInTheSample();
    
    if(nEvents == 0) nEvents = dm->GetEventsInTheSample();

    cout << endl;
    cout << " #==============================================="   << endl;
    cout << " #      sampleName = " << sampleName		  << endl;
    cout << " #       x-section = " << dm->GetCrossSection()	  << endl;
    cout << " #      	nevents = " << dm->GetEventsInTheSample() << endl;
    cout << " #  base file name = " << dm->GetBaseFileName()	  << endl;
    cout << " #      	 weight = " << G_Event_Weight		  << endl;
    cout << " #      	 isData = " << G_IsData 		  << endl;
    cout << " #==============================================="   << endl;
    cout << endl;
  }
  
	// Output file name
  //----------------------------------------------------------------------------
  Bool_t G_Use_CSVM = true;
  TString outputDir = "./temp";

  gSystem->mkdir(outputDir, kTRUE);

  std::ostringstream oss;      
  oss << G_Total_Lumi;

  TString LumiString = oss.str();
  TString outputFile = outputDir;
  if     (Selection == 0)  outputFile += "/Tree_" + sampleName          + ".root";
  else if(Selection == 1)  outputFile += "/Tree_" + sampleName + "Semi" + ".root";
  else if(Selection == 4)  outputFile += "/Tree_" + sampleName + "Fidu" + ".root";
  else if(Selection == 5)  outputFile += "/Tree_" + sampleName + "Part" + ".root";
  else                     outputFile += "/Tree_" + sampleName          + ".root";

  PAF_INFO("RunTree_ReReco5TeV", Form("Output file = %s", outputFile.Data()));
  myProject->SetOutputFile(outputFile);


  // Parameters for the analysis
  //----------------------------------------------------------------------------
  myProject->SetInputParam("sampleName",    sampleName       );
  myProject->SetInputParam("IsData",        G_IsData         );
  myProject->SetInputParam("UseCSVM",       G_Use_CSVM       );
  myProject->SetInputParam("weight",        G_Event_Weight   );
  myProject->SetInputParam("LumiForPU",     G_LumiForPUData  );
  myProject->SetInputParam("TotalLumi",     G_Total_Lumi     );
  myProject->SetInputParam("DoSystStudies", DoSystStudies    );
  myProject->SetInputParam("Selection"    , Selection        );
  myProject->SetInputParam("IsMCatNLO"    , G_IsMCatNLO      );  

 
  if(nEvents != 0) myProject->SetNEvents(nEvents);

  // Name of analysis class
  //----------------------------------------------------------------------------
  myProject->AddSelectorPackage("TOP5TeVAnalyzer");

  // Additional packages
  //----------------------------------------------------------------------------
  //myProject->AddPackage("PUWeight");
  //myProject->AddPackage("BTagSFUtil");
  //myProject->AddPackage("LeptonSF");


  // Let's rock!
  //----------------------------------------------------------------------------
  myProject->Run();
}
