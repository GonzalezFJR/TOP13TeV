// Load DatasetManager in ROOT 6
R__LOAD_LIBRARY(DatasetManager/DatasetManager.C+)
#include "DatasetManager/DatasetManager.h"

/********************************************************
 * Main function
 ********************************************************/
void RunTree_ReReco(TString  sampleName     = "TTbar_Madgraph",
			Int_t    nSlots         =  1,
			Int_t    Selection      = 0,
			Bool_t   DoSystStudies  =  false,
			Long64_t nEvents        = 0) {

  // Selection:
  // 0: dilepton
  // 1: ngenLep < 2
  // 2: emu
  // 3: ee + mumu
  // 4: emu + fiducial
  // 5: particle level
  
  // VARIABLES TO BE USED AS PARAMETERS...
  Float_t G_Total_Lumi    = 19664.225;   
  Float_t G_Event_Weight  = 1.0;         
  Bool_t  G_IsData        = false;       
  Float_t G_LumiForPUData = 19468.3;     // luminosity in http://www.hep.uniovi.es/jfernan/PUhistos
  Bool_t  G_IsMCatNLO     = false;

  // PAF mode
  //----------------------------------------------------------------------------
  cout << endl; 
  PAFIExecutionEnvironment* pafmode = 0;
  if (nSlots <=1 ) {
    PAF_INFO("RunTree_ReReco", "Sequential mode chosen");
    pafmode = new PAFSequentialEnvironment();
  }
  else if (nSlots <=64) {
    PAF_INFO("RunTree_ReReco", "PROOF Lite mode chosen");
    pafmode = new PAFPROOFLiteEnvironment(nSlots);
  }
  else {
    PAF_INFO("RunTree_ReReco", "PoD mode chosen");
    pafmode = new PAFPoDEnvironment(nSlots);
  }

  // Create PAF Project whith that environment
  //----------------------------------------------------------------------------
  PAFProject* myProject = new PAFProject(pafmode);

  // INPUT DATA SAMPLE
  //----------------------------------------------------------------------------
  TString userhome = "/mnt_pool/fanae105/user/$USER/";
  DatasetManager* dm = DatasetManager::GetInstance();
  //dm->SetTab("DR80XasymptoticMiniAODv2");
  dm->SetTab("DR80XSummer16asymptoticMiniAODv2_2");
  //dm->RedownloadFiles();

  // Deal with data samples
  if ((sampleName == "DoubleEG"   ||
       sampleName == "DoubleMuon" ||
       sampleName == "MuonEG"     ||
       sampleName.Contains("Single"))){
    cout << "   + Data..." << endl;
    
    TString datasuffix[] = {
      "16B_03Feb2017",
      "16C_03Feb2017",
      "16D_03Feb2017",
      "16E_03Feb2017",
      "16F_03Feb2017",
      "16G_03Feb2017",
      "16H_03Feb2017_v2",
      "16H_03Feb2017_v3",
    };
    const unsigned int nDataSamples = 8;
    for(unsigned int i = 0; i < nDataSamples; i++) {
      TString asample = Form("Tree_%s_%s",sampleName.Data(), datasuffix[i].Data());
      cout << "   + Looking for " << asample << " trees..." << endl;
      myProject->AddDataFiles(dm->GetRealDataFiles(asample));
    }
    G_Event_Weight = 1.;
    G_IsData = true;
  }
  else{ // Deal with MC samples
    G_IsData = false; //true;  // only for pseudodata
    dm->LoadDataset(sampleName);
    if(sampleName != "TestHeppy")   myProject->AddDataFiles(dm->GetFiles());
    if(
	    sampleName.Contains("TTW")  || 
	    sampleName.Contains("TTZ")  || 
			sampleName.Contains("aMCatNLO") ||
			sampleName.Contains("amcatnlo") ){
			G_Event_Weight = dm->GetCrossSection() / dm->GetSumWeights();

			cout << endl;
      cout << " weightSum(MC@NLO) = " << dm->GetSumWeights()     << endl;
    }
    else if(sampleName == "TestHeppy"){
			TString localpath="/pool/ciencias/users/user/palencia/";
			TString sample = "treeTtbar_jan19.root";
			myProject->AddDataFile(localpath + sample);
			G_Event_Weight = 1;
		}
		else {       
			G_Event_Weight = dm->GetCrossSection() / dm->GetEventsInTheSample();
		}
    
    if(nEvents == 0) nEvents = dm->GetEventsInTheSample();

    cout << endl;
    cout << " #===============================================" << endl;
    cout << " #     sampleName = " << sampleName                 << endl;
		cout << " #	     x-section = " << dm->GetCrossSection()      << endl;
		cout << " #	       nevents = " << dm->GetEventsInTheSample() << endl;
    cout << " # base file name = " << dm->GetBaseFileName()      << endl;
		cout << " #	        weight = " << G_Event_Weight	           << endl;
		cout << " #	        isData = " << G_IsData	                 << endl;
    cout << " #===============================================" << endl;
		cout << endl;
  }
  
	// Output file name
  //----------------------------------------------------------------------------
	Bool_t G_Use_CSVM = true;
	TString outputDir = "./temp";
  if(sampleName.BeginsWith("T2tt")) outputDir += "/Susy";

	gSystem->mkdir(outputDir, kTRUE);

	std::ostringstream oss;      
	oss << G_Total_Lumi;

	TString LumiString = oss.str();
	TString outputFile = outputDir;
	if     (Selection == 0)  outputFile += "/Tree_" + sampleName          + ".root";
  else if(Selection == 1)  outputFile += "/Tree_" + sampleName + "Semi" + ".root";
  else if(Selection == 4)  outputFile += "/Tree_" + sampleName + "Fidu" + ".root";
  else                     outputFile += "/Tree_" + sampleName          + ".root";
  if(outputFile.Contains("_ext2")) outputFile.ReplaceAll("_ext2","");
  if(outputFile.Contains("_ext"))  outputFile.ReplaceAll("_ext","");

  PAF_INFO("RunTree_ReReco", Form("Output file = %s", outputFile.Data()));
  myProject->SetOutputFile(outputFile);

  if(sampleName == "WJetsToLNu_aMCatNLO" || sampleName == "DYJetsToLL_M10to50_aMCatNLO_ext" || sampleName == "DYJetsToLL_M50_aMCatNLO" || sampleName == "TTJets_amcatnlo" || sampleName.Contains("aMCatNLO") || sampleName.Contains("amcatnlo" || sampleName.Contains("TTZ") || sampleName.Contains("TTW"))){ //||
        //     sampleName == "TTWToLNu"  || sampleName == "TTWToQQ" || sampleName == "TTZToQQ" || sampleName == "WWZ" || sampleName == "WZZ" || sampleName == "ZZZ"){
    PAF_INFO("RunTree_ReReco", "This is a MC@NLO sample!");
    G_IsMCatNLO = true;
  }

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
  myProject->AddSelectorPackage("TreeAnalysisTop");
  //myProject->AddSelectorPackage("StopAnalyzer");
  //myProject->AddSelectorPackage("WZ13TeV");

  // Additional packages
  //----------------------------------------------------------------------------
  myProject->AddPackage("mt2");
  myProject->AddPackage("PUWeight");
  myProject->AddPackage("BTagSFUtil");
  myProject->AddPackage("LeptonSF");


  // Let's rock!
  //----------------------------------------------------------------------------
  myProject->Run();
}
