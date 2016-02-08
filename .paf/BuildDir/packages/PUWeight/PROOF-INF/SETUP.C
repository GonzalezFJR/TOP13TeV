

// Template for SETUP.C for packages
Int_t SETUP() {

  // Load library
#if 0
  cout << "Loading PUWeight..." << endl;
#endif

  if (gSystem->Load("libPUWeight") == -1)
     return -1;

#if 0
  // Enlarge include path (for PUWeight.h)
  TString pwdpath = "-I";
  pwdpath+=gSystem->pwd(); //path to local dir
  TString curpath = gSystem->GetIncludePath(); //current path

  //Some info
  if (gProofDebugLevel > 0) {
    TString infomess="Current PWD path is \"";
    infomess+=pwdpath;
    infomess+="\"";
    Info("SETUP.C",infomess.Data());
    infomess = "Current include path is \"";
    infomess+=curpath;
    infomess+="\"";
    Info("SETUP.C",infomess.Data());
    //
  }

  if (!curpath.Contains(pwdpath)) {
    gSystem->AddIncludePath(pwdpath);
    
    if (gProofDebugLevel > 0)
      cout << ">> Include path set to: \"" << gSystem->GetIncludePath() 
       << "\"" << endl;
  }
#endif 

  return 0;
}



