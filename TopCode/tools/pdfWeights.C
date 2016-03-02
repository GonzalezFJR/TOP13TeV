using namespace std; 

// Constants
//const TString path_sum = "/mnt_pool/fanae105/user/juanr/top13TeV_12jun/TopTrees/";
//const TString path_sum = "/nfs/fanae/user/jfernan/www/PUhistos/DR74X/50ns/";
// const TString path = "/mnt_pool/fanae105/user/palencia/april14Run2/TOP/TopTrees/aug10/";
//const TString path = "/mnt_pool/fanae105/user/palencia/april14Run2/TOP/TopTrees/sept25/";
// const TString path = "/mnt_pool/fanae105/user/juanr/top13TeV_12jun/30jul/TopTrees/temp/";
const TString path = "/mnt_pool/fanae105/user/palencia/testHeppy/TOP/TopTrees/mar1/";
float lumi = 2223; // pb-1 
const TString pathsum = "/pool/ciencias/heppyTreesDR76X/v1/";
//"/pool/ciencias/TreesDR74X/heppyTrees/v1/";
const TString sampleName = "Tree_13TeV_EA_TTbar_Powheg.root";
//const TString sampleNameSum = "Histo_TTbar_Powheg.root";

const int nWeights = 248;
const int nChannels = 3;
const int nLevels = 5;
float xsec = 831.8;

// Channels
const TString chan[nChannels] = {
	"Elec",	// 0
	"Muon",	// 1
	"ElMu"	// 2
};

// Levels
const TString level[nLevels] = {
	"dilepton",	// 0
	"ZVeto",	// 1
	"MET",		// 2
	"2jets",	// 3
	"1btag"	// 4
};

TH1* h;
float fac = 0;
TFile* inputfile0 = TFile::Open(pathsum + "Tree_TTbar_Powheg_0.root")->GetObject("CountLHE",h);
bool sumnom = 0;

void GetSumNom(){
  TH1* htemp; 
  for(int i = 1; i < 12; i++){
    TFile* inputfile = TFile::Open(pathsum + Form("Tree_TTbar_Powheg_%i.root", i))->GetObject("CountLHE",htemp);
    h->Add(htemp);
    htemp->SetDirectory(0);
    //cout << "i = " << i << ",  h->GetEntries() = " << h->GetEntries() << endl;
    delete htemp; delete inputfile;
  }
  fac = h->GetBinContent(1002);
}
  
double GetWeightSum(Int_t i){ 
  if(!sumnom){
    GetSumNom();
    sumnom = true;
  }
  //TH1* hist_w;
  //TFile* inputfile = TFile::Open(pathsum);
  //inputfile->GetObject("CountLHE",hist_w);
  double weight = 0.; int j = 0;
  if      (i<10  ) j = i + 1001;   // 1002-1010: muRmuF
  else if (i<112) j = i + 1992;   // 2002-2103: NNPDF
  else if (i<167) j = i + 2890;   // 3002-3056: CT10
  else if (i<223) j = i + 3835;   // 4000-4057: MMHT2014
  else if (i<250) j = i + 4779;   // 5002-5028: muRmuF, hdamp 
  //(CountLHE->GetBinContent(1002)/CountLHE->GetBinContent(j)
  weight = h->GetBinContent(j);
	return weight;// /eventCount;
}

float Yield(Int_t ip_chan, Int_t ip_level){ 
  TH1F* hist; 
  TFile* inputfile = TFile::Open(path + sampleName);
  inputfile->GetObject("H_Yields_" + chan[ip_chan],hist);
  float yield = hist->GetBinContent(ip_level+1);
  delete inputfile;
  return yield*lumi;
}
/*
float Yield0(Int_t ip_level){
  TH1F* hist;
  TFile* inputfile = TFile::Open(path + sampleName);
//  inputfile->GetObject("H_Yields_" + chan[ip_chan],hist);
  inputfile->GetObject("H_nominal_" + level[ip_level],hist);
  float yield = hist->GetBinContent(1);
	delete inputfile;
	return yield;//fac;
}*/


float GetWeight(TString chan, TString lev, Int_t index){ 
  // index from 1 to 249
  TH1F* hist; 
  TFile* inputfile = TFile::Open(path + sampleName);
  inputfile->GetObject("H_LHEweights"+ chan + "_" + lev,hist);
  float weight = hist->GetBinContent(index)*lumi;
	inputfile -> Close();
	delete inputfile;
	return weight/GetWeightSum(index)*fac;
}

void line() { cout << "--------------------------------------------------------------------------------------------------------" << endl;}
void line2(){ cout << "========================================================================================================" << endl;}
void printi(int i){ cout << Form(" i = %3i || ", i);}
void printchan(){ cout << "         Elec        ||         Muon        ||        ElMu        " << endl;}
void print3(TString lev, int i, int inom){ 
  float v1 = GetWeight("Elec", lev, i);
  float v2 = GetWeight("Muon", lev, i);
  float v3 = GetWeight("ElMu", lev, i);
  float e1 = TMath::Abs(GetWeight("Elec", lev, inom)-v1)/GetWeight("Elec", lev, inom)*100;
  float e2 = TMath::Abs(GetWeight("Muon", lev, inom)-v2)/GetWeight("Muon", lev, inom)*100;
  float e3 = TMath::Abs(GetWeight("ElMu", lev, inom)-v3)/GetWeight("ElMu", lev, inom)*100;
  cout << Form("  %6.2f (%2.2f %)  ||  %6.2f (%2.2f %)  ||  %6.2f (%2.2f %)  ", v1, e1, v2, e2, v3, e3);
  if(i == inom) cout << " <---- NOMINAL ";
  cout << endl;
}



void muRmuFvariations(TString lev){
  line2();
  cout << " ### muR/muF ME variations (hdamp = mtop), nominal + 8" << endl;
  line();
  cout << "                                      "; printchan();
  line();
  printi(1);  cout << " muR = 1  , muF = 1    || " ; print3(lev, 1 ,1);
  printi(2);  cout << " muR = 1  , muF = 2    || " ; print3(lev, 2 ,1);
  printi(3);  cout << " muR = 1  , muF = 0.5  || " ; print3(lev, 3 ,1);
  printi(4);  cout << " muR = 2  , muF = 1    || " ; print3(lev, 4 ,1);
  printi(5);  cout << " muR = 2  , muF = 2    || " ; print3(lev, 5 ,1);
  printi(6);  cout << " muR = 2  , muF = 0.5  || " ; print3(lev, 6 ,1);
  printi(7);  cout << " muR = 0.5, muF = 1    || " ; print3(lev, 7 ,1);
  printi(8);  cout << " muR = 0.5, muF = 2    || " ; print3(lev, 8 ,1);
  printi(9);  cout << " muR = 0.5, muF = 0.5  || " ; print3(lev, 9 ,1);
  line();
  float nom = 0; float y = 0; float vmax = 0;
  TString o = " Maximum variation                  ";
  for(int i = 0; i < nChannels; i++){
    nom = GetWeight(chan[i], lev, 1); vmax = nom;
    for(int k = 0; k<9; k++) {
       if (TMath::Abs(GetWeight(chan[i], lev, k+1)-nom) > TMath::Abs(vmax-nom)) vmax = GetWeight(chan[i], lev, k+1);
    }
    o += Form("||  %6.2f (%2.2f %)  ", TMath::Abs(vmax-nom), TMath::Abs(vmax-nom)/nom*100);
  }
  //cout << Form(" Maximum variation   %6.2f (%2.3f %)  ||  %6.2f (%2.3f %)  ||  %6.2f (%2.3f %)  ", v1, e1, v2, e2, v3, e3)" << endl;
  cout << o << endl;
  line2();
}

void NNPDFvariations(TString lev){
  line2();
  cout << " ### NNPDF variations: 100 + 2 (alpha_s) " << endl;
  line();
  cout << "            "; printchan();
  line();
  for(int i = 10; i<110; i++){
    printi(i); print3(lev, i, 1);
  }
  line();
  cout << " - NNPDF alpha_s variations " << endl; 
  line();
  cout << "           "; printchan();
  line();
  printi(110); print3(lev, 110, 1);
  printi(111); print3(lev, 111, 1);
  line2();
}

void CT10variations(TString lev){
  line2();
  cout << " ### CT10 variations: 1 (nominal) + 2x26 (up/down) + 2 (alpha_s) " << endl;
  line();
  cout << "            "; printchan();
  line();
  for(int i = 112; i< 112 + (2*26) + 1; i++){
    printi(i); print3(lev, i, 112);
  }
  line();
  cout << " - CT10 alpha_s variations " << endl;
  line();
  cout << "            "; printchan();
  line();
  printi(165); print3(lev, 165, 112);
  printi(166); print3(lev, 166, 112);
  line2();
}

void MMHTvariations(TString lev){ 
  line2();
  cout << " ### MMHT2014 variations: 1 (nominal) + 25x2 (up/down) + 5 (alpha_s)" << endl;
  line();
  cout << "            "; printchan();
  line();
  for(int i = 167; i< 167 + 50 + 1; i++){
    printi(i); print3(lev, i, 167);
  }
  line();
  cout << " - MMHT2014 alpha_s variations " << endl;
  line();
  cout << "             "; printchan();
  line();
  printi(218); print3(lev, 217, 167);
  printi(219); print3(lev, 218, 167);
  printi(220); print3(lev, 220, 167);
  printi(221); print3(lev, 221, 167);
  printi(222); print3(lev, 222, 167);
  line2();
}

void hdampvariations(TString lev){
  line2();
  cout << " ### hdamp + ME variations: 3x9 (hdamp = 0, mtop/2, mtop*2)" << endl;
  line();
  cout << " hdamp = 0                            "; printchan();
  line();
  printi(223);  cout << " muR = 1  , muF = 1    || " ; print3(lev, 223 ,223);
  printi(224);  cout << " muR = 1  , muF = 2    || " ; print3(lev, 224 ,223);
  printi(225);  cout << " muR = 1  , muF = 0.5  || " ; print3(lev, 225 ,223);
  printi(226);  cout << " muR = 2  , muF = 1    || " ; print3(lev, 226 ,223);
  printi(227);  cout << " muR = 2  , muF = 2    || " ; print3(lev, 227 ,223);
  printi(228);  cout << " muR = 2  , muF = 0.5  || " ; print3(lev, 228 ,223);
  printi(229);  cout << " muR = 0.5, muF = 1    || " ; print3(lev, 229 ,223);
  printi(230);  cout << " muR = 0.5, muF = 2    || " ; print3(lev, 230 ,223);
  printi(231);  cout << " muR = 0.5, muF = 0.5  || " ; print3(lev, 231 ,223);
  line();
  cout << " hdamp = mtop/2                       "; printchan();
  line();
  printi(232);  cout << " muR = 1  , muF = 1    || " ; print3(lev, 232 ,232);
  printi(233);  cout << " muR = 1  , muF = 2    || " ; print3(lev, 233 ,232);
  printi(234);  cout << " muR = 1  , muF = 0.5  || " ; print3(lev, 234 ,232);
  printi(235);  cout << " muR = 2  , muF = 1    || " ; print3(lev, 235 ,232);
  printi(236);  cout << " muR = 2  , muF = 2    || " ; print3(lev, 236 ,232);
  printi(237);  cout << " muR = 2  , muF = 0.5  || " ; print3(lev, 237 ,232);
  printi(238);  cout << " muR = 0.5, muF = 1    || " ; print3(lev, 238 ,232);
  printi(239);  cout << " muR = 0.5, muF = 2    || " ; print3(lev, 239 ,232);
  printi(240);  cout << " muR = 0.5, muF = 0.5  || " ; print3(lev, 240 ,232);
  line();
  cout << " hdamp = mtop*2                       "; printchan();
  line();
  printi(241);  cout << " muR = 1  , muF = 1    || " ; print3(lev, 241 ,241);
  printi(242);  cout << " muR = 1  , muF = 2    || " ; print3(lev, 242 ,241);
  printi(243);  cout << " muR = 1  , muF = 0.5  || " ; print3(lev, 243 ,241);
  printi(244);  cout << " muR = 2  , muF = 1    || " ; print3(lev, 244 ,241);
  printi(245);  cout << " muR = 2  , muF = 2    || " ; print3(lev, 245 ,241);
  printi(246);  cout << " muR = 2  , muF = 0.5  || " ; print3(lev, 246 ,241);
  printi(247);  cout << " muR = 0.5, muF = 1    || " ; print3(lev, 247 ,241);
  printi(248);  cout << " muR = 0.5, muF = 2    || " ; print3(lev, 248 ,241);
  printi(249);  cout << " muR = 0.5, muF = 0.5  || " ; print3(lev, 249 ,241);
  line2();
}


void pdfWeights(TString lev){
  cout << endl;
  cout << "Sample : " << path + sampleName << endl;
  cout << "Level: " << lev << endl;
  cout << "Lumi : " << lumi << endl;
  cout << endl;
  muRmuFvariations(lev);
  cout << endl;
  NNPDFvariations(lev);
  cout << endl;
  printNNPDFsyst(lev);
  cout << endl;
  CT10variations(lev);
  cout << endl;
  printCT10syst(lev); 
  cout << endl;
  MMHTvariations(lev);
  cout << endl;
  hdampvariations(lev);
  cout << endl;
}

float NNPDFsyst(TString ch, TString lev){
  float e = 0; float n = GetWeight(ch, lev, 1);
  float y = 0;
  for(int i = 10; i<110; i++){
    y = GetWeight(ch, lev, i);
    e += (y-n)*(y-n);
  }
  float v110 = TMath::Abs(GetWeight(ch, lev, 110)-n);
  float v111 = TMath::Abs(GetWeight(ch, lev, 111)-n);
  float rms = TMath::Sqrt(e/100);
  return TMath::Sqrt(rms*rms + ((v110-v111)*0.75/2)*((v110-v111)*0.75/2));
}

void printNNPDFsyst(TString lev){
  float nee = GetWeight("Elec", lev, 1); float eee = NNPDFsyst("Elec", lev);
  float nmm = GetWeight("Muon", lev, 1); float emm = NNPDFsyst("Muon", lev);
  float nem = GetWeight("ElMu", lev, 1); float eem = NNPDFsyst("ElMu", lev);
  line2();
  cout << " >>>> NNPDF systematic uncertainty" << endl;
  cout << " Evaluated by taking the RMS under the 100 weights" << endl;
  cout << " Alpha_s variations are added in quadrature after rescaling by 0.75" << endl;
  cout << " The formula is: sqrt(RMS^2 + ((alphas var 1 - alphas var 2)*0.75/2)^2 )" << endl;
  line(); printchan(); line();
  cout << Form("     %4.2f (%2.3f %)  ||  %4.2f (%2.3f %)  ||  %4.2f (%2.3f %)  ", eee, eee/nee*100, emm, emm/nmm*100, eem, eem/nem*100) << endl; 
  line2();
}

TString CT10syst(TString ch, TString lev){
  float Dp = 0; float Dm = 0;
  float n = GetWeight(ch, lev, 112);
  float max(float a, float b, float c){
    float d = (a<b)? b:a;
    return(c<d)? d:c;
  }

  float y1 = 0; float y2 = 0; 
  for(int i = 113; i < 113+2*26; i+=2){
    y1 = GetWeight(ch, lev, i);
    y2 = GetWeight(ch, lev, i+1);
    Dp += max(y1-n, y2-n, 0)**2; 
    Dm += max(n-y1, n-y2, 0)**2; 
  } 
  return TString( Form("DX^+_max = %4.2f (%2.3f %),   DX^-_max = %4.2f (%2.3f %)", TMath::Sqrt(Dp), TMath::Sqrt(Dp)/n*100, TMath::Sqrt(Dm), TMath::Sqrt(Dm)/n*100) );  
}

void printCT10syst(TString lev){
  line2();
  cout << " >>>> CT10 systematic uncertainty" << endl;
  cout << " DX_max defined in eq 3 and 4 in: http://arxiv.org/pdf/hep-ph/0605240v2.pdf" << endl;
  line(); 
  cout << "  Elec:    " << CT10syst("Elec", lev) << endl;
  cout << "  Muon:    " << CT10syst("Muon", lev) << endl;
  cout << "  ElMu:    " << CT10syst("ElMu", lev) << endl;
  line2();
}

