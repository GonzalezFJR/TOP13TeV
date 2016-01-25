void BTagSFUtil::GetBTagPayload(TString BTagAlgorithm, TString DataPeriod) {

  TaggerName = BTagAlgorithm;

  TaggerCut = -1.;

  if (TaggerName==     "TCHPT") TaggerCut = 3.41;
  //if (TaggerName==       "JPL") TaggerCut = 0.275;
  //if (TaggerName==       "JPM") TaggerCut = 0.545;
  //if (TaggerName==       "JPT") TaggerCut = 0.790;
  if (TaggerName==      "CSVL") TaggerCut = 0.244;
  //if (TaggerName==      "CSVM") TaggerCut = 0.679;
  //if (TaggerName==      "CSVM") TaggerCut = 0.423;  // this is the cut for the loose WP of the new discriminator iCSVv2
  if (TaggerName==      "CSVM") TaggerCut = 0.89;  // this is the cut for the medium WP of the new discriminator iCSVv2
   if (TaggerName==      "CSVT") TaggerCut = 0.898;
  if (TaggerName==    "CSVV1L") TaggerCut = 0.405;
  if (TaggerName==    "CSVV1M") TaggerCut = 0.783;
  if (TaggerName==    "CSVV1T") TaggerCut = 0.920;
  if (TaggerName==  "CSVSLV1L") TaggerCut = 0.527;
  if (TaggerName==  "CSVSLV1M") TaggerCut = 0.756;
  if (TaggerName==  "CSVSLV1T") TaggerCut = 0.859;

  if (TaggerCut<0.)
    cout << "BTagSFUtil: " << BTagAlgorithm << " is not a supported b-tagging algorithm for Winter13" << endl;

  nBTagPtBins = 16; if (TaggerName.Contains("CSVSLV1")) nBTagPtBins = 13;
  float ptmin[17] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};

  for (int ptbin = 0; ptbin<=nBTagPtBins; ptbin++) 
    BTagPtBinEdge[ptbin] = ptmin[ptbin];
			
  if (TaggerName=="TCHPT") {

    funSFb = new TF1("funSFb", "0.703389*((1.+(0.088358*x))/(1.+(0.0660291*x)))", 20., 800.);
    float Tagger_SFb_error[] = {
      0.0624031,
      0.034023,
      0.0362764,
      0.0341996,
      0.031248,
      0.0281222,
      0.0316684,
      0.0276272,
      0.0208828,
      0.0223511,
      0.0224121,
      0.0261939,
      0.0268247,
      0.0421413,
      0.0532897,
      0.0506714 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];
    
  }

  if (TaggerName=="CSVL") {

    funSFb = new TF1("funSFb", "0.997942*((1.+(0.00923753*x))/(1.+(0.0096119*x)))", 20., 800.);
    float Tagger_SFb_error[] = {
      0.033299,
      0.0146768,
      0.013803,
      0.0170145,
      0.0166976,
      0.0137879,
      0.0149072,
      0.0153068,
      0.0133077,
      0.0123737,
      0.0157152,
      0.0175161,
      0.0209241,
      0.0278605,
      0.0346928,
      0.0350099 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVM") {

    funSFb = new TF1("funSFb", "(0.938887+(0.00017124*x))+(-2.76366e-07*(x*x))", 20., 800.);
    float Tagger_SFb_error[] = {
     /* 0.0415707,
      0.0204209,
      0.0223227,
      0.0206655,
      0.0199325,
      0.0174121,
      0.0202332,
      0.0182446,
      0.0159777,
      0.0218531,
      0.0204688,
      0.0265191,
      0.0313175,
      0.0415417,
      0.0740446,
      0.0596716 };*/
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVT") {

    funSFb = new TF1("funSFb", "(0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x))", 20., 800.);
    float Tagger_SFb_error[] = {
      0.0515703,
      0.0264008,
      0.0272757,
      0.0275565,
      0.0248745,
      0.0218456,
      0.0253845,
      0.0239588,
      0.0271791,
      0.0273912,
      0.0379822,
      0.0411624,
      0.0786307,
      0.0866832,
      0.0942053,
      0.102403 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];
    
  }

  if (TaggerName=="CSVV1L") {

    funSFb = new TF1("funSFb", "1.7586*((1.+(0.799078*x))/(1.+(1.44245*x)))", 20., 800.);
    float Tagger_SFb_error[] = {
      0.0345802,
      0.0152688,
      0.0149101,
      0.0167145,
      0.0167098,
      0.013472,
      0.0146024,
      0.0156735,
      0.0142592,
      0.0147227,
      0.0167101,
      0.0191159,
      0.0360389,
      0.0331342,
      0.0336916,
      0.0298064 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVV1M") {

    funSFb = new TF1("funSFb", "0.952067+(-2.00037e-05*x)", 20., 800.);
    float Tagger_SFb_error[] = {
      0.0376303,
      0.0187774,
      0.019884,
      0.0215849,
      0.0207925,
      0.0180289,
      0.0178674,
      0.0159339,
      0.019042,
      0.020975,
      0.0189178,
      0.0246477,
      0.0291784,
      0.0428437,
      0.0674624,
      0.0479834 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVV1T") {

    funSFb = new TF1("funSFb", "(0.912578+(0.000115164*x))+(-2.24429e-07*(x*x))", 20., 800.);
    float Tagger_SFb_error[] = {
      0.0564014,
      0.0293159,
      0.0315288,
      0.0301526,
      0.0266047,
      0.0240973,
      0.0254404,
      0.0241548,
      0.0233434,
      0.0303961,
      0.040912,
      0.042942,
      0.0440911,
      0.0555312,
      0.105762,
      0.0886457 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVSLV1L") {

    funSFb = new TF1("funSFb", "0.970168*((1.+(0.00266812*x))/(1.+(0.00250852*x)))", 20., 400.);
    float Tagger_SFb_error[] = {
      0.135344,
      0.0288656,
      0.0259088,
      0.0199242,
      0.0189792,
      0.0178341,
      0.0187104,
      0.0239028,
      0.0211104,
      0.017689,
      0.02823,
      0.0259654,
      0.0614497 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVSLV1M") {

    funSFb = new TF1("funSFb", "((0.939238+(0.000278928*x))+(-7.49693e-07*(x*x)))+(2.04822e-10*(x*(x*x)))", 20., 400.);
    float Tagger_SFb_error[] = {
      0.0918443,
      0.0282557,
      0.0264246,
      0.0242536,
      0.0218046,
      0.0207568,
      0.0207962,
      0.0208919,
      0.0200894,
      0.0258879,
      0.0270699,
      0.0256006,
      0.0438219 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  if (TaggerName=="CSVSLV1T") {

    funSFb = new TF1("funSFb", "(0.928257+(9.3526e-05*x))+(-4.1568e-07*(x*x))", 20., 400.);
    float Tagger_SFb_error[] = {
      0.10761,
      0.0333696,
      0.0339123,
      0.0302699,
      0.0261626,
      0.0274243,
      0.0224287,
      0.0239842,
      0.0267866,
      0.0254787,
      0.0317589,
      0.0365968,
      0.0481259 };

    for (int ptbin = 0; ptbin<nBTagPtBins; ptbin++)
      SFb_error[ptbin] = Tagger_SFb_error[ptbin];

  }

  float MaxJetPtLight = 1000.;

  if (TaggerName=="CSVL" || TaggerName=="JPL") {
    
    nBTagEtaBins = 4;
    BTagEtaBinEdge[0] = 0.; BTagEtaBinEdge[1] = 0.5; BTagEtaBinEdge[2] = 1.0; BTagEtaBinEdge[3] = 1.5; 

  } else if (TaggerName=="CSVM" || TaggerName=="JPM") {
    
    nBTagEtaBins = 3;
    BTagEtaBinEdge[0] = 0.; BTagEtaBinEdge[1] = 0.8; BTagEtaBinEdge[2] = 1.6;

  } else {
    
    nBTagEtaBins = 1;
    BTagEtaBinEdge[0] = 0.; 

  }

  if (DataPeriod!="ABCD")
    cout << "BTagSFUtil: " << DataPeriod << " period is not supported for Winter13, use ABCD instead" << endl;

  for (int etabin = 0; etabin<nBTagEtaBins; etabin++) {

    if (BTagEtaBinEdge[etabin]==1.5 || BTagEtaBinEdge[etabin]==1.6) MaxJetPtLight = 850.;

    if( TaggerName == "CSVL" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVL" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVL" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVL" && etabin==3)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVM" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVM" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVM" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVT" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1L" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.03599+(0.00187708*x))+(-3.73001e-06*(x*x)))+(2.09649e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.995735+(0.00146811*x))+(-2.83906e-06*(x*x)))+(1.5717e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.0763+(0.00228243*x))+(-4.61169e-06*(x*x)))+(2.61601e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1L" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.987393+(0.00162718*x))+(-3.21869e-06*(x*x)))+(1.84615e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.947416+(0.00130297*x))+(-2.50427e-06*(x*x)))+(1.41682e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.02741+(0.00194855*x))+(-3.92587e-06*(x*x)))+(2.27149e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1L" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.950146+(0.00150932*x))+(-3.28136e-06*(x*x)))+(2.06196e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.91407+(0.00123525*x))+(-2.61966e-06*(x*x)))+(1.63016e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((0.986259+(0.00178067*x))+(-3.93596e-06*(x*x)))+(2.49014e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1L" && etabin==3)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.01923+(0.000898874*x))+(-2.57986e-06*(x*x)))+(1.8149e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.979782+(0.000743807*x))+(-2.14927e-06*(x*x)))+(1.49486e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.05868+(0.00105264*x))+(-3.00767e-06*(x*x)))+(2.13498e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1M" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.06383+(0.00279657*x))+(-5.75405e-06*(x*x)))+(3.4302e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.971686+(0.00195242*x))+(-3.98756e-06*(x*x)))+(2.38991e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.15605+(0.00363538*x))+(-7.50634e-06*(x*x)))+(4.4624e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1M" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.03709+(0.00169762*x))+(-3.52511e-06*(x*x)))+(2.25975e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.947328+(0.00117422*x))+(-2.32363e-06*(x*x)))+(1.46136e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.12687+(0.00221834*x))+(-4.71949e-06*(x*x)))+(3.05456e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1M" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.01679+(0.00211998*x))+(-6.26097e-06*(x*x)))+(4.53843e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.922527+(0.00176245*x))+(-5.14169e-06*(x*x)))+(3.61532e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.11102+(0.00247531*x))+(-7.37745e-06*(x*x)))+(5.46589e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVV1T" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.15047+(0.00220948*x))+(-5.17912e-06*(x*x)))+(3.39216e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.936862+(0.00149618*x))+(-3.64924e-06*(x*x)))+(2.43883e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.36418+(0.00291794*x))+(-6.6956e-06*(x*x)))+(4.33793e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1L" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.06344+(0.0014539*x))+(-2.72328e-06*(x*x)))+(1.47643e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((1.01168+(0.000950951*x))+(-1.58947e-06*(x*x)))+(7.96543e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.11523+(0.00195443*x))+(-3.85115e-06*(x*x)))+(2.15307e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1L" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.0123+(0.00151734*x))+(-2.99087e-06*(x*x)))+(1.73428e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.960377+(0.00109821*x))+(-2.01652e-06*(x*x)))+(1.13076e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.06426+(0.0019339*x))+(-3.95863e-06*(x*x)))+(2.3342e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1L" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.975277+(0.00146932*x))+(-3.17563e-06*(x*x)))+(2.03698e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.931687+(0.00110971*x))+(-2.29681e-06*(x*x)))+(1.45867e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.0189+(0.00182641*x))+(-4.04782e-06*(x*x)))+(2.61199e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1L" && etabin==3)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.04201+(0.000827388*x))+(-2.31261e-06*(x*x)))+(1.62629e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.992838+(0.000660673*x))+(-1.84971e-06*(x*x)))+(1.2758e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.09118+(0.000992959*x))+(-2.77313e-06*(x*x)))+(1.9769e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1M" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.06212+(0.00223614*x))+(-4.25167e-06*(x*x)))+(2.42728e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.903956+(0.00121678*x))+(-2.04383e-06*(x*x)))+(1.10727e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.22035+(0.00325183*x))+(-6.45023e-06*(x*x)))+(3.74225e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1M" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.04547+(0.00216995*x))+(-4.579e-06*(x*x)))+(2.91791e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.900637+(0.00120088*x))+(-2.27069e-06*(x*x)))+(1.40609e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.19034+(0.00313562*x))+(-6.87854e-06*(x*x)))+(4.42546e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1M" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.991865+(0.00324957*x))+(-9.65897e-06*(x*x)))+(7.13694e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.868875+(0.00222761*x))+(-6.44897e-06*(x*x)))+(4.53261e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.11481+(0.00426745*x))+(-1.28612e-05*(x*x)))+(9.74425e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "CSVSLV1T" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.09494+(0.00193966*x))+(-4.35021e-06*(x*x)))+(2.8973e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.813331+(0.00139561*x))+(-3.15313e-06*(x*x)))+(2.12173e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.37663+(0.00247963*x))+(-5.53583e-06*(x*x)))+(3.66635e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPL" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.991991+(0.000898777*x))+(-1.88002e-06*(x*x)))+(1.11276e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.930838+(0.000687929*x))+(-1.36976e-06*(x*x)))+(7.94486e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.05319+(0.00110776*x))+(-2.38542e-06*(x*x)))+(1.42826e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPL" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.96633+(0.000419215*x))+(-9.8654e-07*(x*x)))+(6.30396e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.904781+(0.000324913*x))+(-7.2229e-07*(x*x)))+(4.52185e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.0279+(0.00051255*x))+(-1.24815e-06*(x*x)))+(8.07098e-10*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPL" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.968008+(0.000482491*x))+(-1.2496e-06*(x*x)))+(9.02736e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.914619+(0.000330357*x))+(-8.41216e-07*(x*x)))+(6.14504e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.02142+(0.000633484*x))+(-1.6547e-06*(x*x)))+(1.18921e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPL" && etabin==3)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.991448+(0.000765746*x))+(-2.26144e-06*(x*x)))+(1.65233e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.933947+(0.000668609*x))+(-1.94474e-06*(x*x)))+(1.39774e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.04894+(0.000861785*x))+(-2.57573e-06*(x*x)))+(1.90702e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPM" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.991457+(0.00130778*x))+(-2.98875e-06*(x*x)))+(1.81499e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.822012+(0.000908344*x))+(-1.89516e-06*(x*x)))+(1.1163e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.16098+(0.00170403*x))+(-4.07382e-06*(x*x)))+(2.50873e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPM" && etabin==1)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.00576+(0.00121353*x))+(-3.20601e-06*(x*x)))+(2.15905e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.845597+(0.000734909*x))+(-1.76311e-06*(x*x)))+(1.16104e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.16598+(0.00168902*x))+(-4.64013e-06*(x*x)))+(3.15214e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPM" && etabin==2)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.939038+(0.00226026*x))+(-7.38544e-06*(x*x)))+(5.77162e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.803867+(0.00165886*x))+(-5.19532e-06*(x*x)))+(3.88441e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.07417+(0.00285862*x))+(-9.56945e-06*(x*x)))+(7.66167e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "JPT" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((0.953235+(0.00206692*x))+(-5.21754e-06*(x*x)))+(3.44893e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.642947+(0.00180129*x))+(-4.16373e-06*(x*x)))+(2.68061e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.26372+(0.0023265*x))+(-6.2548e-06*(x*x)))+(4.20761e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }
    if( TaggerName == "TCHPT" && etabin==0)
      {
	funSFlight[etabin][1] = new TF1("SFlight"   ,"((1.20175+(0.000858187*x))+(-1.98726e-06*(x*x)))+(1.31057e-09*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][0] = new TF1("SFlightMin","((0.968557+(0.000586877*x))+(-1.34624e-06*(x*x)))+(9.09724e-10*(x*(x*x)))", 20.,MaxJetPtLight);
	funSFlight[etabin][2] = new TF1("SFlightMax","((1.43508+(0.00112666*x))+(-2.62078e-06*(x*x)))+(1.70697e-09*(x*(x*x)))", 20.,MaxJetPtLight);
      }  

  }

}
