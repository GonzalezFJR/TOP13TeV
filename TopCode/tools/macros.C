#include <iomanip>
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLine.h"
#include "TChain.h"
#include "TLatex.h"
#include "stdio.h"
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <TMath.h>
#include <TMatrix.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TRandom3.h>
#include <fstream>
#include <string>


// root -l -b -q makePlots.C

void calcOverflow(TH1F* h){

  Int_t nBins             = h->GetNbinsX();
  
  Double_t lastBinContent = h->GetBinContent(nBins);
  Double_t overflow       = h->GetBinContent(nBins+1);
  h->SetBinContent(nBins, lastBinContent + overflow);

  Double_t lastBinError  = h->GetBinError(nBins);
  Double_t overflowError = h->GetBinError(nBins+1);
  h->SetBinError(nBins, sqrt(lastBinError*lastBinError + overflowError*overflowError));
  
  //cout << " last: " << lastBinContent << " +- " << lastBinError << endl; 
  //cout << " over: " << overflow       << " +- " << overflowError << endl;
  //cout << " new last " << lastBinContent + overflow       << " +- " << sqrt(lastBinError*lastBinError + overflowError*overflowError) << endl;
  
  Double_t  underflow       = h->GetBinContent(0);
  Double_t  firstBinContent = h->GetBinContent(1);
  h->SetBinContent(1,  firstBinContent+underflow);

  Double_t underflowError = h->GetBinError(0);
  Double_t firstBinError  = h->GetBinError(1);
  h->SetBinError(1, sqrt(firstBinError*firstBinError + underflowError*underflowError));
}


void doStackPlots()
{
   gStyle->SetOptStat(0);
   //gStyle->SetErrorX(0);

   TString dir = "./TopTrees/feb24/";
   //TString dir = "./TopTrees/dec14/";
     
   // $ROOTSYS/bin/hadd Tree_13TeV_EA_MuonEGsum.root     Tree_13TeV_EA_MuonEG.root Tree_13TeV_EA_SingleEle.root Tree_13TeV_EA_SingleMu.root
   // $ROOTSYS/bin/hadd Tree_13TeV_EA_DoubleMuonSum.root Tree_13TeV_EA_DoubleMuon.root Tree_13TeV_EA_SingleMu.root
   // $ROOTSYS/bin/hadd Tree_13TeV_EA_DoubleEGsum.root   Tree_13TeV_EA_DoubleEG.root Tree_13TeV_EA_SingleEle.root
   
   TFile *ttdil  = new TFile(dir+"Tree_13TeV_EA_TTbar_Powheg.root"); 
   //TFile *ttdil  = new TFile(dir+"Tree_13TeV_EA_TTJets_aMCatNLO.root"); 
   //TFile *ttdil  = new TFile(dir+"Tree_13TeV_EA_TTbar_Powheg_Herwig.root"); 
   
   TFile *ttsem  = new TFile(dir+"Tree_13TeV_EA_TTbarSemi_Powheg.root"); 
   TFile *ww     = new TFile(dir+"Tree_13TeV_EA_WW.root");
   TFile *wz     = new TFile(dir+"Tree_13TeV_EA_WZ.root");
   TFile *zz     = new TFile(dir+"Tree_13TeV_EA_ZZ.root");
   TFile *tW     = new TFile(dir+"Tree_13TeV_EA_TW.root");
   TFile *tbarW  = new TFile(dir+"Tree_13TeV_EA_TbarW.root");
   TFile *DY     = new TFile(dir+"Tree_13TeV_EA_DYJetsToLL_M50_aMCatNLO.root");
   TFile *DYlow  = new TFile(dir+"Tree_13TeV_EA_DYJetsToLL_M10to50_aMCatNLO_ext.root");
   TFile *wjets  = new TFile(dir+"Tree_13TeV_EA_WJetsToLNu_aMCatNLO.root");
   TFile *data   = new TFile(dir+"Tree_13TeV_EA_MuonEGsum.root"); 
   /*
   TFile *ttdil  = new TFile(dir+"Tree_13TeV_EA_TTJets.root"); 
   TFile *ttsem  = new TFile(dir+"Tree_13TeV_EA_TTJetsSemi.root"); 
   TFile *ww     = new TFile(dir+"Tree_13TeV_EA_WW.root");
   TFile *wz     = new TFile(dir+"Tree_13TeV_EA_WZ.root");
   TFile *zz     = new TFile(dir+"Tree_13TeV_EA_ZZ.root");
   TFile *tW     = new TFile(dir+"Tree_13TeV_EA_TW.root");
   TFile *tbarW  = new TFile(dir+"Tree_13TeV_EA_TbarW.root");
   TFile *DY     = new TFile(dir+"Tree_13TeV_EA_DYJetsToLL_M50_aMCatNLO.root");
   TFile *DYlow  = new TFile(dir+"Tree_13TeV_EA_DYJetsToLL_M10to50_aMCatNLO.root");
   TFile *wjets  = new TFile(dir+"Tree_13TeV_EA_WJetsToLNu_aMCatNLO.root");
   TFile *data   = new TFile(dir+"Tree_13TeV_EA_MuonEG.root"); 
*/
   float fixBR = (0.108*9.)*(0.108*9.); // only for aMC@NLO

   const int nVars = 14, nCuts = 3;


   //TString vars[nVars] = {"Yields", "Vtx", "goodVtx", "Lep0Iso", "Lep1Iso", "Lep0Pt", "Lep1Pt", "Lep0Eta", "Lep1Eta", "DiLepPt", "InvMass", 
   //                    "Jet0Pt", "Jet1Pt", "Jet0Eta", "Jet1Eta", "NJets", "NBtagJets", "NBtagsNJets", "CSVTag", "MET", "HT", 
   //	                 "AbsDelPhiLeps"};
   
   TString vars[nVars] = {"NJets", "HT", "Lep0Pt", "Lep1Pt", "Lep0Eta", "Lep1Eta", "Jet0Pt", "Jet1Pt",
   "Jet0Eta", "Jet1Eta", "NBtagJets", "MET", "AbsDelPhiLeps", "InvMass"};
   
   
   //TString myCut[nCuts] = {"dilepton", "2jets","1btag"};   // _dilepton, _ZVeto, _MET, _2jets, _1btag

   TString cut = "", xtitle = "";

   float lumi = 2173.0*1.023;

   //const unsigned int n=5;   // number of mass points

   TCanvas *c= new TCanvas("c","c",10,100,800,600);
   c->Divide(1,2);
   
   //Plot Pad
   TPad *plot = (TPad*)c->GetPad(1); 
   //plot->SetPad(0.01, 0.23, 0.99, 0.99);
   plot->SetPad(0.0, 0.23, 1.0, 1.0);
   //plot->SetPad(0.0, 0.23, 1.0, 1.0);
   plot->SetTopMargin(0.07);
   plot->SetRightMargin(0.025);
   plot->SetLeftMargin(0.12);
  
   TPad *ratio = (TPad*)c->GetPad(2);
   //ratio->SetPad(0.01, 0.02, 0.99, 0.3);
   ratio->SetPad(0.0, 0.0, 1.0, 0.29);
   //ratio->SetGridx();
   ratio->SetGridy();
   ratio->SetTopMargin(0.03);
   ratio->SetBottomMargin(0.4);
   ratio->SetRightMargin(0.025);
   ratio->SetLeftMargin(0.12);
   
   c->Print("figs/stackPlots.ps["); 

   plot->cd();
   
   for(int i=0; i<nVars; i++){
   
      //if(!(vars[i] == "NJets")) continue;
      
      if     (vars[i] == "NJets") 
         cut="_dilepton";
      else if(vars[i] == "Lep0Pt" || vars[i] == "Lep1Pt" || vars[i] == "Lep0Eta" || vars[i] == "Lep1Eta" || 
     	      vars[i] == "Jet0Pt" || vars[i] == "Jet1Pt" || vars[i] == "Jet0Eta" || vars[i] == "Jet1Eta" ||
     	      vars[i] == "NBtagJets" || vars[i] == "HT" )
	 cut="_2jets";
      else if( vars[i] == "MET" || vars[i] == "InvMass" || vars[i] == "AbsDelPhiLeps")
         cut="_1btag";


      TString histodata = "H_"+vars[i]+"_ElMu"+cut;
      TString histo	= histodata;//+syst[s];
      cout << "\nPainting " << histo << " ... " << histodata<< endl;

      TH1F* h_ttdil = (TH1F*)ttdil->Get(histo); h_ttdil->Scale(lumi);
      TH1F* h_ttsem = (TH1F*)ttsem->Get(histo); h_ttsem->Scale(lumi); 
      TH1F* h_ww    = (TH1F*)ww   ->Get(histo); h_ww   ->Scale(lumi);	      
      TH1F* h_wz    = (TH1F*)wz   ->Get(histo); h_wz   ->Scale(lumi);	      
      TH1F* h_zz    = (TH1F*)zz   ->Get(histo); h_zz   ->Scale(lumi);	      
      TH1F* h_tW    = (TH1F*)tW   ->Get(histo); h_tW   ->Scale(lumi);	      
      TH1F* h_tbarW = (TH1F*)tbarW->Get(histo); h_tbarW->Scale(lumi); 
      TH1F* h_DY    = (TH1F*)DY   ->Get(histo); h_DY   ->Scale(lumi);	      
      TH1F* h_DYlow = (TH1F*)DYlow->Get(histo); h_DYlow->Scale(lumi);	      
      TH1F* h_wjets = (TH1F*)wjets->Get(histo); h_wjets->Scale(lumi); 
      TH1F* h_data  = (TH1F*)data ->Get(histodata); 
      
      h_tW   ->Add(h_tbarW);
      h_ww   ->Add(h_wz)   ; h_ww->Add(h_zz);
      h_ttsem->Add(h_wjets);
      h_DY   ->Add(h_DYlow);

      // normalize nonWZ to DD prediction 
      // if     (cut=="_dilepton") h_ttsem->Scale(8.83/h_ttsem->Integral());
      // else if(cut=="_2jets"   ) h_ttsem->Scale(8.83/h_ttsem->Integral());
      // else if(cut=="_1btag"   ) h_ttsem->Scale(8.83/h_ttsem->Integral());
      
      int rebin = 1;
      if     (vars[i] == "Lep0Pt")	 {rebin =  50; xtitle = "Muon p_{T} (GeV)"           ;}
      else if(vars[i] == "Lep0Eta")	 {rebin =   2; xtitle = "Muon |#eta|"                ;}
      else if(vars[i] == "Lep1Pt")	 {rebin =  50; xtitle = "Electron p_{T} (GeV)"       ;}
      else if(vars[i] == "Lep1Eta")	 {rebin =   2; xtitle = "Electron |#eta|"	     ;}
      else if(vars[i] == "Jet0Pt")	 {rebin = 100; xtitle = "Leading jet p_{T} (GeV)"    ;}
      else if(vars[i] == "Jet0Eta")	 {rebin =   2; xtitle = "Leading jet |#eta|"	     ;}
      else if(vars[i] == "Jet1Pt")	 {rebin = 100; xtitle = "Sub-Leading jet p_{T} (GeV)";}
      else if(vars[i] == "Jet1Eta")	 {rebin =   2; xtitle = "Sub-Leading jet |#eta|"     ;}
      else if(vars[i] == "InvMass")	 {rebin =  10; xtitle = "m_{e#mu} (GeV)"	     ;}
      else if(vars[i] == "HT")  	 {rebin = 100; xtitle = "H_{T} (GeV)"  	             ;}
      else if(vars[i] == "AbsDelPhiLeps"){rebin =   1; xtitle = "|#Delta#phi(e, #mu)| (rad) / #pi";}
      else if(vars[i] == "NJets")	 {rebin =   1; xtitle = "Number of jets"                  ;}
      else if(vars[i] == "NBtagJets")	 {rebin =   1; xtitle = "Number of b jets"                ;}
      else if(vars[i] == "MET") 	 {rebin = 100; xtitle = "Missing Transverse Energy (GeV)" ;}
       
      if(vars[i] == "CSVTag" ||vars[i] == "Lep0Iso" || vars[i] == "Lep1Iso" )  plot->SetLogy();
      else plot->SetLogy(0);
       
      h_ttdil->SetLineColor(633);     h_ttdil->SetFillColor(633);   calcOverflow(h_ttdil);   h_ttdil->Rebin(rebin); 
      h_DY   ->SetLineColor(852);     h_DY   ->SetFillColor(852);   calcOverflow(h_DY	);   h_DY   ->Rebin(rebin);
      h_tW   ->SetLineColor(616);     h_tW   ->SetFillColor(616);   calcOverflow(h_tW	);   h_tW   ->Rebin(rebin);  
      h_ww   ->SetLineColor(390);     h_ww   ->SetFillColor(390);   calcOverflow(h_ww	);   h_ww   ->Rebin(rebin);  
      h_ttsem->SetLineColor(413);     h_ttsem->SetFillColor(413);   calcOverflow(h_ttsem);   h_ttsem->Rebin(rebin);
      
      //h_data->SetMarkerSize(0);
      //h_data->SetLineWidth(0);
      h_data->SetMarkerColor(kBlack);
      h_data->SetLineColor(kBlack);
      h_data->SetMarkerStyle(20);
      h_data ->Rebin(rebin);
      calcOverflow(h_data);
 
      cout << setprecision (2) << fixed;
      cout << "Entries:" << endl;
      cout << " ttbar: " << h_ttdil->Integral() << endl;
      cout << " tW   : " << h_tW   ->Integral() << endl;
      cout << " VV   : " << h_ww   ->Integral() << endl;
      cout << " DY   : " << h_DY   ->Integral() << endl;
      cout << " nonW : " << h_ttsem->Integral() << endl;
      //cout << " total: " << hs->GetStack()->Last()	 ->Integral() << endl;
      cout << " data : " << h_data ->Integral() << endl;      
      //cout << " data : " << h_data ->GetEntries() << endl;  

      THStack *hs = new THStack("hs","Stacked 1D histograms");
      hs->SetTitle("");
      
      hs->Add(h_DY);
      hs->Add(h_tW);
      hs->Add(h_ww);
      hs->Add(h_ttsem);
      hs->Add(h_ttdil);


      double max = -999.;
      if( hs->GetMaximum()>h_data->GetMaximum() ) max = hs->GetMaximum();
      else					  max = h_data->GetMaximum();
      hs->SetMaximum(max*1.35);
      
      if(vars[i] == "AbsDelPhiLeps") hs->SetMaximum(max*1.45);
      else if(vars[i] == "Jet0Eta" || vars[i] == "Jet1Eta") hs->SetMaximum(max*1.48);
      else if(vars[i] == "NBtagJets") hs->SetMaximum(max*1.65);
      
      plot->cd();


      /////////////////////////
      // For systematic unc. //        
      /////////////////////////
      bool drawUnc = false;
      bool printInfo = false;
      
      TH1F* my2hf1 = (TH1F*)hs->GetStack()->Last()->Clone();
      const int nbins = my2hf1->GetNbinsX();

      double bkgYield = 0., ttYield 0., bkgStat = 0., bkgSyst = 0., ttStat = 0., ttSyst = 0., bkgUnc = 0., ttUnc = 0.;

      // yields for signal and background
      ttYield = h_ttdil->Integral();
      bkgYield = my2hf1->Integral() - h_ttdil->Integral();
      
      // total uncertainties for signal and background (from yields table 9 in AN-15-022)
      if (cut == "_dilepton"){
     	 bkgStat =  2.94;
     	 bkgSyst = 14.40;  // this is too small, not using DD here... !!!
     	 ttStat  =  0.70;
     	 ttSyst  = 21.35;
      }else if (cut == "_2jets"){
     	 bkgStat =  3.97;
     	 bkgSyst =  4.10;
     	 ttStat  =  0.60;
     	 ttSyst  = 15.94;
      }
      
      // calculate total relative (%) uncertainties
      bkgUnc = sqrt(bkgStat*bkgStat + bkgSyst*bkgSyst)*100/bkgYield;
      ttUnc  = sqrt(ttStat*ttStat   + ttSyst*ttSyst  )*100/ttYield;
      cout << ttYield << " +- " << ttUnc << "%, bkg: " << bkgYield << " +- " << bkgUnc << "%" << endl;
      
      for(int b=1; b<nbins+1; b++){
     	 double binttUnc = 0., binBkgUnc = 0., totSyst = 0.;
     	 
   	 // calculate total uncertainties for signal and bkg in each bin
     	 binttUnc  = ttUnc  * h_ttdil->GetBinContent(b) 			       / 100.;
     	 binBkgUnc = bkgUnc * ( my2hf1->GetBinContent(b) - h_ttdil->GetBinContent(b) ) / 100.;  	 
     	 
   	 // calculate total uncertainties for MC in each bin
     	 totSyst = sqrt(binttUnc*binttUnc + binBkgUnc*binBkgUnc);
     	 
   	 // replace the error by the total uncertainty
     	 my2hf1->SetBinError(b, totSyst);
     	 
     	 if(printInfo)
     	    cout << b << ": tt = "  << h_ttdil->GetBinContent(b) << " +- " << binttUnc 
     		      << ", bkg = " << my2hf1 ->GetBinContent(b) - h_ttdil->GetBinContent(b) << " +- " << binBkgUnc 
     		      << ", tot = " << my2hf1 ->GetBinContent(b) << " +- " << my2hf1 ->GetBinError(b)<< endl;
      }

      TGraphErrors *systUnc = new TGraphErrors(my2hf1);
      systUnc->SetFillStyle(3444);
      systUnc->SetFillColor(kBlack);
      systUnc->SetLineColor(kWhite);
    
      ///////////////////////////////
      // Done with systematic unc. //	      
      ///////////////////////////////
      
      
      hs->Draw("hist");
      h_data ->Draw("psameE1X0");	       
      if(drawUnc) systUnc->Draw("sameE2");
      //plot->RedrawAxis("same");

      hs->GetYaxis()->SetTitle("Number of events");
      hs->GetYaxis()->SetTitleSize(0.08);
      hs->GetYaxis()->SetTitleOffset(0.82);
      hs->GetYaxis()->SetLabelSize(0.06);
      hs->GetYaxis()->SetNdivisions(505);
      hs->GetXaxis()->SetLabelSize(0.0);
      hs->GetYaxis()->SetTickLength(0.02);    
      hs->GetXaxis()->SetTickLength(0.03);    
      if(vars[i] == "NJets") {  // this is done to have one tick in the center of the bin
     	 hs->GetXaxis()->SetBinLabel(1,"0");  
      }



      // Draw Legend
      if(vars[i] == "AbsDelPhiLeps"){
     	 TLegend *tleg2, *tleg3;
   	 tleg2 = new TLegend(0.14,0.56,0.36,0.84);
   	 tleg3 = new TLegend(0.42,0.56,0.64,0.84);
   	 tleg2->SetTextSize(0.08);
   	 tleg3->SetTextSize(0.08);
   	 tleg2->SetBorderSize(0);
   	 tleg3->SetBorderSize(0);
   	 tleg2->SetFillColor(10);
   	 tleg3->SetFillColor(10);

         tleg2->AddEntry(h_data , " Data"  	       , "p");
         tleg2->AddEntry(h_ttdil, " t#bar{t}"	       , "f");
         tleg2->AddEntry(h_ttsem, " Non W/Z"	       , "f");
         tleg3->AddEntry(h_ww   , " VV"		       , "f");
         tleg3->AddEntry(h_tW   , " tW"		       , "f");
         tleg3->AddEntry(h_DY   , " Z/#gamma* #rightarrow e^{#pm}#mu^{#mp}", "f");
         tleg2->Draw("same");
         tleg3->Draw("same");
      }else{  
         TLegend *tleg;
         if (vars[i] == "HT") tleg = new TLegend(0.62,0.42,0.94,0.92);
         else                 tleg = new TLegend(0.66,0.42,0.98,0.92);
         tleg->SetTextSize(0.08);
         tleg->SetBorderSize(0);
         tleg->SetFillColor(10);

         tleg->AddEntry(h_data , " Data"		      , "p");
         tleg->AddEntry(h_ttdil, " t#bar{t}"	      , "f");
         tleg->AddEntry(h_ttsem, " Non W/Z"	      , "f");
         tleg->AddEntry(h_ww   , " VV"		      , "f");
         tleg->AddEntry(h_tW   , " tW"		      , "f");
         tleg->AddEntry(h_DY   , " Z/#gamma* #rightarrow e^{#pm}#mu^{#mp}", "f");
         tleg->Draw("same");
      }
      
      //hs->Draw("histsame");		       
      h_data ->Draw("psameE1X0");	       

      // Draw Luminosity
      TString titlelabel = Form("%4.1f fb^{-1} (13 TeV) ",lumi/1000.);
      TLatex *title  = new TLatex(-20.,50.,titlelabel);
      title->SetNDC();
      title->SetTextAlign(12);
      title->SetX(0.74);
      title->SetY(0.965);
      title->SetTextFont(42);
      title->SetTextSize(0.06);
      title->SetTextSizePixels(22);
      title->Draw("SAME");
  
      // Draw CMS Preliminary
      TLatex *chtitle;
      chtitle  = new TLatex(-20.,50., "CMS"); 
      chtitle->SetNDC();
      chtitle->SetTextAlign(12);
      chtitle->SetX(0.14);
      if(vars[i] == "AbsDelPhiLeps") chtitle->SetX(0.8);
      chtitle->SetY(0.885);
      chtitle->SetTextFont(61);
      chtitle->SetTextSize(0.08);
      chtitle->SetTextSizePixels(22);
      chtitle->Draw("SAME");

      TLatex *chtitle2;
      chtitle2  = new TLatex(-20.,50., "Preliminary"); //+" "+syst[s]);
      chtitle2->SetNDC();
      chtitle2->SetTextAlign(12);
      chtitle2->SetX(0.14);
      if(vars[i] == "AbsDelPhiLeps") chtitle2->SetX(0.8);
      chtitle2->SetY(0.81);
      chtitle2->SetTextFont(52);
      chtitle2->SetTextSize(0.055);
      chtitle2->SetTextSizePixels(22);
      chtitle2->Draw("SAME");


      // Draw a vertical line on top of the vertical right axis
      TLine t1(plot->GetUxmax(), plot->GetUymin(), plot->GetUxmax(), plot->GetUymax());
      t1.Draw();


      // Draw text for selection
      TString chlabel = Form("");
      if(cut == "_2jets")	  chlabel = Form(" e^{#pm}#mu^{#mp} + #geq 2 jets");
      else if(cut == "_1btag")    chlabel = Form(" e^{#pm}#mu^{#mp} + #geq 2 jets + #geq 1 b tag");
      else if(cut == "_dilepton") chlabel = Form(" e^{#pm}#mu^{#mp}");
     	 
      TLatex *chtitle2 = new TLatex(-20.,50., '*'+chlabel); // the star is needed to have the signs at the same height...
      //TLatex *chtitle2 = new TLatex(-20.,50., chlabel); 
      chtitle2->SetNDC();
      chtitle2->SetX(0.31);
      chtitle2->SetY(0.86);
      if(cut == "_1btag") chtitle2->SetX(0.24);
      chtitle2->SetTextFont(42);
      chtitle2->SetTextSize(0.07);
      chtitle2->SetTextSizePixels(22);
      chtitle2->Draw("SAME");

      // Draw white bullet to cover the star 
      TLatex *chtitle3 = new TLatex(-20.,50., "#bullet");
      chtitle3->SetTextColor(0);
      chtitle3->SetNDC();
      chtitle3->SetTextAlign(12);
      chtitle3->SetX(0.3);
      chtitle3->SetY(0.88);
      if(cut == "_1btag") chtitle3->SetX(0.235);
      chtitle3->SetTextFont(42);
      chtitle3->SetTextSize(0.13);
      chtitle3->SetTextSizePixels(22);
      chtitle3->Draw("SAME"); 	 


      ratio->cd();

      TH1F* myhf1 = (TH1F*)hs->GetStack()->Last();
      
      TH1F *H_Ratio = (TH1F*)h_data->Clone();
      H_Ratio->Divide(myhf1);

      TH1F *H_RatioSyst = (TH1F*)my2hf1->Clone();
      //H_RatioSyst->Sumw2();
      H_RatioSyst->Divide(myhf1);
      H_RatioSyst->SetFillStyle(3444);
      H_RatioSyst->SetFillColor(kBlack);

      // Set the error to MC(up, down)/MC
      for(int b=1; b<H_Ratio->GetNbinsX()+1; b++){
     	 if(my2hf1 ->GetBinContent(b)!=0 ){
     	    H_RatioSyst->SetBinError(b, (my2hf1 ->GetBinContent(b) + my2hf1->GetBinError(b))/my2hf1 ->GetBinContent(b) - 1.);
     	 }
      }

      for(int b=1; b<H_Ratio->GetNbinsX()+1; b++){
     	 if(my2hf1 ->GetBinContent(b)!=0 && printInfo){
     	    cout << "\n" << b << endl;
     	    cout << " tot: "  << my2hf1 ->GetBinContent(b) << " +- " << my2hf1 ->GetBinError(b)<< endl; 
     	    cout << " data: " << h_data ->GetBinContent(b)<< endl; 
     	    cout << " MCup/MC: " << (my2hf1 ->GetBinContent(b) + my2hf1->GetBinError(b))/my2hf1 ->GetBinContent(b)<< ", "<< endl; 
     	    cout << " MCdo/MC: " << (my2hf1 ->GetBinContent(b) - my2hf1->GetBinError(b))/my2hf1 ->GetBinContent(b)<< ", "<< endl; 
     	    cout << " data/tot: " << h_data ->GetBinContent(b)/my2hf1 ->GetBinContent(b) << " vs ratio: " << H_Ratio->GetBinContent(b) << endl; 
     	    cout << " Ratio: " << H_RatioSyst->GetBinContent(b) << endl;
     	    cout << " Ratio err: " << H_RatioSyst->GetBinError(b) << endl;
     	 }
      }

      H_Ratio->SetTitle("");
      H_Ratio->GetYaxis()->SetTitle("Data/MC");
      H_Ratio->GetYaxis()->CenterTitle();
      H_Ratio->GetYaxis()->SetTitleOffset(0.25);
      H_Ratio->GetYaxis()->SetTitleSize(0.17);
      H_Ratio->GetYaxis()->SetLabelSize(0.1);
      H_Ratio->GetYaxis()->SetNdivisions(202);
      H_Ratio->GetXaxis()->SetTitleOffset(0.8);
      H_Ratio->GetXaxis()->SetLabelSize(0.14);
      H_Ratio->GetXaxis()->SetTitleSize(0.2);
      H_Ratio->GetXaxis()->SetTitle(xtitle);	 
      H_Ratio->GetXaxis()->SetTickLength(0.05);    
      
      if(vars[i] == "NJets") {
   	 //H_Ratio->GetXaxis()->SetLabelOffset(0.012);
   	 H_Ratio->GetXaxis()->SetLabelSize(0.22);
     	 H_Ratio->GetXaxis()->SetBinLabel(1,"0");    
     	 H_Ratio->GetXaxis()->SetBinLabel(2,"1");    
     	 H_Ratio->GetXaxis()->SetBinLabel(3,"2");    
     	 H_Ratio->GetXaxis()->SetBinLabel(4,"3");    
     	 H_Ratio->GetXaxis()->SetBinLabel(5,"4");    
     	 H_Ratio->GetXaxis()->SetBinLabel(6,"5");    
     	 H_Ratio->GetXaxis()->SetBinLabel(7,"6");    
     	 H_Ratio->GetXaxis()->SetBinLabel(8,"#geq 7");    
      }
      
      if(vars[i] == "NBtagJets") {
   	 //H_Ratio->GetXaxis()->SetLabelOffset(0.012);
   	 H_Ratio->GetXaxis()->SetLabelSize(0.22);
     	 H_Ratio->GetXaxis()->SetBinLabel(1,"0");    
     	 H_Ratio->GetXaxis()->SetBinLabel(2,"1");    
     	 H_Ratio->GetXaxis()->SetBinLabel(3,"2");    
     	 H_Ratio->GetXaxis()->SetBinLabel(4,"#geq 3");    
      }
     	
      float range = 0.4;	
      H_Ratio->SetMinimum(1. - range);
      H_Ratio->SetMaximum(1. + range);
      
      H_Ratio->Draw("pE1X0");
      if(drawUnc) H_RatioSyst->Draw("sameE2");
      H_Ratio->Draw("psameE1X0");

      c->Print("figs/stackPlots.ps");
      //c->SaveAs("figs/"+histo+".ps");
      //c->SaveAs("figs/"+histo+".eps");
      c->SaveAs("figs/"+histo+".pdf");
      c->SaveAs("figs/"+histo+".png");
      //c->SaveAs("figs/"+histo+".C");
      //c->SaveAs("figs/"+histo+".root");
   }
   
   c->Print("figs/stackPlots.ps]"); 

}
