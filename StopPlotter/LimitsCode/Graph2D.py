'''
 Read the datasets and get values of masses for the neutralino, stop and limits

 Author: Juan Rodrigo Gonzalez Fernandez
'''
import re
from ROOT import *
import numpy as np
from array import array

# input amd output dir
inputdir = "/mnt_pool/fanae105/user/juanr/CMSSW_7_4_7_patch1/StopLimits2016/" 
outputdir = inputdir + "plots/" 

def SetupColors():
  num = 5;
  bands = 255;
  stops = [0.00, 0.34, 0.61, 0.84, 1.00];
  red = [0.50, 0.50, 1.00, 1.00, 1.00];
  green = [0.50, 1.00, 1.00, 0.60, 0.50];
  blue = [1.00, 1.00, 0.50, 0.40, 0.50];
  stopsArray = array('d', stops)
  redArray   = array('d', red)
  greenArray = array('d', green)
  blueArray  = array('d', blue)
  fi = TColor.CreateGradientColorTable(num,stopsArray,redArray,greenArray,blueArray,bands);
  colors=[];
  for i in range(bands):
    colors.append(fi+i);
  colorArray = array('i', colors)
  gStyle.SetNumberContours(bands);
  gStyle.SetPalette(bands, colorArray);


SetupColors();

def getTGraph(gr):
  L = gr.GetContourList(1.);
  out = TGraph(); max_points = 0;
  for i in range(L.GetSize()):
    g = L.At(i)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      out = g
      max_points = n_points
  out.SetLineWidth(2);
  return out

def Graph2D(lumi):
  fname = "limitSR_" + str(lumi) + "invfb.txt" #"invfb.txt.bak"
  f = open(inputdir + "/datacards/" + fname, "r")
  text = f.read()
  textlines = text.split("\n")

  nlin = -1
  vlum = []; y2min = []; y1min = []; y1max = []; y2max = []; e = [];
  mStop = []; mLsp = []; xsec = [];

  for line in textlines:
    nlin += 1
    if("mStop" in line):
      mS = float((re.search("mStop = (.*)", line.strip())).group(1))
      mL = float((re.search("mLsp = (.*)", textlines[nlin+1])).group(1))
      exp2sd = float((re.search("r < (.*)" , textlines[nlin+6])).group(1)) 
      exp1sd = float((re.search("r < (.*)" , textlines[nlin+7])).group(1)) 
      expc = float((re.search("r < (.*)" , textlines[nlin+8])).group(1)) 
      exp1su = float((re.search("r < (.*)" , textlines[nlin+9])).group(1)) 
      exp2su = float((re.search("r < (.*)" , textlines[nlin+10])).group(1))
      y2min.append(exp2sd); y1min.append(exp1sd); e.append(expc); y1max.append(exp1su); y2max.append(exp2su);
      mStop.append(mS); mLsp.append(mL); xsec.append(expc*getxsec(mS))

  n = len(mStop)
  print "n = ", n
  xseclim = TGraph2D(); xseclim.Set(n);
  expected = TGraph2D(); expected.Set(n);
  expectedU = TGraph2D(); expectedU.Set(n);
  expectedD = TGraph2D(); expectedD.Set(n);
  expected2 = TGraph2D(); expected2.Set(n);
  expectedU2 = TGraph2D(); expectedU2.Set(n);
  expectedD2 = TGraph2D(); expectedD2.Set(n);
  for k in range(len(mStop)):
    xseclim.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(xsec[k]))
    if ( abs( (mStop[k] - mLsp[k]) - 175) < 5 ):
      expected.SetPoint(k, float(mStop[k]),float(mLsp[k]), 1.2)
      expectedU.SetPoint(k, float(mStop[k]),float(mLsp[k]), 1.4)
      expectedD.SetPoint(k, float(mStop[k]),float(mLsp[k]), 1.1)
    elif (mStop[k] - mLsp[k] < 175 and float(e[k]) < 1.):
      expected.SetPoint(k, float(mStop[k]),float(mLsp[k]), 2.2)
      expectedU.SetPoint(k, float(mStop[k]),float(mLsp[k]), 2.4)
      expectedD.SetPoint(k, float(mStop[k]),float(mLsp[k]), 2.1)
    else:
      expected.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(e[k]))
      expectedU.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(y1max[k]))
      expectedD.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(y1min[k]))
    if ( abs( (mStop[k] - mLsp[k]) - 175) < 5 ):
      expected2.SetPoint(k, float(mStop[k]),float(mLsp[k]), 1.4)
      expectedU2.SetPoint(k, float(mStop[k]),float(mLsp[k]), 1.7)
      expectedD2.SetPoint(k, float(mStop[k]),float(mLsp[k]), 1.2)
    elif (mStop[k] - mLsp[k] > 150 and float(e[k]) < 1.):
      expected2.SetPoint(k, float(mStop[k]),float(mLsp[k]), 2.4)
      expectedU2.SetPoint(k, float(mStop[k]),float(mLsp[k]), 2.7)
      expectedD2.SetPoint(k, float(mStop[k]),float(mLsp[k]), 2.2)
    else:
      expected2.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(e[k]))
      expectedU2.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(y1max[k]))
      expectedD2.SetPoint(k, float(mStop[k]),float(mLsp[k]), float(y1min[k]))

  hist = xseclim.GetHistogram()
  hist2 = expected.GetHistogram()
  histu = expectedU.GetHistogram()
  histd = expectedD.GetHistogram()
  hist22 = expected2.GetHistogram()
  histu2 = expectedU2.GetHistogram()
  histd2 = expectedD2.GetHistogram()
  L = expected.GetContourList(1.);
  LU = expectedU.GetContourList(1.);
  LD = expectedD.GetContourList(1.);
  L2 = expected2.GetContourList(1.);
  LU2 = expectedU2.GetContourList(1.);
  LD2 = expectedD2.GetContourList(1.);
  eGr = TGraph(); euGr = TGraph(); edGr = TGraph(); 
  eGr2 = TGraph(); euGr2 = TGraph(); edGr2 = TGraph(); 

  max_points = 0;
  #print "L.GetSize() = ", L.GetSize()
  for i in range(L.GetSize()):
    g = L.At(0)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      eGr = g
      max_points = n_points

  max_points = 0;
  for i in range(LU.GetSize()):
    g = LU.At(i)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      euGr = g
      max_points = n_points

  max_points = 0;
  for i in range(LD.GetSize()):
    g = LD.At(i)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      edGr = g
      max_points = n_points

  max_points = 0;
  for i in range(L2.GetSize()):
    g = L2.At(i)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      eGr2 = g
      max_points = n_points

  max_points = 0;
  for i in range(LU2.GetSize()):
    g = LU2.At(i)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      euGr2 = g
      max_points = n_points

  max_points = 0;
  for i in range(LD2.GetSize()):
    g = LD2.At(i)
    if(g == nullptr): continue
    n_points = g.GetN()
    if(n_points > max_points):
      edGr2 = g
      max_points = n_points


  eGr.SetLineWidth(2); euGr.SetLineWidth(2); edGr.SetLineWidth(2);
  euGr.SetLineStyle(3); edGr.SetLineStyle(3)
  eGr2.SetLineWidth(2); euGr2.SetLineWidth(2); edGr2.SetLineWidth(2);
  euGr2.SetLineStyle(3); edGr2.SetLineStyle(3)



  lumin = 10;
  E = 13;
  c = TCanvas("c1","c1",600,600);
  p = c1.GetPad(0);
  p.SetPad(0.01, 0.01, 0.95, 0.95)
  p.cd()
  Yaxis = TLatex(0., 0., "95% CL upper limit on cross section [pb]");
  Yaxis.SetTextAngle(90); 
  Yaxis.SetNDC();
  #Yaxis.SetTextAlign(12);
  Yaxis.SetX(0.88);
  Yaxis.SetY(0.15);
  Yaxis.SetTextFont(42);
  Yaxis.SetTextSize(0.042);
  Yaxis.SetTextSizePixels(22);
  texlumi = TLatex(-20, 50, "L = " + str(lumin) + " fb^{-1}, #sqrt{s} = " + str(E) + " TeV")
  texlumi.SetNDC();
  texlumi.SetTextAlign(12);
  texlumi.SetX(0.55);
  texlumi.SetY(0.92);
  #texlumi.SetTextFont(51);
  texlumi.SetTextSize(0.032);
  texlumi.SetTextSizePixels(22);
  texcms = TLatex(0., 0., "CMS Preliminary")
  texcms.SetNDC();
  texcms.SetTextAlign(12);
  texcms.SetX(0.10);
  texcms.SetY(0.92);
  texcms.SetTextFont(61);
  texcms.SetTextSize(0.032);
  texcms.SetTextSizePixels(23);


  #c1.SetGrid();
  c1.SetLogz();
  xseclim.SetMinimum(0.01);
  xseclim.SetMaximum(200);
  xseclim.SetNpx(120);
  xseclim.SetNpy(60);
  #xseclim.SetTitle("95% CL limit on cross section; M_{stop} [GeV]; M_{LSP} [GeV]");
  xseclim.SetTitle("")
  xseclim.Draw("colz")
  
  texlumi.Draw();
  texcms.Draw();
  Yaxis.Draw();
  
  eGr.Draw("c,same")
  euGr.Draw("c,same")
  edGr.Draw("c,same")
  eGr2.Draw("c,same")
  euGr2.Draw("c,same")
  edGr2.Draw("c,same")
  c1.Print(outputdir + "Limits2D.pdf", "pdf");
  c1.Print(outputdir + "Limits2D.png", "png");


def getxsec(StopMass):
  if StopMass == 125: return 574.981;
  elif StopMass == 150: return 249.409;
  elif StopMass == 175: return 121.416;
  elif StopMass == 200: return 64.5085;
  elif StopMass == 225: return 36.3818;
  elif StopMass == 250: return 21.5949;
  elif StopMass == 275: return 13.3231;
  elif StopMass == 300: return 8.51615;
  elif StopMass == 325: return 5.60471;
  elif StopMass == 350: return 3.78661;
  elif StopMass == 375: return 2.61162;
  elif StopMass == 400: return 1.83537;
  elif StopMass == 425: return 1.31169;
  elif StopMass == 450: return 0.948333;
  elif StopMass == 475: return 0.697075;
  elif StopMass == 500: return 0.51848;
  elif StopMass == 525: return 0.390303;
  elif StopMass == 550: return 0.296128;
  elif StopMass == 575: return 0.226118;
  elif StopMass == 600: return 0.174599;
  elif StopMass == 625: return 0.136372;
  elif StopMass == 650: return 0.107045;
  elif StopMass == 675: return 0.0844877;
  elif StopMass == 700: return 0.0670476;
  elif StopMass == 725: return 0.0536438;
  elif StopMass == 750: return 0.0431418;
  elif StopMass == 775: return 0.0348796;
  elif StopMass == 800: return 0.0283338;
  elif StopMass == 825: return 0.0241099;
  elif StopMass == 850: return 0.0189612;
  elif StopMass == 875: return 0.015625;
  elif StopMass == 900: return 0.0128895;
  elif StopMass == 925: return 0.0106631;
  elif StopMass == 950: return 0.00883465;
  elif StopMass == 975: return 0.00735655;
  else: print "No Cross Section for that mass!!"


Graph2D("0010")
