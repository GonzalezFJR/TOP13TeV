from ROOT import *
import os
import sys

inputpath = "/pool/ciencias/TreesDR74X/heppyTrees/v2/";
nSamples = 9;
samplename = [
"T2tt_150to175LSP1to100", "T2tt_250LSP1to175", "T2tt_275LSP75to200",
"T2tt_200LSP1to125", "T2tt_300to375LSP1to300",     "T2tt_400to475LSP1to400",
"T2tt_225LSP25to150", "T2tt_500to550LSP1to475", "T2tt_600to950LSP1to450"
];

NeutralinoMass = []; StopMass = []; Events = [];
Counts = TH3D();

mstop = 0.; mchi = 0.; nEvents = 0.;
def GetValues(sample):
	inputfile = TFile.Open(inputpath + "Tree_" + sample + "_0.root");
	inputfile.GetObject("CountSMS", Counts);
	val = 0;
	for i in range(Counts.GetNbinsX()):
		for j in range(Counts.GetNbinsY()):
			val = Counts.GetBinContent(i,j,1);
			if (val != 0):
				mstop = Counts.GetXaxis().GetBinCenter(i);
				mchi  = Counts.GetYaxis().GetBinCenter(j); 
				nEvents = Counts.GetBinContent(Counts.FindBin(mstop, mchi, 0));
				if(mchi == 0):
					mchi = 1.;
				NeutralinoMass.append(mchi);
				StopMass.append(mstop);
				Events.append(nEvents); 

def GetAllValues():
	for k in range(len(samplename)):
		GetValues(samplename[k])

def PrintInfo():
	for i in range(len(StopMass)):
		print "Mstop = ", StopMass[i], ", Mchi = ", NeutralinoMass[i], ", nEvents = ", Events[i]

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

def Analyze(sample = "All"):
	if (sample == "All"):
		for k in range(len(samplename)):
			print ">>> Progress: sample ", k+1, "/", len(samplename)
			Analyze(samplename[k])
	else:
		#StopMass = []; NeutralinoMass = []; Events = [];
		print ">>> Sample: ", sample
		print ">>> Getting masses of Stop and Neutralino... "
		del NeutralinoMass[:]; del StopMass[:]; del Events[:];
		GetValues(sample)
		print ">>> Done! \n"
		for m in range(len(StopMass)):
			xsec = getxsec(StopMass[m])
			print ">>> Analyzing sample: ", sample 
			print ">>> MStop   = ", StopMass[m]
			print ">>> MLSP    = ", NeutralinoMass[m]		
			print ">>> nEvents = ", Events[m]
			print ">>> xsec    = ", xsec
			print ">>> Progress: point", m+1, "/", len(StopMass), "..."
			print " "
			weight = xsec/Events[m]
			os.system("root -l -q -b 'RunTree_ReReco.C(\"" + sample + "\", 4, false, 0, true, " + str(StopMass[m]) + ", " + str(NeutralinoMass[m]) + ", " + str(weight) + ")'")

if __name__ == "__main__":
  if   len(sys.argv) == 1:
    Analyze()
  elif len(sys.argv) == 2:
    Analyze(str(sys.argv[1]))
