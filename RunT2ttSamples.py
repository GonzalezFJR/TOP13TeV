"""
 @ Author; Juan R Gonzalez Fernandez
   # Run all samples: 
      >> python RunT2ttSamples
   # Run all mases in one signle sample: 
      >> python RunT2ttSamples "T2tt_150to250"
   # Run all samples sending jobs:
      >> python RunT2ttSamples 0 0 0 0 0 0 
   # Run one single point
      >> python RunT2ttSamples 800 50

"""

from ROOT import *
import os
import sys

inputpath = "/pool/ciencias/HeppyTreesDR80X/v1/";
nSamples = 4;
samplename = [
"T2tt_150to250", 
"T2tt_250to350", 
"T2tt_350to400",
"T2tt_400to1200"
];

NeutralinoMass = []; StopMass = []; Events = [];
Counts = TH3D();
jobscounter = 0;

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
  else: 
    print "No Cross Section for that mass!! Extrapolating..."
    v0 = getxsec(StopMass - StopMass%25)
    vf = getxsec(StopMass - StopMass%25 + 25)
    x  = float(StopMass%25)/25
    return v0 + (vf-v0)*x

def Analyze(sample = "All", sendjobs = False, mstopo = 0, mlsp = 0):
	if (sample == "All"):
		for k in range(len(samplename)):
			print ">>> Progress: sample ", k+1, "/", len(samplename)
			Analyze(samplename[k], False, mstopo, mlsp)
	elif (sample == "jobs"):
		for k in range(len(samplename)):
			print ">>> Sending jobs for processing sample ", k+1, "/", len(samplename)
			Analyze(samplename[k], True) 
	else:
		#StopMass = []; NeutralinoMass = []; Events = [];
		if (sendjobs):
			del NeutralinoMass[:]; del StopMass[:]; del Events[:];
			GetValues(sample)
			for m in range(len(StopMass)):
				#if((StopMass[m] - NeutralinoMass[m]) > 100): continue
				if(StopMass[m] != 250 or NeutralinoMass[m] != 50): continue
				#if(StopMass[m] < 750): continue
				#if(StopMass[m] != 150 or NeutralinoMass[m] != 1): continue
				xsec = getxsec(StopMass[m])
				weight = xsec/Events[m]
        #jobdir = "./temp/jobs/tempjobs" + str(StopMass[m]) + str(NeutralinoMass[m]) + "/"
        #os.system("mkdir " + jobdir)
				#tempfile = jobdir + "tempjobs" + str(StopMass[m]) + str(NeutralinoMass[m]) + ".txt"
				tempfile = "tempjobs" + str(StopMass[m]) + str(NeutralinoMass[m]) + ".txt"
				runroot = "root -l -q -b " + "\'RunStopAnalysis.C(\\\"" + sample + "\\\", 1, true, 0, true, " + str(StopMass[m]) + ", " + str(NeutralinoMass[m]) + ", " + str(weight) + ")\'"

				l1 = "#!/bin/bash"
				l2 = "#PBS -e log/error"  + str(StopMass[m]) + "_" + str(NeutralinoMass[m]) + ".log"
				l3 = "#PBS -o log/output" + str(StopMass[m]) + "_" + str(NeutralinoMass[m]) + ".log"
				l4 = "MAINDIR=/nfs/fanae/user/juanr/Stop_76X/"
				l5 = "WDIR=temp/jobs/job" + str(StopMass[m]) + "_" + str(NeutralinoMass[m])

				os.system("echo \"" + l1 + "\" > " + tempfile)
				os.system("echo \"" + l2 + "\" >> " + tempfile)
				os.system("echo \"" + l3 + "\" >> " + tempfile)
				os.system("echo \"" + l4 + "\" >> " + tempfile)
				os.system("echo \"" + l5 + "\" >> " + tempfile)


				os.system("cat job_template.txt >> " + tempfile)
				os.system("echo \"" + runroot + "\" >> " + tempfile)
				os.system("qsub " + tempfile)	
				print "Submited: " + str(StopMass[m]) + " " + str(NeutralinoMass[m])
				os.remove(tempfile)
		else:
			print ">>> Sample: ", sample
			print ">>> Getting masses of Stop and Neutralino... "
			del NeutralinoMass[:]; del StopMass[:]; del Events[:];
			GetValues(sample)
			#PrintInfo()
			print ">>> Done! \n"
			for m in range(len(StopMass)):
				#if(StopMass[m] != 600 or NeutralinoMass[m] != 1): continue;
				print "xxxxxxxxxxxxx Checking sample, mstopo = ", mstopo, " mlsp = ", mlsp
				if(mstopo != 0 and StopMass[m] != mstopo): continue;
				if(mlsp != 0 and NeutralinoMass[m] != mlsp): continue;
				xsec = getxsec(StopMass[m])
				print ">>> Analyzing sample: ", sample 
				print ">>> MStop   = ", StopMass[m]
				print ">>> MLSP    = ", NeutralinoMass[m]		
				print ">>> nEvents = ", Events[m]
				print ">>> xsec    = ", xsec
				print ">>> Progress: point", m+1, "/", len(StopMass), "..."
				print " "
				weight = xsec/Events[m]
				os.system("root -l -q -b 'RunStopAnalysis.C(\"" + sample + "\", 1, true, 0, true, " + str(StopMass[m]) + ", " + str(NeutralinoMass[m]) + ", " + str(weight) + ")'")

if __name__ == "__main__":
  if   len(sys.argv) == 1:
    print " ___ Analyzing all samples ___"
    Analyze()
  elif len(sys.argv) == 2:
    sample = str(sys.argv[1])
    print " ___ Analyzing one sample: ", sample, " ___"
    Analyze(sample)
  elif len(sys.argv) == 3:
    mstopo = float(sys.argv[1]); mlsp = float(sys.argv[2]);
    print " ___ Analyzing one point [", mstopo, ", ", mlsp,"]  ___"
    Analyze("All", False, mstopo, mlsp)  
  else:
    print " ___ Analyzing all samples, sending jobs ___"
    Analyze(str(sys.argv[1]), True)

