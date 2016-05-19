"""
@ Author: Juan R. Gonzalez

Usage: 
>> 

"""

from ROOT import *
import sys
import os
import subprocess
from numpy import sqrt
from numpy import linspace
import re

# outputpath = "/nfs/fanae/user/juanr/stop/bdt/datacards/";
processes = ["tt", "tW", "DY", "VV", "ttV", "NonWZ"] # Data, signal
path = "/mnt_pool/fanae105/user/juanr/CMSSW_7_4_7_patch1/StopLimits2016/"
macrospath = "/nfs/fanae/user/juanr/StopTOP/Plotter/"
treespath = "/nfs/fanae/user/juanr/StopTOP/feb01/Susy/" 

def getLumi(dataname):
  f = TFile.Open(path + "/rootfiles/" + dataname + ".root")
  lumihist = TH1F(); f.GetObject("thelumi", lumihist);
  thelumi = lumihist.Integral()
  return thelumi

def getAsymptoticLimit(mStop, mLsp, thelumi):
  datacard = "datacard_MT2shape_" + str(mStop) + "_" + str(mLsp) + "_comb.txt"
  filename = "limits_" + str(thelumi/1000) + "invfb.txt"
  if "temp.txt" in os.listdir(path):
    os.remove(path + "temp.txt");
  os.system("echo \" Lumi = " + str(thelumi) +  "\" >> " + path + "datacards/"+ filename)
  os.system("echo \" mStop = " + str(mStop) +  "\" >> " + path + "datacards/"+ filename)
  os.system("echo \" mLsp = " + str(mLsp) +  "\" >> " + path + "datacards/"+ filename)
  os.system("combine -M Asymptotic " + path + "datacards/" + datacard + " >> " + path + "temp.txt") 
  os.system("cat " + path + "temp.txt >> " + path + "datacards/" + filename)

def CalculateLimit(mStop, mLsp, lumi, useSystfac):
  print "\n#########################################"
  print "## mStop = ", mStop
  print "## mLsp  = ", mLsp
  chan = "comb";
  os.system("root -l -b -q " + macrospath + "\'CreateDatacardMT2.C(" + str(mStop) + ", " + str(mLsp) + ",\"" + chan + "\"," + str(lumi) + ", " + str(useSystfac) + ", 0 " + ")\'")
  print "## Datacard Created: datacard_MT2shape_" + str(mStop) + "_" + str(mLsp) + "_comb.txt"
  print "## Calculating asymptotic limit..."
  getAsymptoticLimit(mStop, mLsp, lumi)
  print "## Done!! "
  print "#########################################\n"
 

def getAllLimits(thelumi):
  filename = "limits_" + str(thelumi) + "invfb.txt" 
  count = 0; tot = len(os.listdir(path + "datacards/"))
  if filename in os.listdir(path + "datacards/"): os.remove(path + "datacards/" + filename); 
  for tree in os.listdir(treespath):
    if not ("T2tt" in tree): continue
    mStop = int((re.search("mStop(.+?)_" , tree)).group(1)) 
    mLsp = int((re.search("mLsp(.+?).root" , tree)).group(1)) 
    print "-------> ", count, "/", tot, " (", float(count)/tot*100, "%)  ||"
    CalculateLimit(mStop, mLsp, thelumi*1000, 0)
    count += 1;
  print "\n------------ DONE! ------------"
  print " Results in file: ", path + "datacards/" + filename 
  print "-------------------------------\n"


if __name__ == "__main__":
  if   len(sys.argv) == 1:
    print "Usage: python getAllLimits.py [lumi (fb-1)]"
    print "or 3 args to get one limit: mStop  mLsp  thelumi"
  elif len(sys.argv) == 2:
    getAllLimits(float(sys.argv[1]))
  elif len(sys.argv) == 4:
    CalculateLimit(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), 0)

