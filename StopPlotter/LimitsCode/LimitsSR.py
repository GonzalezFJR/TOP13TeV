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
treespath = "/nfs/fanae/user/juanr/StopTOP/mar19/Susy/" 

def getLumi(dataname):
  f = TFile.Open(path + "/rootfiles/" + dataname + ".root")
  lumihist = TH1F(); f.GetObject("thelumi", lumihist);
  thelumi = lumihist.Integral()
  return thelumi

def getAsymptoticLimit(mStop, mLsp):
  datacard = "datacard_SR_MT2_" + str(mStop) + "_" + str(mLsp) + ".txt"
  dataname = "MT2ADH_"+ str(mStop) + "_" + str(mLsp) + "_ElMu" #+".root"
  thelumi = getLumi(dataname)
  filename = "limitSR_" + str(thelumi/1000) + "invfb_temp.txt"
  #if "temp.txt" in os.listdir(path):
  #  os.remove(path + "temp.txt");
  tempname  = "temp"  + str(mStop) + str(mLsp) + ".txt"
  tempname2 = "temp2" + str(mStop) + str(mLsp) + ".txt"
  os.system("echo \" Lumi = " + str(thelumi) +  "\" > " + path + tempname2)
  os.system("echo \" mStop = " + str(mStop) +  "\" >> " + path + tempname2)
  os.system("echo \" mLsp = " + str(mLsp) +  "\" >> " + path + tempname2)
  os.system("combine -M Asymptotic " + path + "datacards/" + datacard + " >> " + path + tempname) #"temp.txt") 
  os.system("cat " + path + tempname  + " >> " + path + tempname2)
  os.system("cat " + path + tempname2 + " >> " + path + "datacards/" + filename)
  os.remove(path + tempname2)
  os.remove(path + tempname )

def CalculateLimit(mStop, mLsp, lumi, useSystfac):
  print "\n#########################################"
  print "## mStop = ", mStop
  print "## mLsp  = ", mLsp
  os.system("root -l -b -q " + macrospath + "\'CreateDatacard.C(" + str(mStop) + ", " + str(mLsp) + ", " + str(lumi) + ", " + str(useSystfac) + ", 0 " + ")\'")
  print "## Datacard Created: datacard_SR_" + str(mStop) + "_" + str(mLsp) + ".txt"
  print "## Calculating asymptotic limit..."
  getAsymptoticLimit(mStop, mLsp)
  print "## Done!! "
  print "#########################################\n"
 

def getAllLimits(thelumi):
  filename = "limits_" + str(thelumi) + "invfb_temp.txt" 
  count = 0; tot = len(os.listdir(treespath))
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

def jobs(thelumi):
  filename = 'limits_' + str(thelumi) + 'invfb.txt' 
  count = 0; tot = len(os.listdir(treespath))
  if filename in os.listdir(path + "datacards/"): os.remove(path + "datacards/" + filename); 
  for tree in os.listdir(treespath):
    if not ("T2tt" in tree): continue
    mStop = int((re.search("mStop(.+?)_" , tree)).group(1)) 
    mLsp = int((re.search("mLsp(.+?).root" , tree)).group(1)) 
    print "-------> ", count, "/", tot, " (", float(count)/tot*100, "%)  ||"
    command = 'python LimitsSR.py  ' + str(mStop) + ' ' + str(mLsp) + " " + str(thelumi*1000)
    count += 1;
    os.system('cat jobstemplate.txt > job.sh')
    os.system('echo ' + '"' + command + '" >> job.sh')
    os.system('qsub job.sh')
  print "\n------------ DONE! ------------"
  print " Results in file: ", path + "datacards/" + filename 
  print "-------------------------------\n"
  


if __name__ == "__main__":
  if   len(sys.argv) == 1:
    print "Usage: python getAllLimits.py [lumi (fb-1)]                 --> Get all limits secuentialy"
    print "       python getAllLimits.py [lumi (fb-1)] 1               --> Send jobs"
    print "       python getAllLimits.py [mStop] [mLsp] [lumi (pb-1)]  --> Calculate one limit"
  elif len(sys.argv) == 2:
    getAllLimits(float(sys.argv[1]))
  elif len(sys.argv) == 3:
    jobs(float(sys.argv[1]))
  elif len(sys.argv) == 4:
    CalculateLimit(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), 0)
