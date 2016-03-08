from ROOT import *
import sys
import os

pathistos = '/mnt_pool/fanae105/user/juanr/TOP13TeV/fiducial'
patheppy = '/pool/ciencias/heppyTreesDR76X/v2/' #2/noSkim/'

def get(sample):
  s = 0.; s1 = 0.; s2 = 0.;
  for tree in os.listdir(patheppy):
    if ("Herwig" in tree and "Herwig" not in sample): continue
    if ("aMCatNLO" in tree and "aMCatNLO" not in sample): continue
    if ("ext" in tree and "ext" not in sample): continue
    if sample in tree:
      #print "Opening ", patheppy + tree, "..."
      f = TFile.Open(patheppy + tree)
      h = TH1D(); t = TTree();
      f.GetObject('Count', h)
      f.GetObject('tree', t)
      s1 += t.GetEntries();
      s += h.GetBinContent(1)
      if 'aMCatNLO' in sample:
        h2 = TH1D()
        f.GetObject('SumGenWeights', h2)
        s2 += h2.GetBinContent(1)
  print 'Sample: ', sample
  print 'Count = ', s
  print 'nEvents = ', s1
  print 'SumGenWeights = ', s2
  print ''

get('TTbar_Powheg')
get('TTJets_aMCatNLO')
get('TTbar_Powheg_Herwig')
