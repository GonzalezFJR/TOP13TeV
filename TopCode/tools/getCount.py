from ROOT import *
import sys
import os

pathistos = '/mnt_pool/fanae105/user/juanr/TOP13TeV/fiducial'
patheppy = '/pool/ciencias/HeppyTreesDR80X/v1/' #2/noSkim/'

def get(sample):
  s = 0.; s1 = 0.; s2 = 0.;
  for tree in os.listdir(patheppy):
    if ("Herwig" in tree and "Herwig" not in sample): continue
    if ("aMCatNLO" in tree and "aMCatNLO" not in sample): continue
    if ("ext" in tree and "ext" not in sample): continue
    if ("TTW" in tree and "TTW" not in sample): continue
    if ("TTZ" in tree and "TTZ" not in sample): continue
    if ("new" in tree and "new" not in sample): continue
    #if ("T2tt_150to250noSkim_9" in tree): continue
    if sample in tree:
      #print "Opening ", patheppy + tree, "..."
      f = TFile.Open(patheppy + tree)
      h = TH1D(); t = TTree();
      f.GetObject('Count', h)
      f.GetObject('tree', t)
      s1 += t.GetEntries();
      s += h.GetBinContent(1)
      #if 'aMCatNLO' in sample:
      h2 = TH1D()
      if(('Double' not in sample) and ('Single' not in sample) and ('MuonEG' not in sample)):
        f.GetObject('SumGenWeights', h2)
        s2 += h2.GetBinContent(1)
  print 'Sample: ', sample
  print 'Count = ', s
  print 'nEvents = ', s1
  print 'SumGenWeights = ', s2
  print ''

#get('TTbar_Powheg_ext')
#get('DYJetsToLL_M50_aMCatNLO')
#get('DYJetsToLL_M10to50_aMCatNLO')
#get('DYJetsToLL_M50_MLM')
#get('DYJetsToLL_M5to50_MLM')
#get('WJetsToLNu_aMCatNLO')
#get('TW')
#get('TbarW')
#get('WW')
#get('WZ')
#get('ZZ')
#get('TTWToLNu')
#get('TTWToQQ')
#get('TTZToLLNuNu')
#get('TTZToQQ')
#get('TTGJets')
#get('TTJets_aMCatNLO')
#get('TTbar_Powheg_Herwig')
#get('MuonEG_Run2016B_PromptReco_v2')
#get('DoubleMuon_Run2016B_PromptReco_v2')
#get('DoubleEG_Run2016B_PromptReco_v2')
#get('DoubleEG_Run2016B_PromptReco_v2')
#get('T2tt_150to250noSkim')
get('T2tt_150to250')
#get('T2tt_250to350')
#get('T2tt_350to400')
#get('T2tt_400to1200')
#get('T2tt_1200to1500')
#get('T2tt_425_325_FS')
#get('T2tt_500_325_FS')
#get('T2tt_850_100_FS')
#get('T2tt_mStop160to210mLSP1to20')
#get('SingleElectron_Run2016B_PromptReco_v2')
#get('SingleElectron_Run2016C_PromptReco_v2')
#get('SingleElectron_Run2016D_PromptReco_v2')
#get('SingleMuon_Run2016B_PromptReco_v2')
#get('SingleMuon_Run2016C_PromptReco_v2')
#get('SingleMuon_Run2016D_PromptReco_v2')
#get('MuonEG_Run2016G_PromptReco_v1')
#get('DoubleMuon_Run2016G_PromptReco_v1')
#get('DoubleEG_Run2016G_PromptReco_v1')
#get('SingleElectron_Run2016G_PromptReco_v1')
#get('SingleMuon_Run2016G_PromptReco_v1')
#get('TTbar_Powheg_new')

