Download the code
====

Download all the code from github.

    git clone https://github.com/GonzalezFJR/TOP13TeV
    cd TOP13TeV



Now you can set the enviroment use one analyzer: TreeAnalysisTop (top x-section at 13 TeV), TOP5TeVAnalyzer (top x-section at 5 TeV) ot StopAnalyzer (Stop dileptonic OS).


Set the enviroment
====

Load Root 6, start PAF and PoD.

If you run in Oviedo, load Root 6:

    source /cms/slc6_amd64_gcc493/external/gcc/4.9.3/etc/profile.d/init.sh
    source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/init.sh
    source /cms/slc6_amd64_gcc493/external/python/2.7.6/etc/profile.d/dependencies-setup.sh
    source /opt/root6/bin/thisroot.sh
    
Start Pod:

    source /opt/PoD/PoD_env.sh
  
Load PAF:

    source /opt/PAF/PAF_setup.sh
  
  
Top cross section at 13 TeV
====

You have to choose the name of the sample, number of slots, the gen selection (0 = dileptonic, 1 = other, 4 = fiducial...), if systematic studies are done (false by default) and the number of events (all in the sample by default).

    root -l -b -q 'RunTree_ReReco.C(sample, nSlots, dileptonic?, doSyst?, nEvents)'

As an example:

    root -l -b -q 'RunTree_ReReco.C("TTJets", 1, 0, true, 1000)'

Or run all the samples (Root 6, PAF and PoD are loaded authomatically):

    source runAll_TopMC.sh


If you'd like to run the analysis on a local tree you can directly write the 'RunTree_ReReco.C' file and change the local path and name of the sample repectively:

    TString localpath = "/your/local/path/";
    TString sample    = "yousample.root";

Then execute the program using for the name of the sample "TestHeppy":

    root -l -b -q 'RunTree_ReReco.C("TestHeppy", 1)'
    
To extract the cross section, get the systematics and draw the plots you can use the code inside the folder 'TopCode'.

    
Running the top cross section at 5 TeV analysis
====

This code is a simplified version of the analysis at 13 TeV. Follow the previous instructions but changing the analyzer and using the files 'run5TeVAll_TopMC.sh' and 'RunTree_ReReco5TeV.C'. To extract the cross section, yields and systematic uncertainties use the code in the folder 'Plotter5TeV'.


Running the Stop Analysis
====

Use the macro 'RunStopAnalysis.C'. To analyze one stop sample you must specify the masses of stop and neutralino and the weight for normalization. Example:

    root -l -b -q 'RunStopAnalysis.C("T2tt_150to175LSP1to100", 1, true, 0, true, 150, 50, 0.0134)'
    root -l -b -q 'RunStopAnalysis.C("TTbar_Powheg", 10)'
    
To analyze all masses in one file or just evey point use the RunT2ttSamples.py script as follows:

    python RunT2ttSamples.py T2tt_150to175LSP1to100
    python RunT2ttSamples.py

    
To draw the plots use the code in the 'StopPlotter' folder. Change the path in SetPlotter.C and the "plotfolder" path defined in DrawPlots.C. To create plots for a given stop signal execute the T2ttStackPlots function with the required values of: variable to plot ("MET", "MT2", "InvMass", "NJets"...), channel ("ElMu", "Muon", "All"), level ("1btag", "MET", "2jets", "DYVetop"...), StopMass and NeutralinoMass. 

    root -l
    .L DrawPlots.C
    T2ttStackPlots(var, chan, level, StopMass, NeutralinoMass)


Upload dataset Manager
====
    
If there was a new release of Dataset Manager you could download it from gitlab.

    git clone https://gitlab.cern.ch/IFCA-UO-CMS/Utils.git
    mv Utils/DatasetManager/ TOP13TeV

