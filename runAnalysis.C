void runAnalysis()
{
    // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t local = kTRUE;
   
    
    // since we will compile a class, tell root where to look for headers
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -g");
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    
    Int_t year = 2015;
    TString prod = "LHC15n";
    TString ppass = "pass3";

    Int_t runNmin=0;
    //Int_t runNmax=22;
    Int_t runNmax=1;
    //Int_t runList[30]={244628, 244627, 244626, 244542, 244540, 244531, 244484, 244483, 244482, 244481, 244480, 244456, 244453, 244421, 244416, 244377, 244364, 244359, 244355, 244351, 244343, 244340};
    Int_t runList[1]={244628}; //for test

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
   // AliESDInputHandler *ESDH = new AliESDInputHandler();
    AliVEventHandler* esdH = new AliESDInputHandler();
    ((AliESDInputHandler *) esdH)->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);

    // compile the class (locally)
    gROOT->LoadMacro("AliAnalysisTaskXic.cxx++g");
    // load the addtask macro
    gROOT->LoadMacro("AddXic.C");
    // create an instance of your analysis task
    AliAnalysisTaskXic *task = AddXic();

    if(!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("esdTree");
        // add a few files to the chain (change this so that your local files are added)
       // chain->Add("/Volumes/Transcend/PbPb/tmp/170593/AliESDs.root");
        chain->Add("AliESDs.root");

        // start the analysis locally, reading the events from the tchain
        mgr->StartAnalysis("local", chain);
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *plugin = new AliAnalysisAlien();
        
        // plugin->SetUser("Your alien username here");
        plugin->SetUser("blim");
        
        // Specify working disk
        gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
        
        // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
       // plugin->SetRunMode("test");  // VERY IMPORTANT
        // set the Alien API version
        plugin->SetAPIVersion("V1.1x");

        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        plugin->SetAliROOTVersion("v5-06-15");
        plugin->SetAliPhysicsVersion("v5-06-15-01");
       
        // make sure your source files get copied to grid
        plugin->SetAdditionalLibs("AliAnalysisTaskXic.cxx AliAnalysisTaskXic.h");
        plugin->SetAnalysisSource("AliAnalysisTaskXic.cxx");
   
        // select the input data
        plugin->SetGridDataDir("/alice/data/2015/LHC15n");
        plugin->SetDataPattern("/pass4/*AliESDs.root");
        // MC has no prefix, data has prefix 000
        plugin->SetRunPrefix("000");
        // runnumber
        //plugin->AddRunNumber(167813);
        for (Int_t irun=runNmin;irun<runNmax;irun++){
            plugin->AddRunNumber((Int_t )runList[irun]);
        }
        // number of files per subjob
        plugin->SetSplitMaxInputFileNumber(10);
        
        plugin->SetExecutable("Xic.sh");
        // specify how many seconds your job may take
        plugin->SetTTL(40000);
        plugin->SetJDLName("Xic.jdl");

        plugin->SetOutputToRunNo(kTRUE);
        plugin->SetKeepLogs(kTRUE);
               plugin->SetMaxMergeStages(1);
        plugin->SetMergeViaJDL(kTRUE);

        // define the output folders
        plugin->SetGridWorkingDir("Xi_c_Test");
        plugin->SetGridOutputDir("20170513_01");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(plugin);
        // and launch the analysis
        plugin->SetRunMode("test");
        // merging: run with "terminate" to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        
        mgr->StartAnalysis("grid");
        
    }
}
