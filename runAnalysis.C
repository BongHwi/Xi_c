void LoadLibraries();
void LoadMacros(Bool_t isMC=kFALSE);

void runAnalysis(TString pluginmode = "test", const char *dataset = "data.txt")
{
    Bool_t isMC = kFALSE;

    LoadLibraries();
    
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
    AliVEventHandler* esdH = new AliESDInputHandler();
    ((AliESDInputHandler *) esdH)->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);

    LoadMacros(isMC);

    if(!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(pluginmode=="local") {
	printf("LOCAL MODE");
        TChain* chain = new TChain("esdTree");
        // chain->Add("/home/blim/data/AliESDs.root");
	chain = CreateESDChain(dataset);
        mgr->StartAnalysis("local", chain);
    } else {
	printf("GRID MODE");
        AliAnalysisAlien *plugin = new AliAnalysisAlien();
        
        plugin->SetUser("blim"); 
        gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
        plugin->SetAPIVersion("V1.1x");
        plugin->SetAliROOTVersion("v5-06-15");
        plugin->SetAliPhysicsVersion("v5-06-15-01");
        plugin->SetAdditionalLibs("AliAnalysisTaskXic.cxx AliAnalysisTaskXic.h");
        plugin->SetAnalysisSource("AliAnalysisTaskXic.cxx");
   
        plugin->SetGridDataDir("/alice/data/2015/LHC15n");
        plugin->SetDataPattern("/pass4/*AliESDs.root");
        // MC has no prefix, data has prefix 000
        plugin->SetRunPrefix("000");
        // runnumber
        //plugin->AddRunNumber(167813);
        for (Int_t irun=runNmin;irun<runNmax;irun++){
            plugin->AddRunNumber((Int_t )runList[irun]);
        }
        plugin->SetSplitMaxInputFileNumber(10);
        
        plugin->SetExecutable("Xic.sh");
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
        plugin->SetRunMode(pluginmode);
        mgr->StartAnalysis("grid");
        
    }
}
void LoadLibraries()
{
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
    printf("Library Loading Complete");
}
void LoadMacros(Bool_t isMC)
{
    // compile the class (locally)
    gROOT->LoadMacro("AliAnalysisTaskXic.cxx++g");

    // Load macro for event selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection();

    // Load macro for PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(isMC);

    // Load macro for centrality selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

    // load the addtask macro
    gROOT->LoadMacro("AddXic.C");

    // Load Create ESD chain macro
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    // create an instance of your analysis task
    AliAnalysisTaskXic *task = AddXic();
    printf("Macro Loading Complete");
}
