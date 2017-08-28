void LoadLibraries();
void LoadMacros(Bool_t fMCcase=kFALSE);

void runAnalysis(const char* pluginmode = "local")
{
    Bool_t fMCcase = kTRUE;
    Bool_t fAODcase = kFALSE;
    Int_t CutListOption = 0;

    LoadLibraries();

    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -g");
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    Int_t year = 2017;
    TString prod = "LHC17c";
    TString MCprod = "LHC15a2a";
    TString pass = "muon_calo_pass1";

    Int_t runNmin=0;
    //Int_t runNmax=22;
    Int_t runNmax=1;
    //Int_t runList[30]={244628, 244627, 244626, 244542, 244540, 244531, 244484, 244483, 244482, 244481, 244480, 244456, 244453, 244421, 244416, 244377, 244364, 244359, 244355, 244351, 244343, 244340};
    Int_t runList[1]={270667}; //for test
    Int_t MCrunList[1] = {130360};

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskXi_c");
    AliVEventHandler* esdH = new AliESDInputHandler();
    ((AliESDInputHandler *) esdH)->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);

    AliMCEventHandler *mcHandler = NULL;
    if(fMCcase){
      ::Info("AnalysisSetup", "Creating MC handler");
           mcHandler  = new AliMCEventHandler();
           mgr->SetMCtruthEventHandler(mcHandler);
    }

    LoadMacros(fMCcase);

    if(!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(pluginmode=="local") {
	     printf("LOCAL MODE");
        TChain* chain = new TChain("esdTree");
        //chain->Add("/home/blim/data/AliESDs.root"); // real data
        chain->Add("/home/blim/data/MC/LHC15a2a/130360/AliESDs.root"); // MC data
        mgr->StartAnalysis("local", chain);
    } else {
	     printf("GRID MODE");
        AliAnalysisAlien *plugin = new AliAnalysisAlien();

        plugin->SetRunMode(pluginmode);
        plugin->SetUser("blim");
        gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
        plugin->SetAPIVersion("V1.1x");
        plugin->SetAliROOTVersion("v5-06-15");
        plugin->SetAliPhysicsVersion("v5-06-15-01");
        plugin->SetNtestFiles(1);
        if (!fMCcase){
          plugin->SetAdditionalLibs("AliAnalysisTaskXic.cxx AliAnalysisTaskXic.h");
          plugin->SetAnalysisSource("AliAnalysisTaskXic.cxx");
        }
        else {
          plugin->SetAdditionalLibs("AliAnalysisTaskXicMC.cxx AliAnalysisTaskXicMC.h");
          plugin->SetAnalysisSource("AliAnalysisTaskXicMC.cxx");
        }

        if (!fMCcase) plugin->SetGridDataDir(Form("/alice/data/%i/%s",year,prod.Data())); // Real data path
        else plugin->SetGridDataDir(Form("/alice/sim/%i/%s",year,MCprod.Data()));      // MC data path
        //plugin->SetDataPattern(Form("%s/*AliESDs.root",pass.Data()));
        if (!fMCcase) plugin->SetDataPattern(Form("/%s/*ESDs.root",pass.Data())); // Real data has production
        else plugin->SetDataPattern(Form("/*ESDs.root"));                      // MC data doesn't have a production
        // MC has no prefix, data has prefix 000
        if (!fMCcase) plugin->SetRunPrefix("000");
        // runnumber
        //plugin->AddRunNumber(270667);
        Int_t nruns = 0;
        for (Int_t irun=runNmin;irun<runNmax;irun++){
            plugin->AddRunNumber(runList[irun]);
            nruns++;
        }
        nruns++;
        plugin->SetNrunsPerMaster(nruns);
        plugin->SetSplitMaxInputFileNumber(20);
        plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER -I$ALICE_ROOT/TOF -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include");

        plugin->SetExecutable("Xic.sh");
        plugin->SetTTL(72000);
        plugin->SetSplitMode("se");
        plugin->SetJDLName("Xic.jdl");
        plugin->SetAnalysisMacro("Xi_c_blim.C");
        plugin->SetValidationScript("Xic_validation.sh");

        plugin->SetDefaultOutputs(kFALSE);
        plugin->SetOutputToRunNo(kTRUE);
        plugin->SetOneStageMerging(kFALSE);
        plugin->SetKeepLogs(kTRUE);
        plugin->SetMaxMergeStages(1);
        plugin->SetMergeViaJDL(kTRUE);
        // define the output folders
        plugin->SetGridWorkingDir("Xi_c_Test");
        plugin->SetGridOutputDir("20170609_01");
	      plugin->SetOutputFiles("AnalysisResults.root");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(plugin);
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
void LoadMacros(Bool_t fMCcase)
{
    // compile the class (locally)
    if(!fMCcase) gROOT->LoadMacro("AliAnalysisTaskXic.cxx++g");
    else gROOT->LoadMacro('AliAnalysisTaskXicMC.cxx++g("test",fAODcase,fMCcase,CutListOption)');

    // Load macro for event selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    Bool_t applyPileupCuts = kTRUE;
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(fMCcase, applyPileupCuts);

    // Load macro for PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    Bool_t tuneOnData = kTRUE;
    TString recoPass = "2";
    if(fMCcase) AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(fMCcase, kTRUE, tuneOnData, recoPass);
    else AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(fMCcase);;

    // Load macro for centrality selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

    // load the addtask macro
    gROOT->LoadMacro("AddXic.C");
    // create an instance of your analysis task
    AliAnalysisTaskXic *task = AddXic();

    // Load Create ESD chain macro
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    printf("Macro Loading Complete");
}
