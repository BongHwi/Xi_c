void LoadLibraries();
void LoadMacros();

void runAnalysis_KIAF(const char *dataset = "test1.list")
{
    TString name = "Xi_c test";
    Bool_t fMCcase = kTRUE;
    Bool_t fAODcase = kFALSE;
    Int_t CutListOption = 0;
    
    LoadLibraries();
    LoadMacros();
    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisXic");
    AliESDInputHandler *esdH =  AddESDHandler();
    mgr->SetInputEventHandler(esdH);

    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection();
    AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(kFALSE);
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();

    // create an instance of your analysis task
    AliAnalysisTaskXic *task = AddXic(name, fAODcase, fMCcase, CutListOption);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    // if you want to run locally , we need to define some input
    if (!mgr->InitAnalysis()) printf("asdfasdf");
    TChain *chain = new TChain("ESDTree");
    //for list of input files
    chain = CreateESDChain(dataset);
    mgr -> StartAnalysis("local", chain);
}

void LoadLibraries()
{
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -g");

    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");

    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    printf("Library Loading Complete");
}
void LoadMacros()
{
    gROOT->LoadMacro("AliAnalysisTaskXic.cxx++g");
    gROOT->LoadMacro("AddXic.C");

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");

    printf("Macro Loading Complete");
}
