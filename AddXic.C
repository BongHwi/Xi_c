///////////////////////////////////////////////////////////////////
//                                                               //
// AddXic                                                        //
// Author: Redmer A. Bertens, Utrecht University, 2012           //
// Modified: Bong-Hwi Lim, Pusan National University, 2017       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskXic* AddXic(TString name = "name", Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption=0)
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    // get the input event handler, again via a static method.
    // this handler is part of the managing system and feeds events
    // to your task
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }
    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":Xic";      // create a subfolder in the file
    // now we create an instance of your task
    AliAnalysisTaskXic* task = new AliAnalysisTaskXic(name.Data(),AODdecision,MCdecision,CutListOption);
    if(!task) return 0x0;
    // add your task to the manager
    mgr->AddTask(task);
    // your task needs input: here we connect the manager to your task
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    // same for the output
    mgr->ConnectOutput(task,1,mgr->CreateContainer("QA", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    mgr->ConnectOutput(task,2,mgr->CreateContainer("Histograms", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
    // in the end, this macro returns a pointer to your task. this will be convenient later on
    // when you will run your analysis in an analysis train on grid
    return task;
}
