void getResult(){
    TFile *filein = TFile::Open("AnalysisResults.root","READ");
    //TFile *filein = TFile::Open("001/AnalysisResults.root","READ");
    if(filein==NULL) continue;
    TDirectoryFile *dir = (TDirectoryFile*)filein->Get("Xic");
    TList *list_qa = (TList*)dir->Get("QA");
    TList *list_hist = (TList*)dir->Get("Histograms");
    (TH1F*)list_qa->FindObject("fInvLambda")->Draw();
    c1->SaveAs("./Figs/Lambda_Invmass.png");
}
