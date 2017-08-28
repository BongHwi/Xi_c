/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskXic
 * test for commit4
 *
 */

 #include <iostream>
 #include <math.h>
 #include "TChain.h"
 #include "TFile.h"
 #include "TKey.h"
 #include "TObject.h"
 #include "TObjString.h"
 #include "TList.h"
 #include "TTree.h"
 #include "TH1F.h"
 #include "TH1D.h"
 #include "TH2D.h"
 #include "TH3D.h"
 #include "TProfile.h"
 #include "TCanvas.h"

 #include "AliAnalysisTask.h"
 #include "AliAnalysisManager.h"
 #include "AliLog.h"

 #include "AliESDEvent.h"
 #include "AliESDtrackCuts.h"
 #include "AliESDInputHandler.h"

 #include "AliAODEvent.h"
 #include "AliAODInputHandler.h"

 #include "AliVEvent.h"
 #include "AliVEventHandler.h"
 #include "AliCentrality.h"
 #include "AliAODcascade.h"
 #include "AliESDcascade.h"
 #include "AliV0vertexer.h"
 #include "AliCascadeVertexer.h"

 #include "AliMCEventHandler.h" // MC
 #include "AliMCEvent.h"        // MC
 #include "AliMCParticle.h"     // MC
 #include "AliStack.h"          // MC
 #include "AliAODMCParticle.h"  // MC

 #include "AliESDv0.h"
 #include "AliKFParticle.h"
 #include "AliKFVertex.h"
 #include "AliAnalysisTaskXic.h"


class AliAnalysisTaskXic;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskXic) // classimp: necessary for root
double getEnergy(Double_t trueMass, Double_t Px, Double_t Py, Double_t Pz);
double getAngle(Double_t Px1, Double_t Py1, Double_t Pz1, Double_t Px2, Double_t Py2, Double_t Pz2);
void CheckChargeV0(AliESDv0 *v0);

AliAnalysisTaskXic::AliAnalysisTaskXic() : AliAnalysisTaskSE(),
    fESD(0x0),
    fOutputList(0x0),
    fTrackCuts(0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0),
    fEventsToMix(0),
    fMCcase(0),
    fAODcase(0),
    fEventCounter(0),
    fCutList(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskXic::AliAnalysisTaskXic(const char* name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption) : AliAnalysisTaskSE(name),
    fESD(0x0),
    fOutputList(0x0),
    fTrackCuts(0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0),
    fEventsToMix(0),
    fMCcase(MCdecision),
    fAODcase(AODdecision),
    fEventCounter(0),
    fCutList(CutListOption)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
    // this chain is created by the analysis manager, so no need to worry about it,
    // it does its work automatically
    DefineOutput(1, TList::Class());    // QA histograms
    DefineOutput(2, TList::Class());    // Outputs
}
//_____________________________________________________________________________
AliAnalysisTaskXic::~AliAnalysisTaskXic()
{
    // destructor
    if (fESD) delete fESD;
    // if(fPIDResponse) delete fPIDResponse;
    // if(fCentrality) delete fCentrality;
    if (fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
    if (fOutputList2) {
        delete fOutputList2;     // at the end of your task, it is deleted from memory by calling this function
    }
    if (fTrackCuts) delete fTrackCuts;

}
//_____________________________________________________________________________
void AliAnalysisTaskXic::XicInit()
{
  //
  //Inits cuts and analysis settings
  //
  fEventCounter=0;// event counter initialization
  cout<<"AliAnalysisTaskXic XicInit() call"<<endl;
  fZvertexBins = 20;

  if(fMCcase) fEventsToMix = 0;
  else fEventsToMix = 40;
}
//_____________________________________________________________________________
void AliAnalysisTaskXic::UserCreateOutputObjects()
{
    // create output objects
    //
    // this function is called ONCE at the start of your analysis (RUNTIME)
    // here you ceate the histograms that you want to use
    //
    // the histograms are in this case added to a tlist, this list is in the end saved
    // to an output file
    //
    fOutputList = new TList();          // this is a list which will contain all of your histograms
    fOutputList2 = new TList();
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
    fOutputList2->SetOwner(kTRUE);

    fTrackCuts = new AliESDtrackCuts();
    //------------------------------------------------
    // QA histograms
    //------------------------------------------------
    // Deafult Analysis setup
    TH1F *hNofV0 = new TH1F("hNofV0", "hNofV0", 40, 0, 40);
      hNofV0->GetXaxis()->SetTitle("# of v0 through each cut");
    TH1F *hNofV0_2 = new TH1F("hNofV0_2", "hNofV0_2", 40, 0, 40);
      hNofV0_2->GetXaxis()->SetTitle("# of v0 through each cut");
    TH1F *hEventSelecInfo = new TH1F("hEventSelecInfo", "hEventSelecInfo", 10, 0, 10);
      hEventSelecInfo->GetXaxis()->SetBinLabel(1, "NONE");
      hEventSelecInfo->GetXaxis()->SetBinLabel(2, "kMB");
      hEventSelecInfo->GetXaxis()->SetBinLabel(3, "kCentral");
      hEventSelecInfo->GetXaxis()->SetBinLabel(4, "kSemiCentral");
      hEventSelecInfo->GetXaxis()->SetBinLabel(5, "kINT7");
      hEventSelecInfo->GetXaxis()->SetBinLabel(6, "kAny");
      hEventSelecInfo->GetXaxis()->SetBinLabel(7, "kPhysicsALL");
    TH1F *hCentrality = new TH1F("hCentrality", "Centrality", 100, 0, 100);
      hCentrality->GetXaxis()->SetTitle("Centrality");
    TH1F *fMultDist = new TH1F("fMultDist", "Multiplicity Distribution", 200, 0, 20000);
      fMultDist->GetXaxis()->SetTitle("Multiplicity");
    TH3F *fVertexDistXYZ = new TH3F("fVertexDistXYZ", "Vertex Distribution", 20, -1, 1, 20, -1, 1, 60, -30, 30);
      fVertexDistXYZ->GetXaxis()->SetTitle("X Vertex (cm)");
      fVertexDistXYZ->GetYaxis()->SetTitle("Y Vertex (cm)");
      fVertexDistXYZ->GetZaxis()->SetTitle("Z Vertex (cm)");
    // Labmda QA
    TH1F *fNTPCcls_lam = new TH1F("fNTPCcls_lam","fNTPCcls_lam",200,0,200);
      fNTPCcls_lam->GetXaxis()->SetTitle("fNTPCcls_lam");
    TH1F *fHistCosPA_lam = new	TH1F("fHistCosPA_lam", "Cosine of Pointing Angle of V0s; Cos PA; N(v0s)",202,0.8,1.01);
    TH1F *fHistDCAV0Daughters_lam = new	TH1F("fHistDCAV0Daughters_lam", "DCA between V0 daughters; DCA (cm); N V0s", 100, 0, 2);
  	TH1F *fHistDecayL_lam = new	TH1F("fHistDecayL_lam", "Distance between V0 and PV; Distance(cm); N(v0s)",200,-0.1,30);
    TH1F *fHistTauLa = new	TH1F("fHistTauLa", "Lifetime under Lambda mass hypothesis; Lifetime(s); N(v0s)",200,0,100);
    TH2F *fHistBetheBlochTPCNeg_lam = new	TH2F("fHistBetheBlochTPCNeg_lam","-dE/dX against Momentum for negative daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
  	TH2F *fHistBetheBlochTPCPos_lam = new	TH2F("fHistBetheBlochTPCPos_lam","-dE/dX against Momentum for positive daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
    // Pion QA
    TH1F *fNTPCcls_Pi = new TH1F("fNTPCcls_Pi","fNTPCcls_Pi",200,0,200);
      fNTPCcls_Pi->GetXaxis()->SetTitle("fNTPCcls_Pi");
    TH2F *fHistBetheBlochTPC_Pi = new	TH2F("fHistBetheBlochTPC_Pi","-dE/dX against Momentum from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);

    // Armenteros-Podolanski Plot
    TH2F *fArmPod_kaon = new TH2F("fArmPod_kaon", "Armenteros-Podolanski Plot", 800, -1.0, 1.0, 100, 0, 0.25);
    TH2F *fArmPod_kaon_cut = new TH2F("fArmPod_kaon_cut", "Armenteros-Podolanski Plot after cut", 800, -1.0, 1.0, 100, 0, 0.25);
    TH2F *fArmPod_lambda = new TH2F("fArmPod_lambda", "Armenteros-Podolanski Plot", 800, -1.0, 1.0, 100, 0, 0.25);
    TH2F *fArmPod_lambda_cut = new TH2F("fArmPod_lambda_cut", "Armenteros-Podolanski Plot after cut", 800, -1.0, 1.0, 100, 0, 0.25);

    // K0s
    TH1F *fInvPion = new TH1F("fInvPion", "Invariant mass distribution of K0s", 200, 0.1, 0.2);
      fInvPion->GetXaxis()->SetTitle("fInvPion");
    TH1F *fInvPionCut = new TH1F("fInvPionCut", "Invariant mass distribution of K0s after mass window cut", 200, 0.1, 0.2);
      fInvPion->GetXaxis()->SetTitle("fInvPionCut");
    // Lambda0
    TH1F *fInvLambda = new TH1F("fInvLambda", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
      fInvLambda->GetXaxis()->SetTitle("fInvLambda");
    TH1F *fInvLambda_beforePID = new TH1F("fInvLambda_beforePID", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
      fInvLambda_beforePID->GetXaxis()->SetTitle("fInvLambda_beforePID");
    TH1F *fInvLambdaCut = new TH1F("fInvLambdaCut", "Invariant mass distribution of Lambda after mass window cut", 400, 1.0, 1.2);
      fInvLambdaCut->GetXaxis()->SetTitle("fInvLambdaCut");

    if (! fHistSwappedV0Counter) {
        fHistSwappedV0Counter = new TH1F("fHistSwappedV0Counter",
                                         "Swap or not histo;Swapped (1) or not (0); count",
                                         2, 0, 2);
        fOutputList->Add(fHistSwappedV0Counter);
    }
    TH1F *fMultDistMC = new TH1F("fMultDistMC", "Multiplicity Distribution of MC", 200, 0, 20000);
      fMultDistMC->GetXaxis()->SetTitle("Multiplicity");
    TH3F *fVertexDistMC = new TH3F("fVertexDistMC", "Vertex Distribution", 20, -1, 1, 20, -1, 1, 60, -30, 30);
      fVertexDistMC->GetXaxis()->SetTitle("X Vertex (cm)");
      fVertexDistMC->GetYaxis()->SetTitle("Y Vertex (cm)");
      fVertexDistMC->GetZaxis()->SetTitle("Z Vertex (cm)");
    TH3F *fMCinputTotalXi1 = new TH3F("fMCinputTotalXi1","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);
    TH3F *fMCinputTotalXibar1 = new TH3F("fMCinputTotalXibar1","Invariant Mass Distribution", 100,0,10, 40,-2,2, 200,1.2,1.4);

    // Analysis Results
    TH2F *hInvMassWithPt = new TH2F("hInvMassWithPt", "Invariant mass distribution vs Pt", 1000, 2.0, 3.0, 100, 0, 10);
    TH1F *hInvMass = new TH1F("hInvMass", "Invariant mass distribution", 1000, 2.0, 3.0);

    fOutputList->Add(hNofV0);
    fOutputList->Add(hNofV0_2);
    fOutputList->Add(hEventSelecInfo);
    fOutputList->Add(hCentrality);
    fOutputList->Add(fMultDist);
    fOutputList->Add(fVertexDistXYZ);
    fOutputList->Add(fNTPCcls_lam);
    fOutputList->Add(fHistCosPA_lam);
    fOutputList->Add(fHistDCAV0Daughters_lam);
  	fOutputList->Add(fHistDecayL_lam);
  	fOutputList->Add(fHistTauLa);
	  fOutputList->Add(fHistBetheBlochTPCNeg_lam);
    fOutputList->Add(fHistBetheBlochTPCPos_lam);
    fOutputList->Add(fNTPCcls_Pi);
	  fOutputList->Add(fHistBetheBlochTPC_Pi);

    fOutputList->Add(fArmPod_kaon);
    fOutputList->Add(fArmPod_kaon_cut);
    fOutputList->Add(fArmPod_lambda);
    fOutputList->Add(fArmPod_lambda_cut);

    fOutputList->Add(fInvPion);
    //fOutputList->Add(fInvPion_beforePID);
    fOutputList->Add(fInvPionCut);
    fOutputList->Add(fInvLambda);
    //fOutputList->Add(fInvLambda_beforePID);
    fOutputList->Add(fInvLambdaCut);

    fOutputList2->Add(hInvMassWithPt);
    fOutputList2->Add(hInvMass);

    if(fMCcase){
      fOutputList->Add(fMultDistMC);
      fOutputList->Add(fVertexDistMC);
      fOutputList->Add(fMCinputTotalXi1);
      fOutputList->Add(fMCinputTotalXibar1);
    }
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if (!fPIDResponse) AliError ("No PID");

    PostData(1, fOutputList);
    PostData(2, fOutputList2);
}
//_____________________________________________________________________________
void AliAnalysisTaskXic::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    cout<<"===========  Event # "<<fEventCounter+1<<"  ==========="<<endl;
    fEventCounter++;

    if(fAODcase) {cout<<"AODs not fully supported! Exiting event."<<endl; return;}
    if(fAODcase) fAOD = dynamic_cast<AliAODEvent*> (InputEvent());
    else fESD = dynamic_cast<AliESDEvent*> (InputEvent());

    if(fAODcase) {if (!fAOD) {Printf("ERROR: fAOD not available"); return;}}
    else {if (!fESD) {Printf("ERROR: fESD not available"); return;}}

    // ESD Trigger Cut
    if(!fAODcase){
      if(!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected())) {
        cout<<"Event Rejected"<<endl; return;
      }
    }


    ///////////////////////////////////////////////////////////
    const AliAODVertex *PrimaryVertexAOD;
    const AliESDVertex *PrimaryVertexESD;

    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    TClonesArray *mcArray       = 0x0;
    //AliAODMCParticle *mcXi;
    //AliAODMCParticle *mcXi_c;
    AliMCEvent  *mcEvent        = 0x0;
    AliStack    *mcstack        = 0x0;
    TParticle   *MCLamD1esd     = 0x0;
    TParticle   *MCLamD2esd     = 0x0;
    TParticle   *MCLamesd       = 0x0;
    TParticle   *MCXiesd        = 0x0;
    TParticle   *MCXiStaresd    = 0x0;

    Double_t px1,py1,pz1, px2,py2,pz2;
    Double_t p1sq,p2sq,e1,e2,angle;
    Double_t dca3d;
    Float_t dca2[2];
    Double_t xiVtx[3];//, xiStarVtx[3];
    Double_t xiP[3], xiStarP[3];
    Double_t xiCMom;
    Double_t xiMass, xiStarMass;
    Double_t xiPt, xiStarPt;
    Double_t xiY, xiStarY;
    Double_t xiCharge;
    Double_t decayLengthXY;
    Double_t pDaughter1[3];
    Double_t pDaughter2[3];
    Double_t xDaughter1[3];
    Double_t xDaughter2[3];
    //
    Double_t bField=0;
    UInt_t status=0;
    Int_t positiveTracks=0, negativeTracks=0;
    Int_t myTracks=0;
    //
    Double_t primaryVtx[3]={0};
    Int_t mBin=0;
    Int_t zBin=0;
    Double_t zStep=2*10/Double_t(fZvertexBins), zStart=-10.;
    //
    Bool_t mcXiFilled=kFALSE;// So that mctracks are never used uninitialized

    if(fMCcase){
      if(fAODcase){
        mcArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
        if(!mcArray){
  	cout<<"No MC particle branch found"<<endl;
  	return;
        }
      }else {
        mcEvent = MCEvent();
        if (!mcEvent) {Printf("ERROR: Could not retrieve MC event"); return;}

        mcstack = mcEvent->Stack();
        if (!mcstack) {Printf("ERROR: Could not retrieve the stack"); return;}
      }
    }

    if(fAODcase){// AOD case
      cout<<"Currently AOD Case is not supported."<<endl;
    }else {// ESDs
      ((TH1F*)fOutputList->FindObject("fMultDistMC"))->Fill(fESD->GetNumberOfTracks());
      PrimaryVertexESD = fESD->GetPrimaryVertex();
      if(!PrimaryVertexESD) {cout<<"ERROR: Could not retrieve Vertex"<<endl; return;}

      primaryVtx[0]=PrimaryVertexESD->GetX();
      primaryVtx[1]=PrimaryVertexESD->GetY();
      primaryVtx[2]=PrimaryVertexESD->GetZ();
      ((TH3F*)fOutputList->FindObject("fVertexDistMC"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);

      if(fMCcase){
        /////////////////////////////////////////////////
        // Lam mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
          TParticle *mcInputTrack = (TParticle*)mcstack->Particle(it);
          if (!mcInputTrack) {
        	  Error("UserExec", "Could not receive track %d", it);
        	  continue;
        	}
        
        // Xi
      	if(mcInputTrack->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXi1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
      	if(mcInputTrack->GetPdgCode() == -kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXibar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
      }
    }
    if(0){

    Int_t debugmode = 0; // for debuging, 101 for general debuging, 51 for specific debuging, 11 for only check v0 survived
    // Parameters used for cuts.
    double cutCosPa(0.998), cutcTau(2);
	  double cutNImpact(-999), cutDCA(0.9);
  	double cutBetheBloch(3);
  	double cutMinNClustersTPC(70), cutMaxChi2PerClusterTPC(-999);
  	double cutEta(0.8);

    //Track Cuts set here
    if(cutMinNClustersTPC != -999) (fTrackCuts->SetMinNClustersTPC(int(cutMinNClustersTPC)));
    if(cutMaxChi2PerClusterTPC != -999) fTrackCuts->SetMaxChi2PerClusterTPC(cutMaxChi2PerClusterTPC);
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);

    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {Printf("ERROR: fESD not available"); return;}
    if(debugmode > 100) AliInfo("test!");
    //------------------------------------------------
    //Step 1: Check for selected Trigger
    //------------------------------------------------
    Bool_t isSelectedMB = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if (isSelectedMB) ((TH1F*)fOutputList->FindObject("hEventSelecInfo"))->Fill(1);
    Bool_t isSelectedkCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kCentral);
    if (isSelectedkCentral) ((TH1F*)fOutputList->FindObject("hEventSelecInfo"))->Fill(2);
    Bool_t isSelectedkSemiCentral = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kSemiCentral));
    if (isSelectedkSemiCentral) ((TH1F*)fOutputList->FindObject("hEventSelecInfo"))->Fill(3);
    Bool_t isSelectedINT7 = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() && AliVEvent::kINT7);
    if (isSelectedINT7) ((TH1F*)fOutputList->FindObject("hEventSelecInfo"))->Fill(4);
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() && AliVEvent::kAny);
    if (isSelected) ((TH1F*)fOutputList->FindObject("hEventSelecInfo"))->Fill(5);
    //if (!isSelected)Printf("There is events in kANY");
    ////////////******* Do Event selecction *******////////////
    if (!(isSelectedINT7 | isSelectedMB | isSelectedkCentral | isSelectedkSemiCentral)) {cout << "Event Rejected" << endl; return;}
    //cout << "Event Accepted" << endl;
    ((TH1F*)fOutputList->FindObject("fMultDist"))->Fill(fESD->GetNumberOfTracks());
    if(debugmode > 100) AliInfo("after trigger selecction!");
    //------------------------------------------------
    //Step 2: Check for centrality for Pb-Pb
    //------------------------------------------------
    Float_t  centralityV0M = -100;
    fCentrality = fESD->GetCentrality();
    centralityV0M = fCentrality->GetCentralityPercentile("V0M");
    ((TH1F*)fOutputList->FindObject("hCentrality"))->Fill(centralityV0M);

    if(debugmode > 100) AliInfo("after centrality check!");
    //------------------------------------------------
    //Step 3: Check for Vertex-Z position
    //------------------------------------------------
    const AliESDVertex *PrimaryVertexESD;
    PrimaryVertexESD = fESD->GetPrimaryVertex();
    //if(!PrimaryVertexESD) return;
    //if(PrimaryVertexESD->GetNContributors() < 1) return;

    Double_t primaryVtx[3] = {0};
    primaryVtx[0] = PrimaryVertexESD->GetX(); // call primary vertex position of X
    primaryVtx[1] = PrimaryVertexESD->GetY(); // call primary vertex position of Y
    primaryVtx[2] = PrimaryVertexESD->GetZ(); // call primary vertex position of Z
    ((TH3F*)fOutputList->FindObject("fVertexDistXYZ"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
    ////////////******* DO Vertex-Z selecction *******////////////
    if (fabs(primaryVtx[2]) > 10.) return;
    if(debugmode > 100) AliInfo("after Vertex-Z position check!");
    //------------------------------------------------
    //Step 4: Check for SPD Pileup
    //------------------------------------------------
    if (fESD->IsPileupFromSPD()) return; // Reject Pile-up events

    //------------------------------------------------
    //Step 5: V0 Loop
    //------------------------------------------------
    Double_t  lMagneticField = fESD->GetMagneticField();
    // mass constant
    static Double_t k0Mass = 0.497611;
    static Double_t l0Mass = 1.115683;
    static Double_t piMass = 0.13957;
    static Double_t prMass = 0.93827;

    // cut values
    Bool_t fkUseOnTheFly = kFALSE;
    Double_t fMinV0Pt = 0.15;
    Double_t fMaxV0Pt = 1.E10;

    // loop for Lambda
    Int_t nv0s = 0;
    nv0s = fESD->GetNumberOfV0s();
    const Int_t nv0sc = nv0s;
    Int_t v0checklam[nv0sc];
    Int_t v0checkk0s[nv0sc];
    // loop for pion track
    Int_t iTracks(fESD->GetNumberOfTracks());

    if(debugmode > 100) AliInfo("Starting V0 loop!");
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++){
        v0checklam[iV0] = 0;
        bool lambdaCandidate = true;
        bool antilambdaCandidate = true;
        // keep only events of interest for fHistMLa plots
        // from PWGLF/STRANGENESS/LambdaK0PbPb/AliAnalysisTaskLukeV0.cxx

        AliESDv0 *v0i = ((AliESDEvent*)fESD)->GetV0(iV0);
        if (!v0i) continue;
        if(debugmode > 51) ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(1);
        if(debugmode > 100) AliInfo("01");
        //---> Fix On-the-Fly candidates, count how many swapped
        if ( v0i->GetParamN()->Charge() > 0 && v0i->GetParamP()->Charge() < 0 ) {
            fHistSwappedV0Counter -> Fill( 1 );
        } else {
            fHistSwappedV0Counter -> Fill( 0 );
        }
        if ( fkUseOnTheFly ) CheckChargeV0(v0i);

        Int_t    lOnFlyStatus = 0; // nv0sOn = 0, nv0sOff = 0;
        lOnFlyStatus = v0i->GetOnFlyStatus();
        if (lOnFlyStatus == 0) continue;
        if(debugmode > 100) AliInfo("02");
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(3);

        //// Get V0 informations for the cuts
        Double_t lPt = 0;
        lPt = v0i->Pt();
        if (v0i->GetEffMass(4,2) < 1.08 || v0i->GetEffMass(4,2) > 1.2 || TMath::Abs(v0i->Y(3122))>0.5 ) lambdaCandidate = false;
        if (v0i->GetEffMass(2,4) < 1.08 || v0i->GetEffMass(2,4) > 1.2 || TMath::Abs(v0i->Y(-3122))>0.5) antilambdaCandidate = false;
        if(debugmode > 50 && lambdaCandidate) AliInfo("Lambda0 Case");
        if(debugmode > 50 && antilambdaCandidate) AliInfo("Anti-Lambda0 Case");
        // get daughter particle
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0i->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0i->GetNindex());
        AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
        AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
        const AliExternalTrackParam * paramPosl = v0i->GetParamP();
        const AliExternalTrackParam * paramNegl = v0i->GetParamN();
        // TPC nCluster
        Int_t fTPCNcls = -100;
        fTPCNcls = pTrack->GetTPCNcls();
        // Cosine Pointing Angle and DCA Values
        Double_t lV0cosPointAngle = v0i->GetV0CosineOfPointingAngle(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
        Double_t lV0Position[3];
        Double_t lV0Radius = 0;
        v0i->GetXYZ(lV0Position[0], lV0Position[1], lV0Position[2]);
        lV0Radius = TMath::Sqrt(lV0Position[0] * lV0Position[0] + lV0Position[1] * lV0Position[1]);
        Double_t lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(primaryVtx[0], primaryVtx[1], lMagneticField));
        Double_t lDcaNegToPrimVertex = TMath::Abs(pTrack->GetD(primaryVtx[0], primaryVtx[1], lMagneticField));
        Double_t lDcaV0Daughters = v0i->GetDcaV0Daughters();
        Double_t tV0momi[3], tV0momj[3], tV0mom_result[3];
        v0i->GetPxPyPz(tV0momi[0], tV0momi[1], tV0momi[2]);
        // Decay length
        double decayLength = (sqrt((primaryVtx[0]-primaryVtx[0])*(primaryVtx[0]-primaryVtx[0])+(primaryVtx[1]-primaryVtx[1])*(primaryVtx[1]-primaryVtx[1])+(primaryVtx[2]-primaryVtx[2])*(primaryVtx[2]-primaryVtx[2])));
        double cTauLa = decayLength*(v0i->GetEffMass(4,2))/(v0i->P());
        double cTauLb = decayLength*(v0i->GetEffMass(2,4))/(v0i->P());
        // momentums
        double pTrackMomentum[3];
        double nTrackMomentum[3];
        pTrack->GetConstrainedPxPyPz(pTrackMomentum);
        nTrack->GetConstrainedPxPyPz(nTrackMomentum);
        double pPos2 = sqrt(pTrackMomentum[0]*pTrackMomentum[0]+pTrackMomentum[1]*pTrackMomentum[1]+pTrackMomentum[2]*pTrackMomentum[2]);
        double pNeg2 = sqrt(nTrackMomentum[0]*nTrackMomentum[0]+nTrackMomentum[1]*nTrackMomentum[1]+nTrackMomentum[2]*nTrackMomentum[2]);

        if(debugmode > 100) AliInfo("03");
        //// Cuts
        if(!(fTrackCuts->IsSelected(pTrack)) || !(fTrackCuts->IsSelected(nTrack))) {
          lambdaCandidate = false;
          antilambdaCandidate = false;
        }
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(!(lambdaCandidate == false && antilambdaCandidate == false) && debugmode > 50) AliInfo("Track cut!");
        ((TH1F*)fOutputList->FindObject("fNTPCcls_lam"))->Fill(fTPCNcls);
        // Pt cut for mother particle
        if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;
        if(debugmode > 51) ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(5);
        // is daughter particle okay?
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(7);
        // Like sign cut
        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(9);
        //psuedorapidity cut
        if(cutEta != -999) {
          if(TMath::Abs(pTrack->Eta()) > cutEta || TMath::Abs(nTrack->Eta())  >cutEta) {
            lambdaCandidate = false;
            antilambdaCandidate = false;
          }
        }
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(!(lambdaCandidate == false && antilambdaCandidate == false) && debugmode > 50) AliInfo("Eta cut!"); ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(11);
        // CPA cut
        if(cutCosPa != -999) {
          if (lV0cosPointAngle < cutCosPa){
            lambdaCandidate = false;
            antilambdaCandidate = false;
          }
        }
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(!(lambdaCandidate == false && antilambdaCandidate == false) && debugmode > 50) AliInfo("CPA cut!"); ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(13);
        // DCA between daughterscut
        if(cutDCA != -999) {
          if(v0i->GetDcaV0Daughters() > cutDCA) {
            lambdaCandidate = false;
            antilambdaCandidate = false;
          }
        }
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(!(lambdaCandidate == false && antilambdaCandidate == false) && debugmode > 50) AliInfo("DCA cut!"); ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(15);
        //lifetime cut
        /*
        if(cutcTau != -999){
          if(cTauLa < cutcTau){
            lambdaCandidate = false;
          }
          if(cTauLb < cutcTau){
            antilambdaCandidate = false;
          }
        }*/
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(!(lambdaCandidate == false && antilambdaCandidate == false) && debugmode > 50) AliInfo("Lifetime cut!"); ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(17);
        // Bethe Bloch cut. Made sightly complicated as options for crude cuts still included. Should probably reduce to just 'official' cuts
        if(cutBetheBloch != -999) {
          if(pTrack->GetTPCsignal() <0 || nTrack->GetTPCsignal()<0) continue;
          if(lambdaCandidate) {
            if(cutBetheBloch > 0) {
            if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kProton)) > cutBetheBloch ) lambdaCandidate = false;
            if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > cutBetheBloch ) lambdaCandidate = false;
            }
            if(cutBetheBloch == -4) {
              double beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+prMass*prMass))),2);
              double gamma2 = 1.0/(1.0-beta2);
              if(pTrack->GetTPCsignal() < (2.3/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2)) lambdaCandidate = false;
            }
          }
          if(antilambdaCandidate) {
            if(cutBetheBloch > 0) {
              if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kProton)) > cutBetheBloch )
              {antilambdaCandidate = false;}
              if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > cutBetheBloch )
              {antilambdaCandidate = false;}
            }

            if(cutBetheBloch == -4) {
                double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+0.9*0.9))),2);
                double gamma2 = 1.0/(1.0-beta2);
                if(nTrack->GetTPCsignal() < (2.3/beta2)*(TMath::Log(1e6*beta2*gamma2)-beta2)) antilambdaCandidate = false;
            }
          }
        }
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(!(lambdaCandidate == false && antilambdaCandidate == false) && debugmode > 50) AliInfo("PID cut!"); ((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(19);
        //if ((lDcaPosToPrimVertex < 0.1) || (lDcaNegToPrimVertex < 0.1) || (lV0cosPointAngle < 0.998) || (lV0Radius < 0.0) || (lV0Radius > 1000) ) continue;
        if(debugmode > 100) AliInfo("04");

        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(25);
        //remove all non-candidates
        if(lambdaCandidate == false && antilambdaCandidate == false) continue;
        if(antilambdaCandidate) continue; // for the Labmda_c
        if(debugmode > 10) AliInfo("============v0 survived!============");
        if(pPos2 == 0 || pNeg2 ==0) continue;
        ((TH1F*)fOutputList->FindObject("fHistCosPA_lam"))->Fill(lV0cosPointAngle);
        ((TH1F*)fOutputList->FindObject("fHistDecayL_lam"))->Fill(decayLength);
        ((TH1F*)fOutputList->FindObject("fHistTauLa"))->Fill(cTauLa);
        ((TH2F*)fOutputList->FindObject("fHistBetheBlochTPCPos_lam"))->Fill(TMath::Log10(pPos2),pTrack->GetTPCsignal());
        ((TH2F*)fOutputList->FindObject("fHistBetheBlochTPCNeg_lam"))->Fill(TMath::Log10(pNeg2),nTrack->GetTPCsignal());

        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(27);
        // Mass Hypothesis for Lambda
        //v0i->ChangeMassHypothesis(3122);
        //sets assumed particle type of pos/neg daughters.
        // 0 = electron, 1 = Muon, 2 = pion, 3 = kaon, 4 = proton.
        int dPos = 4;
        int dNeg = 2;
        if(!(v0i->GetEffMass(dPos,dNeg) > 1.11 && v0i->GetEffMass(dPos,dNeg) < 1.13)) continue;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(29);
        double lInvMassLambda = 0.;
        if(lambdaCandidate) lInvMassLambda = v0i->GetEffMass(dPos,dNeg);
        if(antilambdaCandidate) lInvMassLambda = v0i->GetEffMass(dPos,dNeg);
        if(debugmode > 100) AliInfo("04-1");

        if (!((pTrack->GetMass() > 0.9 && nTrack->GetMass() < 0.2)||(pTrack->GetMass() < 0.2 && nTrack->GetMass() > 0.9))) continue;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(31);
        if(debugmode > 50) AliInfo("daughter mass cut pass");
        ((TH2F*)fOutputList->FindObject("fArmPod_lambda"))->Fill(v0i->AlphaV0(),v0i->PtArmV0());

        if(debugmode > 100) AliInfo("04-5");
        // Armenteros-Podolansiki Cut
        if (TMath::Abs(0.2 * v0i->AlphaV0()) < v0i->PtArmV0()) continue;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0"))->Fill(33);
        ((TH2F*)fOutputList->FindObject("fArmPod_lambda_cut"))->Fill(v0i->AlphaV0(),v0i->PtArmV0());
        if(debugmode > 100) AliInfo("05");
        if(lambdaCandidate) v0checklam[iV0] = 1;
        //cout << "Lambda : " << v0checklam[iV0] << endl;
        //if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 0 && fkUseOnTheFly == kTRUE ) ){
        ((TH1F*)fOutputList->FindObject("fInvLambda"))->Fill(lInvMassLambda);
        //if (lInvMassLambda > l0Mass + 0.0008 || lInvMassLambda < l0Mass - 0.008) continue; // Mass window
        ((TH1F*)fOutputList->FindObject("fInvLambdaCut"))->Fill(lInvMassLambda); // After Cut
        //}
        for(Int_t i(0); i < iTracks; i++) { // new loop for pion+
          AliESDtrack* track = fESD->GetTrack(i);
          if(!track) continue;
          Int_t fSign = 0;
          fSign = track->GetSign();
          if(fSign < 0) continue;  // only pion+
          if(!fTrackCuts->AcceptTrack(track)) continue;
          Int_t fTPCNcls_Pi = 0;
          track->GetTPCNcls();
          ((TH1F*)fOutputList->FindObject("fNTPCcls_Pi"))->Fill(fTPCNcls_Pi);
          double fPt = 0.;
          double fx = 0.;
          double fy = 0.;
          double fz = 0.;
          double fPxpi = 0.;
          double fPypi = 0.;
          double fPzpi = 0.;
          double fMpi = 0.;
          fPt = track->Pt();
          fx = track->GetX(); // will be used for DCA cut
          fy = track->GetY();
          fz = track->GetZ();
          fPxpi = track->Px();
          fPypi = track->Py();
          fPzpi = track->Pz();
          fMpi = track->M();
          Float_t nsigpi= fabs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
          if(TMath::Abs(nsigpi) > 3) continue;
          double Ppi = sqrt(fPxpi*fPxpi+fPypi*fPypi+fPzpi*fPzpi);
          ((TH2F*)fOutputList->FindObject("fHistBetheBlochTPC_Pi"))->Fill(TMath::Log10(Ppi),track->GetTPCsignal());
          ((TH1F*)fOutputList->FindObject("fInvPion"))->Fill(fMpi);
          if(TMath::Abs(fMpi-piMass) > 0.02) continue;
          ((TH1F*)fOutputList->FindObject("fInvPionCut"))->Fill(fMpi);
          double fDCAr = sqrt( pow(fx - lV0Position[0],2) + pow(fy - lV0Position[1],2) );
          double fDCAz = fabs(fz - lV0Position[2]);
          //if(fDCAr > 0.1 || fDCAz > 0.1) continue; // DCA-r,z cut

          Double_t ei = 0.;
          Double_t ej = 0.;
          Double_t angle = 0.;
          Double_t fMass = 0.;
          Double_t fPt_result = 0.;

          ei = getEnergy(fMpi, fPxpi, fPypi, fPzpi); // Energy of first particle(Pi+)
          ej = getEnergy(lInvMassLambda, tV0momi[0], tV0momi[1], tV0momi[2]); // Energy of first particle(Lambda)
          //cout << "Energy" << endl;
          angle = getAngle(fPxpi, fPypi, fPzpi, tV0momi[0], tV0momi[1], tV0momi[2]);

          //cout << "ei: " << ei << ", ej: " << ej << ", angle: "<< angle << endl;

          fMass = fMpi * fMpi + lInvMassLambda * lInvMassLambda + 2.*ei * ej - 2.*angle;
          if (fMass <= 0) continue;
          fMass = sqrt(fMass);

          tV0mom_result[0] = fPxpi + tV0momj[0];
          tV0mom_result[1] = fPypi + tV0momj[1];
          tV0mom_result[2] = fPzpi + tV0momj[2];
          fPt_result = sqrt(pow(tV0mom_result[0], 2) + pow(tV0mom_result[1], 2));

          ((TH2F*)fOutputList2->FindObject("hInvMassWithPt"))->Fill(fMass, fPt_result); // with Pt
          ((TH1F*)fOutputList2->FindObject("hInvMass"))->Fill(fMass); // Cumulated
        }

    }/*
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++){
      for (Int_t jV0 = iV0; jV0 < nv0s; jV0++){
        //if(v0checklam[iV0]) cout << iV0 << " Lambda: "<< v0checklam[iV0] << endl;
        //if(v0checkk0s[iV0]) cout << jV0 << " K0s: " << v0checkk0s[iV0] << endl;
        //if(v0checklam[iV0]) cout << "Lambda: "<< v0checklam[iV0] << " K0s: "<< v0checkk0s[iV0] << endl;
        if(!(v0checklam[iV0] && v0checkk0s[iV0])) continue;
        cout << "Lambda: "<< v0checklam[iV0] << " K0s: "<< v0checkk0s[iV0] << endl;
        AliESDv0 *v0i = ((AliESDEvent*)fESD)->GetV0(iV0);
        AliESDv0 *v0j = ((AliESDEvent*)fESD)->GetV0(jV0);
        Double_t tV0momi[3], tV0momj[3], tV0mom_result[3];
        v0i->GetPxPyPz(tV0momi[0], tV0momi[1], tV0momi[2]);
        v0j->GetPxPyPz(tV0momj[0], tV0momj[1], tV0momj[2]);
        //// ---- Calculate inv. mass for Xi_c ---- ////
        Double_t ei = 0.;
        Double_t ej = 0.;
        Double_t angle = 0.;
        Double_t fMass = 0.;
        Double_t fPt_result = 0.;

        ei = getEnergy(k0Mass, tV0momi[0], tV0momi[1], tV0momi[2]); // Energy of first particle(K0Short)
        ej = getEnergy(l0Mass, tV0momj[0], tV0momj[1], tV0momj[2]); // Energy of first particle(K0Short)

        angle = getAngle(tV0momi[0], tV0momi[1], tV0momi[2], tV0momj[0], tV0momj[1], tV0momj[2]);

        cout << "ei: " << ei << ", ej: " << ej << ", angle: "<< angle << endl;

        fMass = k0Mass * k0Mass + l0Mass * l0Mass + 2.*ei * ej - 2.*angle;
        if (fMass <= 0) continue;
        fMass = sqrt(fMass);

        tV0mom_result[0] = tV0momi[0] + tV0momj[0];
        tV0mom_result[1] = tV0momi[1] + tV0momj[1];
        tV0mom_result[2] = tV0momi[2] + tV0momj[2];
        fPt_result = sqrt(pow(tV0mom_result[0], 2) + pow(tV0mom_result[1], 2));

        ((TH2F*)fOutputList2->FindObject("hInvMassWithPt"))->Fill(fMass, fPt_result); // with Pt
        ((TH1F*)fOutputList2->FindObject("hInvMass"))->Fill(fMass); // Cumulated
      }
    }*/
  }
    PostData(1, fOutputList);                           // stream the results the analysis of this event to
    PostData(2, fOutputList2);
}
//_____________________________________________________________________________
void AliAnalysisTaskXic::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________

double getEnergy(Double_t trueMass, Double_t Px, Double_t Py, Double_t Pz)
{
    return TMath::Sqrt(trueMass * trueMass + Px * Px + Py * Py + Pz * Pz);
}
double getAngle(Double_t Px1, Double_t Py1, Double_t Pz1, Double_t Px2, Double_t Py2, Double_t Pz2)
{
    return Px1 * Px2 + Py1 * Py2 + Pz1 * Pz2;
}
void CheckChargeV0(AliESDv0 *v0)
{
    // This function is copied from PWGLF/STRANGENESS/LambdaK0/AliAnalysisTaskExtractV0.cxx
    // This function checks charge of negative and positive daughter tracks.
    // If incorrectly defined (onfly vertexer), swaps out.
    if ( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
        //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
        Long_t lCorrectNidx = v0->GetPindex();
        Long_t lCorrectPidx = v0->GetNindex();
        Double32_t    lCorrectNmom[3];
        Double32_t    lCorrectPmom[3];
        v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
        v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );
        AliExternalTrackParam lCorrectParamN(
            v0->GetParamP()->GetX() ,
            v0->GetParamP()->GetAlpha() ,
            v0->GetParamP()->GetParameter() ,
            v0->GetParamP()->GetCovariance()
        );
        AliExternalTrackParam lCorrectParamP(
            v0->GetParamN()->GetX() ,
            v0->GetParamN()->GetAlpha() ,
            v0->GetParamN()->GetParameter() ,
            v0->GetParamN()->GetCovariance()
        );
        lCorrectParamN.SetMostProbablePt( v0->GetParamP()->GetMostProbablePt() );
        lCorrectParamP.SetMostProbablePt( v0->GetParamN()->GetMostProbablePt() );
        //Get Variables___________________________________________________
        Double_t lDcaV0Daughters = v0 -> GetDcaV0Daughters();
        Double_t lCosPALocal     = v0 -> GetV0CosineOfPointingAngle();
        Bool_t lOnFlyStatusLocal = v0 -> GetOnFlyStatus();
        //Create Replacement Object_______________________________________
        AliESDv0 *v0correct = new AliESDv0(lCorrectParamN, lCorrectNidx, lCorrectParamP, lCorrectPidx);
        v0correct->SetDcaV0Daughters          ( lDcaV0Daughters   );
        v0correct->SetV0CosineOfPointingAngle ( lCosPALocal       );
        v0correct->ChangeMassHypothesis       ( kK0Short          );
        v0correct->SetOnFlyStatus             ( lOnFlyStatusLocal );
        //Reverse Cluster info..._________________________________________
        v0correct->SetClusters( v0->GetClusters( 1 ), v0->GetClusters ( 0 ) );
        *v0 = *v0correct;
        //Proper cleanup..._______________________________________________
        v0correct->Delete();
        v0correct = 0x0;
        //Just another cross-check and output_____________________________
        if ( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
            //AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
        } else {
            //AliWarning("Found Swapped Charges and fixed.");
        }
        //________________________________________________________________
    } else {
        //Don't touch it! ---
        //Printf("Ah, nice. Charges are already ordered...");
    }
    return;
  }
/*
    for (Int_t jV0 = 0; jV0 < nv0s; jV0++){
        v0checkk0s[jV0] = 0;
        bool kshortCandidate = true;
        // keep only events of interest for fHistMLa plots
        // from PWGLF/STRANGENESS/LambdaK0PbPb/AliAnalysisTaskLukeV0.cxx

        AliESDv0 *v0j = ((AliESDEvent*)fESD)->GetV0(jV0);
        if (!v0j) continue;
        if(debugmode > 51) ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(1);
        if(debugmode > 100) AliInfo("01");
        //---> Fix On-the-Fly candidates, count how many swapped
        if ( v0j->GetParamN()->Charge() > 0 && v0j->GetParamP()->Charge() < 0 ) {
            fHistSwappedV0Counter -> Fill( 1 );
        } else {
            fHistSwappedV0Counter -> Fill( 0 );
        }
        if ( fkUseOnTheFly ) CheckChargeV0(v0j);

        Int_t    lOnFlyStatus = 0; // nv0sOn = 0, nv0sOff = 0;
        lOnFlyStatus = v0j->GetOnFlyStatus();
        if (lOnFlyStatus == 0) continue;
        if(debugmode > 100) AliInfo("02");
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(3);

        //// Get V0 informations for the cuts
        Double_t lPt = 0;
        lPt = v0j->Pt();
        if (v0j->GetEffMass(2,2) < 0.414 || v0j->GetEffMass(2,2) > 0.582 || TMath::Abs(v0j->Y(310))>0.5) kshortCandidate = false;

        if(debugmode > 50 && kshortCandidate) AliInfo("K0s Case");
        // get daughter particle
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0j->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0j->GetNindex());
        AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
        AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
        const AliExternalTrackParam * paramPosl = v0j->GetParamP();
        const AliExternalTrackParam * paramNegl = v0j->GetParamN();
        // Cosine Pointing Angle and DCA Values
        Double_t kV0cosPointAngle = v0j->GetV0CosineOfPointingAngle(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
        Double_t kV0Position[3];
        Double_t kV0Radius = 0;
        v0j->GetXYZ(kV0Position[0], kV0Position[1], kV0Position[2]);
        kV0Radius = TMath::Sqrt(kV0Position[0] * kV0Position[0] + kV0Position[1] * kV0Position[1]);
        Double_t lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(primaryVtx[0], primaryVtx[1], lMagneticField));
        Double_t lDcaNegToPrimVertex = TMath::Abs(pTrack->GetD(primaryVtx[0], primaryVtx[1], lMagneticField));
        Double_t lDcaV0Daughters = v0j->GetDcaV0Daughters();
        Double_t tV0momi[3], tV0momj[3], tV0mom_result[3];
        v0j->GetPxPyPz(tV0momi[0], tV0momi[1], tV0momi[2]);
        // Decay length
        double decayLength2 = (sqrt((primaryVtx[0]-primaryVtx[0])*(primaryVtx[0]-primaryVtx[0])+(primaryVtx[1]-primaryVtx[1])*(primaryVtx[1]-primaryVtx[1])+(primaryVtx[2]-primaryVtx[2])*(primaryVtx[2]-primaryVtx[2])));
        double cTauK0 = decayLength2*(v0j->GetEffMass(2,2))/(v0j->P());
        // momentums
        double pTrackMomentum[3];
        double nTrackMomentum[3];
        pTrack->GetConstrainedPxPyPz(pTrackMomentum);
            nTrack->GetConstrainedPxPyPz(nTrackMomentum);
        double pPos2 = sqrt(pTrackMomentum[0]*pTrackMomentum[0]+pTrackMomentum[1]*pTrackMomentum[1]+pTrackMomentum[2]*pTrackMomentum[2]);
        double pNeg2 = sqrt(nTrackMomentum[0]*nTrackMomentum[0]+nTrackMomentum[1]*nTrackMomentum[1]+nTrackMomentum[2]*nTrackMomentum[2]);

        if(debugmode > 100) AliInfo("03");
        //// Cuts
        if(!(fTrackCuts->IsSelected(pTrack)) || !(fTrackCuts->IsSelected(nTrack))) {
          kshortCandidate = false;
        }
        if(kshortCandidate == false) continue;
        if(!(kshortCandidate == false) && debugmode > 50) AliInfo("Track cut!");
        // Pt cut for mother particle
        if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;
        if(debugmode > 51) ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(5);
        // is daughter particle okay?
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(7);
        // Like sign cut
        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(9);
        //psuedorapidity cut
            if(cutEta != -999) {
                if(TMath::Abs(pTrack->Eta()) > cutEta || TMath::Abs(nTrack->Eta())  >cutEta) {
                    kshortCandidate = false;

                }
            }
        if(kshortCandidate == false) continue;
        if(!(kshortCandidate == false) && debugmode > 50) AliInfo("Eta cut!"); ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(11);
        // CPA cut
        if(cutCosPa != -999) {
          if (kV0cosPointAngle < cutCosPa){
            kshortCandidate = false;
          }
        }
        if(kshortCandidate == false) continue;
        if(!(kshortCandidate == false) && debugmode > 50) AliInfo("CPA cut!"); ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(13);
        // DCA between daughterscut
            if(cutDCA != -999) {
          if(v0j->GetDcaV0Daughters() > cutDCA) {
            kshortCandidate = false;

          }
        }
        if(kshortCandidate == false) continue;
        if(!(kshortCandidate == false) && debugmode > 50) AliInfo("DCA cut!"); ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(15);
        //lifetime cut
            if(cutcTau != -999){
                if(cTauLa < cutcTau){
                    kshortCandidate = false;
                }
                if(cTauLb < cutcTau){

                }
            }
        if(kshortCandidate == false) continue;
        if(!(kshortCandidate == false) && debugmode > 50) AliInfo("Lifetime cut!"); ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(17);
        // Bethe Bloch cut. Made sightly complicated as options for crude cuts still included. Should probably reduce to just 'official' cuts
        if(cutBetheBloch != -999) {
            if(pTrack->GetTPCsignal() <0 || nTrack->GetTPCsignal()<0) continue;

            if(cutBetheBloch > 0) {
                    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(nTrack, AliPID::kPion)) > cutBetheBloch ) kshortCandidate = false;
                    if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(pTrack, AliPID::kPion)) > cutBetheBloch ) kshortCandidate = false;
                }

            if(cutBetheBloch == -4) {
                double par0 = 0.20;
                double par1 = 4.2;
                double par2 = 1000000;
                double beta2 = TMath::Power((pNeg2/TMath::Sqrt((pNeg2*pNeg2+par0*par0))),2);
                    double gamma2 = 1.0/(1.0-beta2);
                    if(nTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pNeg2) > -0.6) kshortCandidate = false;

                    beta2 = TMath::Power((pPos2/TMath::Sqrt((pPos2*pPos2+par0*par0))),2);
                    gamma2 = 1.0/(1.0-beta2);
                    if(pTrack->GetTPCsignal() > (par1/beta2)*(TMath::Log(par2*beta2*gamma2)-beta2) && TMath::Log10(pPos2) > -0.6) kshortCandidate = false;
            }
        }
        if(kshortCandidate == false) continue;
        if(!(kshortCandidate == false) && debugmode > 50) AliInfo("PID cut!"); ((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(19);
        //if ((lDcaPosToPrimVertex < 0.1) || (lDcaNegToPrimVertex < 0.1) || (kV0cosPointAngle < 0.998) || (kV0Radius < 0.0) || (kV0Radius > 1000) ) continue;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(23);
        if(debugmode > 100) AliInfo("04");

        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(25);
        //remove all non-candidates
        if(kshortCandidate == false) continue;
        if(debugmode > 10) AliInfo("============v0 survived!============");
        if(pPos2 == 0 || pNeg2 ==0) continue;

        ((TH1F*)fOutputList->FindObject("fHistCosPA_k0s"))->Fill(kV0cosPointAngle);
        if(debugmode > 10) AliInfo("============fHistCosPA_k0s!============");
        ((TH1F*)fOutputList->FindObject("fHistDecayL_k0s"))->Fill(decayLength2);
        if(debugmode > 10) AliInfo("============fHistDecayL_k0s!============");
        //((TH1F*)fOutputList->FindObject("fHistTauk0s"))->Fill(cTauK0);
        if(debugmode > 10) AliInfo("============fHistTauk0s!============");
        ((TH2F*)fOutputList->FindObject("fHistBetheBlochTPCPos_k0s"))->Fill(TMath::Log10(pPos2),pTrack->GetTPCsignal());
        if(debugmode > 10) AliInfo("============fHistBetheBlochTPCPos_k0s!============");
        ((TH2F*)fOutputList->FindObject("fHistBetheBlochTPCNeg_k0s"))->Fill(TMath::Log10(pNeg2),nTrack->GetTPCsignal());
        if(debugmode > 10) AliInfo("============k0 qa histograms============");
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(27);
        // Mass Hypothesis for Lambda
        //v0j->ChangeMassHypothesis(3122);
        //sets assumed particle type of pos/neg daughters.
            // 0 = electron, 1 = Muon, 2 = pion, 3 = kaon, 4 = proton.
            int dPos = 2;
            int dNeg = 2;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(29);
        double lInvMassK0s = 0.;
        lInvMassK0s = v0j->GetEffMass(dPos,dNeg);
        if(debugmode > 100) AliInfo("04-1");
        if (pTrack->GetMass() > 0.2 || nTrack->GetMass() > 0.2) continue;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(31);
        if(debugmode > 50) AliInfo("daughter mass cut pass");
        ((TH2F*)fOutputList->FindObject("fArmPod_kaon"))->Fill(v0j->AlphaV0(),v0j->PtArmV0());
        if(debugmode > 10) AliInfo("============ArmPod Kaon============");
        if(debugmode > 100) AliInfo("04-5");
        // Armenteros-Podolansiki Cut
        if (TMath::Abs(0.2 * v0j->AlphaV0()) > v0j->PtArmV0()) continue;
        if(debugmode > 51)((TH1F*)fOutputList->FindObject("hNofV0_2"))->Fill(33);
        ((TH2F*)fOutputList->FindObject("fArmPod_kaon_cut"))->Fill(v0j->AlphaV0(),v0j->PtArmV0());
        if(debugmode > 10) AliInfo("============ArmPod Kaon 2!============");
        if(debugmode > 100) AliInfo("05");

        if(kshortCandidate) v0checkk0s[jV0] = 1;
        //cout << "K0s : " << v0checkk0s[jV0] << endl;
        //if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 0 && fkUseOnTheFly == kTRUE ) ){
        ((TH1F*)fOutputList->FindObject("fInvK0Short"))->Fill(lInvMassK0s);
        //if (lInvMassK0s > l0Mass + 0.0008 || lInvMassK0s < l0Mass - 0.008) continue; // Mass window
        ((TH1F*)fOutputList->FindObject("fInvK0ShortCut"))->Fill(lInvMassK0s); // After Cut
        if(debugmode > 10) AliInfo("============fill invmass!============");
        //}
    }*/
