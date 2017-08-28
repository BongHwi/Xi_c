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

/* AliAnaysisTaskXicMC
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

#include "AliMCEventHandler.h" // MC include "AliMCEvent.h"        // MC include
#"AliMCParticle.h"     // MC include "AliStack.h"          // MC include
#"AliAODMCParticle.h"  // MC

#include "AliESDv0.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliAnalysisTaskXic.h"

#define PI 3.1415927

class AliAnalysisTaskXic;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskXic) // classimp: necessary for root

AliAnalysisTaskXic::AliAnalysisTaskXicMC() : AliAnalysisTaskSE(),
    fESD(0x0),
    fOutputList(0x0),
    fTrackCuts(0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0),
    fMCcase(0),
    fAODcase(0),
    fEventCounter(0),
    fCutList(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskXic::AliAnalysisTaskXicMC(const char* name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption) : AliAnalysisTaskSE(name),
    fESD(0x0),
    fOutputList(0x0),
    fTrackCuts(0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0),
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
AliAnalysisTaskXic::~AliAnalysisTaskXicMC()
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
    TH1F *fHistCosPA_lam = new  TH1F("fHistCosPA_lam", "Cosine of Pointing Angle of V0s; Cos PA; N(v0s)",202,0.8,1.01);
    TH1F *fHistDCAV0Daughters_lam = new TH1F("fHistDCAV0Daughters_lam", "DCA between V0 daughters; DCA (cm); N V0s", 100, 0, 2);
    TH1F *fHistDecayL_lam = new TH1F("fHistDecayL_lam", "Distance between V0 and PV; Distance(cm); N(v0s)",200,-0.1,30);
    TH1F *fHistTauLa = new  TH1F("fHistTauLa", "Lifetime under Lambda mass hypothesis; Lifetime(s); N(v0s)",200,0,100);
    TH2F *fHistBetheBlochTPCNeg_lam = new   TH2F("fHistBetheBlochTPCNeg_lam","-dE/dX against Momentum for negative daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
    TH2F *fHistBetheBlochTPCPos_lam = new   TH2F("fHistBetheBlochTPCPos_lam","-dE/dX against Momentum for positive daughter from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);
    // Pion QA
    TH1F *fNTPCcls_Pi = new TH1F("fNTPCcls_Pi","fNTPCcls_Pi",200,0,200);
      fNTPCcls_Pi->GetXaxis()->SetTitle("fNTPCcls_Pi");
    TH2F *fHistBetheBlochTPC_Pi = new   TH2F("fHistBetheBlochTPC_Pi","-dE/dX against Momentum from TPC; Log10 P (GeV); -dE/dx (keV/cm ?)",1000,-1,1,1000,0,200);

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
        //mcEvent = MCEvent();
        //if (!mcEvent) {cout<<"ERROR: Could not retrieve MC event"<<endl; return;
        if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
            if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()) mcstack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
            if (!mcstack) {cout<<"ERROR: Could not retrieve the stack"<<endl; return;}
        }
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
          TParticle *mcInputTrack = ((AliMCParticle*)mcstack->Particle(it);
          if (!mcInputTrack) {
        	  Error("UserExec", "Could not receive track %d", it);
        	  continue;
        	}
        }

        // Xi
      	if(mcInputTrack->GetPdgCode() == +kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXi1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
      	if(mcInputTrack->GetPdgCode() == -kXiCode) ((TH3F*)fOutputList->FindObject("fMCinputTotalXibar1"))->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
      }
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
