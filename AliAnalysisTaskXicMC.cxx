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
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TList.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliLog.h"

#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliCentrality.h"

#include "AliMCEventHandler.h" // MC
#include "AliMCEvent.h"        // MC
#include "AliMCParticle.h"     // MC
#include "AliStack.h"          // MC

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
TString isMC = kTURE;

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
    fIsMC(kFALSE)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskXic::AliAnalysisTaskXic(const char* name) : AliAnalysisTaskSE(name),
    fESD(0x0),
    fOutputList(0x0),
    fTrackCuts(0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0),
    fIsMC(kFALSE)
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
    if(isMC){
      TH1F *fMultDistMC = new TH1F("fMultDistMC", "Multiplicity Distribution of MC", 200, 0, 20000);
        fMultDistMC->GetXaxis()->SetTitle("Multiplicity");
    }

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

    if(isMC){
      fOutputList->Add(fMultDistMC);
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

    if(isMC) {
      AliMCEvent* fMCE = MCEvent(); // fMCE: MC event
      if (!fMCE) {
        Printf("ERROR: Could not retrieve MC event");
        return;
      }
      AliStack* stack = fMCE->Stack();
      if(!stack){
        Printf("ERROR: Can't load stack from MC Event");
        return;
      }
    }
    else{
      fESD = dynamic_cast<AliESDEvent*>(InputEvent());
      if (!fESD) {Printf("ERROR: fESD not available"); return;}
    }
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
    //Step 5: MC Loop
    //------------------------------------------------

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
