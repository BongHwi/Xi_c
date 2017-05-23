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
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
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

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliAnalysisTaskXic.h"

class AliAnalysisTaskXic;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskXic) // classimp: necessary for root
double getEnergy(Double_t trueMass, Double_t Px, Double_t Py, Double_t Pz);
double getAngle(Double_t Px1, Double_t Py1, Double_t Pz1, Double_t Px2, Double_t Py2, Double_t Pz2);
void GetArPod(Double_t pos[3], Double_t neg[3], Double_t moth[3],  Double_t arpod[2]);
void CheckChargeV0(AliESDv0 *v0);

AliAnalysisTaskXic::AliAnalysisTaskXic() : AliAnalysisTaskSE(),
    fESD(0x0),
    fOutputList(0x0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskXic::AliAnalysisTaskXic(const char* name) : AliAnalysisTaskSE(name),
    fESD(0x0),
    fOutputList(0x0),
    fOutputList2(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0),
    fHistSwappedV0Counter(0)
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

    // example of a histogram
    TH1F *fMultDist = new TH1F("fMultDist", "Multiplicity Distribution", 200, 0, 20000);
    fMultDist->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist);

    TH2F *fArmPod_kaon = new TH2F("fArmPod_kaon", "Armenteros-Podolanski Plot", 800,-1.0, 1.0, 100, 0, 0.25);
    fOutputList->Add(fArmPod_kaon);
    TH2F *fArmPod_lambda = new TH2F("fArmPod_lambda", "Armenteros-Podolanski Plot", 800,-1.0, 1.0, 100, 0, 0.25);
    fOutputList->Add(fArmPod_lambda);


    TH1F *fInvLambda_before = new TH1F("fInvLambda_before", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
    fInvLambda_before->GetXaxis()->SetTitle("fInvLambda_before");
    fOutputList->Add(fInvLambda_before);

    TH1F *fInvLambda = new TH1F("fInvLambda", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
    fInvLambda->GetXaxis()->SetTitle("fInvLambda");
    fOutputList->Add(fInvLambda);
    TH1F *fInvK0Short = new TH1F("fInvK0Short", "Invariant mass distribution of K0s", 400, 0.3, 0.7);
    fInvK0Short->GetXaxis()->SetTitle("fInvK0Short");
    fOutputList->Add(fInvK0Short);
    TH1F *fInvK0Short_beforePID = new TH1F("fInvK0Short_beforePID", "Invariant mass distribution of K0s", 400, 0.3, 0.7);
    fInvK0Short->GetXaxis()->SetTitle("fInvK0Short_beforePID");
    fOutputList->Add(fInvK0Short_beforePID);
    TH1F *fInvLambda_beforePID = new TH1F("fInvLambda_beforePID", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
    fInvLambda_beforePID->GetXaxis()->SetTitle("fInvLambda_beforePID");
    fOutputList->Add(fInvLambda_beforePID);

    TH1F *fInvLambdaCut = new TH1F("fInvLambdaCut", "Invariant mass distribution of Lambda after mass window cut", 400, 1.0, 1.2);
    fInvLambdaCut->GetXaxis()->SetTitle("fInvLambdaCut");
    fOutputList->Add(fInvLambdaCut);
    TH1F *fInvK0ShortCut = new TH1F("fInvK0ShortCut", "Invariant mass distribution of K0s after mass window cut", 400, 0.3, 0.7);
    fInvK0Short->GetXaxis()->SetTitle("fInvK0ShortCut");
    fOutputList->Add(fInvK0ShortCut);
    
    if(! fHistSwappedV0Counter) {
      fHistSwappedV0Counter = new TH1F("fHistSwappedV0Counter", 
         "Swap or not histo;Swapped (1) or not (0); count", 
         2, 0, 2); 		
      fOutputList->Add(fHistSwappedV0Counter);
    }

    TH1F *hEventSelecInfo = new TH1F("hEventSelecInfo", "hEventSelecInfo", 10, 0, 10);
    fOutputList->Add(hEventSelecInfo);
    hEventSelecInfo->GetXaxis()->SetBinLabel(1, "NONE");
    hEventSelecInfo->GetXaxis()->SetBinLabel(2, "kMB");
    hEventSelecInfo->GetXaxis()->SetBinLabel(3, "kCentral");
    hEventSelecInfo->GetXaxis()->SetBinLabel(4, "kSemiCentral");
    hEventSelecInfo->GetXaxis()->SetBinLabel(5, "kINT7");
    hEventSelecInfo->GetXaxis()->SetBinLabel(6, "kAny");
    hEventSelecInfo->GetXaxis()->SetBinLabel(7, "kPhysicsALL");


    TH1F *hCentrality = new TH1F("hCentrality", "Centrality", 100, 0, 100);
    hCentrality->GetXaxis()->SetTitle("Centrality");
    fOutputList->Add(hCentrality);


    TH3F *fVertexDistXYZ = new TH3F("fVertexDistXYZ", "Vertex Distribution", 20, -1, 1, 20, -1, 1, 60, -30, 30);
    fVertexDistXYZ->GetXaxis()->SetTitle("X Vertex (cm)");
    fVertexDistXYZ->GetYaxis()->SetTitle("Y Vertex (cm)");
    fVertexDistXYZ->GetZaxis()->SetTitle("Z Vertex (cm)");
    fOutputList->Add(fVertexDistXYZ);

    TH2F *hInvMassWithPt = new TH2F("hInvMassWithPt", "Invariant mass distribution vs Pt", 1000, 2.0, 3.0, 100, 0, 10);
    fOutputList2->Add(hInvMassWithPt);
    TH1F *hInvMass = new TH1F("hInvMass", "Invariant mass distribution", 1000, 2.0, 3.0);
    fOutputList2->Add(hInvMass);

    TH2F *hTPCPID_K0s = new TH2F("hTPCPID_K0s","PID via TPC",500,0,20,500,0,200);
    fOutputList->Add(hTPCPID_K0s);
    TH2F *hTPCPID_K0s_after = new TH2F("hTPCPID_K0s_after","PID via TPC",500,0,20,500,0,200);
    fOutputList->Add(hTPCPID_K0s_after);
    TH2F *hTPCPID_lam = new TH2F("hTPCPID_lam","PID via TPC",500,0,20,500,0,200);
    fOutputList->Add(hTPCPID_lam);
    TH2F *hTPCPID_lam_after = new TH2F("hTPCPID_lam_after","PID via TPC",500,0,20,500,0,200);
    fOutputList->Add(hTPCPID_lam_after);

    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    if (!fPIDResponse) AliError ("No PID");



    PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
    // fOutputList object. the manager will in the end take care of writing your output to file
    // so it needs to know what's in the output
    PostData(2, fOutputList2);
}
//_____________________________________________________________________________
void AliAnalysisTaskXic::UserExec(Option_t *)
{
    // user exec
    // this function is called once for each event
    // the manager will take care of reading the events from file, and with the static function InputEvent() you
    // have access to the current event.
    // once you return from the UserExec function, the manager will retrieve the next event from the chain
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());    // get an event (called fESD) from the input file
    // there's another event format (ESD) which works in a similar wya
    // but is more cpu/memory unfriendly. for now, we'll stick with ESD's
    if (!fESD) {Printf("ERROR: fESD not available"); return;} // if the pointer to the event is empty (getting it failed) skip this event
    ((TH1F*)fOutputList->FindObject("fMultDist"))->Fill(fESD->GetNumberOfTracks()); // # of tracks --> multiplicity distribution


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

    if (!isSelected)Printf("There is events in kANY");
    ////////////******* Do Event selecction *******////////////
    if (!isSelectedINT7) {cout << "Event Rejected" << endl; return;}


    //------------------------------------------------
    //Step 2: Check for centrality for Pb-Pb
    //------------------------------------------------
    Float_t  centralityV0M = -100;
    fCentrality = fESD->GetCentrality();
    centralityV0M = fCentrality->GetCentralityPercentile("V0M");
    ((TH1F*)fOutputList->FindObject("hCentrality"))->Fill(centralityV0M);


    // Example of GetStandardITSTPCTrackCuts
    fTrackCut = new AliESDtrackCuts();
    fTrackCut->SetPtRange(.15, 100000000); // min. pT cut //1.
    fTrackCut->SetAcceptKinkDaughters(kFALSE); //2.
    fTrackCut->SetRequireTPCRefit(kTRUE);     //3.
    fTrackCut->SetMaxChi2PerClusterTPC(4);    //4.
    fTrackCut->SetMinNClustersTPC(70);        //5. TPC nclus will be changed in order to performed systematic.
    //fTrackCut->SetEtaRange(-0.8, +0.8);       //6.
    //fTrackCut->SetMaxDCAToVertexZ(2);           //7. DCAz cut
    //fTrackCut->SetMaxDCAToVertexXYPtDep("1.");  //8. DCAr cut


    //------------------------------------------------
    //Step 3: Check for Vertex-Z position
    //------------------------------------------------
    const AliESDVertex *PrimaryVertexESD;
    PrimaryVertexESD = fESD->GetPrimaryVertex(); // call primary vertex from ESD
    //  if(!PrimaryVertexESD) return;
// if(PrimaryVertexESD->GetNContributors() < 1) return;

    Double_t primaryVtx[3] = {0};
    primaryVtx[0] = PrimaryVertexESD->GetX(); // call primary vertex position of X
    primaryVtx[1] = PrimaryVertexESD->GetY(); // call primary vertex position of Y
    primaryVtx[2] = PrimaryVertexESD->GetZ(); // call primary vertex position of Z
    ((TH3F*)fOutputList->FindObject("fVertexDistXYZ"))->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);
    ////////////******* DO Vertex-Z selecction *******////////////
    if (fabs(primaryVtx[2]) > 10.) return;


    //------------------------------------------------
    //Step 4: Check for SPD Pileup
    //------------------------------------------------
    if (fESD->IsPileupFromSPD()) return; // Reject Pile-up events

    double lInvMassLambda = 0.;
    double lInvMassK0Short = 0.0;

    static Double_t k0Mass = 0.497611;      static Double_t l0Mass = 1.115683;
    static Double_t piMass = 0.13957;       static Double_t protonMass = 0.93827;

    Double_t fMinV0Pt = 0.15;  Double_t fMaxV0Pt = 1.E10;    Double_t lV0Radius = 0, lPt = 0;
    Double_t fQt = 0;       Double_t falpha = 0;
    Double_t thetaV0pos, thetaV0neg;
    TVector3 momentumVectorPositiveKF, momentumVectorNegativeKF, vecV0;

    Double_t V0momK0ShortDaughterN[3] = {.0, .0, .0};
    Double_t V0momK0ShortDaughterP[3] = {.0, .0, .0};
    Double_t V0momLambda0DaughterN[3] = {.0, .0, .0};
    Double_t V0momLambda0DaughterP[3] = {.0, .0, .0};

    Double_t tV0momi[3];    Double_t tV0momj[3];    Double_t tV0mom_result[3];

    Double_t ei = 0.;       Double_t ej = 0.;
    Double_t angle = 0.;    Double_t fPt_result = 0.;
    Double_t  fMass = 0.;
    Bool_t fkUseOnTheFly = kFALSE;
    Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;

    // loop for Lambda
    Int_t nv0s = 0;
    nv0s = fESD->GetNumberOfV0s();
    /* for test
    //for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
    {   // This is the begining of the V0 loop for first V0(K0Short)
        AliESDv0 *v0i = ((AliESDEvent*)fESD)->GetV0(iV0);
        if (!v0i) continue;

        lPt = v0i->Pt();
        if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;

        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0i->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0i->GetNindex());

        AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
        AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
        const AliExternalTrackParam * paramPos = v0i->GetParamP();
	const AliExternalTrackParam * paramNeg = v0i->GetParamN();

	if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
	// Draw Armenteros-Podolanski Plot
	// from PWGGA/Hyperon/AliAnalysisTaskSigma0.cxx by Alexander Borissov.
	
	// Pion -> pi+ + pi-  ---------------
    	AliKFParticle negKFKpim(*paramNeg,211);
    	AliKFParticle posKFKprot(*paramPos,211);
    	AliKFParticle kaonKF(negKFKpim,posKFKprot);

	Double_t posp[3]= { pTrack->Px(),  pTrack->Py(),  pTrack->Pz() };
	Double_t negp[3]= { nTrack->Px(),  nTrack->Py(),  nTrack->Pz() };
	Double_t moth[3]= { kaonKF.GetPx(), kaonKF.GetPy(), kaonKF.GetPz() };
	Double_t arpod[2]= {0,0};
	GetArPod( posp, negp, moth, arpod );

	((TH2F*)fOutputList->FindObject("fArmPod_kaon"))->Fill(arpod[1],arpod[0]);

        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }

        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if ( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if ( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

        if ( ( ( ( pTrack->GetTPCClusterInfo(2, 1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2, 1) ) < 70 ) )) continue;

        //GetKinkIndex condition
        if ( pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0 ) continue;

        //Findable clusters > 0 condition
        //if ( pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0 ) continue;
	//Float_t nsigmapiP = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
	//Float_t nsigmapiN = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
	
	//Float_t fTPCPIDmom = pTrack->GetTPCmomentum();
        //Float_t sigTPC = pTrack->GetTPCsignal();
        //Float_t nsigpip= fabs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion));
	//Float_t nsigpin= fabs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kPion));

	//((TH2F*)fOutputList->FindObject("hTPCPID_K0s"))->Fill(fTPCPIDmom,sigTPC);
	//if (nsigpip > 3.0 && nsigpin > 3.0 ) continue;
	//((TH2F*)fOutputList->FindObject("hTPCPID_K0s_after"))->Fill(fTPCPIDmom,sigTPC);
	
        // find new v0 for K0Short
        v0i->ChangeMassHypothesis(310); //kK0Short
        v0i->GetPxPyPz(tV0momi[0], tV0momi[1], tV0momi[2]);
        lInvMassK0Short = v0i->GetEffMass();
        ((TH1F*)fOutputList->FindObject("fInvK0Short_beforePID"))->Fill(lInvMassK0Short); // Before PID
	v0i->ChangeMassHypothesis(3122); //kLambda0
        lInvMassLambda = v0i->GetEffMass();
       	((TH1F*)fOutputList->FindObject("fInvLambda_before"))->Fill(lInvMassLambda); // Before Cut
	if (nsigmapiP > 3.0 || nsigmapiN > 3.0) continue;
	((TH1F*)fOutputList->FindObject("fInvK0Short"))->Fill(lInvMassK0Short); // Before Cut
	if (lInvMassK0Short > k0Mass + 0.008 || lInvMassK0Short < k0Mass - 0.008) continue; // Mass window
        ((TH1F*)fOutputList->FindObject("fInvK0ShortCut"))->Fill(lInvMassK0Short); // After Cut

        for (Int_t jV0 = iV0; jV0 < nv0s; jV0++)
        {   // This is the begining of the V0 loop for second V0(Lambda)
            AliESDv0 *v0j = ((AliESDEvent*)fESD)->GetV0(jV0);
              if (!v0j) continue;
    
            lPt = v0j->Pt();
            if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;
    
            UInt_t lKeyPos = (UInt_t)TMath::Abs(v0j->GetPindex());
            UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0j->GetNindex());
    
    
            AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
            AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
	    const AliExternalTrackParam * paramPosl = v0j->GetParamP();
	    const AliExternalTrackParam * paramNegl = v0j->GetParamN();
            if (!pTrack || !nTrack) {
                Printf("ERROR: Could not retreive one of the daughter track");
                continue;
            }
	    // Draw Armenteros-Podolanski Plot
            // from PWGGA/Hyperon/AliAnalysisTaskSigma0.cxx by Alexander Borissov.
	
            // Lambda -> P+ pi-  ---------------
            AliKFParticle negKFKpim(*paramNegl,211);
            AliKFParticle posKFKprot(*paramPosl,2212);
            AliKFParticle lamKF(negKFKpim,posKFKprot);    
	    
	    if (pTrack->GetMass() > 0.5){
	    //printf("this v0 is antilambda");
	    AliKFParticle negKFKpim(*paramNegl,2212);
            AliKFParticle posKFKprot(*paramPosl,211);
            AliKFParticle lamKF(negKFKpim,posKFKprot);
	    }

	    Double_t posp[3]= { pTrack->Px(),  pTrack->Py(),  pTrack->Pz() };
            Double_t negp[3]= { nTrack->Px(),  nTrack->Py(),  nTrack->Pz() };
            Double_t moth[3]= { lamKF.GetPx(), lamKF.GetPy(), lamKF.GetPz() };
            Double_t arpod[2]= {0,0};
            GetArPod( posp, negp, moth, arpod );
	    
            ((TH2F*)fOutputList->FindObject("fArmPod_lambda"))->Fill(arpod[1],arpod[0]);
	    
            if ( pTrack->GetSign() == nTrack->GetSign()) {
                continue;
            }
    
            // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
            if ( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            if ( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
    
            if ( ( ( ( pTrack->GetTPCClusterInfo(2, 1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2, 1) ) < 70 ) )) continue;
    
            //GetKinkIndex condition
            if ( pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0 ) continue;
    
            //Findable clusters > 0 condition
            if ( pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0 ) continue;
     	    if ( pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0 ) continue;

            Float_t nsigmaprP = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
            Float_t nsigmapiN = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
            if (pTrack->GetMass() > 0.5) if (nsigmaprP > 3.0 || nsigmapiN > 3.0) continue;
	    if (pTrack->GetMass() < 0.5){
	    Float_t nsigmaprN = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
    	    Float_t nsigmapiP = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
	    if (nsigmaprN > 3.0 || nsigmapiP > 3.0) continue;
	    }
	    
	    //Float_t fTPCPIDmom = pTrack->GetTPCmomentum();
            //Float_t sigTPC = pTrack->GetTPCsignal();
            //Float_t nsigpi= fabs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kPion));
            //Float_t nsigk= fabs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
            //Float_t nsigprP = fabs(fPIDResponse->NumberOfSigmasTPC(pTrack,AliPID::kProton));
	    //Float_t nsigprN = fabs(fPIDResponse->NumberOfSigmasTPC(nTrack,AliPID::kProton));
            //((TH2F*)fOutputList->FindObject("hTPCPID_lam"))->Fill(fTPCPIDmom,sigTPC);
            //if (nsigprP > 3.0 || nsigprN > 3.0) continue;
            //((TH2F*)fOutputList->FindObject("hTPCPID_lam_after"))->Fill(fTPCPIDmom,sigTPC);
	    
            // find new v0 for Lambda0
            v0j->ChangeMassHypothesis(3122); //kLambda0
            v0j->GetPxPyPz(tV0momj[0], tV0momj[1], tV0momj[2]);
            lInvMassLambda = v0j->GetEffMass();
            ((TH1F*)fOutputList->FindObject("fInvLambda"))->Fill(lInvMassLambda); // Before Cut
            if (lInvMassLambda > l0Mass + 0.0008 || lInvMassLambda < l0Mass - 0.008) continue; // Mass window
            ((TH1F*)fOutputList->FindObject("fInvLambdaCut"))->Fill(lInvMassLambda); // After Cut

            //// ---- Calculate inv. mass for Xi_c ---- ////
            ei = getEnergy(k0Mass, tV0momi[0], tV0momi[1], tV0momi[2]); // Energy of first particle(K0Short)
            ej = getEnergy(l0Mass, tV0momj[0], tV0momj[1], tV0momj[2]); // Energy of first particle(K0Short)

            angle = getAngle(tV0momi[0], tV0momi[1], tV0momi[2], tV0momj[0], tV0momj[1], tV0momj[2]);
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
    }
    */
       for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
    {   // This is the begining of the V0 loop for first V0(K0Short)
        AliESDv0 *v0i = ((AliESDEvent*)fESD)->GetV0(iV0);
        if (!v0i) continue;
	//---> Fix On-the-Fly candidates, count how many swapped
        if( v0i->GetParamN()->Charge() > 0 && v0i->GetParamP()->Charge() < 0 ){
          fHistSwappedV0Counter -> Fill( 1 );
        }else{
          fHistSwappedV0Counter -> Fill( 0 ); 
        }
        if ( fkUseOnTheFly ) CheckChargeV0(v0i); 
	lOnFlyStatus = v0i->GetOnFlyStatus();

        lPt = v0i->Pt();
        if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;

        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0i->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0i->GetNindex());

        AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
        AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
        const AliExternalTrackParam * paramPosl = v0i->GetParamP();
	const AliExternalTrackParam * paramNegl = v0i->GetParamN();

	if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
	// Draw Armenteros-Podolanski Plot
	// from PWGGA/Hyperon/AliAnalysisTaskSigma0.cxx by Alexander Borissov.
	
	// Lambda -> P+ pi-  ---------------
            AliKFParticle negKFKpim(*paramNegl,211);
            AliKFParticle posKFKprot(*paramPosl,2212);
            AliKFParticle lamKF(negKFKpim,posKFKprot);    
	    
	    if (pTrack->GetMass() > 0.5){
	    //printf("this v0 is antilambda");
	    AliKFParticle negKFKpim(*paramNegl,2212);
            AliKFParticle posKFKprot(*paramPosl,211);
            AliKFParticle lamKF(negKFKpim,posKFKprot);
	    }

	    Double_t posp[3]= { pTrack->Px(),  pTrack->Py(),  pTrack->Pz() };
            Double_t negp[3]= { nTrack->Px(),  nTrack->Py(),  nTrack->Pz() };
            Double_t moth[3]= { lamKF.GetPx(), lamKF.GetPy(), lamKF.GetPz() };
            Double_t arpod[2]= {0,0};
            GetArPod( posp, negp, moth, arpod );
	    
            ((TH2F*)fOutputList->FindObject("fArmPod_lambda"))->Fill(arpod[1],arpod[0]);

        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }

        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if ( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if ( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

        if ( ( ( ( pTrack->GetTPCClusterInfo(2, 1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2, 1) ) < 70 ) )) continue;

        //GetKinkIndex condition
        if ( pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0 ) continue;

        //Findable clusters > 0 condition
        if ( pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0 ) continue;
        // find new v0 for Lambda
	v0i->ChangeMassHypothesis(3122);
        v0i->GetPxPyPz(tV0momi[0], tV0momi[1], tV0momi[2]);
        lInvMassLambda = v0i->GetEffMass();
        ((TH1F*)fOutputList->FindObject("fInvLambda_beforePID"))->Fill(lInvMassLambda); // Before PID
	Float_t nsigmaprP = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
        Float_t nsigmapiN = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
        if (pTrack->GetMass() > 0.5) if (nsigmaprP > 3.0 || nsigmapiN > 3.0) continue;
        if (pTrack->GetMass() < 0.5){
        Float_t nsigmaprN = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
        Float_t nsigmapiP = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
        if (nsigmaprN > 3.0 || nsigmapiP > 3.0) continue;
        }
	if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 0 && fkUseOnTheFly == kTRUE ) ){
	((TH1F*)fOutputList->FindObject("fInvLambda"))->Fill(lInvMassLambda);
	if (lInvMassLambda > l0Mass + 0.0008 || lInvMassLambda < l0Mass - 0.008) continue; // Mass window
	((TH1F*)fOutputList->FindObject("fInvLambdaCut"))->Fill(lInvMassLambda); // After Cut
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

double getEnergy(Double_t trueMass, Double_t Px, Double_t Py, Double_t Pz)
{
    return TMath::Sqrt(trueMass * trueMass + Px * Px + Py * Py + Pz * Pz);
}
double getAngle(Double_t Px1, Double_t Py1, Double_t Pz1, Double_t Px2, Double_t Py2, Double_t Pz2)
{
    return Px1 * Px2 + Py1 * Py2 + Pz1 * Pz2;
}
void GetArPod( Double_t pos[3], Double_t neg[3], Double_t moth[3],  Double_t arpod[2] ){
    
    //from PWGGA/Hyperon/AliAnalysisTaskSigma0.cxx by Alexander Borissov
    
    TVector3 momentumVectorPositiveKF(pos[0],pos[1],pos[2]);
    TVector3 momentumVectorNegativeKF(neg[0],neg[1],neg[2]);
    TVector3 vecV0(moth[0],moth[1],moth[2]);
    
    Float_t thetaV0pos=TMath::ACos(( momentumVectorPositiveKF* vecV0)/(momentumVectorPositiveKF.Mag() * vecV0.Mag()));
    Float_t thetaV0neg=TMath::ACos(( momentumVectorNegativeKF* vecV0)/(momentumVectorNegativeKF.Mag() * vecV0.Mag()));
    
    Float_t alfa =((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)-(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg))/
    ((momentumVectorPositiveKF.Mag())*TMath::Cos(thetaV0pos)+(momentumVectorNegativeKF.Mag())*TMath::Cos(thetaV0neg)) ;
    
    Float_t qt = momentumVectorPositiveKF.Mag()*TMath::Sin(thetaV0pos);
    
    arpod[0]=qt;
    arpod[1]=alfa;
    
}
void CheckChargeV0(AliESDv0 *v0)
{
   // This function is copied from PWGLF/STRANGENESS/LambdaK0/AliAnalysisTaskExtractV0.cxx
   // This function checks charge of negative and positive daughter tracks. 
   // If incorrectly defined (onfly vertexer), swaps out. 
   if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
      //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
      Long_t lCorrectNidx = v0->GetPindex();
      Long_t lCorrectPidx = v0->GetNindex();
      Double32_t	lCorrectNmom[3];
      Double32_t	lCorrectPmom[3];
      v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
      v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );
      AliExternalTrackParam	lCorrectParamN(
        v0->GetParamP()->GetX() , 
        v0->GetParamP()->GetAlpha() , 
        v0->GetParamP()->GetParameter() , 
        v0->GetParamP()->GetCovariance() 
      );
      AliExternalTrackParam	lCorrectParamP(
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
      AliESDv0 *v0correct = new AliESDv0(lCorrectParamN,lCorrectNidx,lCorrectParamP,lCorrectPidx);
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
      if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
        //AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
      }else{
        //AliWarning("Found Swapped Charges and fixed.");
      }
      //________________________________________________________________
   }else{
      //Don't touch it! ---
      //Printf("Ah, nice. Charges are already ordered...");
   }
   return;
} 
