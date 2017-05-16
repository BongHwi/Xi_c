/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliCentrality.h"

#include "AliAnalysisTaskXic.h"

class AliAnalysisTaskXic;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskXic) // classimp: necessary for root
double getEnergy(Double_t trueMass, Double_t Px, Double_t Py, Double_t Pz);
double getAngle(Double_t Px1, Double_t Py1, Double_t Pz1, Double_t Px2, Double_t Py2, Double_t Pz2);

AliAnalysisTaskXic::AliAnalysisTaskXic() : AliAnalysisTaskSE(),
    fESD(0x0),
    fOutputList(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskXic::AliAnalysisTaskXic(const char* name) : AliAnalysisTaskSE(name),
    fESD(0x0),
    fOutputList(0x0),
    fPIDResponse(0x0),
    fCentrality(0),
    fTrackCut(0x0),
    fHistPt(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
    // this chain is created by the analysis manager, so no need to worry about it,
    // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
    // you can add more output objects by calling DefineOutput(2, classname::Class())
    // if you add more output objects, make sure to call PostData for all of them, and to
    // make changes to your AddTask macro!
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
    // at the end of the analysis, the contents of this list are written
    // to the output file
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
    // if requested (dont worry about this now)

    // example of a histogram
    TH1F *fMultDist = new TH1F("fMultDist", "Multiplicity Distribution", 200, 0, 20000);
    fMultDist->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist);

    TH2F *fArmPod = new TH2F("fArmPod", "Armenteros-Podolski Plot", 100, 0, 0.25, 200, -1, 1);
    fOutputList->Add(fArmPod);


    TH1F *fInvLambdaCheck = new TH1F("fInvLambdaCheck", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
    fInvLambdaCheck->GetXaxis()->SetTitle("fInvLambdaCheck");
    fOutputList->Add(fInvLambdaCheck);

    TH1F *fInvLambda = new TH1F("fInvLambda", "Invariant mass distribution of Lambda", 400, 1.0, 1.2);
    fInvLambda->GetXaxis()->SetTitle("fInvLambda");
    fOutputList->Add(fInvLambda);
    TH1F *fInvK0Short = new TH1F("fInvK0Short", "Invariant mass distribution of K0s", 400, 0.3, 0.7);
    fInvK0Short->GetXaxis()->SetTitle("fInvK0Short");
    fOutputList->Add(fInvK0Short);


    TH1F *fInvLambdaCut = new TH1F("fInvLambdaCut", "Invariant mass distribution of Lambda after mass window cut", 400, 1.0, 1.2);
    fInvLambdaCut->GetXaxis()->SetTitle("fInvLambdaCut");
    fOutputList->Add(fInvLambdaCut);
    TH1F *fInvK0ShortCut = new TH1F("fInvK0ShortCut", "Invariant mass distribution of K0s after mass window cut", 400, 0.3, 0.7);
    fInvK0Short->GetXaxis()->SetTitle("fInvK0ShortCut");
    fOutputList->Add(fInvK0ShortCut);


    TH1F *hEventSelecInfo = new TH1F("hEventSelecInfo", "hEventSelecInfo", 10, 0, 10);
    fOutputList->Add(hEventSelecInfo);
    hEventSelecInfo->GetXaxis()->SetBinLabel(1, "NONE");
    hEventSelecInfo->GetXaxis()->SetBinLabel(2, "kMB");
    hEventSelecInfo->GetXaxis()->SetBinLabel(3, "kCentral");
    hEventSelecInfo->GetXaxis()->SetBinLabel(4, "kSemiCentral");
    hEventSelecInfo->GetXaxis()->SetBinLabel(5, "kINT7");
    hEventSelecInfo->GetXaxis()->SetBinLabel(6, "kAny");


    TH1F *hCentrality = new TH1F("hCentrality", "Centrality", 100, 0, 100);
    hCentrality->GetXaxis()->SetTitle("Centrality");
    fOutputList->Add(hCentrality);


    TH3F *fVertexDistXYZ = new TH3F("fVertexDistXYZ", "Vertex Distribution", 20, -1, 1, 20, -1, 1, 60, -30, 30);
    fVertexDistXYZ->GetXaxis()->SetTitle("X Vertex (cm)");
    fVertexDistXYZ->GetYaxis()->SetTitle("Y Vertex (cm)");
    fVertexDistXYZ->GetZaxis()->SetTitle("Z Vertex (cm)");
    fOutputList->Add(fVertexDistXYZ);

    TH2F *hInvMassWithPt = new TH2F("hInvMassWithPt", "Invariant mass distribution vs Pt", 1000, 2.0, 3.0, 100, 0, 10);
    fOutputList->Add(hInvMassWithPt);
    TH1F *hInvMass = new TH1F("hInvMass", "Invariant mass distribution", 1000, 2.0, 3.0);
    fOutputList->Add(hInvMass);


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

    /*
        double fx = -100;
        double fy = -100;
        double fz = -100;


        int fTPCNcls = -100;
        double fDCA = -100; // sqrt (xx+yy+zz)
        double fDCAr = -100; // sqrt(xx+yy)
        double fDCAZ = -100; // z infor

        double  fTrueMassPr = .93827, fTrueMassPi = .13957;
        double  fMass = 0.;

        //------------------------------------------------
        //Step 5: Loops on tracks
        //------------------------------------------------
            Int_t iTracks(fESD->GetNumberOfTracks());           // see how many tracks there are in the event
            for(Int_t i(0); i < iTracks; i++) {                 // loop ove rall these tracks
                AliESDtrack* track = fESD->GetTrack(i);         // get a track (type AliESDTrack) from the event
                if(!track) continue;                            // if we failed, skip this track
                fHistPt->Fill(track->Pt());                     // plot the pt value of the track in a histogram

                if(!fTrackCut->AcceptTrack(track)) continue;


                ((TH1F*)fOutputList->FindObject("fPtDist"))->Fill(track->Pt());
                ((TH1F*)fOutputList->FindObject("fPhiDist"))->Fill(track->Phi());
                ((TH1F*)fOutputList->FindObject("fEtaDist"))->Fill(track->Eta());



                fTPCNcls = track->GetTPCNcls();
                if(fDCAZ > 2.)continue;
                // fDCAr --> Homework
                if(fTPCNcls < 70.) continue;

                //Home work
                fx = track->GetX();
                fy = track->GetY();
                fz = track->GetZ();

                fDCAr = sqrt((fx-primaryVtx[0])*(fx-primaryVtx[0])+(fy-primaryVtx[1])*(fy-primaryVtx[1]));
                fDCAZ = fabs(fz-primaryVtx[2]);

                ((TH1F*)fOutputList->FindObject("fDCArDist"))->Fill(fDCAr);




                ((TH1F*)fOutputList->FindObject("fNTPCcls"))->Fill(fTPCNcls);
                ((TH1F*)fOutputList->FindObject("fPtDistTight"))->Fill(track->Pt());


                //------------------------------------------------
                //Step 6: Particle IDentification
                //------------------------------------------------

                Float_t fTPCPIDmom = track->GetTPCmomentum();
                Float_t sigTPC = track->GetTPCsignal();
                Float_t nsigpi= fabs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion));
                Float_t nsigk= fabs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon));
                Float_t nsigpr= fabs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));

                ((TH2F*)fOutputList->FindObject("hTPCPID"))->Fill(fTPCPIDmom,sigTPC);


                Double_t tDecayVertexV0[3];
                track->GetXYZ(tDecayVertexV0);


                Double_t tV0mom[3];
                track->GetPxPyPz( tV0mom);

                Int_t fCharge = -100.;
                fCharge = track->Charge();

                Double_t pV0mom[3];
                Double_t nV0mom[3];

                    if(nsigpi<3.)track->GetPxPyPz(nV0mom);
                    if(nsigpr<3.)track->GetPxPyPz(pV0mom);
                    if((nV0mom[0]*nV0mom[0] + nV0mom[1]*nV0mom[1] + nV0mom[2]*nV0mom[2])<=0) continue;
                    if((pV0mom[0]*pV0mom[0] + pV0mom[1]*pV0mom[1] + pV0mom[2]*pV0mom[2])<=0) continue;

                    Double_t ep=TMath::Sqrt(fTrueMassPr*fTrueMassPr + pV0mom[0]*pV0mom[0] + pV0mom[1]*pV0mom[1] + pV0mom[2]*pV0mom[2]);

                    Double_t en=TMath::Sqrt(fTrueMassPi*fTrueMassPi + nV0mom[0]*nV0mom[0] + nV0mom[1]*nV0mom[1] + nV0mom[2]*nV0mom[2]);

                    //cout<< "en === " << en <<" ep == " << ep <<endl;

                    // Double_t pl=TMath::Sqrt(tV0mom[0]*tV0mom[0] + tV0mom[1]*tV0mom[1] + tV0mom[2]*tV0mom[2]);
               // Double_t pl=TMath::Sqrt((nV0mom[0]+pV0mom[0])*(nV0mom[0]+pV0mom[0]) + (nV0mom[1]+pV0mom[1])*(nV0mom[1]+pV0mom[1]) + (nV0mom[2]+pV0mom[2])*(nV0mom[2]+pV0mom[2]));


                Double_t angle=nV0mom[0]*pV0mom[0]+nV0mom[1]*pV0mom[1]+nV0mom[2]*pV0mom[2];
                fMass = fTrueMassPr*fTrueMassPr+fTrueMassPi*fTrueMassPi+2.*en*ep-2.*angle;
                if(fMass<=0) continue;
                fMass=sqrt(fMass);
                //cout<< " Lambda Invariant mass " << fMass <<endl;

                ((TH2F*)fOutputList->FindObject("fInvLambdaTrack"))->Fill(fMass);
                }
        */ //No Track use for this analysis
    double lInvMassLambda = 0.;
    double lInvMassK0Short = 0.0;

    static Double_t k0Mass = 0.497611;      static Double_t l0Mass = 1.115683;
    static Double_t piMass = 0.13957;       static Double_t protonMass = 0.93827;

    Double_t fMinV0Pt = 0.15;  Double_t fMaxV0Pt = 1.E10;    Double_t lV0Radius = 0, lPt = 0;
    Double_t fQt = 0;       Double_t falpha = 0;

    Double_t V0momK0ShortDaughterN[3] = {.0, .0, .0};
    Double_t V0momK0ShortDaughterP[3] = {.0, .0, .0};
    Double_t V0momLambda0DaughterN[3] = {.0, .0, .0};
    Double_t V0momLambda0DaughterP[3] = {.0, .0, .0};

    Double_t tV0momi[3];    Double_t tV0momj[3];    Double_t tV0mom_result[3];

    Double_t ei = 0.;       Double_t ej = 0.;
    Double_t angle = 0.;    Double_t fPt_result = 0.;
    Double_t  fMass = 0.;

    // loop for Lambda
    Int_t nv0s = 0;
    nv0s = fESD->GetNumberOfV0s();

    //for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
    for (Int_t iV0 = 0; iV0 < nv0s; iV0++)
    {   // This is the begining of the V0 loop for first V0(K0Short)
        AliESDv0 *v0i = ((AliESDEvent*)fESD)->GetV0(iV0);
        if (!v0i) continue;

        //// ---- Default Cut ---- ////
        lPt = v0i->Pt();
        if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;
        if (v0i->GetDcaV0Daughters() > 1.5) continue;
        if (v0i->GetV0CosineOfPointingAngle() > 0.97) continue;
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0i->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0i->GetNindex());
        AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
        AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        // Armenteros-Podolski Plot
	//fQt = pTrack->Get

        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }
        // Pt range for tracks
	// 
	// TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        //if ( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        //if ( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        //if ( ( ( ( pTrack->GetTPCClusterInfo(2, 1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2, 1) ) < 70 ) )) continue;
        //GetKinkIndex condition
        //if ( pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0 ) continue;
        //Findable clusters > 0 condition
        //if ( pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0 ) continue;
        if(!fTrackCut->AcceptTrack(pTrack)) continue;
	if(!fTrackCut->AcceptTrack(nTrack)) continue;
	//// -------------------- ////


        // find new v0 for K0Short
        v0i->ChangeMassHypothesis(310); //kK0Short
        v0i->GetPxPyPz(tV0momi[0], tV0momi[1], tV0momi[2]);
        lInvMassK0Short = v0i->GetEffMass();
        ((TH1F*)fOutputList->FindObject("fInvK0Short"))->Fill(lInvMassK0Short); // Before Cut
        if (lInvMassK0Short > k0Mass + 0.008 || lInvMassK0Short < k0Mass - 0.008) continue; // Mass window
        ((TH1F*)fOutputList->FindObject("fInvK0ShortCut"))->Fill(lInvMassK0Short); // After Cut

        for (Int_t jV0 = iV0; jV0 < nv0s; jV0++)
        {   // This is the begining of the V0 loop for second V0(Lambda)
            if (iV0 == jV0) continue;
            AliESDv0 *v0j = ((AliESDEvent*)fESD)->GetV0(jV0);

            //// ---- Default Cut ---- ////
            if (!v0j) continue;
            lPt = v0j->Pt();
            if ((lPt < fMinV0Pt) || (fMaxV0Pt < lPt)) continue;
            if (v0j->GetDcaV0Daughters() > 1.5) continue;
            if (v0j->GetV0CosineOfPointingAngle() > 0.97) continue;
            UInt_t lKeyPos = (UInt_t)TMath::Abs(v0j->GetPindex());
            UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0j->GetNindex());
            AliESDtrack *pTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyPos);
            AliESDtrack *nTrack = ((AliESDEvent*)fESD)->GetTrack(lKeyNeg);
            if (!pTrack || !nTrack) {
                Printf("ERROR: Could not retreive one of the daughter track");
                continue;
            }
            if ( pTrack->GetSign() == nTrack->GetSign()) {
                continue;
            }
            // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
            //if ( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            //if ( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            //if ( ( ( ( pTrack->GetTPCClusterInfo(2, 1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2, 1) ) < 70 ) )) continue;
            //GetKinkIndex condition
            //if ( pTrack->GetKinkIndex(0) > 0 || nTrack->GetKinkIndex(0) > 0 ) continue;
            //Findable clusters > 0 condition
            //if ( pTrack->GetTPCNclsF() <= 0 || nTrack->GetTPCNclsF() <= 0 ) continue;
            if(!fTrackCut->AcceptTrack(pTrack)) continue;
	    if(!fTrackCut->AcceptTrack(nTrack)) continue;
	    //// -------------------- ////

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

            ((TH2F*)fOutputList->FindObject("hInvMassWithPt"))->Fill(fMass, fPt_result); // with Pt
            ((TH1F*)fOutputList->FindObject("hInvMass"))->Fill(fMass); // Cumulated
        }
    }

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
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }

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

        // find new v0 for K0Short
        v0i->ChangeMassHypothesis(kLambda0);
        lInvMassLambda = v0i->GetEffMass();
        ((TH1F*)fOutputList->FindObject("fInvLambdaCheck"))->Fill(lInvMassLambda);

    }

    PostData(1, fOutputList);                           // stream the results the analysis of this event to
    // the output manager which will take care of writing

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
