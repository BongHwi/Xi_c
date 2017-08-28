/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskXic_H
#define AliAnalysisTaskXic_H

class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TProfile;

class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODPid.h"
#include "AliESDpid.h"

class AliAnalysisTaskXic : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskXicMC();
                                AliAnalysisTaskXicMC(const char *name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption=0);
        virtual                 ~AliAnalysisTaskXicMC();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        void                    XicInit();      //! initialization of fixed values
        AliESDEvent*            fESD;         //! input event
        TList*                  fOutputList;  //! output list
        AliESDtrackCuts*        fTrackCuts;			// Track cuts
	      TList*                  fOutputList2; //! output list
        AliPIDResponse*         fPIDResponse; //! PID object
        AliCentrality*          fCentrality;  //! Centrality object
        AliESDtrackCuts        *fTrackCut;   //! ESD track cuts


        TH1F*                   fHistPt;        //! dummy histogram
        Bool_t                  fMCcase;        //! switch for MC data or real data
        Bool_t                  fAODcase;       //! switch for AODs or ESDs
        Int_t                   fEventCounter;  //! The event counter
        Int_t                   fCutList;       //! Cut List option (mean values or systematic variations)

	TH1F    *fHistSwappedV0Counter;     					        //! Swapped V0 Counter
        AliAnalysisTaskXic(const AliAnalysisTaskXic&); // not implemented
        AliAnalysisTaskXic& operator=(const AliAnalysisTaskXic&); // not implemented

        ClassDef(AliAnalysisTaskXic, 1);
};

#endif
