/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskXic_H
#define AliAnalysisTaskXic_H

#include "AliAnalysisTaskSE.h"
#include "AliESDpid.h"

class AliAnalysisTaskXic : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskXic();
                                AliAnalysisTaskXic(const char *name);
        virtual                 ~AliAnalysisTaskXic();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

    private:
        AliESDEvent*            fESD;         //! input event
        TList*                  fOutputList;  //! output list
        AliESDtrackCuts*        fTrackCuts;			// Track cuts
	      TList*                  fOutputList2; //! output list
        AliPIDResponse*         fPIDResponse; //! PID object
        AliCentrality*          fCentrality;  //! Centrality object
        AliESDtrackCuts        *fTrackCut;   //! ESD track cuts


        TH1F*                   fHistPt;        //! dummy histogram
	TH1F    *fHistSwappedV0Counter;     					        //! Swapped V0 Counter
        AliAnalysisTaskXic(const AliAnalysisTaskXic&); // not implemented
        AliAnalysisTaskXic& operator=(const AliAnalysisTaskXic&); // not implemented

        ClassDef(AliAnalysisTaskXic, 1);
};

#endif
