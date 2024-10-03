#pragma once

#include <map>

#include "HitTube.h"
#include "PMTNoise.h"
#include "HitTubeCollection.h"
#include "HitDigitizer.h"
#include "TriggerInfo.h"
#include "TriggerAlgo.h"
#include "Configuration.h"
#include "PMTAfterpulse.h"
#include "MTRandom.h"

using std::map;

class MDTManager
{
    public:
        MDTManager(int seed=78923);
        virtual ~MDTManager();

        HitTubeCollection* GetHitTubeCollection(const string &s="Def") { return fPHC[s]; }
        TriggerInfo* GetTriggerInfo(const string &s="Def") { return fTrigInfo[s]; }
        
        void RegisterPMTType(const string &s="Def", PMTResponse *p=0);

        void DoInitialize();
        void DoAddDark(const string &s="Def");
        void DoDigitize(const string &s="Def");   
        void DoTrigger(const string &s="Def");
        void DoAddAfterpulse(const string &s="Def");
        
        void SetHitTubeCollection(HitTubeCollection*, const string &s="Def");
        bool HasThisPMTType(const string&);

        PMTResponse* GetPMTResponse(const string &s="Def") { return fPMTResp[s]; }

    private:
        map<string, HitDigitizer*> fDgtzr;
        MTRandom *fRndm;

        map<string, TriggerAlgo*> fTrigAlgo;
        map<string, TriggerInfo*> fTrigInfo;
        map<string, PMTResponse*> fPMTResp;
        map<string, PMTNoise*> fDark;
        map<string, HitTubeCollection*> fPHC; 
        string fDefName;
};
