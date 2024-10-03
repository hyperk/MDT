#pragma once

#include "PMTResponse.h"

class PMTResponse3inchR12199_02 : public GenericPMTResponse
{
    public: 
        PMTResponse3inchR12199_02(int, const string &s="");
        PMTResponse3inchR12199_02();
        virtual ~PMTResponse3inchR12199_02();
        float HitTimeSmearing(float);
        float HitTimeSmearing(float, int);
        void Initialize(int, const string &s="");
        bool ApplyDE(const TrueHit* th=NULL, const HitTube* ht=NULL);

    private:
        string fPMTDEFile;
        int fLoadDE;
        std::vector<double> fDE;
        void LoadPMTDE(const string &s);
        string fPMTTFile;
        int fLoadT;
        std::vector<double> fT;
        void LoadPMTTime(const string &s);
};
