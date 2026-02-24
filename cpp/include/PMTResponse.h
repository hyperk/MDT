#pragma once

#include <string>
#include <TGraph.h>
#include "MTRandom.h"
#include "HitTube.h"

using std::string;

// Manage:
//
//  - generation of SPE
//  - timing resolution
//
class PMTResponse
{
    public:
        virtual ~PMTResponse() {};
        virtual double GetRawSPE(const TrueHit* th=NULL, const HitTube* ht=NULL) = 0; 
        virtual float HitTimeSmearing(float) = 0;
        virtual void Initialize(int, const string &s="") = 0;
        virtual bool ApplyDE(const TrueHit* th=NULL, const HitTube* ht=NULL) = 0;

    protected:
        MTRandom *fRand;
        string fPMTType;
        float fSclFacTTS;
};

//  Based on the implementation of R12199_02 3'' PMT for mPMT in WCSim
//   - WCSimPMTObject.cc
class GenericPMTResponse : public PMTResponse
{
    public:
        GenericPMTResponse(const int, const string &s="");
        GenericPMTResponse();
        virtual ~GenericPMTResponse();
        void Initialize(int, const string &s="");
        virtual double GetRawSPE(const TrueHit* th=NULL, const HitTube* ht=NULL);
        virtual float HitTimeSmearing(float);
        virtual bool ApplyDE(const TrueHit* th=NULL, const HitTube* ht=NULL);

    protected:
        void LoadCDFOfSPE(const string &s);
        float fqpe0[501];
        string fTxtFileSPECDF;

    //private:
    protected:
        float fTResConstant;
        float fTResMinimum;
};


class ResponseBoxandLine20inchHQE : public GenericPMTResponse
{
    public: 
        ResponseBoxandLine20inchHQE(int, const string &s="");
        ResponseBoxandLine20inchHQE();
        virtual ~ResponseBoxandLine20inchHQE();
        float HitTimeSmearing(float);
        void Initialize(int, const string &s="");

    private:
        static const double ksig_param[4];
        static const double klambda_param[2];
        double fhighcharge_param[2];
};

class Response3inchR14374 : public GenericPMTResponse
{
    public: 
        Response3inchR14374(int, const string &s="");
        Response3inchR14374();
        virtual ~Response3inchR14374();
        float HitTimeSmearing(float);
        void Initialize(int, const string &s="");

    private:
        double fTimeResAt1PE;
};

class Response3inchR14374_WCTE : public GenericPMTResponse
{
    public: 
        Response3inchR14374_WCTE(int, const string &s="");
        Response3inchR14374_WCTE();
        virtual ~Response3inchR14374_WCTE();
        float HitTimeSmearing(float);
        void Initialize(int, const string &s="");
        
    private:
        TGraph  *gTResol = nullptr; // timing resolution
};