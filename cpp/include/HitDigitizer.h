#include <iostream>
#include <vector>
#include <algorithm>

#include "HitTubeCollection.h"
#include "PMTResponse.h"
#include "MTRandom.h"

#include "TH1F.h"

using std::cout;
using std::endl;
using std::vector;

class HitDigitizer
{
    public:
        HitDigitizer(int s=67823);
        virtual ~HitDigitizer();
        virtual void Digitize(HitTubeCollection*, PMTResponse*);
        virtual void DigitizeTube(HitTube*, PMTResponse*);

        virtual void ApplyThreshold(double&, bool&);
        virtual double DoTruncate(const double, const double);

    protected:
        float fPrecisionCharge;
        float fPrecisionTiming;
        float fEfficiency;
        float fIntegWindow;
        int fApplyEff;
    
        MTRandom *fRand;
};

// mPMT specific digitizer
class HitDigitizer_mPMT : public HitDigitizer
{
    public:
        HitDigitizer_mPMT(int s=67823);
        virtual ~HitDigitizer_mPMT();
        void LoadWaveform(const string &filename);
        void DigitizeTube(HitTube*, PMTResponse*);
        TH1F BuildWavetrain(const vector<TrueHit*> PEs, double waveform_window);
        void FitWavetrain(TH1F hist, double& digiT, double& digiQ);

    private:
        TH1F* hWF;
        float fDt;
        float fDv;
        float fWaveformOffset;
        float fADCToPE;
        float fSigmaGuess;
};