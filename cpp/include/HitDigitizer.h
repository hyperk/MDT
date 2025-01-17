#include <iostream>
#include <vector>
#include <algorithm>

#include "HitTubeCollection.h"
#include "PMTResponse.h"
#include "MTRandom.h"

using std::cout;
using std::endl;
using std::vector;

class HitDigitizer
{
    public:
        HitDigitizer(int s=67823);
        virtual ~HitDigitizer();
        void Digitize(HitTubeCollection*, PMTResponse*);
        void DigitizeTube(HitTube*, PMTResponse*);

        void ApplyThreshold(double&, bool&);
        double DoTruncate(const double, const double);

    private:
        float fPrecisionCharge;
        float fPrecisionTiming;
        float fEfficiency;
        float fIntegWindow;
        int fApplyEff;
    
        MTRandom *fRand;
};
