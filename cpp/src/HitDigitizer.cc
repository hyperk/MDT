#include "HitDigitizer.h"
#include "Configuration.h"
#include "TF1.h"
#include "TFitResult.h"
        
HitDigitizer::HitDigitizer(int seed) :
fPrecisionCharge(0.), 
fPrecisionTiming(0.1),
fEfficiency( 0.985 ),
fIntegWindow( 200. ),
fApplyEff( 1 )
{
    fRand = new MTRandom( seed );
    Configuration *Conf = Configuration::GetInstance();
    Conf->GetValue<float>("DigiHitIntegrationWindow", fIntegWindow);
    Conf->GetValue<float>("PrecisionTiming", fPrecisionTiming);
    Conf->GetValue<float>("PrecisionCharge", fPrecisionCharge);
    Conf->GetValue<int>("ApplyDAQEfficiency", fApplyEff);
}

HitDigitizer::~HitDigitizer()
{
    if( fRand ){ delete fRand; fRand = NULL; }
}

void HitDigitizer::Digitize(HitTubeCollection *hc, PMTResponse *pr)
{
    for(hc->Begin(); !hc->IsEnd(); hc->Next())
    {
        if( (&(*hc)())->GetNRawPE()>0 )
        {
            this->DigitizeTube(&(*hc)(), pr);
        }
    }
}

// Based on WCSimWCDigitizerSKI::DigitizeHits in WCSimWCDigitizer.cc
void HitDigitizer::DigitizeTube(HitTube *aHT, PMTResponse *pr)
{
    aHT->SortTrueHits();
    const int NPE = aHT->GetNRawPE();
    const vector<TrueHit*> PEs = aHT->GetPhotoElectrons();

    // Taken from WCSimWCDigitizerSKI::DigitizeHits
    double sumSPE = 0.;
    double tSmeared = -9999.;
    bool isAlive = false;
    double digiT = 0.;
    double digiQ = 0.;

    double intgr_srt = (double)PEs[0]->GetTime();
    double intgr_end = intgr_srt+fIntegWindow;
    vector<int> parent_composition;  
    parent_composition.clear();

    for(int iPE=0; iPE<NPE; iPE++)
    {
        if( PEs[iPE]->GetTime()>=intgr_srt && PEs[iPE]->GetTime()<intgr_end )
        {
            sumSPE += pr->GetRawSPE(PEs[iPE], aHT);
            parent_composition.push_back( PEs[iPE]->GetParentId() );
        }
        else
        {
            digiT = intgr_srt;
            digiQ = sumSPE;

            isAlive = true;
            if (fApplyEff ){ this->ApplyThreshold(digiQ, isAlive); }
            if( isAlive ) 
            {
                tSmeared = pr->HitTimeSmearing(digiQ);
                digiQ *= fEfficiency;
                digiT += tSmeared;
                digiQ = this->DoTruncate(digiQ, fPrecisionCharge);
                digiT = this->DoTruncate(digiT, fPrecisionTiming);
                aHT->AddDigiHit(digiT, digiQ, parent_composition);
            }
            sumSPE = 0.;
            parent_composition.clear(); 

            intgr_srt = PEs[iPE]->GetTime();
            intgr_end = intgr_srt+fIntegWindow;
            sumSPE += pr->GetRawSPE(PEs[iPE], aHT);
            parent_composition.push_back( PEs[iPE]->GetParentId() );
        }
    }

    digiT = intgr_srt;
    digiQ = sumSPE;
    tSmeared = pr->HitTimeSmearing(digiQ);
    isAlive = true;
    if (fApplyEff ){ this->ApplyThreshold(digiQ, isAlive); }
    if( isAlive )
    {
        digiQ *= fEfficiency;
        digiT += tSmeared ;
        digiQ = this->DoTruncate(digiQ, fPrecisionCharge);
        digiT = this->DoTruncate(digiT, fPrecisionTiming);
        aHT->AddDigiHit(digiT, digiQ, parent_composition);
    }
}


void HitDigitizer::ApplyThreshold(double& pe, bool& pass)
{
// Apply DAQ efficiency
// Taken from WCSimWCDigitizerSKI::Threshold 
// in WCSimWCDigitizer.hh
//
// NOTE: the input charge, pe, will be altered in this function
    pass=true;
    double x=pe+0.1;
    double thr=0.;
    if( x<1.1 )
    {
      thr = std::min(1.0,
		     -0.06374+x*(3.748+x*(-63.23+x*(452.0+x*(-1449.0+x*(2513.0
									+x*(-2529.+x*(1472.0+x*(-452.2+x*(51.34+x*2.370))))))))));
    } 
    else 
    {
      thr = 1.0;
    }

    if( thr<fRand->Rndm() )
    {
        pe=0.;
        pass=false;
    }
    else
    {
        pe += fRand->Gaus(0., 0.03);
    }
}

double HitDigitizer::DoTruncate(const double x, const double precision)
{
// Based on WCSimWCDigitizerBase::AddNewDigit in WCSimWCDigitizer.hh
// digitised hit information does not have infinite precision
// so need to round the charge and time information
//
// The following is based on WCSimWCDigitizerBas::Truncate
// in WCSimWCDigitizer.hh
    if(precision < 1E-10){ return x; }
    return precision * (int)(x/precision); 
}

// mPMT specific digitizer
HitDigitizer_mPMT::HitDigitizer_mPMT(int seed) : HitDigitizer(seed)
{
    hWF = nullptr;
    fDt = 8;
    fDv = 1;
    fWaveformOffset = 95;
    fADCToPE = 108;
    fSigmaGuess = 6;

    Configuration *Conf = Configuration::GetInstance();
    Conf->GetValue<float>("SamplingInterval", fDt);
    Conf->GetValue<float>("SamplingResolution", fDv);
    Conf->GetValue<float>("WaveformOffset", fWaveformOffset);
    Conf->GetValue<float>("ADCToPE", fADCToPE);
    Conf->GetValue<float>("SigmaGuess", fSigmaGuess);

    this->LoadWaveform(Conf->GetValue<string>("WaveformFile"));
}

HitDigitizer_mPMT::~HitDigitizer_mPMT()
{
    if( fRand ){ delete fRand; fRand = NULL; }
    if(hWF != nullptr) delete hWF;
}

void HitDigitizer_mPMT::LoadWaveform(const string &filename)
{
    if(hWF != nullptr) delete hWF;

    ifstream ifs(filename.c_str());
    if (!ifs.is_open())
    {
        cout<<" [ERROR] HitDigitizer_mPMT::LoadWaveform" <<endl;
        cout<<"  - Cannot open: " << filename <<endl;
        cout<<"  -> EXIT" <<endl;
        exit(-1);
    }
    string aLine;
    //vector<float> t_val;
    float time_val[20000]; // just a large number to avoid overflow
    vector<float> volt_val;
    int nEntries = 0;
    while( std::getline(ifs, aLine) )
    {
        if( aLine[0] == '#' ){ continue; }
        stringstream ssline(aLine);
        float t, v;
        ssline >> t >> v;
        //t_val.push_back(t);
        time_val[nEntries++] = t*1.e9; // ns
        volt_val.push_back(v*1000);    // mV
    }
    ifs.close();

    hWF = new TH1F("hWF","hWF",nEntries-1,time_val);
    hWF->SetDirectory(0);
    for (int i=0;i<nEntries-1;i++)
    {
        hWF->SetBinContent(i+1,volt_val.at(i));
    }
}

void HitDigitizer_mPMT::DigitizeTube(HitTube *aHT, PMTResponse *pr)
{
    aHT->SortTrueHits();
    const int NPE = aHT->GetNRawPE();
    const vector<TrueHit*> PEs = aHT->GetPhotoElectrons();

    // Taken from WCSimWCDigitizerSKI::DigitizeHits
    double sumSPE = 0.;
    bool isAlive = false;
    double digiT = 0.;
    double digiQ = 0.;

    double intgr_srt = (double)PEs[0]->GetTime();
    double intgr_end = intgr_srt+fIntegWindow;
    vector<int> parent_composition;  
    parent_composition.clear();
    vector<TrueHit*> digiPEs;
    digiPEs.clear();

    for(int iPE=0; iPE<NPE; iPE++)
    {
        if( PEs[iPE]->GetTime()>=intgr_srt && PEs[iPE]->GetTime()<intgr_end )
        {
            sumSPE += pr->GetRawSPE(PEs[iPE], aHT);
            parent_composition.push_back( PEs[iPE]->GetParentId() );
            digiPEs.push_back(PEs[iPE]);
            //intgr_end = PEs[iPE]->GetTime()+fIntegWindow;
        }
        else
        {

            TH1F hWT = BuildWavetrain(digiPEs, fIntegWindow);
            this->FitWavetrain(hWT,digiT,digiQ);
            digiT += pr->HitTimeSmearing(digiQ,aHT->GetTubeID());

            isAlive = true;
            if (fApplyEff ){ this->ApplyThreshold(digiQ, isAlive); }
            if( isAlive ) 
            {
                digiQ *= fEfficiency;
                digiQ = this->DoTruncate(digiQ, fPrecisionCharge);
                digiT = this->DoTruncate(digiT, fPrecisionTiming);
                aHT->AddDigiHit(digiT, digiQ, parent_composition);
                if (!aHT->GetDigiWF()) aHT->SetDigiWF(hWT);
                aHT->SetDigiPulls(digiT-digiPEs.front()->GetTime(),digiQ-digiPEs.size());
                aHT->SetTrueTQ(digiPEs.front()->GetTime(),digiPEs.size());
            }
            digiPEs.clear();
            sumSPE = 0.;
            parent_composition.clear(); 

            intgr_srt = PEs[iPE]->GetTime();
            intgr_end = intgr_srt+fIntegWindow;
            sumSPE += pr->GetRawSPE(PEs[iPE], aHT);
            parent_composition.push_back( PEs[iPE]->GetParentId() );
            digiPEs.push_back(PEs[iPE]);
        }
    }

    TH1F hWT = BuildWavetrain(digiPEs, fIntegWindow);
    this->FitWavetrain(hWT,digiT,digiQ);
    digiT += pr->HitTimeSmearing(digiQ,aHT->GetTubeID());
    isAlive = true;
    if (fApplyEff ){ this->ApplyThreshold(digiQ, isAlive); }
    if( isAlive )
    {
        digiQ *= fEfficiency;
        digiQ = this->DoTruncate(digiQ, fPrecisionCharge);
        digiT = this->DoTruncate(digiT, fPrecisionTiming);
        aHT->AddDigiHit(digiT, digiQ, parent_composition);
        if (!aHT->GetDigiWF()) aHT->SetDigiWF(hWT);
        aHT->SetDigiPulls(digiT-digiPEs.front()->GetTime(),digiQ-digiPEs.size());
        aHT->SetTrueTQ(digiPEs.front()->GetTime(),digiPEs.size());
    }
}

TH1F HitDigitizer_mPMT::BuildWavetrain(const vector<TrueHit*> PEs, double waveform_window)
{
    double dt = fDt; // interval, start, end of sampling
    double tmin = floor((PEs.front()->GetTime())/dt)*dt;
    double tmax = ceil((PEs.back()->GetTime()+waveform_window)/dt)*dt; 
    TH1F hWT("","",(int)(tmax-tmin)/dt,tmin,tmax);

    double waveform_offset = fWaveformOffset;

    for(long unsigned int iPE=0; iPE<PEs.size(); iPE++)
    {
        double trueT = (double)PEs[iPE]->GetTime();

        int bin_start = hWT.FindBin(trueT);
        int bin_end = hWT.FindBin(trueT+waveform_window);

        for (int iBin = bin_start; iBin<=bin_end; iBin++)
        {
            double adcT = hWT.GetBinLowEdge(iBin);
            double adcV = hWF->GetBinContent(hWF->FindBin(adcT-trueT+waveform_offset));
            hWT.Fill(adcT+dt/2,adcV);
        }
    }

    for (int i=1;i<=hWT.GetNbinsX();i++)
    {
        double val = (int)(hWT.GetBinContent(i)/fDv);
        val *= fDv;
        hWT.SetBinContent(i,val);
        hWT.SetBinError(i,fDv);
    }

    return hWT;
}

void HitDigitizer_mPMT::FitWavetrain(TH1F hist, double& digiT, double& digiQ)
{
    double dt = fDt; 
    int maxBin = hist.GetMaximumBin();
    double digiT_guess = hist.GetBinCenter(maxBin);
    double digiQ_guess = hist.GetBinContent(maxBin);
    double sigma_guess = fSigmaGuess;
    double range_min = hist.GetBinLowEdge(maxBin);
    double range_max = hist.GetBinLowEdge(maxBin+1);
    double adc_to_pe = fADCToPE;
    for (int i=maxBin-1;i>=1;i--)
    {
        if (hist.GetBinContent(i)>0) range_min = hist.GetBinLowEdge(i);
        else break;
    }
    for (int i=maxBin+1;i<=hist.GetNbinsX();i++)
    {
        if (hist.GetBinContent(i)>0) range_max = hist.GetBinLowEdge(i+1);
        else break;
    }
    TF1 f1("f1", "gaus", range_min,range_max);
    f1.SetParameters(digiQ_guess,digiT_guess,sigma_guess);
    TFitResultPtr r = hist.Fit("f1","S Q N","",range_min,range_max);
    digiQ = r->Parameter(0)*r->Parameter(2)/adc_to_pe;
    digiT = r->Parameter(1) - 3*r->Parameter(2) - dt;
}