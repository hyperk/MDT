#include "PMTResponse3inchR12199_02.h"

#include "Configuration.h"

#include <map>
#include <math.h> 

using std::map;

PMTResponse3inchR12199_02::PMTResponse3inchR12199_02(int seed, const string &pmtname)
{
    this->Initialize(seed, pmtname);
}

PMTResponse3inchR12199_02::PMTResponse3inchR12199_02()
{
}

PMTResponse3inchR12199_02::~PMTResponse3inchR12199_02()
{
}

void PMTResponse3inchR12199_02::Initialize(int seed, const string &pmtname)
{
    fPMTType = pmtname;
    fRand = new MTRandom(seed);

    map<string, string> s;
    s["TimingResConstant"] = "TimingResConstant";
    s["TimingResMinimum"] = "TimingResMinimum"; 
    s["ScalFactorTTS"] = "ScalFactorTTS";
    s["SPECDFFile"] = "SPECDFFile";
    s["PMTDE"] = "PMTDE";
    s["PMTTime"] = "PMTTime";
    if( fPMTType!="" )
    {
        map<string, string>::iterator i;
        for(i=s.begin(); i!=s.end(); i++)
        {
            i->second += "_" + fPMTType;
        }
    }
    Configuration *Conf = Configuration::GetInstance();
    Conf->GetValue<float>(s["TimingResConstant"], fTResConstant);
    Conf->GetValue<float>(s["TimingResMinimum"], fTResMinimum);
    Conf->GetValue<float>(s["ScalFactorTTS"], fSclFacTTS);
    Conf->GetValue<string>(s["SPECDFFile"], fTxtFileSPECDF);
    Conf->GetValue<string>(s["PMTDE"], fPMTDEFile);
    Conf->GetValue<string>(s["PMTTime"], fPMTTFile);
    this->LoadCDFOfSPE(fTxtFileSPECDF);
    this->LoadPMTDE(fPMTDEFile);
    this->LoadPMTTime(fPMTTFile);
}

void PMTResponse3inchR12199_02::LoadPMTTime(const string &filename)
{
    fLoadT = 0;
    fT.clear();
    ifstream ifs(filename.c_str());
    if (!ifs)
    {
        cout<<" PMTResponse3inchR12199_02::LoadPMTTime" <<endl;
        cout<<"  - No PMT Time file: " << filename <<endl;
        cout<<"  - Do not apply individual PMT timing " << endl;
    }
    string aLine;
    while( std::getline(ifs, aLine) )
    {
        if( aLine[0] == '#' ){ continue; }
        stringstream ssline(aLine);
        string item;
        while (getline(ssline, item, ssline.widen(' ')))
        {
            fT.push_back( atof(item.c_str()) );
        }
    }
    ifs.close();

    fLoadT = fT.size();
    if (fLoadT>0)
    {
        cout<<" PMTResponse3inchR12199_02::LoadPMTTime" <<endl;
        cout<<"  - Load PMT Time file: " << filename <<endl;
        cout<<"  - # Entries = " << fLoadT << endl;
    }
}

float PMTResponse3inchR12199_02::HitTimeSmearing(float Q)
{
    Q = (Q > 0.5) ? Q : 0.5;
    float timingResolution = 0.5*fSclFacTTS*(0.33 + sqrt(fTResConstant/Q));
    if( timingResolution<fTResMinimum ){ timingResolution = fTResMinimum; }
    return fRand->Gaus(0.0,timingResolution);
}

float PMTResponse3inchR12199_02::HitTimeSmearing(float Q, int tubeID)
{
    if (fLoadT>0)
    {
        if (tubeID>=fLoadT)
        {
            cout<<" PMTResponse3inchR12199_02::HitTimeSmearing" <<endl;
            cout<<"  - tubeID = " << tubeID << " >= fLoadT = " << fLoadT << endl;
            cout<<"  -> EXIT" <<endl;
            exit(-1);
        }
        return fT[tubeID];
    }

    return 0;
}

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////