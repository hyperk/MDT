#include <iostream>
#include <string>
#include <vector>

#include "MDTManager.h"
#include "WCRootData.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int fNEvtToProc = -1;
int fSeed = 67592;
string fParFileName = "";
string fInFileName = "";
string fOutFileName = "";

bool ParseCmdArguments(int, char**);

int main(int argc, char **argv)
{
    ParseCmdArguments(argc, argv);

    Configuration *Conf = Configuration::GetInstance();
    Conf->ReadParameter(fParFileName);
    Conf->PrintParameters();
    Conf = 0;


    const int NPMTType = 2;
    string fPMTType[NPMTType];
    fPMTType[0] = "BoxandLine20inchHQE";
    fPMTType[1] = "PMT3inchR14374";

    MDTManager *MDT = new MDTManager(fSeed);
    MDT->RegisterPMTType(fPMTType[0], new ResponseBoxandLine20inchHQE());
    MDT->RegisterPMTType(fPMTType[1], new Response3inchR14374());

    const vector<string> listWCRootEvt{"wcsimrootevent", "wcsimrootevent2"};
    const vector<string> listWCRootCopyTree{"wcsimGeoT","Settings","wcsimRootOptionsT"};

    WCRootData *inData = new WCRootData();
    WCRootData *outData = new WCRootData();

    inData->ReadFile(fInFileName.c_str(), listWCRootEvt);
    outData->CreateTree(fOutFileName.c_str(), listWCRootEvt);

    const int nEntries = inData->GetEntries();
    const float toffset = 0.; // for IWCD pile-up event generation 

    cout<<" Start processing " << nEntries <<" entries........ " <<endl;
    for(int iEntry=0; iEntry<nEntries; iEntry++)
    {
        inData->GetEntry(iEntry);

        outData->AddTracks(inData->GetTrigger(0, 0), toffset, 0);
        for(int j=0; j<NPMTType; j++)
        {
            inData->AddTrueHitsToMDT(MDT->GetHitTubeCollection(fPMTType[j]), MDT->GetPMTResponse(fPMTType[j]), toffset, j);
            MDT->DoAddDark(fPMTType[j]);
            MDT->DoDigitize(fPMTType[j]);
        }
        
        // !!!!!!!!!!! Trigger is done only for the 20'' PMT, and the digi hits of both PMT types are added based on that !!!!!!!!!!!
        MDT->DoTrigger(fPMTType[0]);
        TriggerInfo* ti20 = MDT->GetTriggerInfo(fPMTType[0]);
        outData->AddDigiHits(MDT->GetHitTubeCollection(fPMTType[0]), ti20, iEntry, 0);
        outData->AddDigiHits(MDT->GetHitTubeCollection(fPMTType[1]), ti20, iEntry, 1);

        outData->FillTree();

        MDT->DoInitialize();
    }
    outData->WriteTree();
    inData->CloseFile();
}


bool ParseCmdArguments(int argc, char **argv)
{
    cout<<" Parsing command line arguments..... " <<endl;
	for(int i=1; i<argc; i++)
	{
		cout<<"     - argv[" << i <<"] :" <<argv[i] <<endl;
	}
	cout<<endl;

	for(int i=1; i<argc; i++)
	{
		if( string(argv[i])=="-i" )			            { fInFileName   = TString( argv[i+1] );	i++;}
        else if( string(argv[i])=="-p" )			    { fParFileName  = TString( argv[i+1] );	i++;}
		else if( string(argv[i])=="-o" )		        { fOutFileName	= TString( argv[i+1] );	i++;}
		else if( string(argv[i])=="-s" )		        { fSeed	        = atoi( argv[i+1] ); i++;}
		else if( string(argv[i])=="-n" )		        { fNEvtToProc	= atoi( argv[i+1] ); i++;}
		else
		{
            cout<<" i: " << argv[i] <<endl;
			return false;
		}
	}
	return true;
}

