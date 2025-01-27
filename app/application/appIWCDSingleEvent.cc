// This is a MDT application for IWCD single event
// This takes true hits in the input file that is an output of WCSim 
// and digitize those hits with the same method as used in WCSim (i.e. /DAQ/Digitizer SKI)
// Optionally, dark hits can be added to the true hits before the digitization is done.
// In addition, the same triggering algorithm as implemeted in WCSim (i.e. /DAQ/Trigger NDigits)
// The output of this code is a ROOT file which has the same format of the input file
// The digitized hits in the output file should be consistent with the ones in the input file
// wihtin a statistical fluctuation (, provided that the impelementation in MDT is exactly the same as in WCSim)

#include <iostream>
#include <string>
#include <vector>

#include "MDTManager.h" 
#include "WCRootData.h"

// PMT type used for 3-inch PMTs of IWCD/WCTE
#include "PMTResponse3inchR12199_02.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

int fNEvtToProc = -1;
int fSeed = 67592;
string fParFileName = "";
string fInFileName = "";
string fOutFileName = "";
bool fUseOD = false;

bool ParseCmdArguments(int, char**);

int main(int argc, char **argv)
{
    ParseCmdArguments(argc, argv);

    Configuration *Conf = Configuration::GetInstance();
    Conf->ReadParameter(fParFileName);
    Conf->PrintParameters();
    Conf = 0;

	// IWCD PMT type for ID and OD
    int NPMTType = !fUseOD ? 1 : 2 ;
    vector<string> fPMTType(NPMTType);
    fPMTType[0] = "PMT3inchR12199_02";
    if (fUseOD) fPMTType[1] = "PMT3inchR12199_02_OD";

    MDTManager *MDT = new MDTManager(fSeed);
    MDT->RegisterPMTType(fPMTType[0], new PMTResponse3inchR12199_02());
    if (fUseOD) MDT->RegisterPMTType(fPMTType[1], new PMTResponse3inchR12199_02());

    vector<string> listWCRootEvt(NPMTType);
    listWCRootEvt[0] = "wcsimrootevent";
    if (fUseOD) listWCRootEvt[1] = "wcsimrootevent_OD";

	// WCRootData is an interface class between MDT and WCSim root file
    WCRootData *inData = new WCRootData();
    WCRootData *outData = new WCRootData();

	// Read input WCSim file
    inData->ReadFile(fInFileName.c_str(), listWCRootEvt);

	// Create TTrees to be saved in the output file
    outData->CreateTree(fOutFileName.c_str(), listWCRootEvt);

    int nEntries = inData->GetEntries();
    const float toffset = 0.; // for IWCD pile-up event generation -> just ignore

    if (fNEvtToProc>0 && fNEvtToProc<nEntries) nEntries = fNEvtToProc;

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
            MDT->DoTrigger(fPMTType[j]);

            TriggerInfo *ti = MDT->GetTriggerInfo(fPMTType[j]);
            outData->AddDigiHits(MDT->GetHitTubeCollection(fPMTType[j]), ti, iEntry, j);
        }
        outData->FillTree();

        MDT->DoInitialize();
    }
    outData->WriteTree();
    inData->CloseFile();

    // Copy geo trees, etc.
    vector<string> listTrees{"Settings","wcsimGeoT","wcsimRootOptionsT"};
    for (auto s:listTrees) outData->CopyTree(fInFileName.c_str(),s.c_str());
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
        else if( string(argv[i])=="-od" )		        { fUseOD = true;}
		else
		{
            cout<<" i: " << argv[i] <<endl;
			return false;
		}
	}
	return true;
}
