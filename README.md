# MDT
C++ Stand-alone library for Merging, Digitizing, and Triggering hits of photomultiplier tubes (PMTs) based on [WCSim](https://github.com/WCSim/WCSim). A Python interface is optionally provided, requiring pre-installation of [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html).

The library provides C++ classes that manage three tasks:

 - Merging true hits (simulated photoelectrons by WCSim) of an events with those of any number of events
    - Intrinsic dark noise can be added during this stage
 - Digitizing true hits
    - Timing and charge of a digitized hit is simulated by a simple model used in WCSim
    - The library is capable of accomodating alternative models 
    - Parameters describing PMT charcteristics such as timing resolution can be varied
    - Digitized hits by PMT afterpusing can be added
 - Triggering digitized hits
    - A simple algorithm that counts number of digitized hits falling in a sliding time winodw
    - Selectable trigger window and thershold

## WCTE/WCSim usage
Simple runtime procedures. First set up your ROOT and WCSIM:
```
source your_thisroot.sh
source your_this_wcsim.sh # or just export WCSIM_BUILD_DIR
```
Then set up the MDT environment.
```
source envMDT.sh
cd $MDTROOT/cpp; make clean; make all
cd $MDTROOT/app/utilities/WCRootData; make clean; make all
cd $MDTROOT/app/application; make clean; make all
cd $MDTROOT
# edit variables properly in run_test_mdt4wcte.sh
bash run_test_mdt4wcte.sh
```

## IWCD usage
Essentially the same usage but with OD support.
```
# edit variables properly in run_test_mdt4iwcd.sh
bash run_test_mdt4iwcd.sh
```

## IWCD pile-up generation
Example usage
```
$MDTROOT/app/application/appGenPileUpSpill $MDTROOT/example/genPileUpConfig.txt
```
The application generates pile-up events by combining ID neutrino interaction events and beam background events in a spill. 

Input variables are:
- `ListIDNuIntFiles`,`ListBeamBkgFiles`: list of WCSim output files for the ID and background events
- `OutFileNamePrefix`: output file name will be something like "OutFileNamePrefix".00000.root
- `MDTParFile`: MDT config file
- `InitialSeed`: Random seed
- `IDNuIntRate`,`BeamBkgRate`: Mean number of ID and background events per spill. In each spill, the actual number of interactions are drawn from a Poisson distribution, and interaction timing according to the bunch structure (see `BeamTiming` class under `app/utilities/WCRootData/`)
- `UseOD`: Process OD hits
- `NumOfSpillsSavedPerFile`, `TotalNumOfSpills`: output spill setup

## How to simulate and access digitized pulses
For each true hit, a digitized waveform is simulated by sampling the single PE pulse (defined by the `<WaveformFile>` parameter) every 8 ns with 1 mV resolution. If there is another PE arriving within the same pulse window, the waveforms are added. 

To do pulse fitting, the peak of each pulse is found, then a Gaussian fit is done and the fitted parameters are used to calculate the digitized time and charge.

The waveform of the first pulse of each PMT in each event is saved in `TClonesArray`. To read the pulses,
```
// open the file and get the digitzed waveform tree
TTree* wcsimDigiWFTree = (TTree*)f->Get("wcsimDigiWFTree");
TClonesArray *arr = new TClonesArray("TH1F");
wcsimDigiWFTree->GetBranch("wcsimrootevent_waveform")->SetAutoDelete(kFALSE);
wcsimDigiWFTree->SetBranchAddress("wcsimrootevent_waveform",&arr);
// In each event, each array index corresponds to PMT id (from 0 to nPMTs-1)
wcsimDigiWFTree->GetEntry(0); // event id
TH1F* h = (TH1F*)arr->At(0); // PMT id
```
