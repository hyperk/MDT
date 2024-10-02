#include "WCRootData.h"
#include "Configuration.h"

#include "TBranch.h"

WCRootData::WCRootData()
{
    fDetCentreY = 0.;
    fWCGeom = 0;
    fWCSimC = 0;
    fWCSimT = 0;
    fSpEvt.clear();
    isOD.clear();
    fOutFileName = "";

    int mult_flag = 1;
    fHitTimeOffset = 0.;

    Configuration *Conf = Configuration::GetInstance();
    Conf->GetValue<float>("TimeOffset", fHitTimeOffset);
    Conf->GetValue<int>("FlagMultDigits", mult_flag);

    fMultDigiHits = bool(mult_flag);
}

WCRootData::~WCRootData()
{
    if( fWCGeom ){ delete fWCGeom; fWCGeom = 0; }
    if( fWCSimT ){ delete fWCSimT; fWCSimT = 0; }
    if( fWCSimC ){ delete fWCSimC; fWCSimC = 0; }
    fSpEvt.clear();
    fSpEvt.shrink_to_fit();
    isOD.clear();
    isOD.shrink_to_fit();
}


void WCRootData::GetGeometryInfo(const char *filename)
{
    TFile *f=TFile::Open(filename);
    TTree *t=(TTree*)f->Get("wcsimGeoT");
    t->SetBranchAddress("wcsimrootgeom", &fWCGeom);
    t->GetEntry( 0 );
    fDetCentreY = fWCGeom->GetWCOffset(1);
    cout<<" DetCentreY: " << fDetCentreY <<endl;
    f->Close();
}


void WCRootData::AddTrueHitsToMDT(MDTManager *mdt)
{
    HitTubeCollection *hc = mdt->GetHitTubeCollection();
    this->AddTrueHitsToMDT(hc, mdt->GetPMTResponse(), 0., 0);
    hc = 0;
}


void WCRootData::AddTrueHitsToMDT(HitTubeCollection *hc, PMTResponse *pr, float intTime, int iPMT)
{
	// Set all tubes no matter how many of true hits, anyway
	// Position and orientation are also set to each tube
	this->SetTubes(hc, iPMT);

    const WCSimRootTrigger *aEvt = fSpEvt[iPMT]->GetTrigger(0);
    const int nCkovHits = aEvt->GetNcherenkovhits();
	//std::cout<<" nCkovHits: " << nCkovHits <<std::endl;
    TClonesArray *hitTimeArray = aEvt->GetCherenkovHitTimes();
    for(int iHit=0; iHit<nCkovHits; iHit++)
    {
        WCSimRootCherenkovHit *aHit=
                dynamic_cast<WCSimRootCherenkovHit*>((aEvt->GetCherenkovHits())->At(iHit));
    
        int timeArrayIndex = aHit->GetTotalPe(0);
        int peForTube = aHit->GetTotalPe(1);   
        int tubeID = aHit->GetTubeID();
        const int maxPE = timeArrayIndex + peForTube;
        for(int jPE=timeArrayIndex; jPE<maxPE; jPE++)
        {
            WCSimRootCherenkovHitTime *aHitTime=
                dynamic_cast<WCSimRootCherenkovHitTime*>(hitTimeArray->At(jPE));
    
            float truetime = aHitTime->GetTruetime();
            if( truetime<0. ){ continue; }
            if( aHitTime->GetParentID()<0 ){ continue; }

            TrueHit *th = new TrueHit(truetime+intTime, aHitTime->GetParentID());

            for(int k=0; k<3; k++){ th->SetPosition(k, aHitTime->GetPhotonEndPos(k)); }
            for(int k=0; k<3; k++){ th->SetDirection(k, aHitTime->GetPhotonEndDir(k)); }
            for(int k=0; k<3; k++){ th->SetStartDirection(k, aHitTime->GetPhotonStartDir(k)); }

            th->SetStartTime(aHitTime->GetPhotonStartTime()+intTime);
            for(int k=0; k<3; k++){ th->SetStartPosition(k, aHitTime->GetPhotonStartPos(k)); }
            th->SetCreatorProcess((int)(aHitTime->GetPhotonCreatorProcess()));
            if( !pr->ApplyDE(th) ){ continue; }

            (&(*hc)[tubeID])->AddRawPE(th);
        }
    }
}

void WCRootData::ReadFile(const char *filename, const vector<string> &list)
{
    fWCSimC = new TChain("wcsimT");
    fWCSimC->Add(filename);

	const char *filename0 = NULL;
	filename0 = fWCSimC->GetFile()->GetName();
	if( filename0!=NULL )
	{
		this->GetGeometryInfo(filename0);
	}

    if( list.size()==0 ) // default
    {
        fSpEvt.push_back( 0 );
        isOD.push_back(false);
        fWCSimC->SetBranchAddress("wcsimrootevent", &fSpEvt[0]);
    }
    else
    {
        fSpEvt.clear();
        fSpEvt = vector<WCSimRootEvent*>(list.size(),0);
        isOD.clear();
        isOD = vector<bool>(list.size(),false);
        for(unsigned int i=0; i<list.size(); i++)
        {
            fWCSimC->SetBranchAddress(list[i].c_str(), &fSpEvt[i]);
            if ( list[i].find("OD")!=std::string::npos ) isOD[i] = true;
        }
    }
    fWCSimC->SetAutoDelete();
}

void WCRootData::CloseFile()
{
    delete fWCSimC; fWCSimC = 0;
    fSpEvt.clear();
    fSpEvt.shrink_to_fit();
    isOD.clear();
    isOD.shrink_to_fit();
}

int WCRootData::GetEntries() 
{
    return fWCSimC->GetEntries();
}

void WCRootData::GetEntry(int i)
{
    fWCSimC->GetEntry(i);
}


void WCRootData::CreateTree(const char *filename, const vector<string> &list)
{
    fOutFileName = TString(filename);
    TFile *fout = new TFile(fOutFileName, "recreate");
    fout->SetCompressionLevel(2);

    Int_t branchStyle = 1;
    Int_t bufferSize = 64000;
    TString bName = "wcsimrootevent";
    TString bAddress = "WCSimRootEvent";

    fWCSimT = new TTree("wcsimT", "Processed WCSim Tree for pile-up");
    if( list.size()==0 )
    {
        fSpEvt.push_back( new WCSimRootEvent() );
        fSpEvt[0]->Initialize();
        TTree::SetBranchStyle(branchStyle);
        fWCSimT->Branch(bName, bAddress, &fSpEvt[0], bufferSize, 2);
    }
    else  
    {
        fSpEvt.clear();
        fSpEvt = vector<WCSimRootEvent*>(list.size(), 0);
        isOD.clear();
        isOD = vector<bool>(list.size(),false);
        for(unsigned int i=0; i<list.size(); i++)
        {
//            fSpEvt.push_back( new WCSimRootEvent() );
            fSpEvt[i] = new WCSimRootEvent();
            fSpEvt[i]->Initialize();
            TTree::SetBranchStyle(branchStyle);
            fWCSimT->Branch(list[i].c_str(), bAddress, &fSpEvt[i], bufferSize, 2);
            if ( list[i].find("OD")!=std::string::npos ) isOD[i] = true;
        }
    }
}

void WCRootData::AddDigiHits(MDTManager *mdt, int eventID, int iPMT)
{
    HitTubeCollection *hc = mdt->GetHitTubeCollection();
    TriggerInfo *ti = mdt->GetTriggerInfo();
    this->AddDigiHits(hc, ti, eventID, iPMT);
}

void WCRootData::AddDigiHits(HitTubeCollection *hc, TriggerInfo *ti, int eventID, int iPMT)
{
    WCSimRootTrigger* anEvent = fSpEvt[iPMT]->GetTrigger(0);
    // Save raw hits
    // container for photon info
    std::vector<double> truetime;
    std::vector<int>   primaryParentID;
    std::vector<float> photonStartTime;
    std::vector<TVector3> photonStartPos;
    std::vector<TVector3> photonEndPos;
    std::vector<TVector3> photonStartDir;
    std::vector<TVector3> photonEndDir;
    std::vector<ProcessType_t> photonCreatorProcess;
    for(hc->Begin(); !hc->IsEnd(); hc->Next())
    {
        // Get tube ID
        HitTube *aPH = &(*hc)();
        int tubeID = aPH->GetTubeID();
        int mPMTID = aPH->GetmPMTID();
        int mPMT_PMTID = aPH->GetmPMT_PMTID();


        const int NPE = aPH->GetNRawPE();
        const vector<TrueHit*> PEs = aPH->GetPhotoElectrons();
        for(int iPE=0; iPE<NPE; iPE++)
        {
            truetime.push_back(PEs[iPE]->GetTime());
            primaryParentID.push_back(PEs[iPE]->GetParentId());
            photonStartTime.push_back(PEs[iPE]->GetStartTime());
            photonStartPos.push_back(TVector3(PEs[iPE]->GetStartPosition(0),PEs[iPE]->GetStartPosition(1),PEs[iPE]->GetStartPosition(2)));
            photonEndPos.push_back(TVector3(PEs[iPE]->GetPosition(0),PEs[iPE]->GetPosition(1),PEs[iPE]->GetPosition(2)));
            photonStartDir.push_back(TVector3(PEs[iPE]->GetStartDirection(0),PEs[iPE]->GetStartDirection(1),PEs[iPE]->GetStartDirection(2)));
            photonEndDir.push_back(TVector3(PEs[iPE]->GetDirection(0),PEs[iPE]->GetDirection(1),PEs[iPE]->GetDirection(2)));
            photonCreatorProcess.push_back((ProcessType_t)(PEs[iPE]->GetCreatorProcess()));
        }

        anEvent->AddCherenkovHit(tubeID,
                                mPMTID,
                                mPMT_PMTID,
                                truetime,
                                primaryParentID,
                                photonStartTime,
                                photonStartPos,
                                photonEndPos,
                                photonStartDir,
                                photonEndDir,
                                photonCreatorProcess);

        truetime.clear();
        primaryParentID.clear();
        photonStartTime.clear();
        photonStartPos.clear();
        photonEndPos.clear();
        photonStartDir.clear();
        photonEndDir.clear();
        photonCreatorProcess.clear();
    }
    const int nTriggers = ti->GetNumOfTrigger(); 
    for(int iTrig=0; iTrig<nTriggers; iTrig++) 
    {
        TriggerType_t trigType = kTriggerNDigits;
        float hitTimeOffset = fHitTimeOffset;
        if( ti->GetType(iTrig)==(int)TriggerType::eFailure )
        {
            trigType = kTriggerFailure; 
        }
        else if ( ti->GetType(iTrig)==(int)TriggerType::eNoTrig )
        {
            trigType = kTriggerNoTrig; 
            hitTimeOffset = 0;
        }

        if( iTrig>=1 )
        {
            fSpEvt[iPMT]->AddSubEvent();
            anEvent = fSpEvt[iPMT]->GetTrigger(iTrig);
            anEvent->SetHeader(eventID, 0, 0, iTrig+1);
            anEvent->SetMode(0);
        }

        vector<Double_t> info(1, ti->GetNHits(iTrig));
        info.push_back(hitTimeOffset);
        info.push_back(ti->GetTriggerTime(iTrig));
        anEvent->SetTriggerInfo(trigType, info);

        const float triggerTime = ti->GetTriggerTime(iTrig);
        const float tWinLowEdge = ti->GetLowEdge(iTrig);
        const float tWinUpEdge = ti->GetUpEdge(iTrig);

        int nHits = 0;
        float sumQ = 0.;
        map<int, bool> tubeCount;
        for(hc->Begin(); !hc->IsEnd(); hc->Next())
        {
            HitTube *aPH = &(*hc)();
            int tubeID = aPH->GetTubeID();
            int mPMTID = aPH->GetmPMTID();
            int mPMT_PMTID = aPH->GetmPMT_PMTID();
            int nCount = 0;
            for(int i=0; i<aPH->GetNDigiHits(); i++)
            {
                float t = aPH->GetTimeDigi(i);

                if( t>=tWinLowEdge && t<=tWinUpEdge )
                {
                    bool doFill = false;
                 
                    // The first hit is always allowed to be saved.
                    // Later hits depend on your setting
                    if( tubeCount.count(tubeID)==0 ){ doFill = true; } 
                    else { doFill = fMultDigiHits ? true : false; }  

                    if( doFill ) 
                    {
                        vector<int> true_pe_comp = aPH->GetParentCompositionDigi(i);
                        t = t +  hitTimeOffset - triggerTime;
                        float q = aPH->GetChargeDigi(i);
                        anEvent->AddCherenkovDigiHit(q,
                                                     t,
                                                     tubeID,
                                                     mPMTID,
                                                     mPMT_PMTID,
                                                     true_pe_comp);
                        nHits += 1;
                        sumQ += q;
                        nCount+=1;
                        tubeCount[tubeID] = 1;
                    }
                }
            } // digi hit loop
            aPH = NULL;
        } // PMT loop
        anEvent->SetNumDigitizedTubes(nHits); 
        anEvent->SetSumQ(sumQ);
        WCSimRootEventHeader *eh = anEvent->GetHeader();
        eh->SetDate( int(triggerTime) );
    }
}

void WCRootData::FillTree()
{
    TFile *f = fWCSimT->GetCurrentFile();
    f->cd();
    fWCSimT->Fill();     
    fWCSimT->Write("",TObject::kOverwrite);
    for(unsigned int i=0; i<fSpEvt.size(); i++)
    {
        fSpEvt[i]->ReInitialize();
    }
}


void WCRootData::AddTracks(const WCSimRootTrigger *aEvtIn, float offset_time, int iPMT)
{
    WCSimRootTrigger *aEvtOut = fSpEvt[iPMT]->GetTrigger(0);
    TClonesArray *tracks = aEvtIn->GetTracks();
    const int ntrack = tracks->GetEntries();
    for(int i=0; i<ntrack; i++)
    {
        WCSimRootTrack *aTrack = (WCSimRootTrack*)tracks->At(i);
        Int_t     ipnu = aTrack->GetIpnu();
        Int_t     flag = aTrack->GetFlag();
        Int_t     startvol = aTrack->GetStartvol();
        Int_t     stopvol = aTrack->GetStopvol();
        Int_t     parenttype = aTrack->GetParenttype();
        Int_t     id = aTrack->GetId();
        Int_t     idPrnt = aTrack->GetParentId();
        ProcessType_t creatorProcess = aTrack->GetCreatorProcess();

        Double_t   dir[3];
        Double_t   pdir[3];
        Double_t   stop[3];
        Double_t   start[3];
        Double_t   m = aTrack->GetM();
        Double_t   p = aTrack->GetP();
        Double_t   E = aTrack->GetE();
        Double_t   time = aTrack->GetTime() + offset_time;
        for(int j=0; j<3; j++)
        {
            dir[j] = aTrack->GetDir(j); 
            pdir[j] = aTrack->GetPdir(j);
            stop[j] = aTrack->GetStop(j);
            start[j] = aTrack->GetStart(j);
        }

        std::vector<std::vector<float>> bPs = aTrack->GetBoundaryPoints();
        std::vector<float> bKEs = aTrack->GetBoundaryKEs();
        std::vector<double> bTimes = aTrack->GetBoundaryTimes();
        std::vector<int> bTypes = aTrack->GetBoundaryTypes();

        aEvtOut->AddTrack(ipnu, 
			  flag, 
			  m, 
			  p, 
			  E, 
			  startvol, 
			  stopvol, 
			  dir, 
			  pdir, 
			  stop,
			  start,
			  parenttype,
              creatorProcess,
			  time,
			  id,
			  idPrnt,
              bPs,
              bKEs,
              bTimes,
              bTypes);

    }
}

void WCRootData::WriteTree()
{
    TFile *f = fWCSimT->GetCurrentFile();
    f->cd();
    fWCSimT->Write();
    f->Close();
}


void WCRootData::CopyTree(const char *filename,
                          const char *treename,
                          const vector<bool> &savelist)
{

	TFile *fin = TFile::Open(filename);
	if( !fin->FindKey(treename) )
	{
		fin->Close();
		std::cout<<" [ERROR] WCRootData::CopyTree " <<std::endl;
		std::cout<<"  - No " << treename <<" tree is saved in : " << filename <<std::endl;
        exit(-1);
	}

    // Copy all the entries
    if( savelist.size()==0 )
    {
        TTree *tin = (TTree*)fin->Get(treename);

        TFile *fout = TFile::Open(fOutFileName, "update");
        fout->cd();
        TTree *tout = tin->CloneTree(-1, "fast");
        tout->Write();
        fout->Close();

        // Not sure why the special treatment is needed at the first place
        // if( strcmp(treename,"Settings")!=0 )
        // {
        //     TTree *tin = (TTree*)fin->Get(treename);

        //     TFile *fout = TFile::Open(fOutFileName, "update");
        //     fout->cd();
        //     TTree *tout = tin->CloneTree(-1, "fast");
        //     tout->Write();
        //     fout->Close();
        // }
        // else // This is a special treatment
        // {
        //     TTree *tin = (TTree*)fin->Get(treename);

        //     TBranch *bWCXRotation = tin->GetBranch("WCXRotation");
        //     TBranch *bWCYRotation = tin->GetBranch("WCYRotation");
        //     TBranch *bWCZRotation = tin->GetBranch("WCZRotation");
        //     TBranch *bWCDetCentre = tin->GetBranch("WCDetCentre");
        //     TBranch *bWCDetRadius = tin->GetBranch("WCDetRadius");
        //     TBranch *bWCDetHeight = tin->GetBranch("WCDetHeight");

        //     Float_t WCXRotation[3] = {0.};
        //     Float_t WCYRotation[3] = {0.};
        //     Float_t WCZRotation[3] = {0.};
        //     Float_t WCDetCentre[3] = {0.};
        //     Float_t WCDetRadius = 0.;
        //     Float_t WCDetHeight = 0.;

        //     bWCXRotation->SetAddress(WCXRotation);
        //     bWCYRotation->SetAddress(WCYRotation);
        //     bWCZRotation->SetAddress(WCZRotation);
        //     bWCDetCentre->SetAddress(WCDetCentre);
        //     bWCDetRadius->SetAddress(&WCDetRadius);
        //     bWCDetHeight->SetAddress(&WCDetHeight);

        //     bWCXRotation->GetEntry(0);
        //     bWCYRotation->GetEntry(0);
        //     bWCZRotation->GetEntry(0);
        //     bWCDetCentre->GetEntry(0);
        //     bWCDetRadius->GetEntry(0);
        //     bWCDetHeight->GetEntry(0);

        //     TFile *fout = TFile::Open(fOutFileName, "update");
        //     fout->cd();
        //     TTree *tout = new TTree("Settings", "");
        //     tout->Branch("WCXRotation", WCXRotation, "WCXRotation[3]/F");
        //     tout->Branch("WCYRotation", WCYRotation, "WCYRotation[3]/F");
        //     tout->Branch("WCZRotation", WCZRotation, "WCZRotation[3]/F");
        //     tout->Branch("WCDetCentre", WCDetCentre, "WCDetCentre[3]/F");
        //     tout->Branch("WCDetRadius", &WCDetRadius, "WCDetRadius/F");
        //     tout->Branch("WCDetHeight", &WCDetHeight, "WCDetHeight/F");
        //     tout->Fill();
        //     tout->Write();
        //     fout->Close();
        // }
    }
    else
    {
        TTree *tin = (TTree*)fin->Get(treename);

        TFile *fout = TFile::Open(fOutFileName, "update");
        fout->cd();
        TTree *tout = tin->CloneTree(0);

        const int nEntries = tin->GetEntries();
        for(int iEntry=0; iEntry<nEntries; iEntry++)
        {
            tin->GetEntry(iEntry);
            if( !savelist[iEntry] ){ continue; }
            fout->cd();
            tout->Fill();
        }
        fout->cd();
        tout->Write();
        fout->Close();
    }
    fin->Close();
}

void WCRootData::SetTubes(HitTubeCollection *hc, const int iPMT)
{
	const int nTubes = !isOD.at(iPMT) ? fWCGeom->GetWCNumPMT(bool(iPMT)) : fWCGeom->GetODWCNumPMT() ;
	for(int i=0; i<nTubes; i++)
	{
		const WCSimRootPMT *tube = !isOD.at(iPMT) ? fWCGeom->GetPMTPtr(i, bool(iPMT)) : fWCGeom->GetODPMTPtr(i) ;

		const int tubeID = tube->GetTubeNo();
        const int mPMTID = tube->GetmPMTNo();
        const int mPMT_PMTID = tube->GetmPMT_PMTNo();
		hc->AddHitTube(tubeID);
		for(int j=0; j<3; j++)
		{
			(&(*hc)[tubeID])->SetPosition(j, tube->GetPosition(j));
			(&(*hc)[tubeID])->SetOrientation(j, tube->GetOrientation(j));
		}
        (&(*hc)[tubeID])->SetmPMTID(mPMTID);
        (&(*hc)[tubeID])->SetmPMT_PMTID(mPMT_PMTID);
	}
}
