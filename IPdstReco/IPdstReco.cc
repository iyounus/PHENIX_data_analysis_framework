#include "IPdstReco.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include "PHCompositeNode.h"
#include "PHIODataNode.h"
#include "RunHeader.h"
#include "EventHeader.h"
#include "PHGlobal.h"
#include "PHCentralTrack.h"
#include "PHMuoTracksOut.h"
#include "emcClusterContainer.h"
#include "emcClusterContent.h"
#include "ErtOut.h"
#include "SmdOut.h"
#include "SpinDataEventOut.h"
#include "TrigLvl1.h"
#include "MWGVertex.h"
#include "lpcRaw.h"
#include "Fun4AllReturnCodes.h"

#include "IHeader.hh"
#include "IHeaderList.hh"
#include "IPhoton.hh"
#include "IPhotonList.hh"
#include "ITrack.hh"
#include "ITrackList.hh"
#include "IElectron.hh"
#include "IElectronList.hh"
#include "IMuon.hh"
#include "IMuonList.hh"
#include "IDiMuon.hh"
#include "IDiMuonList.hh"
#include "IBbc.hh"
#include "IBbcList.hh"
#include "IZdcSmd.hh"
#include "IZdcSmdList.hh"


#include "TString.h"
#include "TFile.h"
#include "TTree.h"


#define NTOWERS        24768  // total number of towers in EMCal
#define MAXTREESIZE    700*1048576 // 700MB

#define DUMMYPOSINT     9999
#define DUMMYNEGINT    -9999

#define DUMMYPOSFLOAT   9999.0
#define DUMMYNEGFLOAT  -9999.0

#define ZVTXERR        2.0 // for p+p the z vtx error is 2 cm.
//#define ZVTXERR        0.5 // for heavy ions the z vtx error is 0.5 cm.

#define MaxNoMuons     10 // maximum number of muons after cuts. I don't
                          // expect to be more than 10 in p+p.


using namespace std;


typedef PHIODataNode<TObject>              TObjectNode_t;
typedef PHIODataNode<PHGlobal>             PHGlobalNode_t;
typedef PHIODataNode<PHCentralTrack>       PHCntNode_t;
typedef PHIODataNode<PHMuoTracksOut>       PHMuoNode_t;
typedef PHIODataNode<RunHeader>            RunHeaderNode_t;
typedef PHIODataNode<EventHeader>          EventHeaderNode_t;
typedef PHIODataNode<emcClusterContainer>  PHEmcClusterNode_t;
typedef PHIODataNode<ErtOut>               ErtOutNode_t;
typedef PHIODataNode<SmdOut>               SmdOutNode_t;
typedef PHIODataNode<SpinDataEventOut>     SpinDataNode_t;
typedef PHIODataNode<TrigLvl1>             TrigLvl1Node_t;
typedef PHIODataNode<lpcRaw>               lpcRaw_t;


IPdstReco::IPdstReco(char* outFileName):
  _cZVtxCentral(35.),
  _cZVtxMuon(999.),
  _cMaxZedDC(80.),
  _cMinTrackPt(0.2),
  _cMinPhotonEn(0.1),
  _cMinPhotonProb(0.02),
  _cMinTrackHighPt(0.2),
  _cMinPhotonHighEn(0.8), // this cut does not effect the number of pi0s.
  _bSpin(true),
  _bErt(false),       // by default, ert trig data is assumed, see Init()
  _bMuTrig(false),    // muon trig flag is flase
  _bMinBias(false),   // min bias trig flag is false
  _bPhoton(true),
  _bTrack(true),
  _bMuon(true),
  _bDiMuon(true),
  _bBbc(true),
  _bZdcSmd(true)
{
  cout << "IPdstReco::IPdstReco  ---------------------------------------------"
       << endl;
  ThisName = "IPdstReco";

  _outFileName = outFileName;

  d_evthdr   = NULL;
  d_cnttrack = NULL;
  d_muotrack = NULL;
  d_emcclust = NULL;
  d_ertout   = NULL;
  d_smdout   = NULL;
  d_spindata = NULL;
  d_triglvl1 = NULL;
  d_lpcRaw   = NULL;

  d_photon = NULL;
  d_track  = NULL;
  d_elect  = NULL;
  d_muon   = NULL;
  d_dimuon = NULL;
  d_bbc    = NULL;
  d_zdcsmd = NULL;

  _runNumber   = 0;
  _eventNumber = DUMMYNEGINT;
  _nTotalEvts  = DUMMYNEGINT;
  _ival = 0;
  _ErrCount = 0; // this counts how many times a node wasn't found. It stops
                 // the error to be printed for every event;

  d_header = new IHeaderList();

  _outFile = new TFile(outFileName, "RECREATE");

  _pdst = new TTree("pDST","pDST for UNM LANL analyses");
  _pdst->SetMaxTreeSize(MAXTREESIZE);

  // header branch is always created. but it will not be filled if
  // no other particle branch is created.
  _pdst->Branch("header",  "IHeaderList",  &d_header);

  cout << "Maximum tree size .... " << _pdst->GetMaxTreeSize() << endl;
  cout << "IPdstReco::IPdstReco done  ----------------------------------------"
       << endl;
  return;
}
//-----------------------------------------------------------------------------


void IPdstReco::InitTree()
{
  cout << "IPdstReco::InitTree  ----------------------------------------------"
       << endl;

  if (_bPhoton)
    {
      d_photon = new IPhotonList();
      _pdst->Branch("photon",  "IPhotonList",  &d_photon);
      cout << "Adding photon branch." << endl;
    }

  if (_bTrack)
    {
      d_track  = new ITrackList();
      _pdst->Branch("track",   "ITrackList",   &d_track);
      cout << "Adding track branch." << endl;

      d_elect  = new IElectronList();
      _pdst->Branch("electron","IElectronList",&d_elect);
      cout << "Adding electron branch." << endl;
    }

  if (_bMuon)
    {
      d_muon   = new IMuonList();
      _pdst->Branch("muon",    "IMuonList",    &d_muon);
      cout << "Adding muon branch." << endl;
    }


  if (_bDiMuon)
    {
      d_dimuon   = new IDiMuonList();
      _pdst->Branch("dimuon",  "IDiMuonList",  &d_dimuon);
      cout << "Adding dimuon branch." << endl;
    }


  if (_bBbc)
    {
      d_bbc    = new IBbcList();
      _pdst->Branch("bbc",     "IBbcList",     &d_bbc);
      cout << "Adding bbc branch." << endl;
    }

  if (_bZdcSmd)
    {
      d_zdcsmd = new IZdcSmdList();
      _pdst->Branch("zdcsmd",  "IZdcSmdList",  &d_zdcsmd);
      cout << "Adding zdcsmd branch." << endl;
    }


  cout << "IPdstReco::InitTree done  -----------------------------------------"
       << endl;
}
//-----------------------------------------------------------------------------


IPdstReco::~IPdstReco()
{  
  cout << "IPdstReco::~IPdstReco" << endl;
  delete d_header;

  if (_bPhoton)   delete d_photon;
  if (_bTrack)
    {
      delete d_track;
      delete d_elect;
    }
  if (_bMuon)     delete d_muon;
  if (_bDiMuon)   delete d_dimuon;
  if (_bBbc)      delete d_bbc;
  if (_bZdcSmd)   delete d_zdcsmd;

  delete [] _EMC_badChan;
}
//-----------------------------------------------------------------------------


int IPdstReco::Init(PHCompositeNode *topNode)
{
  cout << "IPdstReco::Init  --------------------------------------------------"
       << endl;


  short trig = short(_bErt) + short(_bMuTrig) + short(_bMinBias);
  if (trig > 1)
    { // this is to make sure only one trigger is selected.
      cout << "Confusing trigger assignemtent. exiting .... " << endl;
      exit(EXIT_FAILURE);
    }

  // if no trigger is selected by the user then select ert
  if (trig == 0) _bErt = true;

  if (_bErt)     { _bMuTrig = false; _bMinBias = false;}
  if (_bMuTrig)  { _bErt = false; _bMinBias = false;}
  if (_bMinBias) { _bErt = false; _bMuTrig = false;}


  InitTree();

  _nTotalEvts = 0;

  cout << "Zvtx for central arm " << _cZVtxCentral << endl;
  cout << "Zvtx for muon arm    " << _cZVtxMuon << endl;
  cout << "Max ZedDC            " << _cMaxZedDC << endl;
  cout << "Min track pt         " << _cMinTrackPt << endl;
  cout << "Min photon En        " << _cMinPhotonEn << endl;
  cout << "Min photon prob      " << _cMinPhotonProb << endl;
  cout << "At least one track  should have Pt > " << _cMinTrackHighPt  << endl;
  cout << "At least one photon should have En > " << _cMinPhotonHighEn << endl;
  if (!_bSpin)   cout << "SpinDataNode_t will not be read" << endl;
  if (_bMinBias) cout << "Analysing MinBias data .... " << endl;
  if (_bMuTrig)  cout << "Analysing Muon trigger data .... " << endl;
  if (_bErt)     cout << "Analysing ERT data .... " << endl;
  if (!_bErt)    cout << "ERT will not be required for photons" << endl;


  cout << "IPdstReco::Init done  ---------------------------------------------"
       << endl;

  return 0;
}
//-----------------------------------------------------------------------------


int IPdstReco::InitRun(PHCompositeNode *topNode)
{
  cout << "IPdstReco::InitRun  -----------------------------------------------"
       << endl;
  //Get runnumber, check to see if it is a good run
  RunHeader *d_runhdr=NULL;
  PHTypedNodeIterator<RunHeader> iRUN(topNode);

  RunHeaderNode_t *RUN = iRUN.find("RunHeader");

  if (RUN) d_runhdr = RUN->getData(); 
  if (!d_runhdr) 
    {
      cout << PHWHERE 
	   << "IPdstReco::InitRun ERROR RunHeader not in Node Tree" 
	   << endl;

      return ABORTRUN;
    }

  _runNumber = d_runhdr->get_RunNumber();
  cout << "RunNumber: " << _runNumber << endl;

  _eventNumber = 0;

  TString warnmap =
    "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_EdgeTowers.dat";

  if (_runNumber >= 167415 && _runNumber <= 179846) // Run5 pp ERT/MB
    warnmap = 
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run5.dat";

  if (_runNumber >= 189579 && _runNumber <= 204639) // Run6 pp 200 ERT/MB
    warnmap = 
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run6.dat";

  if (_runNumber >= 205154 && _runNumber <= 206495) // Run6 pp 62 ERT/MB
    warnmap =
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run6.dat";

  if (_runNumber >= 246214 && _runNumber <= 253701) // Run8 dAu ERT/MB
    warnmap = 
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run8.dat";

  if (_runNumber >= 256450 && _runNumber <= 259575) // Run8 pp ERT/MB
    warnmap = 
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run8.dat";

  if (_runNumber >= 276324 && _runNumber <= 280242) // Run9 pp 500 GeV ERT/MB
    warnmap = 
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run9_500GeV.dat";

  if (_runNumber >= 281911 && _runNumber <= 291579) // Run9 pp 200 GeV ERT/MB
    warnmap =
      "/phenix/u/workarea/iyounus/EMCalWarnMaps/EMCal_WarnMap_Run9_200GeV.dat";


  ReadEMCalWarnMap(warnmap.Data());


  cout << "IPdstReco::InitRun done  ------------------------------------------"
       << endl;
  return 0;
}
//-----------------------------------------------------------------------------


int IPdstReco::process_event(PHCompositeNode *topNode)
{
  _nTotalEvts++;

  PHTypedNodeIterator<PHGlobal> iGLOBAL(topNode);
  PHGlobalNode_t *global = iGLOBAL.find("PHGlobal");
  if(global) d_global = (PHGlobal*)global->getData();
  if(!d_global)
    {
      cout << PHWHERE
	   << "IPdstReco::PHGlobal data not in the Node Tree" 
	   << endl;

      return ABORTEVENT;
    }

  PHTypedNodeIterator<EventHeader> iEvtHdr(topNode);
  EventHeaderNode_t *evthdr = iEvtHdr.find("EventHeader");
  if(evthdr) d_evthdr = (EventHeader*)evthdr->getData();
  if(!d_evthdr)
    {
      cout << PHWHERE
	   << "IPdstReco::EventHeader data not in the Node Tree"
	   << endl;

      return ABORTEVENT;
    }

  PHTypedNodeIterator<PHCentralTrack> iCNT(topNode);
  PHCntNode_t *cnt = iCNT.find("PHCentralTrack");
  if(cnt) d_cnttrack = (PHCentralTrack*)cnt->getData();
  if(!d_cnttrack && _ErrCount < 7)
    {
      cout << PHWHERE 
	   << "IPdstReco::PHCentralTrack data not in the Node Tree" 
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<PHMuoTracksOut> iMUO(topNode);
  PHMuoNode_t *muo = iMUO.find("PHMuoTracksOO");
  if(muo) d_muotrack = (PHMuoTracksOut*)muo->getData();
  if(!d_muotrack && _ErrCount < 7)
    {
      cout << PHWHERE
	   << "IPdstReco::PMuoTracksOut data not in the Node Tree"
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<emcClusterContainer> iEMC(topNode);
  PHEmcClusterNode_t *emcCluster = iEMC.find("emcClusterContainer");
  if(emcCluster) d_emcclust = (emcClusterContainer*)emcCluster->getData();
  if(!d_emcclust && _ErrCount < 7)
    {
      cout << PHWHERE 
	   << "IPdstReco::emcClusterContainer data not in the Node Tree" 
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<SpinDataEventOut> iSPIN(topNode);
  SpinDataNode_t *spin = iSPIN.find("SpinDataEventOut");
  if(spin) d_spindata = spin->getData();
  if (_bSpin && !d_spindata && _ErrCount < 7)
    {
      cout << PHWHERE
	   << "IPdstReco::SpinDataEventOut not in Node Tree" 
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<ErtOut> iERT(topNode);
  ErtOutNode_t *ert = iERT.find("ErtOut");
  if(ert) d_ertout = ert->getData();
  if (!d_ertout && _ErrCount < 7)
    {
      cout << PHWHERE
	   << "IPdstReco::ErtOut not in Node Tree" 
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<SmdOut> iSMD(topNode);
  SmdOutNode_t *smd = iSMD.find("SmdOut");
  if(smd) d_smdout = smd->getData();
  if (!d_smdout && _ErrCount < 7)
    {
      cout << PHWHERE
	   << "IPdstReco::SmdOut not in Node Tree" 
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<TrigLvl1> iTRIG(topNode);
  TrigLvl1Node_t *trig = iTRIG.find("TrigLvl1");
  if(trig) d_triglvl1 = trig->getData();
  if (!d_triglvl1 && _ErrCount < 7)
    {
      cout << PHWHERE
	   << "IPdstReco::TriggerHelper not in Node Tree" 
	   << endl;
      _ErrCount++;
    }

  PHTypedNodeIterator<lpcRaw> iLPC(topNode);
  lpcRaw_t *lpc = iLPC.find("lpcRaw");
  if(lpc) d_lpcRaw = lpc->getData();
  if (!d_lpcRaw && _ErrCount < 7)
    {
      cout << PHWHERE
	   << "IPdstReco::lpcRaw not in Node Tree" 
	   << endl;
      _ErrCount++;
    }

  float zVertex = d_global->getBbcZVertex();
  if (zVertex == 0) return ABORTEVENT; // for some historic reasons

  // first fill Muons, Track and Photon, If there are no particles in
  // an event, then, no need to fill the header, Bbc etc.

  short nPhotons=0, nTracks=0, nMuons=0;
  // fill Muons first because zVtx cut is different for muons
  if (_bMuon && d_muotrack) nMuons   = FillMuons();

  // for muon trigger data, every event should have at least one muon
  if (_bMuTrig && nMuons == 0) return ABORTEVENT;

  // wider zvertex cut for muons
  if (nMuons > 0  && fabs(zVertex) > _cZVtxMuon)    return ABORTEVENT;
  if (nMuons == 0 && fabs(zVertex) > _cZVtxCentral) return ABORTEVENT;

  // NOTE FillTracks() must be called before FillPhotons(), because in 
  // FillPhotons(), towerID comparision is made for photons and tracks.

  if (_bTrack  && d_cnttrack) nTracks  = FillTracks();
  if (_bPhoton && d_emcclust) nPhotons = FillPhotons();

  if (_bErt && (nPhotons+nTracks)==0) return ABORTEVENT;


  //create a header and fill it
  IHeader *hdr = d_header->AddHeader();

  hdr->SetRunID(_runNumber);
  hdr->SetEventID(d_evthdr->get_EvtSequence());
  hdr->SetZVertex(zVertex);

  if (d_spindata)
    {
      hdr->SetGL1CrossingID(d_spindata->GetGL1CrossingID());
      hdr->SetSpinGL1CrossingID(d_spindata->GetSpinGL1CrossingID());
    }
  else
    {
      hdr->SetGL1CrossingID(DUMMYPOSINT);
      hdr->SetSpinGL1CrossingID(DUMMYPOSINT);
    }

  if (d_triglvl1)
    {
      hdr->SetTrigScaled(d_triglvl1->get_lvl1_trigscaled());
      hdr->SetTrigLive  (d_triglvl1->get_lvl1_triglive());
    }
  else
    {
      // this is unsigend int. so, if no trigger available, I don't want to
      // set it to 0. but set it to the largest integer
      hdr->SetTrigScaled((unsigned int)(4294967295.));
      hdr->SetTrigLive  ((unsigned int)(4294967295.));
    }

  if (_bBbc)    FillBbc();
  if (_bZdcSmd) FillZdcSmd();

  _pdst->Fill();
  _eventNumber++;

  return EVENT_OK;
}
//-----------------------------------------------------------------------------


short IPdstReco::FillTracks()
{
  float maxPt = DUMMYNEGFLOAT;
  short iTrack = 0;

  for(int itk=0; itk<(int)d_cnttrack->get_npart(); itk++)
    {
      short quality = d_cnttrack->get_quality(itk);
      if (quality <= 3) continue;

      float mom = d_cnttrack->get_mom(itk);
      if (mom > 200.) continue;   // sanity check;

      float the0 = d_cnttrack->get_the0(itk);
      float pt = mom * sin(the0);

      if (pt < _cMinTrackPt) continue;

      float zed = fabs(d_cnttrack->get_zed(itk));
      if (zed > _cMaxZedDC) continue;

      float pc3sdz   = d_cnttrack->get_pc3sdz(itk);
      float pc3sdphi = d_cnttrack->get_pc3sdphi(itk);
      float emcsdz   = d_cnttrack->get_emcsdz(itk);
      float emcsdphi = d_cnttrack->get_emcsdphi(itk);

      if (pc3sdz   == DUMMYNEGINT && emcsdz   == DUMMYNEGINT) continue;
      if (pc3sdphi == DUMMYNEGINT && emcsdphi == DUMMYNEGINT) continue;


      if (maxPt < pt) maxPt = pt;  // this is to chosse the highest pt track

      ITrack *track = d_track->AddTrack();

      track->SetCharge(short(d_cnttrack->get_charge(itk)));
      track->SetQuality(quality);


      short dcarm = short(d_cnttrack->get_dcarm(itk));
      short arm  = (dcarm + 1)%2;
      short sect = short(d_cnttrack->get_sect(itk));
      short iz   = short(d_cnttrack->get_zsect(itk));
      short iy   = short(d_cnttrack->get_ysect(itk));

      track->SetTowerID(TowerID(arm, sect, iy, iz));

      short ertTrig = GetERT(arm, sect, iy, iz);
      track->SetErt(ertTrig);

      track->SetMomentum(mom);
      track->SetTheta0(the0);
      track->SetPhi0(d_cnttrack->get_phi0(itk));
      track->SetPhiDC(d_cnttrack->get_phi(itk));
      track->SetZedDC(d_cnttrack->get_zed(itk));
      track->SetAlpha(d_cnttrack->get_alpha(itk));

      track->SetPc3sdZ(pc3sdz);
      track->SetPc3sdPhi(pc3sdphi);
      track->SetEmcsdZ(emcsdz);
      track->SetEmcsdPhi(emcsdphi);

      track->SetEmcE(d_cnttrack->get_emce(itk));
      track->SetEcore(d_cnttrack->get_ecore(itk));
      track->SetProb(d_cnttrack->get_prob(itk));

      // save electron specific variables.
      short N0 = short(d_cnttrack->get_n0(itk));
      short N1 = short(d_cnttrack->get_n1(itk));

      if ( N0 >= 0 || N1 >=0 )
	{
	  IElectron *elect = d_elect->AddElectron();

	  elect->SetTrack(iTrack);
	  elect->SetN0(N0);
	  elect->SetN1(N1);

	  elect->SetNpe0(d_cnttrack->get_npe0(itk));
	  elect->SetChi2(d_cnttrack->get_chi2(itk));

	  elect->SetRICHdisp(d_cnttrack->get_disp(itk));
	  elect->SetEmcsdPhiE(d_cnttrack->get_emcsdphi_e(itk));
	  elect->SetEmcsdZE(d_cnttrack->get_emcsdz_e(itk));
	}
      iTrack++;
    }//loop over cnt tracks

  // if the highest pt track is less than _cMinTrackHighPt GeV,
  // then then all tracks are discarded.
  if (maxPt < _cMinTrackHighPt) d_track->Reset();

  return d_track->GetNTracks();
}
//-----------------------------------------------------------------------------


short IPdstReco::FillPhotons()
{
  /*
    NOTE about ecore():
    Ecore was introduced to clean up noisy channels in run02 AuAu data,
    and it is reconstructed from the cluster core (max energy towers) only.
    For any "electromagnetic" analysis the Ecore  should be used,
    for "charged hadronic" analysis like pi^\pm or high_pT DC match with 
    EMCal one should use E.
  */

  float maxEn = DUMMYNEGFLOAT;

  float zVtx = d_global->getBbcZVertex();

  short nERT = 0;

  for(unsigned int iclus=0; iclus<d_emcclust->size(); iclus++)
    {
      emcClusterContent *emcContent = d_emcclust->getCluster(iclus);

      // NOTE emcContent->towerid(0)] give the same towerID as the funciton
      // IPdstReco::TowerID. Confirmed with run8 data on 24Jan2010

      int twID = emcContent->towerid(0);

      if (_EMC_badChan[twID]) continue;
      if (emcContent->prob_photon() < _cMinPhotonProb) continue;

      float ecore = emcContent->ecore();
      if (ecore < _cMinPhotonEn) continue;
      if (ecore > 200.) continue;     // sanity chack

      //if (CompareTowerID(twID)) continue;

      IPhoton *pht = d_photon->AddPhoton();

      pht->SetTowerID(short(twID));
      pht->SetNTowers(short(emcContent->multiplicity()));

      pht->SetEnergy(ecore);

      pht->SetProb(emcContent->prob_photon());
      pht->SetTof(emcContent->tof());

      pht->SetX(emcContent->x());
      pht->SetY(emcContent->y());

      float z = emcContent->z() - zVtx;
      pht->SetZ(z);

      short iz   = short(emcContent->izpos());
      short iy   = short(emcContent->iypos());
      short arm  = short(emcContent->arm());
      short sect = short(emcContent->sector());

      short ertTrig = GetERT(arm, sect, iy, iz);
      pht->SetErt(ertTrig);

      nERT += short(ertTrig > 0);
      if (maxEn < ecore) maxEn = ecore;

      /*
      cout << setw(5) << nPhotons << setw(10) << ecore
	   << setw(10) << emcContent->prob_photon()
	   << setw(8) << emcContent->iypos()
	   << setw(8) << emcContent->izpos()
	   << setw(8) << emcContent->towerid(0)
	   << setw(5) << ertTrig
	   << setw(5) << _EMC_badChan[emcContent->towerid(0)]
	   << endl;
      */
    }//loop over the clusters

  // if the highest energy photon is less than _cMinPhotonHighEn GeV,
  // then the event is not useful for prompt photon analysis
  if (maxEn < _cMinPhotonHighEn) d_photon->Reset();

  // If not photon has ERT trigger, no need to record photons.
  // These photons cannot be used to pi0 or prompt photon ana.
  if (_bErt && nERT == 0) d_photon->Reset();

  return d_photon->GetNPhotons();;
}
//-----------------------------------------------------------------------------


short IPdstReco::FillMuons()
{
  float zVtx = d_global->getBbcZVertex();

  int iSelected[MaxNoMuons];

  int count = 0;
  for (unsigned int imuo=0; imuo < d_muotrack->get_npart(); imuo++)
    {
      // this is to calculate dg0
      float Xmut = d_muotrack->get_xpos(4, imuo); // kalman projection
      float Ymut = d_muotrack->get_ypos(4, imuo); // kalman projection
      float Zmut = d_muotrack->get_zpos(4, imuo); // kalman projection


      if (Xmut == 0 && Ymut == 0 && Zmut == 0)
	{
	  // if kalman projection is not available then use these
	  Xmut = d_muotrack->get_xpos(3, imuo); // at station 3
	  Ymut = d_muotrack->get_ypos(3, imuo); // at station 3
	  Zmut = d_muotrack->get_zpos(3, imuo); // at station 3
	}

      float dxdz = d_muotrack->get_px(3, imuo)/d_muotrack->get_pz(3, imuo);
      float dydz = d_muotrack->get_py(3, imuo)/d_muotrack->get_pz(3, imuo);

      // select best muid road
      float min_dg0 = DUMMYPOSFLOAT;
      int bestRoad = -9;
      for (int iroad=0; iroad<3; iroad++)
	{
	  if (d_muotrack->get_muIDOOhits(iroad, imuo) == 0) continue;

	  float Xmui = d_muotrack->get_muIDOO_gap0(0, iroad, imuo);
	  float Ymui = d_muotrack->get_muIDOO_gap0(1, iroad, imuo);
	  float Zmui = d_muotrack->get_muIDOO_gap0(2, iroad, imuo);

	  float Xproj = Xmut + dxdz*(Zmui - Zmut);
	  float Yproj = Ymut + dydz*(Zmui - Zmut);

	  float dg0 = sqrt( (Xproj-Xmui)*(Xproj-Xmui) +
			    (Yproj-Ymui)*(Yproj-Ymui) );

	  if (min_dg0 > dg0)
	    {
	      min_dg0 = dg0;
	      bestRoad = iroad;
	    }
	}

      if (min_dg0 == DUMMYPOSFLOAT)
	continue;              // if no road is found, no need to save muon 


      IMuon *muon = d_muon->AddMuon();

      muon->SetCharge(short(d_muotrack->get_charge(imuo)));
      muon->SetMuTrChi2(d_muotrack->get_chisquare(imuo));
      muon->SetMuTrHitPat((unsigned short)d_muotrack->get_muTRhits(imuo));

      for (int i=0; i<5; i++) // 5 positions: vtx, st1 st2, st3, kalman proj
	{
	  muon->SetX(i, d_muotrack->get_xpos(i, imuo));
	  muon->SetY(i, d_muotrack->get_ypos(i, imuo));
	  muon->SetZ(i, d_muotrack->get_zpos(i, imuo));

	  muon->SetPx(i, d_muotrack->get_px(i, imuo));
	  muon->SetPy(i, d_muotrack->get_py(i, imuo));
	  muon->SetPz(i, d_muotrack->get_pz(i, imuo));
	}

      //set muid parameters for best road.
      muon->SetMuIdHitPat((unsigned short)
			  d_muotrack->get_muIDOOhits(bestRoad, imuo));
      muon->SetMuIdChi2  (d_muotrack->get_muIDOOchi (bestRoad, imuo));

      muon->SetX(IMuon::GP0, d_muotrack->get_muIDOO_gap0(0, bestRoad, imuo));
      muon->SetY(IMuon::GP0, d_muotrack->get_muIDOO_gap0(1, bestRoad, imuo));
      muon->SetZ(IMuon::GP0, d_muotrack->get_muIDOO_gap0(2, bestRoad, imuo));

      muon->SetDx( d_muotrack->get_muIDOO_gap0(3, bestRoad, imuo));
      muon->SetDy( d_muotrack->get_muIDOO_gap0(4, bestRoad, imuo));

      // vertex refit. adopted from MWGpico code SngmuonsRun8pp.C

      iSelected[count++] = imuo; // this is needed for dimuon vtx refit
      MWGVertex vertex;
      vertex.set_verbosity(0);
      vertex.add_track(imuo, d_muotrack);
      vertex.add_vertex(zVtx, ZVTXERR);
      vertex.fit();

      muon->SetPx(IMuon::VTX, vertex.get_px(0));
      muon->SetPy(IMuon::VTX, vertex.get_py(0));
      muon->SetPz(IMuon::VTX, vertex.get_pz(0));
      muon->SetVtxRefitChi2(vertex.get_chisquare()/vertex.get_ndf());
    }

  if (_bDiMuon && count > 1) FillDiMuons(count, iSelected);

  return d_muon->GetNMuons();
}
//-----------------------------------------------------------------------------


void IPdstReco::FillDiMuons(int nMu, int *MuIndex)
{
  float zVtx = d_global->getBbcZVertex();

  int idx[2];
  int iMu[2];

  bool MuSelected = false;

  for (int idimu=0; idimu<d_muotrack->get_ndimu(); idimu++)
    {
      idx[0] = d_muotrack->get_ditrkIndex(0,idimu);
      idx[1] = d_muotrack->get_ditrkIndex(1,idimu);

      // check if the muons have been selected in FillMuon function
      // if both muons are already selected, it implies that they passed
      // single muon cuts
      MuSelected = false;
      for (int i=0; i<nMu; i++)
	if (MuIndex[i] == idx[0]) {MuSelected = true; iMu[0] = i;}
      if (!MuSelected) continue;

      MuSelected = false;
      for (int i=0; i<nMu; i++)
	if (MuIndex[i] == idx[1]) {MuSelected = true; iMu[1] = i;}
      if (!MuSelected) continue;

      IDiMuon *dimu = d_dimuon->AddDiMuon();

      for (int k=0; k<2; k++)
	{
	  IMuon *muon = d_muon->GetMuon(iMu[k]);

	  for (int l=0; l<6; l++)
	    {
	      dimu->SetX(k, l, muon->GetX(l));
	      dimu->SetY(k, l, muon->GetY(l));
	      dimu->SetZ(k, l, muon->GetZ(l));
	    }

	  for (int l=1; l<5; l++)
	    {
	      dimu->SetPx(k, l, muon->GetPx(l));
	      dimu->SetPy(k, l, muon->GetPy(l));
	      dimu->SetPz(k, l, muon->GetPz(l));
	    }

	  dimu->SetDx(k, muon->GetDx());
	  dimu->SetDy(k, muon->GetDy());

	  dimu->SetCharge(k, muon->GetCharge());

	  dimu->SetMuTrHitPat(k, muon->GetMuTrHitPat());
	  dimu->SetMuIdHitPat(k, muon->GetMuIdHitPat());

	  dimu->SetMuTrChi2(k, muon->GetMuTrChi2());
	  dimu->SetMuIdChi2(k, muon->GetMuIdChi2());
	}

      MWGVertex vertex;
      vertex.set_verbosity(0);

      // add tracks
      vertex.add_track(idx[0], d_muotrack);
      vertex.add_track(idx[1], d_muotrack);
      vertex.add_vertex(zVtx, ZVTXERR);

      vertex.fit();

      dimu->SetPx(0, IMuon::VTX, vertex.get_px(0));
      dimu->SetPy(0, IMuon::VTX, vertex.get_py(0));
      dimu->SetPz(0, IMuon::VTX, vertex.get_pz(0));

      dimu->SetPx(1, IMuon::VTX, vertex.get_px(1));
      dimu->SetPy(1, IMuon::VTX, vertex.get_py(1));
      dimu->SetPz(1, IMuon::VTX, vertex.get_pz(1));

      dimu->SetVtxRefitChi2(vertex.get_chisquare()/vertex.get_ndf());
      dimu->SetVtxRefitZ   (vertex.get_vtx_z());
    }
}
//-----------------------------------------------------------------------------


void IPdstReco::FillBbc()
{
  IBbc *bbc = d_bbc->AddBbc();

  bbc->SetBbcT0(d_global->getBbcTimeZero());
  bbc->SetBbcMultN(short(d_global->getBbcMultN()));
  bbc->SetBbcMultS(short(d_global->getBbcMultS()));
  bbc->SetBbcChargeN(d_global->getBbcChargeN());
  bbc->SetBbcChargeS(d_global->getBbcChargeS());
}
//-----------------------------------------------------------------------------


void IPdstReco::FillZdcSmd()
{
  float ZdcN_E = d_global->getZdcEnergyN();
  float ZdcS_E = d_global->getZdcEnergyS();

  // this is to make sure if the information exists
  if (ZdcN_E > 5. || ZdcS_E > 5.)
    {
      IZdcSmd *zdcsmd = d_zdcsmd->AddZdcSmd();

      zdcsmd->SetZdcN_E(ZdcN_E);
      zdcsmd->SetZdcS_E(ZdcS_E);

      zdcsmd->SetSmdN_E(d_global->get_SmdEN());
      zdcsmd->SetSmdS_E(d_global->get_SmdES());

      zdcsmd->SetSmdN_X(d_global->get_SmdXN());
      zdcsmd->SetSmdS_X(d_global->get_SmdXS());
      zdcsmd->SetSmdN_Y(d_global->get_SmdYN());
      zdcsmd->SetSmdS_Y(d_global->get_SmdYS());

      if (d_lpcRaw)
	{
	  zdcsmd->SetScintN_ADC(d_lpcRaw->get_AdcPost(2) -
				d_lpcRaw->get_AdcPre(2));
	  zdcsmd->SetScintS_ADC(d_lpcRaw->get_AdcPost(0) -
				d_lpcRaw->get_AdcPre(0));
	  zdcsmd->SetScintN_TDC(d_lpcRaw->get_Tdc0(2));
	  zdcsmd->SetScintS_TDC(d_lpcRaw->get_Tdc0(0));
	}
      else
	{
	  zdcsmd->SetScintN_ADC(DUMMYPOSINT);
	  zdcsmd->SetScintS_ADC(DUMMYPOSINT);
	  zdcsmd->SetScintN_TDC(DUMMYPOSINT);
	  zdcsmd->SetScintS_TDC(DUMMYPOSINT);
	}

      if (d_smdout)
	{
	  for (int i=0; i<8; i++)
	    {
	      zdcsmd->SetSmdN_chargeY(i, d_smdout->get_Charge(i));
	      zdcsmd->SetSmdS_chargeY(i, d_smdout->get_Charge(i+16));
	    }

	  for (int i=0; i<7; i++)
	    {
	      zdcsmd->SetSmdN_chargeX(i, d_smdout->get_Charge(i+8));
	      zdcsmd->SetSmdS_chargeX(i, d_smdout->get_Charge(i+24));
	    }
	}
      else
	{
	  for (int i=0; i<8; i++)
	    {
	      zdcsmd->SetSmdN_chargeY(i, DUMMYNEGINT);
	      zdcsmd->SetSmdS_chargeY(i, DUMMYNEGINT);
	    }

	  for (int i=0; i<7; i++)
	    {
	      zdcsmd->SetSmdN_chargeX(i, DUMMYNEGINT);
	      zdcsmd->SetSmdS_chargeX(i, DUMMYNEGINT);
	    }
	}
    }
}
//-----------------------------------------------------------------------------


int IPdstReco::ResetEvent(PHCompositeNode *topNode)
{
  //clearing the data for the next event
  d_header->Reset();
  if (_bPhoton)   d_photon->Reset();
  if (_bTrack)
    {
      d_track->Reset();
      d_elect->Reset();
    }

  if (_bMuon)     d_muon->Reset();
  if (_bDiMuon)   d_dimuon->Reset();
  if (_bBbc)      d_bbc->Reset();
  if (_bZdcSmd)   d_zdcsmd->Reset();


  if (_nTotalEvts%50000 == 0)
    {
      cout << "nEvent " << _nTotalEvts << " and nGoodEvents " << _eventNumber
	   << " = " << double(_eventNumber)/double(_nTotalEvts) << endl;

      //cout << "Tree size .... " << _pdst->GetTotBytes() << "\t" 
      //<< _pdst->GetZipBytes() << endl;
    }

  return 0;
}
//-----------------------------------------------------------------------------


int IPdstReco::End(PHCompositeNode *topNode)
{
  //write out the ntuples, etc.
  cout << "IPdstReco::End  ---------------------------------------------------"
       << endl;
  cout << "nEvent " << _nTotalEvts << " and nGoodEvents " << _eventNumber 
       << " = " << (double)_eventNumber/(double)_nTotalEvts << endl;
  cout << "Tree size .... " << _pdst->GetTotBytes() << "\t"
       << _pdst->GetZipBytes() << endl;


  _outFile = _pdst->GetCurrentFile();
  cout << "\nWRITING ...........\n" << _outFile->GetName() << "\n" << endl;

  _outFile->Write();

  _outFile->Close();
  delete _outFile;

  cout << "IPdstReco::End done  ----------------------------------------------"
       << endl;

  return 0;
}
//-----------------------------------------------------------------------------


void IPdstReco::ReadEMCalWarnMap(const char* file)
{
  // read bad towers from the file
  cout << "IPdstReco::ReadEMCalWarnMap" << endl;
  cout << "Reading EMCal warnmap from " << endl;
  cout << file << endl;

  _EMC_badChan = new bool[NTOWERS];
  for (int i=0; i<NTOWERS; i++)
    _EMC_badChan[i] = false;

  ifstream fin(file);
  int towerID = 0, badFlag = 0, nChan = 0;

  while (fin >> towerID >> badFlag)
    if (badFlag > 0)
      {
        _EMC_badChan[towerID] = true;
        nChan++;
      }

  fin.close();

  cout << "Number of masked (hot/dead/edge) EMCal towers: "
       << nChan << endl;
}
//-----------------------------------------------------------------------------


short IPdstReco::GetERT(short arm, short sect, short ypos, short zpos)
{
  if (!d_ertout) return DUMMYNEGINT;

  /*
    trigger bit in ErtOut are 
    0:4x4a, 1:4x4b, 2:4x4c, 3:2x2, 4:RICH

    So, here I'm calculating ert as a 5 digit decimel number as

    RICH*10^4 + 2x2*10^3 + 4x4c*10^2 + 4x4b*10 + 4x4a

    This is how you can get which trigger was fired:
    for example, if you need to check 4x4c:

    (ert%pow(10, 3) - ert%pow(10, 2))/ert%pow(10, 2)
  */

  short sm = GetEMCalSM(arm, sect, ypos, zpos);
  short ert = 0;

  for (int i=0; i<5; i++)
    ert += pow(10, i)*d_ertout->get_ERTbit(i, int(arm), int(sect), int(sm));

  return ert;
}
//-----------------------------------------------------------------------------


bool IPdstReco::CompareTowerID(short towerID)
{
  if (d_track->GetNTracks() < 1) return false;

  bool trackExists = false;

  for (unsigned int i = 0; i<d_track->GetNTracks(); i++)
    if (towerID == d_track->GetTrack(i)->GetTowerID())
      {
	trackExists = true;
	break;
      }

  return trackExists;
}
//-----------------------------------------------------------------------------
