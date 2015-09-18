void run_IPdstMuon(int run=257651)
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIPdstShort.so");


  TStopwatch t;
  t.Start();

  IConsts ic(8);

  //ic.SetInDir("/data3/IPdstReco_Run8pp200_taxi216/Muon/");
  ic.SetInDir("/media/data/IPdstReco_Run8pp200/Muon");
  ic.SetOutDir("/home/iyounus/workdir/IPdstShor_Run8pp_Muon/");
  ic.SetFilePrefix("IPdstReco_MU_");


  ic.Print();

  IPdstShort *pdst = new IPdstShort(run, "IPdstMuonAN_");

  pdst->SetPhotonFlag(false);
  pdst->SetTrackFlag(false);
  pdst->SetZdcSmdFlag(false);
  pdst->SetBbcFlag(false);

  // This means that there must be atleast one muon in the final event
  // but number of photons or tracks are irrelevent
  //pdst->SetMinNPhotons(0);
  //pdst->SetMinNTracks(0);

  pdst->SetMinNMuons(1);

  pdst->Process();


  cout << endl;
  t.Stop();
  t.Print();
}
