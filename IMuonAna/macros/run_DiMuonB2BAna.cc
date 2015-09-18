void run_DiMuonB2BAna()
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIMuonAna.so");


  TStopwatch t;
  t.Start();


  IConsts ic(8);
  ic.SetInDir("/home/iyounus/workdir/Run6CondorOutput/");
  ic.SetOutDir("./");
  ic.SetRunList("Run6_goodRunList_North_AN636.txt");
  ic.SetSpinFile("../IPdstAna/Run6GL1P_data_11Oct09.root");
  ic.SetFilePrefix("IPdstReco_Run6_MU_");


  double ptt[] = {0., 1.4, 6.};
  ic.SetTriggerPtBins(2, ptt);

  IDiMuonB2BAna *dimu = new
    IDiMuonB2BAna("diMuonB2B_run6pp.root");

  cout << endl;
  t.Stop();
  t.Print();
}
