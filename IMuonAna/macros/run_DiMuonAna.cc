void run_DiMuonAna()
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIMuonAna.so");


  TStopwatch t;
  t.Start();


  IConsts ic(8);
  //ic.SetInDir("/data3/IPdstReco_Run8pp200_taxi205/Muon/");
  ic.SetInDir("/home/iyounus/workdir/IPdstShort_Run8pp_Muon/");
  ic.SetOutDir("./");
  ic.SetRunList("Run8_goodRunList_JPsiAN.txt");
  ic.SetSpinFile("../IPdstAna/Run8GL1P_data_11Oct09.root");
  ic.SetFilePrefix("IPdstMuonAN_");


  double ptt[] = {0., 6.};
  ic.SetTriggerPtBins(1, ptt);

  IDiMuonAna *dimu = new
    IDiMuonAna("diMuonAN_run8pp_JpsiMassRange.root");

  cout << endl;
  t.Stop();
  t.Print();
}
