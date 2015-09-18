void run_IGammaHadron()
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIGammaHadronTree.so");
  gSystem->Load("libIGammaHadronAna.so");


  TStopwatch t;
  t.Start();

  double ptt[] = {2.,3.,5, 10.,30.};
  double pta[] = {2.,3.,5.,10.,30.};

  IConsts ic(8);
  ic.SetInDir("/data2/IPdstReco_taxi168/Run8/");
  ic.SetOutDir("./");
  //ic.SetRunList("Run5kt_finalGoodRuns.txt");
  ic.SetRunList("Run8pp_runlist.txt");


  ic.SetTriggerPtBins(4, ptt);
  ic.SetAssociatedPtBins(4, pta);
  //ic.Print();


  IGammaHadron *corr = new IGammaHadron("IGammaHadron_Run8pp.root");

  //corr->Write("IGammaHadron_test.root");
  //delete corr;

  cout << endl;
  t.Stop();
  t.Print();
}
