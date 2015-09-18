void run_IPi0Reco(int run)
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  //gSystem->Load("libIGammaHadronTree.so");
  gSystem->Load("libIPi0Ana.so");

  TStopwatch t;
  t.Start();

  IConsts ic(8);

  double pt[] = {2., 3., 4., 5., 6., 8., 10., 15.,};
  ic.SetTriggerPtBins(7, pt);


  ic.SetInDir("/data1/IPdstPi0_Run8pp200ERT/");
  ic.SetOutDir("/data1/IPdstPi0Histos_Run8pp200ERT/");

  ic.SetFilePrefix("IPdstPi0_");
  ic.Print();

  IPi0Reco *pi0 = new IPi0Reco(run);

  //delete pi0;

  cout << endl;
  t.Stop();
  t.Print();
}
