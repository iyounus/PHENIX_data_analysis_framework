
void run_IHadronAna()
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIHadronAna.so");


  IConsts ic(8);
  ic.SetInDir("/data2/Run9CondorOutput_VtxFit/");
  ic.SetOutDir("./");
  ic.SetOutFile("whatever.root");
  ic.SetRunList("run9_templist.txt");

  double ptt[] = {4.,5.,7.,10.,15.};
  ic.SetTriggerPtBins(4, ptt);


  IHadronAna *hd = new IHadronAna();

}
