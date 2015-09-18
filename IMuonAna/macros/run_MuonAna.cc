
void run_MuonAna()
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIMuonAna.so");


  IConsts ic(8);
  ic.SetInDir("/data2/IPdstReco_Run8pp200/Muon/");
  ic.SetOutDir("./");
  ic.SetOutFile("test_single_muons_run8pp.root");
  ic.SetRunList("test_runlist.txt");

  double ptt[] = {2.,4.,5.,7.,10.,15.};
  ic.SetTriggerPtBins(5, ptt);

  IMuonAna();
}
