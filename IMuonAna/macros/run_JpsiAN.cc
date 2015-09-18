void run_JpsiAN()
{
  gROOT->Reset();

  gSystem->
    SetIncludePath("-I/home/iyounus/workdir/analysis/unm_pdst_ana/IPdstAna/");


  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  //gSystem->Load("libIMuonAna.so");

  gROOT->ProcessLine(".L IJpsiAN.cc+");


  IConsts ic(8);
  ic.SetRunList("Run8_goodRunList_JPsiAN.txt");
  ic.SetSpinFile("../IPdstAna/Run8GL1P_data_21Mar10.root");


  double ptt[] = {0., 6.};
  ic.SetTriggerPtBins(1, ptt);

  IJpsiAN();

}
