void run_IPdstGamma(int run=259568)
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIPdstShort.so");

  gSystem->AddIncludePath("-I/home/iyounus/workdir/analysis/unm_pdst_ana/"
			  "IContainers/");
  gSystem->AddIncludePath("-I/home/iyounus/workdir/analysis/unm_pdst_ana/"
			  "IPdstAna/");


  IConsts ic(8);

  ic.SetInDir("/media/WD/IPdstReco_Run8pp_taxi256/ERT/");
  ic.SetOutDir("/home/iyounus/workdir/IPdstGamma_Run8pp200ERT_taxi256/");
  ic.SetFilePrefix("IPdstReco_ERT_");


  ic.Print();

  IPdstShort *pdst = new IPdstShort(run, "IPdstDirectGamma_");

  pdst->SetMinPhotonHighPt(4.0);

  // pdst->SetEMCalWarnMap("/home/iyounus/workdir/analysis/unm_pdst_ana/IGammaAN/"
  // 			"EMCal_WarnMap_Run8_noGaurdVeto.dat");


  pdst->SetMuonFlag(false);
  pdst->SetZdcSmdFlag(false);
  pdst->SetElectronFlag(false);

  TStopwatch t;
  t.Start();


  pdst->Process();

  cout << endl;
  t.Stop();
  t.Print();
}
