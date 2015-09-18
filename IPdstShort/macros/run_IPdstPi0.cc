void run_IPdstPi0(int run=259568)
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");
  gSystem->Load("libIPdstShort.so");


  TStopwatch t;
  t.Start();

  IConsts ic(8);

  ic.SetInDir("/data3/IPdstReco_Run8pp200_taxi216/ERT/");
  //ic.SetOutDir("/data1/IPdst_Run8pp200_taxi174_GammaAN/");
  //ic.SetInDir("/data2/IPdstReco_Run8dAu_taxi149/ERT/");
  //ic.SetOutDir("/data2/IPdst_Run8dAu_corrTrees/ERT/");
  ic.SetOutDir("/data1/IPdstPi0_Run8pp200ERT/");
  ic.SetFilePrefix("IPdstReco_ERT_");

  ic.Print();

  IPdstShort *pdst = new IPdstShort(run, "IPdstPi0_");

  pdst->SetMuonFlag(false);
  pdst->SetZdcSmdFlag(false);
  //pdst->SetBbcFlag(false);
  pdst->SetTrackFlag(false);
  pdst->SetMinNPhotons(2);

  pdst->SetEMCalWarnMap("/home/iyounus/workdir/analysis/unm_pdst_ana/IGammaAN/"
			"EMCal_WarnMap_Run8_11Oct10.dat");


  pdst->Process();


  cout << endl;
  t.Stop();
  t.Print();
}
