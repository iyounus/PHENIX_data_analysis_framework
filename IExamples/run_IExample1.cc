{
  // load the liberaries
  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");

  // set include paths before compiling your code
  gSystem->AddIncludePath("-I/home/iyounus/workdir/analysis/unm_pdst_ana/"
			  "IContainers/");
  gSystem->AddIncludePath("-I/home/iyounus/workdir/analysis/unm_pdst_ana/"
			  "IPdstAna/");

  // compile your code
  gROOT->ProcessLine(".L IExample1.cc+");

  TStopwatch t;
  t.Start();

  // define IConsts
  IConsts ic(8); // this ctor takes a number but is not used anywhere
                 // this is just to specify which RHIC run is being analysed.

  // specify the location of pdsts
  ic.SetInDir("/media/WD/IPdstReco_Run8pp_taxi256/ERT/");
  ic.SetOutDir("./");

  // for this example, I'm reading ert trigger data
  ic.SetFilePrefix("IPdstReco_ERT_");

  // run your code
  IExample1(259575);

  cout << endl;
  t.Stop();
  t.Print();
}
