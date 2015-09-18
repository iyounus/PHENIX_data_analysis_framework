

void runExample()
{
  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");

  gSystem->
    SetIncludePath("-I/home/iyounus/workdir/analysis/unm_pdst/gammaHadronAna "
		   "-I/home/iyounus/workdir/analysis/unm_pdst/IContainers "
		   "-I/home/iyounus/workdir/analysis/unm_pdst/IPdstAna ");


  gROOT->ProcessLine(".L Example.C+");


  IConsts *iconsts = new IConsts(8);

  iconsts->SetInDir("/data2/IPdstReco_taxi149/ERT/");
  iconsts->SetRunList("testRunList.txt");
  //iconsts->SetPTA(2.,5.);


  Example();

}
