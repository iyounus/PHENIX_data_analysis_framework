
void test_ISpinPattern()
{
  gROOT->Reset();

  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");


  // ISpinPattern *sp = new
  //   ISpinPattern("Run6GL1P_data_11Oct09.root",
  // 		 "../ITextFiles/Run6kt_finalGoodRuns.txt","Fill");

  // for (int i=7622; i<7958; i++) sp->Print(i);


  ISpinPattern *sp = new
    ISpinPattern("Run8GL1P_data_13Jan2011.root",
  		 "../ITextFiles/Run8pp_taxi174_runlist.txt","Fill");

  for (int i=9884; i<10001; i++) sp->Print(i);

  // TGraphErrors *gr = sp->RelLumiGraph(ISpinPattern::R_ANYELLOW);
  // gr->Draw("AP");
  // gr->Draw("P");


  // TCanvas *cc = new TCanvas("cc","cc",1200,900);
  // cc->Divide(1,3);

  // cc->cd(1);
  // sp->DrawBBC(258664);
  // cc->cd(2);
  // sp->DrawSpinB(258664);
  // cc->cd(3);
  // sp->DrawSpinY(258664);


}
