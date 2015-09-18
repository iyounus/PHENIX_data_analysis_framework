
{
  gSystem->Load("libIContainers.so");
  gSystem->Load("libIPdstAna.so");

  ifstream fin("../IGammaAN/Run8pp_GammaAN_runList.txt");


  IConsts c(8);

  IFillLookup *ff = new IFillLookup();

  int run;

  while (fin >> run)
    {
      cout << run << "  " << ff->GetFill(run) << endl;
    }

}
