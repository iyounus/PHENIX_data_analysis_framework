
{
  gROOT->Reset();

  ifstream fin("../IPdstShort/taxi216_runList_ert.txt");

  ofstream fout1("run8pp_1.csh");
  ofstream fout2("run8pp_2.csh");
  ofstream fout3("run8pp_3.csh");
  ofstream fout4("run8pp_4.csh");
  ofstream fout5("run8pp_5.csh");


  int run;

  while (fin >> run)
    {
      fout1 << "root -l -q \'run_IPi0Reco.cc(" << run
	    << ")\' > logFiles/" << run << ".log" << endl;

      fin >> run;
      fout2 << "root -l -q \'run_IPi0Reco.cc(" << run
      	    << ")\' > logFiles/" << run << ".log" << endl;

      fin >> run;
      fout3 << "root -l -q \'run_IPi0Reco.cc(" << run
      	    << ")\' > logFiles/" << run << ".log" << endl;

      fin >> run;
      fout4 << "root -l -q \'run_IPi0Reco.cc(" << run
      	    << ")\' > logFiles/" << run << ".log" << endl;

      fin >> run;
      fout5 << "root -l -q \'run_IPi0Reco.cc(" << run
      	    << ")\' > logFiles/" << run << ".log" << endl;

    }

  fout1.close();
  fout2.close();
  fout3.close();
  fout4.close();
  fout5.close();

}
