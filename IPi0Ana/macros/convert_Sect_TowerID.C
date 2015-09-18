
{
  gROOT->Reset();


  ifstream fin("run5_emctof_tower_by_tower_corrections.txt");
  ofstream fout("Run5_emcTofCorr_tower_by_tower.dat");

  int sect=-1, ypos=-1, zpos=-1;
  float mean=0;

  fin.ignore(1024, '\n');

  int towerID = -99;
  while (fin >> sect >> ypos >> zpos >> mean)
    {
      fin.ignore(500, '\n');

      if (sect < 4)
	towerID = zpos + 72*ypos + sect*2592;

      if (sect > 5)
	towerID = zpos + 72*ypos + (sect-2)*2592;

      if (sect == 4 || sect == 5)
	towerID = zpos + 96*ypos + (sect-4)*4608 + 15552;


      fout << towerID << "  " << mean << "  " 
	   << sect << "  " << ypos << "  " << zpos << endl;
    }
  fout.close();

  fin.close();
}
