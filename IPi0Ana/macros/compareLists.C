
{
  gROOT->Reset();

  int nList1=0, List1[10000];
  ifstream fin1("tempList.txt");
  while (fin1 >> List1[nList1])
    {
      fin1.ignore(80, '\n');
      nList1++;
    }

  cout << nList1 << endl;
  fin1.close();

  int nList2=0, List2[10000];
  ifstream fin2("Run5FinalBadRunList_noRelLumi.txt");
  while (fin2 >> List2[nList2])
    nList2++;

  cout << nList2 << endl;
  fin2.close();


  for (int i=0; i<nList1; i++)
    {
      bool goodrun = false;
      for (int j=0; j<nList2; j++)
	if (List1[i] == List2[j])
	  goodrun = true;

      if (!goodrun) cout << List1[i] << endl;
    }
}
