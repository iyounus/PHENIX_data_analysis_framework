#include "IFillLookup.hh"

#include <fstream>
#include <iostream>
#include <assert.h>

#include "IConsts.hh"

using namespace std;

IFillLookup::IFillLookup()
{
  cout << "IFillLookup::IFillLookup" << endl;
  assert(IConsts::Defined);
  cout << "Reading " << IConsts::FillTable << endl;

  _fill = new map<int, int>;

  ifstream fin(IConsts::FillTable);

  int run, fill;

  while (fin >> run >> fill)
    (*_fill)[run] = fill;


  fin.close();
}
//-----------------------------------------------------------------------------


IFillLookup::IFillLookup(const char* fillTable)
{
  cout << "IFillLookup::IFillLookup" << endl;
  assert(IConsts::Defined);
  cout << "Reading " << IConsts::FillTable << endl;

  _fill = new map<int, int>;

  ifstream fin(fillTable);

  int run, fill;

  while (fin >> run >> fill)
    (*_fill)[run] = fill;


  fin.close();
}
//-----------------------------------------------------------------------------


IFillLookup::~IFillLookup()
{
  delete _fill;
}
//-----------------------------------------------------------------------------


int IFillLookup::GetFill(int run)
{
  return (*_fill)[run];
}
//-----------------------------------------------------------------------------


int IFillLookup::GetFillList(int *fill, const char* runList)
{
  // NOTE: runList should be text file with single column list of runs
  ifstream fin(runList);
  int run;

  int nFill = 0;
  int oldFill=-9, newFill=-9;

  fin >> run;
  oldFill = GetFill(run);

  while (fin >> run)
    {
      newFill = GetFill(run);
      if (oldFill != newFill)
	{
	  fill[nFill++] = oldFill;
	  oldFill = newFill;
	}
    }

  return nFill;
}
//-----------------------------------------------------------------------------


short IFillLookup::GetRuns(int fill, int *runs)
{
  short n=0;
  map<int, int>::iterator itr;

  for (itr = _fill->begin(); itr != _fill->end(); itr++)
    if (itr->second == fill) runs[n++] = itr->first;

  return n;
}
//-----------------------------------------------------------------------------


short IFillLookup::GetNRuns(int fill)
{
  short n=0;
  map<int, int>::iterator itr;

  for (itr = _fill->begin(); itr != _fill->end(); itr++)
    if (itr->second == fill) n++;

  return n;
}
//-----------------------------------------------------------------------------

