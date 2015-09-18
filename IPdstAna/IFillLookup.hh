#ifndef IFILLLOOKUP_HH
#define IFILLLOOKUP_HH

#include <map>

class IFillLookup
{
public:
  IFillLookup();
  IFillLookup(const char* fillTable);
  ~IFillLookup();

  int   GetFill(int run);
  short GetRuns(int fill, int *runs); // return value is # of runs in fill
  short GetNRuns(int fill);           // returns number of runs in the fill


  // this function takes list of runs and gives corresponding list for fill
  // numbers in fillList array. Returs number of fills;
  int GetFillList(int *fillList, const char* runList);

private:
  std::map <int, int> *_fill;
};
#endif
