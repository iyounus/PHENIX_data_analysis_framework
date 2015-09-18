#ifndef IHADRONANA_HH
#define IHADRONANA_HH

#include <vector>

class IPdst;
class IEvent;
class TH1D;
class TFile;


class IHadronAna
{
public:
  IHadronAna();
  ~IHadronAna();

private:
  TFile *_outFile;

  IPdst *_pdst;

  std::vector <IEvent *> *_vEvent;
  std::vector <TH1D *> *_vHist;

  TH1D **_hPt;        // array of histograms for different pt bins
  TH1D **_hEP;        // E/p ratio

  void CreateHistos();
  void ReadPdst(int run);
  void FillHistos();
  void ResetEventVector();
  void ResetHistos();
};
#endif
