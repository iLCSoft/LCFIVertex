#ifndef HistMap_h
#define HistMap_h 1

#include <string>
#include <map>
#include "marlin/Processor.h"

#ifdef MARLIN_USE_AIDA

#include <AIDA/IHistogramFactory.h>
#include <AIDA/ITree.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>

#endif

using namespace std;


class HistMap {
  
 public:
 
  HistMap(marlin::Processor* processor, string const& folderName="");
  ~HistMap();
  HistMap(const HistMap&) = delete;
  HistMap& operator=(const HistMap&) = delete;

  void fill(string histname,double val, double weight=1,
	    string title="", int nbins=0, double min=0, double max=0,
	    string axistitle="");
  void fill(string histname,double xval, double yval, double weight=1,
	    string title="", int nbinsx=0, double xmin=0, double xmax=0,
	    int nbinsy=0, double ymin=0, double ymax=0,
	    string xaxistitle="", string yaxistitle=""); 

 private:

  string _dirName;


#ifdef MARLIN_USE_AIDA

  AIDA::IHistogramFactory* _histFact=nullptr;
  AIDA::ITree* _aidaTree=nullptr;

  map<string,AIDA::IHistogram1D*> _histmap1D{};
  map<string,AIDA::IHistogram2D*> _histmap2D{};

#endif

  map<string,double> _histmap1D_xmin{};
  map<string,double> _histmap2D_xmin{};
  map<string,double> _histmap2D_ymin{};
  map<string,double> _histmap1D_xmax{};
  map<string,double> _histmap2D_xmax{};
  map<string,double> _histmap2D_ymax{};
  map<string,char> _histmap1D_xnan{};
  map<string,char> _histmap2D_xnan{};
  map<string,char> _histmap2D_ynan{};

} ;

#endif
