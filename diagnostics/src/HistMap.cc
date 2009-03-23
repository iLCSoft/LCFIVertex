#include "HistMap.h"
#include "streamlog/streamlog.h"

#include <cstdlib>

#ifdef MARLIN_USE_AIDA
#include "marlin/AIDAProcessor.h"
#endif

HistMap::HistMap(marlin::Processor* processor, string folderName) {
  _dirName="/"+processor->name()+"/"+folderName;

#ifdef MARLIN_USE_AIDA
  _histFact=marlin::AIDAProcessor::histogramFactory(processor);
  _aidaTree=marlin::AIDAProcessor::tree(processor);
  _aidaTree->mkdir("/"+processor->name());
  _aidaTree->mkdir(_dirName);
  _histmap1D.clear();
  _histmap2D.clear();
#endif

  streamlog_out(MESSAGE) << "creating HistMap in "+_dirName << endl;
}

HistMap::~HistMap() {
}


void HistMap::fill(string histname,double val, double weight,
		   string title, int nbins, double min, double max,
		   string axistitle) {


#ifdef MARLIN_USE_AIDA
  if (!_histmap1D[histname]) {
    if (nbins==0) {
      streamlog_out(ERROR) << "must define histogram "+histname+
	" properties on first call of fillhisto(...)";
      exit(1);
    }
    _aidaTree->cd(_dirName);
    _histmap1D[histname]
      = _histFact->createHistogram1D(histname,title,nbins,min,max);
    if (!_histmap1D[histname]) {
      streamlog_out(ERROR) << "failed to create histogram "
	+_dirName+"/"+histname+" (out of memory?)";
      exit(1);
    }
    _histmap1D_xnan[histname]=0;
  }

  // fill histogram and record any unusual entries (under-/overflow, NaN)
  if ((val>=0) || (val<0)) {
    if (_histmap1D[histname]->allEntries()==0) {
      _histmap1D_xmin[histname] = val;
      _histmap1D_xmax[histname] = val;
    } else {
      _histmap1D_xmin[histname] = std::min(val,_histmap1D_xmin[histname]);
      _histmap1D_xmax[histname] = std::max(val,_histmap1D_xmax[histname]);
    }
    _histmap1D[histname]->fill(val,weight);
  } else {
    ++_histmap1D_xnan[histname];
  }
#endif

}


void HistMap::fill(string histname,double xval, double yval,
		   double weight,
		   string title, int nbinsx, double xmin, double xmax,
		   int nbinsy, double ymin, double ymax,
		   string xaxistitle, string yaxistitle) {


#ifdef MARLIN_USE_AIDA
  if (!_histmap2D[histname]) {
    if (nbinsx==0) {
      streamlog_out(ERROR) << "must define histogram "+histname
	+" properties on first call of fillhisto(...)";
      exit(1);
    }

    _aidaTree->cd(_dirName);
    _histmap2D[histname] = _histFact->createHistogram2D( histname, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax ); 
    if (!_histmap2D[histname]) {
      streamlog_out(ERROR) << "failed to create histogram "
	+_dirName+"/"+histname+" (out of memory?)";
      exit(1);
    }
    _histmap1D_xnan[histname]=0;
  }
  _histmap2D[histname]->fill(xval,yval,weight);
#endif

}
