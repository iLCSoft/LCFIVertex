#ifndef ConversionTagger_h
#define ConversionTagger_h 1

#include "HistMap.h"

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>


using namespace lcio;
using namespace marlin;
using namespace std;


class ConversionTagger : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new ConversionTagger ; }
  ConversionTagger() ;
  ConversionTagger(const ConversionTagger&) = delete;
  ConversionTagger& operator=(const ConversionTagger&) = delete;

  virtual void init() ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void end() ;
  
  
 private:

  void tagger( LCEvent *evt, const string collectionName);
  ReconstructedParticle* CreateRecoPart(Track* trk);
  double diParticleMass(const double* mom1, const double* mom2,
			const double mass1, const double mass2,
			const float* vertmom);

  MCParticle* isFromV0(ReconstructedParticle* rp, vector<LCRelationNavigator*>relCols);

  HistMap* histos=nullptr;

  std::vector<std::string> _InputCollections{};
  double _massRangePhoton=0.0;
  double _massRangeKaon=0.0;
  double _massRangeLambda=0.0;
  double _distCut=0.0;
  bool   _cheatMode=false;
  bool   _cheatEvenMore=false;
  std::vector<int> _PdgToTag{};
  std::map<int,bool> _TagPDG{};
  double _minDistFromIP=0.0;

  double _BField=0.0;
} ;

#endif
