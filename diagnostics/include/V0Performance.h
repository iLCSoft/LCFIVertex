#ifndef V0Performance_h
#define V0Performance_h 1

#include "HistMap.h"

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>

#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>

using namespace lcio;
using namespace marlin;
using namespace std;


class V0Performance : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new V0Performance ; }
  V0Performance() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  
  
 private:

  void treeAnalysis( const LCEvent *evt, const string collectionName);
  void hitAnalysis( const LCEvent *evt );
  void trackAnalysis( const LCEvent *evt, const string collectionName);
  void recoAnalysis( const LCEvent *evt, const string collectionName );

  // processor parameters


  // number of hits that at least two charged conv/V0 decay particles
  // need to leave in the detector in order for this conv/V0 to be counted
  int _minHits; 


  // structure to store V0 and conversion candidates
  // - pointer to original MCParticle
  // - pointers to all decay MCParticles
  // - pointers to all tracks corresponding to decay MCParticles
  // - pointers to all V0 candidates corresponding to original MCParticle
  typedef struct {
    const MCParticle* mother;
    std::vector<MCParticle*> daughters;
    std::map<string,vector<Track*> > tracks;
    std::map<string,vector<double> > track_weights;
    std::map<string,vector<ReconstructedParticle*> > recopart;
    int V0Type;
    double radius;
    double z;
  } V0Candidate_type;
  std::vector<V0Candidate_type*> V0Candidates;

  map<MCParticle*,int> numHits;

  void dump_V0Candidate(V0Candidate_type* cand);
  size_t _maxCollNameLength;
  std::vector<std::string> _CollectionsToPrint;
  
  // LCRelation handling
  map<string,LCRelationNavigator*> lcRel;

  // efficiency/purity determination
  enum V0Types { V0Gamma, V0K0, V0Lambda, V0LastType};
  string V0Name[V0LastType];

  int mc_num[V0LastType];
  int mc_tracks[V0LastType];

  map<string,int> num_tracks_total;
  map<string,int> num_composites_total;

  map<string,int> foundNOtracks[V0LastType];
  map<string,int> foundONEtrack[V0LastType];
  map<string,int> foundTWOtracks[V0LastType];
  map<string,int> foundALLtracks[V0LastType];
  map<string,int> foundComposite[V0LastType];

  HistMap *histos;
} ;

#endif



