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
    bool is_conversion;
    double radius;
    double z;
  } V0Candidate_type;
  std::vector<V0Candidate_type*> V0Candidates;

  map<MCParticle*,int> numHits;
  int minHits;

  void dump_V0Candidate(V0Candidate_type* cand);

  // LCRelation handling
  map<string,LCRelationNavigator*> lcRel;

  // efficiency/purity determination
  int total_num_v0;
  int ntrue_tracks;
  map<string,int> ncand_tracks;
  map<string,int> ngood_tracks;
  map<string,int> ngood_tracks_in_composites;
  map<string,int> ntracks_in_composites;

  HistMap *histos;
} ;

#endif



