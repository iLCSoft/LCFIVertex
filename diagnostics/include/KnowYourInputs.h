#ifndef KnowYourInputs_h
#define KnowYourInputs_h 1


#include "HistMap.h"

//Marlin and LCIO includes 
#include "marlin/Processor.h" 
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"

#include <string>
#include <map>
#include <iostream>
#include <fstream>



using namespace lcio;
using namespace marlin;
using namespace std;

/**  Diagnostics processor for (Sim)TrackerHits, (Sim)CalorimeterHits,
 *   Tracks and ReconstructedParticles.
 *
 * @author Kristian Harder, RAL
 * @version $Id: KnowYourInputs.h,v 1.2 2008-04-29 07:30:37 harderk Exp $ 
 */

class KnowYourInputs : public Processor {
  
 public:
 
  virtual Processor*  newProcessor() { return new KnowYourInputs ; }
   
  KnowYourInputs() ;
  
  virtual void init() ;

  virtual void processEvent( LCEvent * evt ) ; 

  virtual void end() ;
  
  
 private:

  void trackerPlots( const LCEvent *evt,
		     const string collName, const string subdet );
  void calPlots( const LCEvent *evt,
		 const string collName, const string subdet );
  void trackerHitPlots( const LCEvent *evt,
			const string collName, const string subdet );
  void simTrackerHitPlots( const LCEvent *evt,
			   const string collName, const string subdet );
 
  bool _trackEnable=false;
  bool _recParEnable=false;
  bool _simTrkHitEnable=false;
  bool _trkHitEnable=false;
  bool _calHitEnable=false;
  bool _simCalHitEnable=false;
  
  double _bField=0.0;

  vector<string> _knownCollections{};
  vector<string> _knownCollectionTypes{};

  map<string,HistMap*> _histograms{};
};

#endif
