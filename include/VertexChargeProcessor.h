#ifndef VertexChargeProcessor_h
#define VertexChargeProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>

#include "inc/algo.h"
#include "inc/decaychain.h"
#include "inc/jet.h"
#include "util/inc/projection.h"


using namespace lcio ;
using namespace marlin ;
using vertex_lcfi::DecayChain;
using vertex_lcfi::Jet;
using vertex_lcfi::util::Projection;

/** Calculates the Vertex Charge. 
 *
 *  @author Erik Devetak(erik.devetak1@physics.ox.ac.uk)
*/
class VertexChargeProcessor : public Processor {
  
 public:
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new VertexChargeProcessor ; }
  VertexChargeProcessor() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  protected:
  std::string _JetRPColName;
  std::string _DecayChainRPColName;
  std::string _RelationColName;
  std::string _IPVertexCollectionName;
  std::string _VertexChargeCollectionName;
  
  std::vector<std::string> _JetVariableNames;
  
  vertex_lcfi::Algo<DecayChain*,double>* _VertexCharge;
  vertex_lcfi::Algo<DecayChain*,DecayChain* >* _Attach;

  bool _ChargeAddAllTracksFromSecondary;
  double _ChargeLoDCutmin;
  double _ChargeLoDCutmax;
  double _ChargeCloseapproachCut;

  int _nRun ;
  int _nEvt ;
} ;

#endif



