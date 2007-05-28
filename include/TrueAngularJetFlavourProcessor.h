#ifndef TrueAngularJetFlavourProcessor_h
#define TrueAngularJetFlavourProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include "inc/algo.h"


using namespace lcio ;
using namespace marlin ;
using EVENT::Track;
using EVENT::MCParticle;

/**  Determine MC Jet Flavour by angular matching of heavy quarks to jets
 * 
 *  The processor looks at all the PDG codes of all MC particles and recognises recognises all particles containing b- and c-quarks.
 *  It then looks at the momentum of the heavy MC particles and at the momentum of the jets.
 *  The association is done by matching matching heavy flavour hadrons to the jet that is angularly closest. 
 *  More than one heavy particle can therefore be associated with the same jet. If this happens the jet flavour is 
 * the flavour of the heaviest particle associated with the jet.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *  Needs the collection of Reconstructed Particles that represent the jets.
 *  
 *
 *  <h4>Output</h4> 
 *  It writes to lcio the calculated flavours of the jets. This is stored in a collection of LCIntVec. By default
 *  the collection is called TrueJetFlavour.
 * 
 * @param MCParticleColName Name of the MCParticle collection.
 * @param JetRPColName Name of the ReconstructedParticle collection that represents jets.
 * @param TrueJetFlavourCollection Name of the output collection where the jet flavours will be stored. 
 * @param MaximumAngle Maximum value allowed between MCParticle and jet momentum expressed in degrees.
 * If the closest jet is at a wider angle than MaximumAngle the MC particle does not get assigned.
 *
 * @author Erik Devetak (e.devetak1@physics.ox.ac.uk), 
 * <br> interface by Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
 * 
 */

class TrueAngularJetFlavourProcessor : public Processor {
  
 public:
  
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new TrueAngularJetFlavourProcessor ; }
  TrueAngularJetFlavourProcessor() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;
  protected:
  std::string _JetRPColName;
  std::string _MCParticleColName;
  std::string _TrueJetFlavourColName;
  double _MaximumAngle;
  int _nRun ;
  int _nEvt ;
} ;

#endif


