#ifndef DSTCollectionProcessor_h
#define DSTCollectionProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"



using namespace lcio ;
using namespace marlin ;



class DSTCollectionProcessor : public Processor {
  
public:
  
  //The usual Marlin processor methods
  virtual Processor*  newProcessor() { return new DSTCollectionProcessor ; }
  DSTCollectionProcessor() ;
  virtual void init() ;
  virtual void processRunHeader( LCRunHeader* run ) ;
  virtual void processEvent( LCEvent * evt ) ; 
  virtual void check( LCEvent * evt ) ; 
  virtual void end() ;

protected:

  std::string _FlavourTagCollectionName{};
  std::string _FlavourTagInputsCollectionName{};
  std::string _TrueJetFlavourColName{};
  
  std::string _JetCollectionName{};

  int _lastRunHeaderProcessed=-1;

  std::map<std::string,unsigned int> _IndexOfForEachTag{};
  std::map<std::string,unsigned int>  _FlavourIndex{};
  std::map<std::string,unsigned int>_InputsIndex{};

  int _nRun=-1;
  int _nEvt=-1;

  int _debug=0;
} ;

#endif



