#include "DSTCollectionProcessor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <EVENT/LCParameters.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/LCIntVec.h>
#include "EVENT/LCFloatVec.h"


#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <vector>
#include <string>
#include <map>
#include <set>

using namespace marlin ;
using namespace lcio;
using namespace std;

using std::vector;
using std::string;
using std::map;
using EVENT::Track;


DSTCollectionProcessor aDSTCollectionProcessor ;

DSTCollectionProcessor::DSTCollectionProcessor() : Processor("DSTCollectionProcessor") {
  
  // modify processor description
  _description = "DSTCollectionProcessor - Takes the flavour tag info, mc truth info, and some NN  input info, and adds it as Particle ID info to the jets, for the DSTs" ;
  

  // register steering parameters: name, description, class-variable, default value



  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" , 
			   "Name of the collection of ReconstructedParticles that is the jet"  ,
			   _JetCollectionName ,
			   std::string("Jets") );


  registerInputCollection( LCIO::LCFLOATVEC,
			    "FlavourTagCollection" , 
			    "Name of the LCFloatVec Collection containing the flavour tags (in same order as jet collection)"  ,
			    _FlavourTagCollectionName,
			      "FlavourTag" ) ;

  registerInputCollection( LCIO::LCFLOATVEC,
  			      "FlavourTagInputsCollection" , 
			      "Name of the LCFloatVec Collection that contains the flavour tag inputs (in same order as jet collection)"  ,
			      _FlavourTagInputsCollectionName,
			      "FlavourTagInputs" ) ;


  
  registerInputCollection( LCIO::LCFLOATVEC,
			   "TrueJetFlavourCollection" , 
			   "Name of the LCIntVec collection containing the true flavour of the jets (same order as jets)"  ,
			   _TrueJetFlavourColName ,
			   std::string("TrueJetFlavour") ) ;

  registerProcessorParameter( "Debug",
                              "Debugging option",
                              _debug,
                              int(0)); 
  

  }

void DSTCollectionProcessor::init() 
{ 
	// usually a good idea to
	printParameters() ;
	
	_nRun = 0 ;
	_nEvt = 0 ;
}

void DSTCollectionProcessor::processRunHeader( LCRunHeader* pRun) { 


  // Marlin doesn't necessarily process the run header, e.g. if you use the "SkipNEvents"
  // parameter in the steering file. The flavour tag variable/tag value names are held in
  // the run header though, so this processor has to have access to it. Set this variable
  // so that "processEvent" can tell if "processRunHeader" has been called for the run
  // it's in.
  _lastRunHeaderProcessed=pRun->getRunNumber();


  //map the parameter names to the indexes, for the flavout tags, NN input vales, and mc truth info

  std::vector<std::string> VarNames;
  (pRun->parameters()).getStringVals(_FlavourTagCollectionName,VarNames);  
  for (size_t i = 0;i < VarNames.size();++i)
    _IndexOfForEachTag[VarNames[i]] = i;
  
    
  //do the same for the true jet values
  std::vector<std::string> trueJetFlavourVarNames;
  (pRun->parameters()).getStringVals(_TrueJetFlavourColName,trueJetFlavourVarNames);
  for (size_t i = 0;i < trueJetFlavourVarNames.size();++i)
    _FlavourIndex[trueJetFlavourVarNames[i]] = i;

  
  //do the same for the FT inputs  
  std::vector<std::string> InputVarNames;
  (pRun->parameters()).getStringVals(_FlavourTagInputsCollectionName,InputVarNames);
  for (size_t i = 0;i < InputVarNames.size();++i)
    _InputsIndex[InputVarNames[i]] = i;

      
  
  _nRun++ ;
  
}


void DSTCollectionProcessor::processEvent(  LCEvent* pEvent) { 

  try{
    LCCollection* JetCollection=pEvent->getCollection( _JetCollectionName );

   try{
     LCCollection* trueJetFlavourCollection = pEvent->getCollection(_TrueJetFlavourColName);   

     try{
       LCCollection* flavourTagInputsCollection = pEvent->getCollection( _FlavourTagInputsCollectionName);
       
       try{
	      LCCollection* TagCollection=pEvent->getCollection( _FlavourTagCollectionName );
  

	 //make a PIDHandler to add the info with
	 PIDHandler jetPID ( JetCollection ) ;

	 //the string vec containing the parameter names for FT info
	 StringVec pNames ;
	 pNames.push_back( "BTag" );
	 pNames.push_back( "CTag" );
	 pNames.push_back( "BCTag" ) ;

	 pNames.push_back("NumVertices");
	 pNames.push_back("JointProbRPhi"); 
	 pNames.push_back("JointProbZ");
	 pNames.push_back("NumTracksInVertices");
	 pNames.push_back("DecayLength");
	 pNames.push_back("DecayLengthSignificance"); 
	 pNames.push_back("RawMomentum");
	 pNames.push_back("PTCorrectedMass");
	 pNames.push_back("SecondaryVertexProbability");


	 //the string vec containing the parameter names for MC info
	 StringVec mcNames ;

	 mcNames.push_back("TruePDGCode");    
	 mcNames.push_back("TruePartonCharge");
	 mcNames.push_back("TrueHadronCharge");
	 mcNames.push_back("TrueJetFlavour");
    
	 //make two new PID algorithms, one for the FT info, one for the MC info
	 int ftagID = jetPID.addAlgorithm( "LCFIFlavourTag",pNames ) ;
	 int mcID = jetPID.addAlgorithm( "MCTruth",mcNames ) ;

	 //fill the float vec that will contain the paramters with dummy values = 9999
	 FloatVec fv(pNames.size());
	 //fill with dummy values
	 for(int i = 0; i < (int)pNames.size(); i++)
	   fv[i] = 9999;

	 FloatVec mcv(mcNames.size());
	 //fill with dummy values
	 for(int i = 0; i < (int)mcNames.size(); i++)
	   mcv[i] = 9999;
      
	 if (_debug>0)
	   cout<<"Number of jets in the event "<<JetCollection->getNumberOfElements()<<endl;
	 
	 //loop over all jets
	 for( int jet=0; jet<JetCollection->getNumberOfElements(); jet++ )
	   {
	     

	     ReconstructedParticle*  Jet  = dynamic_cast<ReconstructedParticle*>( JetCollection->getElementAt(jet) );
	 
	     //get the FT info, add it to the parameter float vec
	     LCFloatVec* flavourTags = dynamic_cast<LCFloatVec*>(TagCollection->getElementAt( jet ));
	     fv[0]= (*flavourTags)[_IndexOfForEachTag["BTag"]];
	     fv[1]= (*flavourTags)[_IndexOfForEachTag["CTag"]];
	     fv[2]=(*flavourTags)[_IndexOfForEachTag["BCTag"]];
	       
	     //get the NN input info, add it to the parameter float vec
	     LCFloatVec* tagInputs = dynamic_cast<LCFloatVec*>(flavourTagInputsCollection->getElementAt( jet ));
	     int  NumVertices = int((*tagInputs)[_InputsIndex["NumVertices"]]);
	     fv[3]= NumVertices;
	     //only add this info if the no of vertices >=2
	     if(NumVertices > 1)
	       {
		 fv[4]=(*tagInputs)[_InputsIndex["JointProbRPhi"]];
		 fv[5]=(*tagInputs)[_InputsIndex["JointProbZ"]];
		 fv[6]=(*tagInputs)[_InputsIndex["NumTracksInVertices"]];
		 fv[7]=(*tagInputs)[_InputsIndex["DecayLength"]];
		 fv[8]=(*tagInputs)[_InputsIndex["DecayLengthSignificance"]];
		 fv[9]=(*tagInputs)[_InputsIndex["RawMomentum"]];
		 fv[10]=(*tagInputs)[_InputsIndex["PTCorrectedMass"]];
		 fv[11]=(*tagInputs)[_InputsIndex["SecondaryVertexProbability"]];
	       }
	     

	       //add the FT particle id to the jet collection
	       jetPID.setParticleID( Jet , 42 , // user type
				     99999 , // PDG
				     0.0, // likelihood
				     ftagID ,
				     fv ) ;
	       
	       //get the mc truth info, add it to the mc parameter float vec
	       LCFloatVec* pJetFlavour = dynamic_cast<LCFloatVec*>(trueJetFlavourCollection->getElementAt( jet ));
	       if(pJetFlavour==0) std::cerr << "The wrong type of true jet flavour collection was found, dynamic cast failed" << std::endl;
	       mcv[0] = (*pJetFlavour)[_FlavourIndex["TruePDGCode"]];
	       mcv[1] =(*pJetFlavour)[_FlavourIndex["TruePartonCharge"]];
	       mcv[2]= (*pJetFlavour)[_FlavourIndex["TrueHadronCharge"]];
	       mcv[3]= (*pJetFlavour)[_FlavourIndex["TrueJetFlavour"]];
	       
	       //add the MC particle id to the jet collection
	       jetPID.setParticleID( Jet , 43 , // user type
				     99998 , // PDG
				     0.0, // likelihood
				     mcID ,
				     mcv ) ;
	     }


       }
       catch(DataNotAvailableException &e){
	 if (_debug == 1) 
	   std::cout << "Collection " <<_FlavourTagCollectionName << " is unavailable in event " << _nEvt << std::endl;
       }
     }
     catch(DataNotAvailableException &e){
       if (_debug == 1) 
	 std::cout << "Collection " <<_FlavourTagInputsCollectionName << " is unavailable in event " << _nEvt << std::endl;
     }
   }
   catch(DataNotAvailableException &e){
    if (_debug == 1) 
      std::cout << "Collection " << _TrueJetFlavourColName << " is unavailable in event " << _nEvt << std::endl;
   }
  }
  catch(DataNotAvailableException &e){
    if (_debug == 1) 
      std::cout << "Collection " <<_JetCollectionName << " is unavailable in event " << _nEvt << std::endl;
  }




  _nEvt ++ ;

}




void DSTCollectionProcessor::check( LCEvent* ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DSTCollectionProcessor::end(){ 
	
	std::cout << "DSTCollectionProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}

