// First of all, make sure that AIDA was enabled by setting the environment
// variable MARLIN_USE_AIDA when Marlin was compiled. The makefile will then
// have done the setup and defined this macro.

#ifndef MARLIN_USE_AIDA

#warning "--------------------------------------------------------------------------------"
#warning "- LCFIAIDAPlotProcessor requires MARLIN_USE_AIDA to be defined. Did you enable -"
#warning "- AIDA when compiling Marlin? LCFIAIDAPlotProcessor will not be compiled.      -"
#warning "--------------------------------------------------------------------------------"

// Can't do anything else.
#else

#include "LCFIAIDAPlotProcessor.h" 
 
// standard library includes 
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <set>
#include <math.h>

// LCIO includes... 
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/Vertex.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Track.h"


// Marlin includes
#include <marlin/Exceptions.h>

// AIDA includes...
#include <marlin/AIDAProcessor.h>
#include <AIDA/IMeasurement.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IAxis.h>
#include <AIDA/ITree.h>
#include <AIDA/ITupleFactory.h>
#include <AIDA/ICloud2D.h>
#include <UTIL/LCRelationNavigator.h>

//change USING_RAIDA to USING_JAIDA if you are using JAIDA/AIDAJNI - you will obtain more functionality!
#define USING_RAIDA

#ifdef USING_RAIDA 
#warning "USING_RAIDA defined"
#else
#define USING_JAIDA
#warning "USING_JAIDA defined"
#endif


#ifdef USING_JAIDA//Data point sets aren't implemented in RAIDA - which is a shame as they have functionality not given by histograms
//such as the facility to set the error
#include <AIDA/IDataPointSet.h>
#include <AIDA/IDataPointSetFactory.h>
#include <AIDA/IDataPoint.h>
#endif

#include "TypesafeCollection.h"

// There needs to be at least one instantiation for the base constructor to register the processor with 
// the Marlin processor manager. This is it. 
LCFIAIDAPlotProcessor aLCFIAIDAPlotProcessor; 

LCFIAIDAPlotProcessor::LCFIAIDAPlotProcessor() : marlin::Processor( "LCFIAIDAPlotProcessor" ) 
{ 

  _description="Creates an AIDA plot of the LCFIVertex tagging efficiency-purity values and various other things.  Make sure that MarlinAIDAProcessor is run before this.";
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "JetCollectionName" , 
			   "Name of the collection of ReconstructedParticles that is the jet"  ,
			   _JetCollectionName ,
			   std::string("FTSelectedJets") );
  
  std::vector<std::string> FlavourTagCollectionNamesDefault;
  FlavourTagCollectionNamesDefault.push_back("FlavourTag");
  registerProcessorParameter("FlavourTagCollections" , 
			     "Names of the LCFloatVec Collections that contain the flavour tags (one purity efficiency plot per tag) (in same order as jet collection)"  ,
			     _FlavourTagCollectionNames,
			     FlavourTagCollectionNamesDefault) ;
  
  //registerInputCollection( LCIO::LCINTVEC,
  //		   "TrueJetFlavourCollection" , 
  //		   "Name of the LCIntVec collection containing the true flavour of the jets (same order as jets)"  ,
  //		   _TrueJetFlavourColName ,
  //		   std::string("TrueJetFlavour") ) ;

  registerOptionalParameter( "TrueJetFlavourCollection" , 
			     "Name of the LCIntVec collection containing the true flavour of the jets (same order as jets)"  ,
			     _TrueJetFlavourColName ,
			     std::string("TrueJetFlavour") ) ;
  
  //registerInputCollection( LCIO::LCFLOATVEC,
  //		   "TrueJetHadronChargeCollection",
  //		   "Name of the LCFloatVec collection containing the true hadron charge of the jets (same order as jets)"  ,
  //		   _TrueJetHadronChargeColName ,
  //		   std::string("TrueJetHadronCharge") ) ;  
  
  registerOptionalParameter("TrueJetHadronChargeCollection",
			    "Name of the LCFloatVec collection containing the true hadron charge of the jets (same order as jets)"  ,
			    _TrueJetHadronChargeColName ,
			    std::string("TrueJetHadronCharge") ) ; 

  //registerInputCollection( LCIO::LCINTVEC,
  //			   "TrueJetPDGCodeCollection" , 
  //			   "Name of the LCIntVec collection containing the true PDG code of the jets (same order as jets)"  ,
  //			   _TrueJetPDGCodeColName,
  //			   std::string("TrueJetPDGCode") ) ;

  registerOptionalParameter( "TrueJetPDGCodeCollection" , 
			     "Name of the LCIntVec collection containing the true PDG code of the jets (same order as jets)"  ,
			     _TrueJetPDGCodeColName,
			     std::string("TrueJetPDGCode") ) ;

  //registerInputCollection( LCIO::LCFLOATVEC,
  //			   "TrueJetPartonChargeCollection",
  //			   "Name of the LCFloatVec collection containing the true parton charge of the jets (same order as jets)"  ,
  //			   _TrueJetPartonChargeColName ,
  //			   std::string("TrueJetPartonCharge") ) ;    
  
  registerOptionalParameter("TrueJetPartonChargeCollection",
			    "Name of the LCFloatVec collection containing the true parton charge of the jets (same order as jets)"  ,
			    _TrueJetPartonChargeColName ,
			    std::string("TrueJetPartonCharge") ) ;    
  
  registerInputCollection( lcio::LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the collection that holds all MC particles. "  ,
			   _MCParticleColName ,
			   std::string("MCParticle") ) ;

  //registerInputCollection( lcio::LCIO::VERTEX,
  //			   "VertexCollection",
  //			   "Name of the collection that holds the Vertices",
  //			   _VertexColName,
  //			   std::string("ZVRESVertices") ) ;

  registerOptionalParameter("VertexCollection",
			    "Name of the collection that holds the Vertices",
			    _VertexColName,
			    std::string("ZVRESVertices") ) ;

  //registerInputCollection(LCIO::LCFLOATVEC,
  //			  "CVertexChargeCollection",
  //			  "Name of collection containing the vertex charge of the jets, assuming they are C-jets",
  //			  _CVertexChargeCollection,
  //			  std::string("CCharge") );
  
  registerOptionalParameter("CVertexChargeCollection",
			    "Name of collection containing the vertex charge of the jets, assuming they are C-jets",
			    _CVertexChargeCollection,
			    std::string("CCharge") );
  
  //registerInputCollection( LCIO::LCFLOATVEC,
  //			   "BVertexChargeCollection",
  //			   "Name of collection containing the vertex charge of the jets, assuming they are B-jets",
  //			   _BVertexChargeCollection,
  //			   std::string("BCharge") ) ;
  
  registerOptionalParameter( "BVertexChargeCollection",
			      "Name of collection containing the vertex charge of the jets, assuming they are B-jets",
			      _BVertexChargeCollection,
			      std::string("BCharge") ) ;
  
  //registerInputCollection( LCIO::LCCOLLECTION,
  //			   "TrueTracksToMCPCollection",
  //			   "Name of collection linking the tracks and the Monte Carlo Particles",
  //			   _TrueTracksToMCPCollection,
  //			   std::string("LDCTracksMCP") ) ;

  registerOptionalParameter("TrueTracksToMCPCollection",
			    "Name of collection linking the tracks and the Monte Carlo Particles",
			    _TrueTracksToMCPCollection,
			    std::string("LDCTracksMCP") ) ;
  
  ///registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  //			   "ZVRESDecayChainCollection" , 
  //			   "Name of the ZVRES DecayChain collection"  ,
  //			   _ZVRESDecayChainCollection ,
  //			   std::string("ZVRESDecayChains") );
  
  registerOptionalParameter( "ZVRESDecayChainCollection" , 
			     "Name of the ZVRES DecayChain collection"  ,
			     _ZVRESDecayChainCollection ,
			     std::string("ZVRESDecayChains") );
  
  //registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  //		   "ZVRESSelectedJetsCollection" , 
  //		   "Name of the ZVRES Selected Jets collection"  ,
  //		   _ZVRESSelectedJetsCollection ,
  //		   std::string("ZVRESSelectedJets") );
  
  registerOptionalParameter( "ZVRESSelectedJetsCollection" , 
			     "Name of the ZVRES Selected Jets collection"  ,
			     _ZVRESSelectedJetsCollection ,
			     std::string("ZVRESSelectedJets") );
  
  //registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
  //		   "ZVRESDecayChainTrackCollection" , 
  //		   "Name of the ZVRES Decay Chain Tracks Collection"  ,
  //		   _ZVRESDecayChainRPTracksCollection ,
  //		   std::string("ZVRESDecayChainRPTracks") );
  
  registerOptionalParameter("ZVRESDecayChainTrackCollection" , 
			    "Name of the ZVRES Decay Chain Tracks Collection"  ,
			    _ZVRESDecayChainRPTracksCollection ,
			    std::string("ZVRESDecayChainRPTracks") );

  FlavourTagCollectionNamesDefault.clear();
  FlavourTagCollectionNamesDefault.push_back("FlavourTagInputs");
  registerProcessorParameter("TagInputsCollections" , 
			     "Names of the LCFloatVec Collections that contain the flavour tag inputs (in same order as jet collection)"  ,
			     _FlavourTagInputsCollectionNames,
			     FlavourTagCollectionNamesDefault) ;

  registerOptionalParameter( "CosThetaJetMax",
			     "Cut determining the maximum cos(theta) of the jet.  Default: |cos(theta)|<0.9"  ,
			     _CosThetaJetMax,
			     double(0.9)) ;
   
  registerOptionalParameter( "CosThetaJetMin",
			     "Cut determining the minimum cos(theta) of the jet.  Default: no lower cut."  ,
			     _CosThetaJetMin,
			     double(0.0)) ;

  registerOptionalParameter("PJetMax",
			    "Cut determining the maximum momentum of the jet.  Default: 10000 GeV/c"  ,
			    _PJetMax,
			    double(10000.)) ;
   
  registerOptionalParameter( "PJetMin",
			     "Cut determining the minimum momentum of the jet.  Default: no lower cut."  ,
			     _PJetMin,
			     double(0.0)) ;

  registerOptionalParameter( "PrintTrackVertexOutput",
			     "Set true if you want a print-out of the track-vertex association purity",
			     _PrintTrackVertexOutput,
			     bool(false));

  registerOptionalParameter( "MakePurityEfficiencyPlots",
			     "Set true if you want to make the purity-efficiency plots, and leakage rates plots for the various flavour tags",
			     _MakePurityEfficiencyPlots,
			     bool(true));
  
  registerOptionalParameter( "PrintPurityEfficiencyValues",
			     "Set true if you want a print-out of the purity-efficiency for the various flavour tags",
			     _PrintPurityEfficiencyValues,
			     bool(true));
  
  registerOptionalParameter( "MakeAdditionalPlots",
			     "Set true if you want to make all other plots (i.e. the non purity-efficiency plots) and the other functionality provided by LCFIAIDAPlotProcessor",
			     _MakeAdditionalPlots,
			     bool(false));

  registerOptionalParameter( "PurityEfficiencyOutputFile" , 
			     "Output filename for the Purity-Efficiency values.  Only used if PrintPurityEfficiencyValues parameter is true.  If left blank, output will be directed to standard out.",
			     _PurityEfficiencyOutputFile,
			     std::string("PurityEfficiencyOutput.txt") ) ;

  registerOptionalParameter( "TrackVertexOutputFile" , 
			     "Output filename for the table of the Track-Vertex association.  Only used if PrintTrackVertexOutput parameter is true.  If left blank, output will be directed to standard out.",
			     _TrackVertexOutputFile,
			     std::string("TrackVertexOutput.txt") ) ;

  registerOptionalParameter( "MakeTuple",
			     "Set true to make a tuple of the flavour tag input variables.  Default is false (only works with jaida).",
			     _MakeTuple,
			     bool(false));

  registerOptionalParameter( "CTagNNCut",
			     "Cut determining the Neural Net cut used to select C-Jets",
			     _CTagNNCut,
			     double(0.7));

  registerOptionalParameter( "BTagNNCut",
			     "Cut determining the Neural Net cut used to select B-Jets",
			     _BTagNNCut,
			     double(0.7));

  registerOptionalParameter( "UseFlavourTagCollectionForVertexCharge",
			     "Integer parameter determing which FlavourTag Collection to use the determine C-Jets and B-Jets in Vertex Charge Plots",
			     _iVertexChargeTagCollection,
			     int(0));

} 

LCFIAIDAPlotProcessor::~LCFIAIDAPlotProcessor() 
{ 
} 

void LCFIAIDAPlotProcessor::init()
{
  
  if ((_iVertexChargeTagCollection >=  int(_FlavourTagCollectionNames.size()) || _iVertexChargeTagCollection < 0) && _FlavourTagCollectionNames.size()!=0) {
    std::cerr << " In " << __FILE__ << "(" << __LINE__ << "): Invalid parameter for UseFlavourTagCollectionForVertexCharge.  Setting to 0." << std::endl;
    _myVertexChargeTagCollection = 0;
  } else if (_FlavourTagCollectionNames.size()==0) {
    _myVertexChargeTagCollection = 0;
  } else {
    _myVertexChargeTagCollection = uint(_iVertexChargeTagCollection);
  }

  _ZoomedVarNames.push_back("D0Significance1"); 
  _ZoomedVarNames.push_back("D0Significance2");
  _ZoomedVarNames.push_back("Z0Significance1");
  _ZoomedVarNames.push_back("Z0Significance2");

  _VertexCatNames.resize(N_VERTEX_CATEGORIES+1);
  _VertexCatNames[0]="AnyNumberOfVertices";
  _VertexCatNames[1]="OneVertex";
  _VertexCatNames[2]="TwoVertices";
  _VertexCatNames[3]="ThreeOrMoreVertices";

  
  _NumVertexCatDir.resize(N_VERTEX_CATEGORIES+1);
  _NumVertexCatDir[1]="OneVertex";
  _NumVertexCatDir[2]="TwoVertices";
  _NumVertexCatDir[3]="ThreeOrMoreVertices";
  _NumVertexCatDir[0]="AnyNumberOfVertices";


  _numberOfPoints=100;

  
 
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  
  if(  pHistogramFactory!=0 )
    {
      if (!(pTree->cd( "/" + name() + "/"))) {
	pTree->mkdir( "/" + name() + "/" );
	pTree->cd( "/" + name() + "/");
      }
      
       
      if (_MakePurityEfficiencyPlots)  CreateTagPlots();
      if (_MakeAdditionalPlots) CreateAdditionalPlots();
      if (_MakeTuple) CreateFlavourTagTuple();
      
    } else {

    std::cerr  << "### " << __FILE__ << "(" << __LINE__ << "): Unable to get the histogram factory! No histograms will be made."<< std::endl;
  }
  
  _lastRunHeaderProcessed=-1;
  _suppressOutputForRun=-1;
  

  InternalVectorInitialisation();

}

void LCFIAIDAPlotProcessor::processRunHeader( LCRunHeader* pRun ) 
{

	// Marlin doesn't necessarily process the run header, e.g. if you use the "SkipNEvents"
	// parameter in the steering file. The flavour tag variable/tag value names are held in
	// the run header though, so this processor has to have access to it. Set this variable
	// so that "processEvent" can tell if "processRunHeader" has been called for the run
	// it's in.
	_lastRunHeaderProcessed=pRun->getRunNumber();

	InitialiseFlavourTagInputs(pRun);
	

	//
	// Perform a check to see if the variable names we need are here
	//
	for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag) // Loop over the different tag collection names given in the steering
	{
	  std::vector<std::string> VarNames;
	  (pRun->parameters()).getStringVals(_FlavourTagCollectionNames[iTag],VarNames);
	  
	  
	  //Fill a map so that we can get the array index from just the string
	  std::set<std::string> AvailableNames;
	  std::map<std::string,unsigned int> IndexOf;
	  
	  for (size_t i = 0;i < VarNames.size();++i)
	    {
	      AvailableNames.insert(VarNames[i]);
	      IndexOf[VarNames[i]] = i;
	    }
	  
	  //Add the index to the list
	  _IndexOfForEachTag.push_back(IndexOf);
	  
	  //Check the required information is in the LCFloatVec
	  std::set<std::string> RequiredNames;
	  RequiredNames.insert("BTag");
	  RequiredNames.insert("CTag");
	  RequiredNames.insert("BCTag");
	  
	  if (!std::includes(AvailableNames.begin(),AvailableNames.end(),RequiredNames.begin(),RequiredNames.end()))
	    {
	      std::cerr << __FILE__ << "(" << __LINE__ << "): The collection \"" << _FlavourTagCollectionNames[iTag]
			<< "\" (if it exists) does not contain the tag values required by " << type() << "." << std::endl;
	      std::cerr <<   __FILE__ << "(" << __LINE__ << "): The collection \"" << _FlavourTagCollectionNames[iTag]
			<< "\" (if it exists) does not contain the tag values required by " << type() << "." << std::endl;
	    }
	}

	if (_MakeTuple) CreateFlavourTagInputPlots(pRun);
	
}

void LCFIAIDAPlotProcessor::CreateFlavourTagTuple()
{
 
  //AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  

#ifdef USING_JAIDA
  //something dosen't work for me with the tuples in RAIDA
  AIDA::ITupleFactory* pTupleFactory=marlin::AIDAProcessor::tupleFactory( this );


  if (_MakeTuple) {
    pTree->cd( "/" + name());
    
    //make the ntuple
    //this breaks the paradigm of reading these in from the flavour tag collections themselves
    std::string columnNames="int TrueJetFlavour=-1,  int NumberOfVertices=-1, int NumberOfTracksInVertices=-1, float D0Significance1 = -999., float D0Significance2 = -999., float DecayLength = -999., float DecayLength_SeedToIP= -999., float DecayLengthSignificance= -999., float JointProbRPhi= -999., float JointProbZ= -999., float Momentum1= -999.,float Momentum2= -999., float PTCorrectedMass= -999., float RawMomentum= -999., float SecondaryVertexProbability= -999., float Z0Significance1= -999., float Z0Significance2= -999., int BQVtx=-10, int CQVtx=-10";

   	
    if (!pTree->cd( "/" + name() + "/" +_FlavourTagCollectionNames[_myVertexChargeTagCollection] + "/")) {
      pTree->cd( "/" + name() + "/"); 
      if (!pTree->cd( _FlavourTagCollectionNames[_myVertexChargeTagCollection])) {
	pTree->mkdir( _FlavourTagCollectionNames[_myVertexChargeTagCollection] + "/");
	pTree->cd( _FlavourTagCollectionNames[_myVertexChargeTagCollection]+ "/");
      }
    }
    
    if (!pTree->cd( "TupleDir/")) {
      pTree->mkdir( "TupleDir/");
      pTree->cd( "TupleDir/");
    }
    
    _pMyTuple=pTupleFactory->create( "FlavourTagInputsTuple","FlavourTagInputsTuple", columnNames);
 }
#endif	
 
 
}

void LCFIAIDAPlotProcessor::InitialiseFlavourTagInputs(LCRunHeader* pRun )
{
  for (unsigned int iInputCollection=0; iInputCollection < _FlavourTagInputsCollectionNames.size(); ++iInputCollection)
    {
      
      std::vector<std::string> VarNames;
      (pRun->parameters()).getStringVals(_FlavourTagInputsCollectionNames[iInputCollection],VarNames);
      //Fill the map relating names and indexes
      std::map<std::string,unsigned int> IndexOf;
      for (size_t i = 0;i < VarNames.size();++i)
	{ 	  
	  IndexOf[VarNames[i]] = i;
	}
	    
      _InputsIndex.push_back(IndexOf);
    }
}



void LCFIAIDAPlotProcessor::CreateFlavourTagInputPlots(LCRunHeader* pRun )
{
	
  _inputsHistogramsBJets.resize( _FlavourTagInputsCollectionNames.size() );
  _inputsHistogramsCJets.resize( _FlavourTagInputsCollectionNames.size() );
  _inputsHistogramsUDSJets.resize( _FlavourTagInputsCollectionNames.size() );

  _zoomedInputsHistogramsBJets.resize( _FlavourTagInputsCollectionNames.size() );
  _zoomedInputsHistogramsCJets.resize( _FlavourTagInputsCollectionNames.size() );
  _zoomedInputsHistogramsUDSJets.resize( _FlavourTagInputsCollectionNames.size() );
  
  	AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
	AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
	
	for (unsigned int iInputCollection=0; iInputCollection < _FlavourTagInputsCollectionNames.size(); ++iInputCollection)
	  {
	    	  
	    std::vector<std::string> VarNames;
	    (pRun->parameters()).getStringVals(_FlavourTagInputsCollectionNames[iInputCollection],VarNames);
	    
	    //Fill the map relating names and indexes
	    for (size_t i = 0;i < VarNames.size();++i)
	      {
		
		// If there is no histogram for this name then create one
		if( _inputsHistogramsBJets[iInputCollection][VarNames[i]]==0 )
		  {
		    pTree->cd( "/" + name() + "/");
		    if( !pTree->cd(_FlavourTagInputsCollectionNames[iInputCollection] ) )
		      {
			pTree->mkdir( _FlavourTagInputsCollectionNames[iInputCollection] ) ; 
			pTree->cd(_FlavourTagInputsCollectionNames[iInputCollection] ) ;
		      }
		    
		    int numberOfPoints=_numberOfPoints/4;
		    double lowerBin=-1;
		    double higerBin=1;
		   
		    //binning variables: if the name is not listed here it will use the default above
		    if( VarNames[i]=="BQVtx" || VarNames[i]=="CQVtx" )
		      {
			numberOfPoints=9;
			lowerBin=-4.5;
			higerBin=4.5;
		      }
		    else if( VarNames[i]=="NumVertices" )
		      {
			numberOfPoints=5;
			lowerBin=0.5;
			higerBin=5.5;
		      }
		    else if( VarNames[i]=="NumTracksInVertices" )
		      {
			numberOfPoints=16;
			lowerBin=-0.5;
			higerBin=15.5;
		      }
		    else if ( VarNames[i]=="D0Significance1" || VarNames[i]=="D0Significance2" )
		      {
			numberOfPoints=120;
			lowerBin=-20.;
			higerBin=100.;
		      }
		    else if ( VarNames[i]=="Z0Significance1" || VarNames[i]=="Z0Significance2")
		      {
			numberOfPoints=100;
			lowerBin=-50.;
			higerBin=50.;
		      }
		    else if (VarNames[i]=="DecayLengthSignificance") 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=100.;
		      }
		    else if (VarNames[i]=="DecayLength" || VarNames[i]=="DecayLength(SeedToIP)" ) 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=10.;
		      }
		    else if (VarNames[i]=="JointProbRPhi" || VarNames[i]=="JointProbZ"|| VarNames[i]=="SecondaryVertexProbability") 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=1.0;
		      }
		    else if (VarNames[i]=="Momentum1" || VarNames[i]=="Momentum2" ||  VarNames[i]=="RawMomentum" ) 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=50.;
		      }
		    else if (VarNames[i]=="PTCorrectedMass" ) 
		      {
			numberOfPoints=100;
			lowerBin=0.;
			higerBin=10.;
		      }
		    
		    pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection]);
		    if( !pTree->cd( "bJets" ) )
		    {
		      pTree->mkdir( "bJets" ) ; 
		      pTree->cd(    "bJets" ) ;
		    }
		    
		    _inputsHistogramsBJets[iInputCollection][VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		    
		    pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection]);
		    if( !pTree->cd( "cJets" ) )
		      {
			pTree->mkdir( "cJets" ) ; 
			pTree->cd(    "cJets" ) ;
		      }

		    _inputsHistogramsCJets[iInputCollection][VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		    
		    pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection]);
		    if( !pTree->cd( "udsJets/" ) )
		      {
			pTree->mkdir( "udsJets/" ) ; 
			pTree->cd(    "udsJets/" ) ;
		      }
		    
		    _inputsHistogramsUDSJets[iInputCollection][VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		  
		  }//end of histogram creation
	      }
	    
	    if (isFirstEvent()) {
	      //We'd like to make zoomed histograms of some of the flavour tag inputs too
	      for (size_t i = 0;i < _ZoomedVarNames.size();++i) {
		
		std::string zoomed_name = _ZoomedVarNames[i] + " (zoomed)";
		
		pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection]);
		  if( !pTree->cd( "bJets" ) )
		    {
		      pTree->mkdir( "bJets" ) ; 
		      pTree->cd(    "bJets" ) ;
		    }
		
		_zoomedInputsHistogramsBJets[iInputCollection][zoomed_name] = pHistogramFactory->createHistogram1D( zoomed_name, 100, -10., 20.);
		
		pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection]);
		if( !pTree->cd( "cJets" ) )
		    {
		      pTree->mkdir( "cJets" ) ; 
		      pTree->cd(    "cJets" ) ;
		    }

		_zoomedInputsHistogramsCJets[iInputCollection][zoomed_name] = pHistogramFactory->createHistogram1D( zoomed_name, 100, -10., 20.);
		
		
		pTree->cd( "/" + name() + "/" + _FlavourTagInputsCollectionNames[iInputCollection]);
		if( !pTree->cd( "udsJets" ) )
		    {
		      pTree->mkdir( "udsJets" ) ; 
		      pTree->cd(    "udsJets" ) ;
		    }

		_zoomedInputsHistogramsUDSJets[iInputCollection][zoomed_name] = pHistogramFactory->createHistogram1D( zoomed_name, 100, -10., 20.);
	      }
	    }
	  }
}

 
void LCFIAIDAPlotProcessor::processEvent( LCEvent* pEvent ) 
{ 

  // Make sure that "processRunHeader" has been called for this run (see the comment in that method).
  if( (_lastRunHeaderProcessed != pEvent->getRunNumber()) && (_suppressOutputForRun != pEvent->getRunNumber()) )
    {
      std::cerr << __FILE__ << "(" << __LINE__ << "): processRunHeader() was not called for run " << pEvent->getRunNumber()
		<< " (did you use \"SkipNEvents\"?). The order of the information in the flavour tag collection(s) is going to be guessed." << std::endl;
      
      //Only want to do this once for this run, so set a marker that this run has been done
      _suppressOutputForRun=pEvent->getRunNumber();
      
      // Just assume that the elements are in the order "BTag", "CTag", "BCTag"
      std::map<std::string,unsigned int> guessedOrder;
      guessedOrder["BTag"]=0; guessedOrder["CTag"]=1; guessedOrder["BCTag"]=2;
      _IndexOfForEachTag.clear();
      for (unsigned int iTag=0; iTag < _FlavourTagCollectionNames.size(); ++iTag) _IndexOfForEachTag.push_back( guessedOrder );
    }
  
  
  // Try and get the jet collection. If unable, show why. No need to worry about quitting
  // because getNumberOfElements will return zero, so flow will never go into the loop.
  // TypesafeCollection is just a wrapper around LCCollection with more error checking and
  // things. Just makes the code a bit easier to read (in my opinion).

  TypesafeCollection<lcio::ReconstructedParticle> jetCollection( pEvent, _JetCollectionName );
   
  
  
  //apply any cuts on the event here
  if( PassesEventCuts(pEvent) )
    {

      ReconstructedParticle* pJet;
      //loop over the jets
      for( int jetNumber=0; jetNumber<jetCollection.getNumberOfElements(); ++jetNumber )
	{
	  pJet=jetCollection.getElementAt(jetNumber);
	  
	  //only do anything if the jet passes the jet cuts
	  if( PassesJetCuts(pJet) )
	    {
	      
	      if (_MakePurityEfficiencyPlots || _PrintPurityEfficiencyValues) FillTagPlots( pEvent, jetNumber );
	      if (_MakeAdditionalPlots) FillVertexChargePlots( pEvent, jetNumber );
	      if (_MakeTuple) FillInputsPlots( pEvent, jetNumber );
	      if (_MakeAdditionalPlots) FillVertexPlots( pEvent, jetNumber );
	      //Look at the track-vertex matching
	      if (_PrintTrackVertexOutput) FillZVRESTable(pEvent);
	      
	    }
	} 
    }
}

 
 
void LCFIAIDAPlotProcessor::check( LCEvent* pEvent ) 
{
}



void LCFIAIDAPlotProcessor::end() 
{
  if (_MakePurityEfficiencyPlots||_PrintPurityEfficiencyValues) CalculateIntegralAndBackgroundPlots();
  if (_MakePurityEfficiencyPlots) CalculateEfficiencyPurityPlots();
  if (_MakeAdditionalPlots) CalculateAdditionalPlots();
  if (_PrintPurityEfficiencyValues) PrintNNOutput();  
  if (_PrintTrackVertexOutput) PrintZVRESTable();
  
}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently applies no cuts at all*/
bool LCFIAIDAPlotProcessor::PassesEventCuts( LCEvent* pEvent )
{
  //
  // No event cuts at present
  //
  
  return true;

}
// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs

bool LCFIAIDAPlotProcessor::PassesJetCuts( ReconstructedParticle* pJet )
{
  //
  // This cut added on the suggestion of Sonja Hillert 12/Jan/07.
  //
  // Selects jets for which the cosine of the jet polar
  // angle theta for all jets is not too large.
  //
  // Make something that's easy to search for to track down erroneous cuts:
  // (too many bad experiences of long forgotten about cuts hiding somewhere)
  // GREPTAG_CUT : Jet cut on abs(cos(theta)) of jet axis
  //
  
  
  const double* jetMomentum=pJet->getMomentum(); 
  
  double momentumMagnitude = sqrt(pow(jetMomentum[0],2)+pow(jetMomentum[1],2)+pow(jetMomentum[2],2));

  double cosTheta = jetMomentum[2] / momentumMagnitude;
  if( fabs(cosTheta) < _CosThetaJetMin || fabs(cosTheta) > _CosThetaJetMax ) return false;
  
  if (momentumMagnitude > _PJetMax ||  momentumMagnitude < _PJetMin) return false;

  
  // If control gets to this point then the jet has passed
  return true;
}


void LCFIAIDAPlotProcessor::CalculateIntegralAndBackgroundPlots() {
  
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  
  _pBJetBTagIntegral.resize( _FlavourTagCollectionNames.size() );     
  _pCJetBTagIntegral.resize( _FlavourTagCollectionNames.size() );	  
  _pLightJetBTagIntegral.resize( _FlavourTagCollectionNames.size() ); 
  _pBJetCTagIntegral.resize( _FlavourTagCollectionNames.size() );  	  
  _pCJetCTagIntegral.resize( _FlavourTagCollectionNames.size() );	  
  _pLightJetCTagIntegral.resize( _FlavourTagCollectionNames.size() ); 
  _pBJetBCTagIntegral.resize( _FlavourTagCollectionNames.size() ); 	  
  _pCJetBCTagIntegral.resize( _FlavourTagCollectionNames.size() );	  
  _pLightJetBCTagIntegral.resize( _FlavourTagCollectionNames.size() );
  

  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    {
      
      for (unsigned int iVertexCat=1;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	//sum over the different vertex catagories, this information goes into the "AnyNumberOfVertices" directory
	
	_pBJetBTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pBJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pBJetCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pBJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pBJetBCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pBJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCJetBTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pCJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCJetCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pCJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCJetBCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pCJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pLightJetBTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pLightJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pLightJetCTag[iTagCollection][_VertexCatNames[0]] -> 	add(*_pLightJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pLightJetBCTag[iTagCollection][_VertexCatNames[0]] -> add(*_pLightJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
      }

       for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	//add up all the background values
	 
	pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection]+ "/" +_NumVertexCatDir[iVertexCat]);
	
	_pBTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-B jets by B-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetBTag[iTagCollection][_VertexCatNames[iVertexCat]],*_pCJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]);
	_pCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-C jets by C-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetCTag[iTagCollection][_VertexCatNames[iVertexCat]],*_pBJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]); 
	_pBCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-C jets by BC-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]],*_pBJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]);
      }
       
    }
}



void LCFIAIDAPlotProcessor::CalculateEfficiencyPurityPlots()
{
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

#ifdef USING_JAIDA
  AIDA::IDataPointSetFactory* pDataPointSetFactory=marlin::AIDAProcessor::dataPointSetFactory(this);
#endif


  //now calculate the efficiencies, leakage rate and purity
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    {
      
      for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {

	pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection]+ "/" +_NumVertexCatDir[iVertexCat] );
	
	std::string nvname = _VertexCatNames[iVertexCat];
	
#ifdef USING_JAIDA
	AIDA::IDataPointSet* _pBJetBTagEfficiency = CreateEfficiencyPlot( _pBJetBTag[iTagCollection][nvname] , pDataPointSetFactory->create("B-Tag efficiency  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetCTagEfficiency = CreateEfficiencyPlot( _pCJetCTag[iTagCollection][nvname] , pDataPointSetFactory->create("C-Tag efficiency  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetBCTagEfficiency = CreateEfficiencyPlot( _pCJetBCTag[iTagCollection][nvname] , pDataPointSetFactory->create("BC-Tag efficiency  ("+ nvname +")",2));
#endif

	
	_pBJetBTagIntegral[iTagCollection][nvname] =   	
	  CreateIntegralHistogram( _pBJetBTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pBJetBTag[iTagCollection][nvname]->axis().bins(),_pBJetBTag[iTagCollection][nvname]->axis().lowerEdge(),_pBJetBTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pCJetBTagIntegral[iTagCollection][nvname] =     
	  CreateIntegralHistogram( _pCJetBTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pCJetBTag[iTagCollection][nvname]->axis().bins(),_pCJetBTag[iTagCollection][nvname]->axis().lowerEdge(),_pCJetBTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pLightJetBTagIntegral[iTagCollection][nvname] = 
	  CreateIntegralHistogram( _pLightJetBTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pLightJetBTag[iTagCollection][nvname]->axis().bins(),_pLightJetBTag[iTagCollection][nvname]->axis().lowerEdge(),_pLightJetBTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pBJetCTagIntegral[iTagCollection][nvname] = 
	CreateIntegralHistogram( _pBJetCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pBJetCTag[iTagCollection][nvname]->axis().bins(),_pBJetCTag[iTagCollection][nvname]->axis().lowerEdge(),_pBJetCTag[iTagCollection][nvname]->axis().upperEdge()));
      
	_pCJetCTagIntegral[iTagCollection][nvname] =     
	CreateIntegralHistogram( _pCJetCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pCJetCTag[iTagCollection][nvname]->axis().bins(),_pCJetCTag[iTagCollection][nvname]->axis().lowerEdge(),_pCJetCTag[iTagCollection][nvname]->axis().upperEdge()));

	_pLightJetCTagIntegral[iTagCollection][nvname] = 
      CreateIntegralHistogram( _pLightJetCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pLightJetCTag[iTagCollection][nvname]->axis().bins(),_pLightJetCTag[iTagCollection][nvname]->axis().lowerEdge(),_pLightJetCTag[iTagCollection][nvname]->axis().upperEdge()));
  
	_pBJetBCTagIntegral[iTagCollection][nvname] =    
	  CreateIntegralHistogram( _pBJetBCTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pBJetBCTag[iTagCollection][nvname]->axis().bins(),_pBJetBCTag[iTagCollection][nvname]->axis().lowerEdge(),_pBJetBCTag[iTagCollection][nvname]->axis().upperEdge()));

	_pCJetBCTagIntegral[iTagCollection][nvname] =   
	  CreateIntegralHistogram( _pCJetBCTag[iTagCollection][nvname], 
				   pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pCJetBCTag[iTagCollection][nvname]->axis().bins(),_pCJetBCTag[iTagCollection][nvname]->axis().lowerEdge(),_pCJetBCTag[iTagCollection][nvname]->axis().upperEdge()));
	
	_pLightJetBCTagIntegral[iTagCollection][nvname] = 
	  CreateIntegralHistogram( _pLightJetBCTag[iTagCollection][nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",
								    _pLightJetBCTag[iTagCollection][nvname]->axis().bins(),_pLightJetBCTag[iTagCollection][nvname]->axis().lowerEdge(),_pLightJetBCTag[iTagCollection][nvname]->axis().upperEdge()));
	
	//Examples of the integral plots - instead of histograms - the histogram calculate the errors wrongly
	
	//integralplots	_pBJetBTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pBJetBTag[iTagCollection][nvname], pDataPointSetFactory->create("B-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pCJetBTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pCJetBTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pLightJetBTagIntegral[iTagCollection][nvname] = 
	//integralplots//	  CreateIntegralPlot( _pLightJetBTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pBJetCTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pBJetCTag[iTagCollection][nvname], pDataPointSetFactory->create("B-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pCJetCTagIntegral[iTagCollection][nvname] =     
	//integralplots//	  CreateIntegralPlot( _pCJetCTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pLightJetCTagIntegral[iTagCollection][nvname] = 
	//integralplots//	  CreateIntegralPlot( _pLightJetCTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pBJetBCTagIntegral[iTagCollection][nvname] =    
	//integralplots//	  CreateIntegralPlot( _pBJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("B-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pCJetBCTagIntegral[iTagCollection][nvname] =    
	//integralplots//	  CreateIntegralPlot( _pCJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",2));
	//integralplots//	_pLightJetBCTagIntegral[iTagCollection][nvname] = 
	//integralplots//	  CreateIntegralPlot( _pLightJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",2));
	
#ifdef USING_JAIDA
	AIDA::IDataPointSet* _pBJetBTagPurity =  CreatePurityPlot( _pBJetBTag[iTagCollection][nvname],  _pBTagBackgroundValues[iTagCollection][nvname] , pDataPointSetFactory->create("B-Jet purity for B-Tag  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetCTagPurity =  CreatePurityPlot( _pCJetCTag[iTagCollection][nvname],  _pCTagBackgroundValues[iTagCollection][nvname] , pDataPointSetFactory->create("C-Jet purity for C-Tag  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pCJetBCTagPurity = CreatePurityPlot( _pCJetBCTag[iTagCollection][nvname], _pBJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("C-Jet purity for BC-Tag  ("+ nvname +")",2));      
	
	AIDA::IDataPointSet* _pCJetBTagLeakage =      CreateLeakageRatePlot( _pCJetBTag[iTagCollection][nvname],      pDataPointSetFactory->create("C-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pLightJetBTagLeakage =  CreateLeakageRatePlot( _pLightJetBTag[iTagCollection][nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pBJetCTagLeakage =      CreateLeakageRatePlot( _pBJetCTag[iTagCollection][nvname],      pDataPointSetFactory->create("B-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pLightJetCTagLeakage =  CreateLeakageRatePlot( _pLightJetCTag[iTagCollection][nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pBJetBCTagLeakage =     CreateLeakageRatePlot( _pBJetBCTag[iTagCollection][nvname],     pDataPointSetFactory->create("B-Jets: Leakage Rate into BC-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pLightJetBCTagLeakage = CreateLeakageRatePlot( _pLightJetBCTag[iTagCollection][nvname], pDataPointSetFactory->create("Light-Jets: Leakage Rate into BC-Tag Sample  ("+ nvname +")",2));     
	AIDA::IDataPointSet* _pNonBJetBTagLeakage =   CreateLeakageRatePlot( _pBTagBackgroundValues[iTagCollection][nvname],      pDataPointSetFactory->create("C-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pNonCJetCTagLeakage =   CreateLeakageRatePlot( _pCTagBackgroundValues[iTagCollection][nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
	AIDA::IDataPointSet* _pNonCJetBCTagLeakage =  CreateLeakageRatePlot( _pBCTagBackgroundValues[iTagCollection][nvname],      pDataPointSetFactory->create("B-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
	
	CreateXYPlot(_pBJetBTagEfficiency, _pBJetBTagPurity, pDataPointSetFactory->create("Purity-Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pCJetCTagPurity, pDataPointSetFactory->create("Purity-Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetBCTagEfficiency, _pCJetBCTagPurity, pDataPointSetFactory->create("Purity-Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	
	CreateXYPlot(_pBJetBTagEfficiency, _pNonBJetBTagLeakage,  pDataPointSetFactory->create("Leakage Rate-Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pNonCJetCTagLeakage,  pDataPointSetFactory->create("Leakage Rate-Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetBCTagEfficiency,_pNonCJetBCTagLeakage, pDataPointSetFactory->create("Leakage Rate-Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	
	CreateXYPlot(_pBJetBTagEfficiency, _pCJetBTagLeakage,  pDataPointSetFactory->create("Leakage Rate (C into B) vs Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pBJetBTagEfficiency, _pLightJetBTagLeakage,  pDataPointSetFactory->create("Leakage Rate (Light into B) Efficiency for B-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pBJetCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (B into C) vs Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pLightJetCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (Light into C) Efficiency for C-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetBCTagEfficiency, _pBJetBCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (B into BC) vs Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
	CreateXYPlot(_pCJetCTagEfficiency, _pLightJetBCTagLeakage,  pDataPointSetFactory->create("Leakage Rate (Light into BC) Efficiency for BC-Tag  ("+ nvname +")",2), 1, 1);
#endif	
	
      }//end of loop of different number of vertices
  
    }//end loop over iTagCollection

}



void LCFIAIDAPlotProcessor::CalculateAdditionalPlots()
{
  
  //AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

#ifdef USING_JAIDA
  AIDA::IDataPointSetFactory* pDataPointSetFactory=marlin::AIDAProcessor::dataPointSetFactory(this);
#endif


  //now calculate the efficiencies, leakage rate and purity
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    {
      
      for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {

	pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection]+ "/" +_NumVertexCatDir[iVertexCat] );
	
	std::string nvname = _VertexCatNames[iVertexCat];
		
	
      }//end of loop of different number of vertices
 
      pTree->cd( "/" + name() + "/" + _FlavourTagCollectionNames[iTagCollection]+ "/");

#ifdef USING_JAIDA
      CreateEfficiencyPlot2( _pCDecayLengthAll[iTagCollection], _pCDecayLengthTwoVertices[iTagCollection], pDataPointSetFactory->create("Efficiency to reconstruct secondary vertex vs true decay length for true C-events",2));
      CreateEfficiencyPlot2( _pBDecayLengthAll[iTagCollection], _pBDecayLengthTwoVertices[iTagCollection], pDataPointSetFactory->create("Efficiency to reconstruct secondary vertex vs true decay length for true B-events",2));
#endif 

    }//end loop over iTagCollection

  //create vertex charge leakage rate plots
  pTree->cd( "/" + name() + "/VertexChargePlots/");


#ifdef USING_JAIDA
  AIDA::IDataPointSet* pBJetVtxChargeDPS = pDataPointSetFactory->create("B-Jets: Vertex Charge Leakage",2);	
  AIDA::IDataPointSet* pCJetVtxChargeDPS = pDataPointSetFactory->create("C-Jets: Vertex Charge Leakage",2);
  CreateVertexChargeLeakagePlot(pBJetVtxChargeDPS, pCJetVtxChargeDPS);
#endif
   
}


void LCFIAIDAPlotProcessor::FillInputsPlots( LCEvent* pEvent, unsigned int jetNumber )
{  
  int jetType=FindJetType( pEvent, jetNumber );
  if( jetType==0 ) return;
  

  
  for (unsigned int iInputsCollection=0; iInputsCollection < _FlavourTagInputsCollectionNames.size(); ++iInputsCollection )
    {
      TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _FlavourTagInputsCollectionNames[iInputsCollection] );
      if( !inputsCollection.is_valid() )
	{
	  std::cerr << _FlavourTagCollectionNames[iInputsCollection] << "is not valid, not filling flavour tag input plots" << std::endl;
	  std::cerr << "set MakeAdditionalPlots to false to remove this warning" << std::endl;
	}
      else
	{
	  //Do stuff...
	  lcio::LCFloatVec* pInputs=inputsCollection.getElementAt( jetNumber );
	  if( !pInputs )
	    {
	    }
	  else
	    {
	      //ok everything is okay with the data
#ifdef USING_JAIDA	      
	      //I have a problem with tuples in RAIDA
	      
	      int CQVtx =  FindCQVtx(pEvent, jetNumber);
	      int BQVtx =  FindBQVtx(pEvent, jetNumber);
	      
	      if (_MakeTuple && iInputsCollection==0) {
		
		//this could probably be done automatically
		
		int  NumVertices = int((*pInputs)[_InputsIndex[iInputsCollection]["NumVertices"]]);
		int  NumTracksInVertices = int((*pInputs)[_InputsIndex[iInputsCollection] ["NumTracksInVertices"]]);
		float D0Significance1=(*pInputs)[_InputsIndex[iInputsCollection]     ["D0Significance1"]];
		float D0Significance2=(*pInputs)[_InputsIndex[iInputsCollection]     ["D0Significance2"]];
		float DecayLength=(*pInputs)[_InputsIndex[iInputsCollection]	       ["DecayLength"]];
		float DecayLength_SeedToIP=(*pInputs)[_InputsIndex[iInputsCollection]["DecayLength(SeedToIP)"]];
		float DecayLengthSignificance=(*pInputs)[_InputsIndex[iInputsCollection]  ["DecayLengthSignificance"]];
		float JointProbRPhi=(*pInputs)[_InputsIndex[iInputsCollection]	    ["JointProbRPhi"]];
		float JointProbZ=(*pInputs)[_InputsIndex[iInputsCollection]	       ["JointProbZ"]];
		float Momentum1=(*pInputs)[_InputsIndex[iInputsCollection]	       ["Momentum1"]];
		float Momentum2=(*pInputs)[_InputsIndex[iInputsCollection]	       ["Momentum2"]];
		float PTCorrectedMass=(*pInputs)[_InputsIndex[iInputsCollection]    ["PTCorrectedMass"]];
		float RawMomentum=(*pInputs)[_InputsIndex[iInputsCollection]	       ["RawMomentum"]];
		float SecondaryVertexProbability=(*pInputs)[_InputsIndex[iInputsCollection]["SecondaryVertexProbability"]];
		float Z0Significance1=(*pInputs)[_InputsIndex[iInputsCollection]     ["Z0Significance1"]];
		float Z0Significance2=(*pInputs)[_InputsIndex[iInputsCollection]     ["Z0Significance2"]];
		
		if (_pMyTuple) {

		  _pMyTuple->fill( 0, jetType );
		  _pMyTuple->fill( 1, NumVertices );
		  _pMyTuple->fill( 2, NumTracksInVertices );
		  _pMyTuple->fill( 3, D0Significance1);
		  _pMyTuple->fill( 4, D0Significance2);
		  _pMyTuple->fill( 5, DecayLength);
		  _pMyTuple->fill( 6, DecayLength_SeedToIP);
		  _pMyTuple->fill( 7, DecayLengthSignificance);
		  _pMyTuple->fill( 8, JointProbRPhi);
		  _pMyTuple->fill( 9, JointProbZ);
		  _pMyTuple->fill( 10, Momentum1);
		  _pMyTuple->fill( 11, Momentum2);
		  _pMyTuple->fill( 12, PTCorrectedMass);
		  _pMyTuple->fill( 13, RawMomentum);
		  _pMyTuple->fill( 14, SecondaryVertexProbability);
		  _pMyTuple->fill( 15, Z0Significance1);
		  _pMyTuple->fill( 16 ,Z0Significance2);
		  
		  _pMyTuple->fill( 17, BQVtx );
		  _pMyTuple->fill( 18, CQVtx );
		  
		  _pMyTuple->addRow();
	
		}//endif _pMyTuple
	      }
#endif
	      
	      for( std::map<std::string,unsigned int>::iterator iTagNames=_InputsIndex[iInputsCollection].begin(); 
		   iTagNames!=_InputsIndex[iInputsCollection].end(); ++iTagNames ) {
		
		double input=(*pInputs)[(*iTagNames).second]; 
		
		//if the quantity relates to the second vertex, and there is no second vertex, then don't plot it		    
		if (! ((*pInputs)[_InputsIndex[iInputsCollection]["NumVertices"]] < 2 && 
		       ((*iTagNames).first == "DecayLength" || (*iTagNames).first == "RawMomentum"  ||
			(*iTagNames).first == "SecondaryVertexProbability" || (*iTagNames).first == "PTCorrectedMass" ||
			(*iTagNames).first == "DecayLength(SeedToIP)" || (*iTagNames).first == "DecayLengthSignificance") )) {
		  
		  if( jetType==B_JET ) _inputsHistogramsBJets[iInputsCollection][(*iTagNames).first]->fill(input);
		  else if( jetType==C_JET ) _inputsHistogramsCJets[iInputsCollection][(*iTagNames).first]->fill(input);
		  else _inputsHistogramsUDSJets[iInputsCollection][(*iTagNames).first]->fill(input);
		}
		
		//fill a few extra histograms created by hand
		for (std::vector<std::string>::const_iterator iter = _ZoomedVarNames.begin() ; iter != _ZoomedVarNames.end(); iter++){
		  if ((*iTagNames).first == *iter) {
		    std::string zoomed_name = (*iTagNames).first + " (zoomed)";
		    if( jetType==B_JET ) _zoomedInputsHistogramsBJets[iInputsCollection][zoomed_name]->fill(input);
		    else if( jetType==C_JET ) _zoomedInputsHistogramsCJets[iInputsCollection][zoomed_name]->fill(input);
		    else _zoomedInputsHistogramsUDSJets[iInputsCollection][zoomed_name]->fill(input);
		    
		  }
		}
		
	      }
	      
	    }
	}
    }
}



void LCFIAIDAPlotProcessor::FillTagPlots( LCEvent* pEvent, unsigned int jetNumber)
{
  int jetType=FindJetType( pEvent, jetNumber );
  if( jetType==0 ) return;
  
  //needs to be tidied up
  TypesafeCollection<lcio::ReconstructedParticle> jetCollection( pEvent, _JetCollectionName );
  ReconstructedParticle* pJet;
  pJet=jetCollection.getElementAt(jetNumber);
  
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection)
    {
      TypesafeCollection<lcio::LCFloatVec> tagCollection( pEvent, _FlavourTagCollectionNames[iTagCollection] );
      if( !tagCollection.is_valid() )
	{
	  std::cerr << _FlavourTagCollectionNames[iTagCollection] << "is not valid, not filling tag plots" << std::endl;
	  std::cerr << "set MakeAdditionalPlots to false to remove this warning" << std::endl;

	}
      else
	{
	  lcio::LCFloatVec* pJetFlavourTags=tagCollection.getElementAt( jetNumber );
	  if( !pJetFlavourTags )
	    {
	    }
	  else
	    {
	      
	      //const double* jetMomentum = pJet->getMomentum();
	      //double cosTheta = jetMomentum[2] / sqrt(pow(jetMomentum[0],2)+pow(jetMomentum[1],2)+pow(jetMomentum[2],2));
	      
	      double bTag= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["BTag"]];
	      double cTag= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["CTag"]];
	      double cTagBBack= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["BCTag"]];
	      unsigned int NumVertices = FindNumVertex(pEvent, jetNumber, iTagCollection);
	      //int CQVtx =  FindCQVtx(pEvent, jetNumber);
	      //int BQVtx =  FindBQVtx(pEvent, jetNumber);
	      //int trueJetCharge = int(FindJetHadronCharge(pEvent,jetNumber));

	      std::string nvname = _VertexCatNames[ (NumVertices>=N_VERTEX_CATEGORIES) ? (N_VERTEX_CATEGORIES) : (NumVertices)];
	      
	      if( jetType==B_JET )  {

		if( bTag<=1 && bTag>=0 )
		  {
		      _pBJetBTag[iTagCollection][nvname]->fill( bTag );
		    } 
		  else 
		    {
		      _pBJetBTag[iTagCollection][nvname]->fill( -0.005 );
		    }
		  
		  if( cTag<=1 && cTag>=0 )
		    {
		      _pBJetCTag[iTagCollection][nvname]->fill( cTag );
		    }
		  else 
		    {
		      _pBJetCTag[iTagCollection][nvname]->fill( -0.005 );
		    }
		  if( cTagBBack<=1 && cTagBBack>=0 ) 
		    {
		      _pBJetBCTag[iTagCollection][nvname]->fill( cTagBBack );
		    }
		  else 
		    {
		      _pBJetBCTag[iTagCollection][nvname]->fill( -0.005 );
		    }
		
	      } else if( jetType==C_JET ) {
		
		if( bTag<=1 && bTag>=0 )
		  {
		    _pCJetBTag[iTagCollection][nvname]->fill( bTag );
		  } 
		else 
		  {
		    _pCJetBTag[iTagCollection][nvname]->fill( -0.005 );
		  }
		
		if( cTag<=1 && cTag>=0 )
		  {
		      _pCJetCTag[iTagCollection][nvname]->fill( cTag );
		  }
		else 
		  {
		    _pCJetCTag[iTagCollection][nvname]->fill( -0.005 );
		  }
		
		if( cTagBBack<=1 && cTagBBack>=0 ) 
		  {
		    _pCJetBCTag[iTagCollection][nvname]->fill( cTagBBack );
		  }
		else 
		  {
		    _pCJetBCTag[iTagCollection][nvname]->fill( -0.005 );
		  }
	      } else {
		if( bTag<=1 && bTag>=0 )
		  {
		    _pLightJetBTag[iTagCollection][nvname]->fill( bTag );
		  } 
		else 
		  {
		    _pLightJetBTag[iTagCollection][nvname]->fill( -0.005 );
		  }
		
		if( cTag<=1 && cTag>=0 )
		  {
		    _pLightJetCTag[iTagCollection][nvname]->fill( cTag );
		  }
		else 
		  {
		    _pLightJetCTag[iTagCollection][nvname]->fill( -0.005 );
		  } 
		
		if( cTagBBack<=1 && cTagBBack>=0 ) 
		  {
		    _pLightJetBCTag[iTagCollection][nvname]->fill( cTagBBack );
		  }
		else 
		  {
		    _pLightJetBCTag[iTagCollection][nvname]->fill( -0.005 );
		  }
	      }
	    }
	}
    }
}

void LCFIAIDAPlotProcessor::CreateAdditionalPlots()
{
  _pBJetCharge.resize( _FlavourTagCollectionNames.size() );
  _pCJetCharge.resize( _FlavourTagCollectionNames.size() );
  _pCDecayLengthAll.resize( _FlavourTagCollectionNames.size() );
  _pBDecayLengthAll.resize( _FlavourTagCollectionNames.size() );
  _pCDecayLengthTwoVertices.resize( _FlavourTagCollectionNames.size() );
  _pBDecayLengthTwoVertices.resize( _FlavourTagCollectionNames.size() );


  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

     if (!pTree->cd("/" + name() + "/ZVTOPPlots/")) {
	pTree->cd( "/" + name() + "/");
	pTree->mkdir("ZVTOPPlots/");
	pTree->cd( "ZVTOPPlots/");
      }
      
      

      _decayLengthBJet2D = pHistogramFactory->createHistogram2D( "B jets: Reconstructed secondary decay length vs MC B decay length (mm)",50,0.,20.,50,0.,20.);
      _decayLengthCJet2D = pHistogramFactory->createHistogram2D( "B jets: Reconstructed sec-ter decay length vs MC D decay length (mm)",50,0.,20.,50,0.,20.);
      _decayLengthBJetCloud2D = pHistogramFactory->createCloud2D( "B jets: Reconstructed secondary decay length vs MC B decay length (mm) (cloud)");
      _decayLengthCJetCloud2D = pHistogramFactory->createCloud2D( "B jets: Reconstructed sec-ter decay length vs MC D decay length (mm) (cloud)");

      _reconstructedSecondaryDecayLength = pHistogramFactory->createHistogram1D( "All jets: Reconstructed secondary decay length (mm)",100,0.,20.);
      _reconstructedSecTerDecayLength = pHistogramFactory->createHistogram1D( "All jets: Reconstructed secondary-tertiary decay length",100,0.,20.);
       
      _recoDecayLengthBJet = pHistogramFactory->createHistogram1D( "B jets: Reconstructed secondary decay length (mm)",50,0.,20.);
      _recoDecayLengthCJet = pHistogramFactory->createHistogram1D( "C jets: Reconstructed secondary decay length (mm)",50,0.,20.);
      _recoDecayLengthBCJet = pHistogramFactory->createHistogram1D( "B jets: Reconstructed secondary-tertiary decay length (mm)",50,0.,20.);
      _recoDecayLengthLightJet = pHistogramFactory->createHistogram1D( "Light jets: Reconstructed secondary decay length (mm)",50,0.,10.);

      _nVerticesBJet = pHistogramFactory->createHistogram1D("B jets: Number of reconstructed vertices",10,-0.5,9.5);
      _nVerticesCJet = pHistogramFactory->createHistogram1D( "C jets: Number of reconstructed vertices",10,-0.5,9.5);
      _nVerticesLightJet = pHistogramFactory->createHistogram1D( "Light jets: Number of reconstructed vertices",10,-0.5,9.5);
      
      _decayLengthBJetTrue = pHistogramFactory->createHistogram1D( "B jets: MC secondary decay length (mm)",50,0.,20.);;
      _decayLengthBCJetTrue = pHistogramFactory->createHistogram1D( "B jets: MC secondary-tertiary decay length (mm)",50,0.,20.);;
      _decayLengthBJetTrue = pHistogramFactory->createHistogram1D( "C jets: MC secondary decay length (mm)",50,0.,20.);;  
      

      //some plots of vertex charge
      if (!pTree->cd("/" + name() + "/VertexChargePlots/")) {
	pTree->cd( "/" + name() + "/");
	pTree->mkdir("VertexChargePlots/");
	pTree->cd( "VertexChargePlots/");
      }
      
      _pBJetCharge2D = pHistogramFactory->createHistogram2D( "B Jets: Reconstructed Vertex Charge vs True Jet Charge",7,-3.5,+3.5,7,-3.5,+3.5);
      _pCJetCharge2D = pHistogramFactory->createHistogram2D( "C Jets: Reconstructed Vertex Charge vs True Jet Charge",7,-3.5,+3.5,7,-3.5,+3.5);
      
      _pBJetVertexCharge = pHistogramFactory->createHistogram1D( "B Jets: Reconstructed Vertex Charge",9,-4.5,+4.5);
      _pCJetVertexCharge = pHistogramFactory->createHistogram1D( "C Jets: Reconstructed Vertex Charge",9,-4.5,+4.5);
      
      _pCJetLeakageRate = pHistogramFactory->createHistogram1D("C Jets: Charged Leakage Rate  (DON'T TRUST ERRORS)", N_JETANGLE_BINS,0.,1.);
      _pBJetLeakageRate = pHistogramFactory->createHistogram1D("B Jets: Charged Leakage Rate  (DON'T TRUST ERRORS)", N_JETANGLE_BINS,0.,1.);

     
      
      for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
	{
	  
	  if (!pTree->cd( "/" + name() + "/" +_FlavourTagCollectionNames[iTagCollection] + "/")) {
	    pTree->cd( "/" + name() + "/"); 
	    if (!pTree->cd( _FlavourTagCollectionNames[iTagCollection])) {
	      pTree->mkdir( _FlavourTagCollectionNames[iTagCollection] + "/");
	      pTree->cd( _FlavourTagCollectionNames[iTagCollection]+ "/");
	    }
	  }

	  if (!pTree->cd( "VertexChargePlots/")) {
	    pTree->mkdir( "VertexChargePlots/");
	    pTree->cd( "VertexChargePlots/");
	  }
	    
	  _pCDecayLengthAll[iTagCollection] = pHistogramFactory->createHistogram1D( "True MC decay length of all reconstructed C-jets", 25, 0, 10.0 );
	  _pBDecayLengthAll[iTagCollection] = pHistogramFactory->createHistogram1D( "True MC decay length of all reconstructed B-jets", 25, 0, 10.5 );
      	  _pCDecayLengthTwoVertices[iTagCollection] = pHistogramFactory->createHistogram1D( "True MC decay length of C-jets with more than one vertex", 25, 0, 10.0 );
	  _pBDecayLengthTwoVertices[iTagCollection] = pHistogramFactory->createHistogram1D( "True MC decay length of B-jets with more than one vertex", 25, 0, 10.5 );
      
	  _pBJetCharge[iTagCollection] = pHistogramFactory->createHistogram2D( "B Jets: Reconstructed Jet Charge vs True Jet Charge",7,-3.5,+3.5,7,-3.5,+3.5);
	  _pCJetCharge[iTagCollection] = pHistogramFactory->createHistogram2D( "C Jets: Reconstructed Jet Charge vs True Jet Charge",7,-3.5,+3.5,7,-3.5,+3.5);
	 
	}

      pTree->cd( "/"  + name());

      if (!pTree->cd( "VertexPlots/")) {
 	pTree->mkdir( "VertexPlots/");
	pTree->cd( "VertexPlots/");
      }
      
      _pVertexDistanceFromIP = pHistogramFactory->createHistogram1D( "Reconstructed Vertex distance from IP",100, 0., 10.);
      _pVertexPositionX = pHistogramFactory->createHistogram1D( "Non-primary vertex: x-position", 100, -10., 10.) ;
      _pVertexPositionY = pHistogramFactory->createHistogram1D( "Non-primary vertex: y-position", 100, -10., 10.);
      _pVertexPositionZ = pHistogramFactory->createHistogram1D( "Non-primary vertex: z-position", 100, -10., 10.);
      
      _pPrimaryVertexPullX = pHistogramFactory->createHistogram1D( "Non-primary vertex: x-pull", 100, -10., 10.);
      _pPrimaryVertexPullY = pHistogramFactory->createHistogram1D( "Non-primary vertex: y-pull", 100, -10., 10.);
      _pPrimaryVertexPullZ = pHistogramFactory->createHistogram1D( "Non-primary vertex: z-pull", 100, -10., 10.);
      _pPrimaryVertexPositionX = pHistogramFactory->createHistogram1D( "Primary vertex: x-position", 100, -10., 10.);
      _pPrimaryVertexPositionY = pHistogramFactory->createHistogram1D( "Primary vertex: y-position", 100, -10., 10.);
      _pPrimaryVertexPositionZ = pHistogramFactory->createHistogram1D( "Primary vertex: z-position", 100, -10., 10.);
    

}


void LCFIAIDAPlotProcessor::CreateTagPlots()
{
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

  _pBJetBTag.resize( _FlavourTagCollectionNames.size() );
  _pBJetCTag.resize( _FlavourTagCollectionNames.size() );
  _pBJetBCTag.resize( _FlavourTagCollectionNames.size() );
  
  _pCJetCTag.resize( _FlavourTagCollectionNames.size() );
  _pCJetBTag.resize( _FlavourTagCollectionNames.size() );
  _pCJetBCTag.resize( _FlavourTagCollectionNames.size() );
  
  _pLightJetBTag.resize( _FlavourTagCollectionNames.size() ); 
  _pLightJetCTag.resize( _FlavourTagCollectionNames.size() ); 
  _pLightJetBCTag.resize( _FlavourTagCollectionNames.size() ); 
  
  _pBTagBackgroundValues.resize( _FlavourTagCollectionNames.size() );
  _pCTagBackgroundValues.resize( _FlavourTagCollectionNames.size() );
  _pBCTagBackgroundValues.resize( _FlavourTagCollectionNames.size() );

  
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
    { 
      for (unsigned int iVertexCat=0;  iVertexCat <  N_VERTEX_CATEGORIES+1; ++iVertexCat ){
	_pLightJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
	_pLightJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	_pLightJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	
	_pBJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;	 
	_pBJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	_pBJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;
	
	_pCJetBTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;	 
	_pCJetCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;  
	_pCJetBCTag[iTagCollection][_VertexCatNames[iVertexCat]]=0;  
	
	_pBTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
	_pCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
	_pBCTagBackgroundValues[iTagCollection][_VertexCatNames[iVertexCat]]=0; 
      }
    }

      
      for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection )
	{
	  
	  if (!pTree->cd( "/" + name() + "/" +_FlavourTagCollectionNames[iTagCollection] + "/")) {
	    pTree->cd( "/" + name() + "/"); 
	    if (!pTree->cd( _FlavourTagCollectionNames[iTagCollection])) {
	      pTree->mkdir( _FlavourTagCollectionNames[iTagCollection] + "/");
	      pTree->cd( _FlavourTagCollectionNames[iTagCollection]+ "/");
	    }
	  }	  

	  for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ){
	    
	    std::string nvname = _VertexCatNames[iVertexCat];
	    
	    
	    pTree->cd( "/" + name() + "/" +_FlavourTagCollectionNames[iTagCollection] + "/");
	    
	    if (!pTree->cd( _NumVertexCatDir[iVertexCat]+"/")) {
	      pTree->mkdir( _NumVertexCatDir[iVertexCat]+"/");
	      pTree->cd( _NumVertexCatDir[iVertexCat]+"/");
	    }	      
	    
	    _pLightJetBTag[iTagCollection][nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by B-tag NN value. ("+ nvname +")",_numberOfPoints , 0, 1.0 );
	    _pLightJetCTag[iTagCollection][nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pLightJetBCTag[iTagCollection][nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pBJetBTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of B jets by B-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
	    _pBJetCTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of B jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
	    _pBJetBCTag[iTagCollection][nvname]    = pHistogramFactory->createHistogram1D( "Numbers of B jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pCJetBTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of C jets by B-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
	    _pCJetCTag[iTagCollection][nvname]     = pHistogramFactory->createHistogram1D( "Numbers of C jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
	    _pCJetBCTag[iTagCollection][nvname]    = pHistogramFactory->createHistogram1D( "Numbers of C jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );    
	
	  }
	}
}

void LCFIAIDAPlotProcessor::FillVertexChargePlots(LCEvent* pEvent, unsigned int jetNumber)
{
  int jetType=FindJetType( pEvent, jetNumber );
  if( jetType==0 ) return;
  
  //needs to be tidied up
  TypesafeCollection<lcio::ReconstructedParticle> jetCollection( pEvent, _JetCollectionName );
  ReconstructedParticle* pJet;
  pJet=jetCollection.getElementAt(jetNumber);
  
  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection)
    {
      TypesafeCollection<lcio::LCFloatVec> tagCollection( pEvent, _FlavourTagCollectionNames[iTagCollection] );
      if( !tagCollection.is_valid() )
	{
	  std::cerr << _FlavourTagCollectionNames[iTagCollection] << "is not valid, not filling vertex charge plots" << std::endl;
	  std::cerr << "set MakeAdditionalPlots to false to remove this warning" << std::endl;
	}
      else
	{
	  lcio::LCFloatVec* pJetFlavourTags=tagCollection.getElementAt( jetNumber );
	  if( !pJetFlavourTags )
	    {
	    }
	  else
	    {
	      const double* jetMomentum = pJet->getMomentum();
	      double cosTheta = jetMomentum[2] / sqrt(pow(jetMomentum[0],2)+pow(jetMomentum[1],2)+pow(jetMomentum[2],2));
	      
	      
	      if (iTagCollection == _myVertexChargeTagCollection) {
		
		
		double bTag= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["BTag"]];
		double cTag= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["CTag"]];
		//double cTagBBack= (*pJetFlavourTags)[_IndexOfForEachTag[iTagCollection]["BCTag"]];
		unsigned int NumVertices = FindNumVertex(pEvent, jetNumber, iTagCollection);
		int CQVtx =  FindCQVtx(pEvent, jetNumber);
		int BQVtx =  FindBQVtx(pEvent, jetNumber);
		int trueJetCharge = int(FindJetHadronCharge(pEvent,jetNumber));
		
		std::string nvname = _VertexCatNames[ (NumVertices>=N_VERTEX_CATEGORIES) ? (N_VERTEX_CATEGORIES) : (NumVertices)];
		
		//vertex charge plots
		if( jetType==C_JET && cTag > _CTagNNCut) {
		  
		  int bin = _pCJetLeakageRate->coordToIndex(fabs(cosTheta));
		  
		  if (trueJetCharge==+2)   _cJet_truePlus2++;
		  if (trueJetCharge==+1)   _cJet_truePlus++;
		  if (trueJetCharge==0)    _cJet_trueNeut++;
		  if (trueJetCharge==-1)   _cJet_trueMinus++;
		  if (trueJetCharge==-2)   _cJet_trueMinus2++;
		  
		  if (trueJetCharge==+2){ 
		    if(CQVtx>0)  _cJet_truePlus2_recoPlus++; 
		    if(CQVtx==0) _cJet_truePlus2_recoNeut++;
		    if(CQVtx<0)  _cJet_truePlus2_recoMinus++;
		  }
		  if (trueJetCharge==+1){ 
		    if(CQVtx>0)  _cJet_truePlus_recoPlus++; 
		    if(CQVtx==0) _cJet_truePlus_recoNeut++;
		    if(CQVtx<0)  _cJet_truePlus_recoMinus++;
		  }
		  if (trueJetCharge==0) { 
		    if(CQVtx>0)  _cJet_trueNeut_recoPlus++; 
		    if(CQVtx==0) _cJet_trueNeut_recoNeut++;
		    if(CQVtx<0)  _cJet_trueNeut_recoMinus++;
		  }
		  if (trueJetCharge==-1)  { 
		    if(CQVtx>0)  _cJet_trueMinus_recoPlus++; 
		    if(CQVtx==0) _cJet_trueMinus_recoNeut++;
		    if(CQVtx<0)  _cJet_trueMinus_recoMinus++;
		  }
		  if (trueJetCharge==-2) { 
		    if(CQVtx>0)  _cJet_trueMinus2_recoPlus++; 
		    if(CQVtx==0) _cJet_trueMinus2_recoNeut++;
		    if(CQVtx<0)  _cJet_trueMinus2_recoMinus++;
		  }
		  
		  _pCJetVertexCharge->fill(CQVtx);
		  _pCJetCharge2D->fill(trueJetCharge,CQVtx);

		  if (trueJetCharge==+2)   _cJet_truePlus2_angle[bin]++;
		  if (trueJetCharge==+1)   _cJet_truePlus_angle[bin]++;
		  if (trueJetCharge==0)    _cJet_trueNeut_angle[bin]++;
		  if (trueJetCharge==-1)   _cJet_trueMinus_angle[bin]++;
		  if (trueJetCharge==-2)   _cJet_trueMinus2_angle[bin]++;
		  
		  if (trueJetCharge==+2){ 
		    if(CQVtx>0)  _cJet_truePlus2_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_truePlus2_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_truePlus2_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==+1){ 
		    if(CQVtx>0)  _cJet_truePlus_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_truePlus_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_truePlus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==0) { 
		    if(CQVtx>0)  _cJet_trueNeut_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_trueNeut_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_trueNeut_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==-1)  { 
		    if(CQVtx>0)  _cJet_trueMinus_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_trueMinus_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_trueMinus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==-2) { 
		    if(CQVtx>0)  _cJet_trueMinus2_recoPlus_angle[bin]++; 
		    if(CQVtx==0) _cJet_trueMinus2_recoNeut_angle[bin]++;
		    if(CQVtx<0)  _cJet_trueMinus2_recoMinus_angle[bin]++;
		  }    
		  
		} else if ( jetType==B_JET && bTag > _BTagNNCut) {
		  
		  int bin = _pBJetLeakageRate->coordToIndex(fabs(cosTheta));
		  
		  if (trueJetCharge==+2)   _bJet_truePlus2++;
		  if (trueJetCharge==+1)   _bJet_truePlus++;
		  if (trueJetCharge==0)    _bJet_trueNeut++;
		  if (trueJetCharge==-1)   _bJet_trueMinus++;
		  if (trueJetCharge==-2)   _bJet_trueMinus2++;
		  
		  if (trueJetCharge==+2){ 
		    if(BQVtx>0)  _bJet_truePlus2_recoPlus++; 
		    if(BQVtx==0) _bJet_truePlus2_recoNeut++;
		    if(BQVtx<0)  _bJet_truePlus2_recoMinus++;
		  }
		  if (trueJetCharge==+1){ 
		    if(BQVtx>0)  _bJet_truePlus_recoPlus++; 
		    if(BQVtx==0) _bJet_truePlus_recoNeut++;
		    if(BQVtx<0)  _bJet_truePlus_recoMinus++;
		  }
		  if (trueJetCharge==0) { 
		    if(BQVtx>0)  _bJet_trueNeut_recoPlus++; 
		    if(BQVtx==0) _bJet_trueNeut_recoNeut++;
		    if(BQVtx<0)  _bJet_trueNeut_recoMinus++;
		  }
		  if (trueJetCharge==-1)  { 
		    if(BQVtx>0)  _bJet_trueMinus_recoPlus++; 
		    if(BQVtx==0) _bJet_trueMinus_recoNeut++;
		    if(BQVtx<0)  _bJet_trueMinus_recoMinus++;
		  }
		  if (trueJetCharge==-2) { 
		    if(BQVtx>0)  _bJet_trueMinus2_recoPlus++; 
		    if(BQVtx==0) _bJet_trueMinus2_recoNeut++;
		    if(BQVtx<0)  _bJet_trueMinus2_recoMinus++;
		  }
		  
		  
		  if (trueJetCharge==+2) _bJet_truePlus2_angle[bin]++;
		  if (trueJetCharge==+1) _bJet_truePlus_angle[bin]++;	
		  if (trueJetCharge==0)  _bJet_trueNeut_angle[bin]++;	
		  if (trueJetCharge==-1) _bJet_trueMinus_angle[bin]++;	
		  if (trueJetCharge==-2) _bJet_trueMinus2_angle[bin]++;
		  
		  if (trueJetCharge==+2){ 
		    if(BQVtx>0)  _bJet_truePlus2_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_truePlus2_recoNeut_angle[bin]++;
		    if(BQVtx<0)  _bJet_truePlus2_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==+1){ 
		    if(BQVtx>0)  _bJet_truePlus_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_truePlus_recoNeut_angle[bin]++;
		    if(BQVtx<0)  _bJet_truePlus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==0) { 
		    if(BQVtx>0) _bJet_trueNeut_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_trueNeut_recoNeut_angle[bin]++;
		    if(BQVtx<0) _bJet_trueNeut_recoMinus_angle[bin]++;
		}
		  if (trueJetCharge==-1)  { 
		    if(BQVtx>0) _bJet_trueMinus_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_trueMinus_recoNeut_angle[bin]++;
		    if(BQVtx<0) _bJet_trueMinus_recoMinus_angle[bin]++;
		  }
		  if (trueJetCharge==-2) { 
		    if(BQVtx>0) _bJet_trueMinus2_recoPlus_angle[bin]++; 
		    if(BQVtx==0) _bJet_trueMinus2_recoNeut_angle[bin]++;
		    if(BQVtx<0) _bJet_trueMinus2_recoMinus_angle[bin]++;
		  }
		  
		  _pBJetVertexCharge->fill(BQVtx);
		  _pBJetCharge2D->fill(trueJetCharge,BQVtx);
		}
	      }
	    }
	}
    }
}





int LCFIAIDAPlotProcessor::FindJetPDGCode( LCEvent* pEvent, unsigned int jetNumber )
{
	TypesafeCollection<lcio::LCIntVec> trueJetPDGCodeCollection( pEvent, _TrueJetPDGCodeColName );
	if( !trueJetPDGCodeCollection.is_valid() )
	{
	  std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " <<  _TrueJetPDGCodeColName << " is not valid " << std::endl;
	  return 0; //can't do anything without this collection
	}

	int pdgCode;
	lcio::LCIntVec* pTruePDGCodeVector=trueJetPDGCodeCollection.getElementAt( jetNumber );
	if( pTruePDGCodeVector )
	{
		if( pTruePDGCodeVector->size()==1 ) pdgCode=pTruePDGCodeVector->back();
		else
		{
		  std::cerr << __FILE__ << "(" << __LINE__ << "): The LCIntVec for jet " << jetNumber << " from the collection "
			    << _TrueJetFlavourColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
			    << " is not of size 1." << std::endl;
		  return 0; //can't fill any plots if we don't know the true flavour
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCIntVec for jet " << jetNumber << " from the collection " << _TrueJetFlavourColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return pdgCode;
}

float LCFIAIDAPlotProcessor::FindJetPartonCharge( LCEvent* pEvent, unsigned int jetNumber )
{
	TypesafeCollection<lcio::LCFloatVec> trueJetPartonChargeCollection( pEvent, _TrueJetPartonChargeColName );
	if( !trueJetPartonChargeCollection.is_valid() )
	{
	  std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " <<  _TrueJetPartonChargeColName << " is not valid " << std::endl;
		return 0; //can't do anything without this collection
	}

	float partonCharge;
	lcio::LCFloatVec* pTruePartonChargeVector=trueJetPartonChargeCollection.getElementAt( jetNumber );
	if( pTruePartonChargeVector )
	{
		if( pTruePartonChargeVector->size()==1 ) partonCharge=pTruePartonChargeVector->back();
		else
		{
			std::cerr << __FILE__ << "(" << __LINE__ << "): The LCFloatVec for jet " << jetNumber << " from the collection "
					<< _TrueJetPartonChargeColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
					<< " is not of size 1." << std::endl;
			return 0; //can't fill any plots if we don't know the true flavour
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCFloatVec for jet " << jetNumber << " from the collection " << _TrueJetPartonChargeColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return partonCharge;
}


int LCFIAIDAPlotProcessor::FindJetType( LCEvent* pEvent, unsigned int jetNumber )
{
	TypesafeCollection<lcio::LCIntVec> trueJetFlavourCollection( pEvent, _TrueJetFlavourColName );
	if( !trueJetFlavourCollection.is_valid() )
	  {
	    std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " <<  _TrueJetFlavourColName << " is not valid " << std::endl;
		return 0; //can't do anything without this collection
	}

	int jetType;
	lcio::LCIntVec* pTrueJetTypeVector=trueJetFlavourCollection.getElementAt( jetNumber );
	if( pTrueJetTypeVector )
	{
		if( pTrueJetTypeVector->size()==1 ) jetType=pTrueJetTypeVector->back();
		else
		{
			std::cerr << __FILE__ << "(" << __LINE__ << "): The LCIntVec for jet " << jetNumber << " from the collection "
					<< _TrueJetFlavourColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
					<< " is not of size 1." << std::endl;
			return 0; //can't fill any plots if we don't know the true flavour
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCIntVec for jet " << jetNumber << " from the collection " << _TrueJetFlavourColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return jetType;
}

float LCFIAIDAPlotProcessor::FindJetHadronCharge( LCEvent* pEvent, unsigned int jetNumber )
{
  TypesafeCollection<lcio::LCFloatVec> trueJetHadronChargeCollection( pEvent, _TrueJetHadronChargeColName );
  if( !trueJetHadronChargeCollection.is_valid() )
    {
      std::cerr << " In " << __FILE__ << "(" << __LINE__ << "):  Collection " << _TrueJetHadronChargeColName << " is not valid " << std::endl;
      return -99; //can't do anything without this collection
    }

	float hadronCharge;
	lcio::LCFloatVec* pTrueJetChargeVector=trueJetHadronChargeCollection.getElementAt( jetNumber );
	if( pTrueJetChargeVector )
	{
		if( pTrueJetChargeVector->size()==1 ) hadronCharge=pTrueJetChargeVector->back();
		else
		{
			std::cerr << __FILE__ << "(" << __LINE__ << "): The LCFloatVec for jet " << jetNumber << " from the collection "
					<< _TrueJetHadronChargeColName << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber()
					<< " is not of size 1." << std::endl;
			return -99; 
		}
	}
	else
	{
		std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCFloatVec for jet " << jetNumber << " from the collection " << _TrueJetHadronChargeColName
				<< " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
		return 0; //can't fill any plots if we don't know the true flavour
	}

	return hadronCharge;
}

int LCFIAIDAPlotProcessor::FindNumVertex( LCEvent* pEvent, unsigned int jetNumber, unsigned int iInputsCollection)
{
  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _FlavourTagInputsCollectionNames[iInputsCollection] );

  if( !inputsCollection.is_valid() ) {
    
    std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCFloatVec from the collection " << _FlavourTagInputsCollectionNames[iInputsCollection]
	      << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
    
    return 0; //can't do anything without this collection 
  }
  else
    {
      //Do stuff...
      lcio::LCFloatVec* pInputs=inputsCollection.getElementAt( jetNumber );
      
      if( !pInputs )
	{
	  std::cerr << __FILE__ << "(" << __LINE__ << "): Unable to get the LCFloatVec for jet " << jetNumber << " from the collection " << _FlavourTagInputsCollectionNames[iInputsCollection]
		    << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "." << std::endl;
	}
      else
	{
	  return  int((*pInputs)[_InputsIndex[iInputsCollection]["NumVertices"]]);
	}
      
    }
  return 0;
}

int LCFIAIDAPlotProcessor::FindBQVtx( LCEvent* pEvent, unsigned int jetNumber) 
{
  
  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _BVertexChargeCollection);
  
  if( !inputsCollection.is_valid() )  {
    
    std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection " << _BVertexChargeCollection << "  for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " BQVtx will be invalid" << std::endl;
    return -99;
    
  } else { //inputsCollection.is_valid() 
    
    float bqvtx;
    lcio::LCFloatVec* pBVtxChargeVector =inputsCollection.getElementAt( jetNumber );
    
    //bool evaluation is done left to right...
    if( pBVtxChargeVector && pBVtxChargeVector->size() == 1) {
      
      bqvtx = pBVtxChargeVector->back();
      
    } else {
      
      std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection element in  " << _BVertexChargeCollection << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " corresponding to jet number: " << jetNumber << " BQVtx will be invalid" << std::endl;
      return -99;
    }
    
    return int(bqvtx);
  } 

  //should never get here
  return -99;
}


int LCFIAIDAPlotProcessor::FindCQVtx( LCEvent* pEvent, unsigned int jetNumber) 
{
  
  TypesafeCollection<lcio::LCFloatVec> inputsCollection( pEvent, _CVertexChargeCollection);
  
  if( !inputsCollection.is_valid() )  {
    
    std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection " << _CVertexChargeCollection << "  for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " CQVtx will be invalid" << std::endl;
    return -99;
    
  } else { //inputsCollection.is_valid() 
    
    float cqvtx;
    lcio::LCFloatVec* pCVtxChargeVector =inputsCollection.getElementAt( jetNumber );
    
    //bool evaluation is done left to right...
    if( pCVtxChargeVector && pCVtxChargeVector->size() == 1) {
      
      cqvtx = pCVtxChargeVector->back();
      
    } else {
      
      std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): Cannot find collection element in  " << _CVertexChargeCollection << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << " corresponding to jet number: " << jetNumber << " CQVtx will be invalid" << std::endl;
      return -99;
    }
    
    return int(cqvtx);
  }
  //should never get here
  return -99;
}


AIDA::IHistogram1D* LCFIAIDAPlotProcessor::CreateIntegralHistogram(const AIDA::IHistogram1D* pNN, AIDA::IHistogram1D* pIntegral)
{
  //the errors on these entries are wrong...
  
  const int numberOfBins=pNN->axis().bins();

  double integral=pNN->binHeight(AIDA::IAxis::OVERFLOW_BIN);
  pIntegral->fill(pNN->axis().binLowerEdge(AIDA::IAxis::OVERFLOW_BIN)+pNN->axis().binWidth(numberOfBins-1)/2.,integral);
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber )
    {
      integral+= pNN->binHeight( binNumber );
      pIntegral->fill( pNN->axis().binLowerEdge(binNumber)+pNN->axis().binWidth(binNumber)/2.,integral);
    }
  
  integral+= pNN->binHeight(AIDA::IAxis::UNDERFLOW_BIN);
  pIntegral->fill(pNN->axis().binUpperEdge(AIDA::IAxis::UNDERFLOW_BIN)-pNN->axis().binWidth(0)/2.,integral);

  return pIntegral;
}



void  LCFIAIDAPlotProcessor::CalculateTagEfficiency(const AIDA::IHistogram1D* pSignal, std::vector<double>& tagEfficiency, std::vector<double>& tagEfficiencyError)
{
  //makes an effiency plot from one histogram - assuming the last bin contains all the events (eg NN distribution)
  double totalSignal=pSignal->sumBinHeights();
  double signalPassedCut=0; 
  const int numberOfBins=pSignal->axis().bins();

  tagEfficiency.clear();
  tagEfficiencyError.clear();

  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber)
    {
      signalPassedCut+=pSignal->binHeight( binNumber );
      
      double efficiency = signalPassedCut/totalSignal;
      
      double efficiencyError = efficiency * (1. - efficiency) / totalSignal;
      if (efficiencyError>0) efficiencyError = sqrt(efficiencyError);
      
      tagEfficiency.push_back(efficiency);
      tagEfficiencyError.push_back(efficiencyError);

    }
  return;
}

void LCFIAIDAPlotProcessor::CalculateTagPurity(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, std::vector<double>& tagPurity, std::vector<double>& tagPurityError)
{
  const int numberOfBins=pSignal->axis().bins();
  
  double signalPassedCut=0;
  double backgroundPassedCut=0;
  
  tagPurity.clear();
  tagPurityError.clear();


  for (int binNumber = numberOfBins-1; binNumber >= 0 ; --binNumber)  
    {
      
      signalPassedCut+=pSignal->binHeight( binNumber );
      backgroundPassedCut+=pBackground->binHeight( binNumber );
      
      double purity = ((signalPassedCut+backgroundPassedCut)>0) ? signalPassedCut/(signalPassedCut+backgroundPassedCut) : 1.;
      double purityError = ((signalPassedCut+backgroundPassedCut)>0) ? purity * (1. - purity) / (signalPassedCut+backgroundPassedCut) : 0.;
      if (purityError>0) purityError = sqrt(purityError);
      else  purityError = -sqrt(-purityError);
 
      tagPurity.push_back(purity);
      tagPurityError.push_back(purityError);

    }  

  return;
}



#ifdef USING_JAIDA
AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateEfficiencyPlot(const AIDA::IHistogram1D* pSignal, AIDA::IDataPointSet* pDataPointSet)
{
  //makes an effiency plot from one histogram - assuming the last bin contains all the events (eg NN distribution)
  double totalSignal=pSignal->sumBinHeights();
  double signalPassedCut=0;
    
  const int numberOfBins=pSignal->axis().bins();
  int iPoint=0;
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber, iPoint++ )
    {
      signalPassedCut+=pSignal->binHeight( binNumber );
      
      double efficiency = signalPassedCut/totalSignal;
      
      double efficiencyError = efficiency * (1. - efficiency) / totalSignal;
      if (efficiencyError>0) efficiencyError = sqrt(efficiencyError);
      
      
      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pSignal->axis().binLowerEdge(binNumber)+pSignal->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue( efficiency );
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(efficiencyError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(efficiencyError);
    }
  
  return pDataPointSet;
}  


AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateEfficiencyPlot2(const AIDA::IHistogram1D* pAllEvents, const AIDA::IHistogram1D *pPassEvents, AIDA::IDataPointSet* pDataPointSet)
{
  //I'm returning a datapoint set as then I can set the errors to the correct values - not possible for a histogram
  //make a "traditional" efficiency plot, divides the 2nd histogram by the 1st one
      
  const int numberOfBins=pAllEvents->axis().bins();
  if (pPassEvents->axis().bins() != numberOfBins)  return pDataPointSet;

  int iPoint=0;
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber, iPoint++ )
    {
      double all = pAllEvents->binHeight(binNumber);
      double pass = pPassEvents->binHeight(binNumber);

      double efficiency = (all>0) ? pass/all : 0.;
      
      double efficiencyError = (all>0) ? efficiency * (1. - efficiency) / all : 0.;
      if (efficiencyError>0.) efficiencyError = sqrt(efficiencyError);
      
      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pAllEvents->axis().binLowerEdge(binNumber)+pAllEvents->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue( efficiency );
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(efficiencyError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(efficiencyError);
    }
  
  return pDataPointSet;
}



AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateIntegralPlot(const AIDA::IHistogram1D* pNN, AIDA::IDataPointSet* pDataPointSet )
{
  //it might make more sense to make this a histogram, but then the error on the entries are wrong

  const int numberOfBins=pNN->axis().bins();

  double integral=0; 
  
  for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ )  
    {
      
      integral+= pNN->binHeight( binNumber );
      pDataPointSet->addPoint();
      pDataPointSet->point(binNumber)->coordinate(0)->setValue(pNN->axis().binLowerEdge(binNumber)+pNN->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(binNumber)->coordinate(1)->setValue(integral);
      pDataPointSet->point(binNumber)->coordinate(1)->setErrorPlus(sqrt(integral));
      pDataPointSet->point(binNumber)->coordinate(1)->setErrorMinus(sqrt(integral));
    }
  return pDataPointSet;
}

AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreatePurityPlot(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet)  
{
  const int numberOfBins=pSignal->axis().bins();
  int iPoint=0;
  
  double signalPassedCut=0;
  double backgroundPassedCut=0;

  for (int binNumber = numberOfBins-1; binNumber >= 0 ; --binNumber, iPoint++ )  
    {
      
      signalPassedCut+=pSignal->binHeight( binNumber );
      backgroundPassedCut+=pBackground->binHeight( binNumber );
      
      double purity = signalPassedCut/(signalPassedCut+backgroundPassedCut);
      double purityError = purity * (1. - purity) / (signalPassedCut+backgroundPassedCut);
      if (purityError>0) purityError = sqrt(purityError);
 
      
      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pSignal->axis().binLowerEdge(binNumber)+pSignal->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue(purity);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(purityError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(purityError);

      
    }  

  return pDataPointSet;
}

AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateLeakageRatePlot(const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet)
{
  
  double totalBackground = pBackground->sumBinHeights();
  double backgroundPassedCut=0;
  
  const int numberOfBins=pBackground->axis().bins();
  int iPoint=0;
  
  for( int binNumber=numberOfBins-1; binNumber>=0; --binNumber , iPoint++ )
    {

      backgroundPassedCut+=pBackground->binHeight( binNumber );
          
      double leakageRate = backgroundPassedCut/totalBackground;
    
      double leakageRateError = leakageRate * (1. - leakageRate) / totalBackground;
      if (leakageRateError>0) leakageRateError = sqrt(leakageRateError);
    

      pDataPointSet->addPoint();
      pDataPointSet->point(iPoint)->coordinate(0)->setValue(pBackground->axis().binLowerEdge(binNumber)+pBackground->axis().binWidth(binNumber)/2.);
      pDataPointSet->point(iPoint)->coordinate(1)->setValue(leakageRate);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorPlus(leakageRateError);
      pDataPointSet->point(iPoint)->coordinate(1)->setErrorMinus(leakageRateError);
    }
  
  return pDataPointSet;
}



AIDA::IDataPointSet* LCFIAIDAPlotProcessor::CreateXYPlot(const AIDA::IDataPointSet* pDataPointSet0, const AIDA::IDataPointSet* pDataPointSet1, AIDA::IDataPointSet* xyPointSet, const int dim0, const int dim1 )
{
  
  //need to do some comparision here
  if (pDataPointSet0->size() == pDataPointSet1->size()) {

    for (int iPoint = 0 ; iPoint != pDataPointSet1->size(); iPoint++) 
      {
	xyPointSet->addPoint();
	xyPointSet->point(iPoint)->coordinate(0)->setValue(pDataPointSet0->point(iPoint)->coordinate(dim0)->value());
	xyPointSet->point(iPoint)->coordinate(1)->setValue(pDataPointSet1->point(iPoint)->coordinate(dim1)->value());
	xyPointSet->point(iPoint)->coordinate(0)->setErrorPlus(pDataPointSet0->point(iPoint)->coordinate(dim0)->errorPlus());
	xyPointSet->point(iPoint)->coordinate(1)->setErrorPlus(pDataPointSet1->point(iPoint)->coordinate(dim1)->errorPlus());
	xyPointSet->point(iPoint)->coordinate(0)->setErrorMinus(pDataPointSet0->point(iPoint)->coordinate(dim0)->errorMinus());
	xyPointSet->point(iPoint)->coordinate(1)->setErrorMinus(pDataPointSet1->point(iPoint)->coordinate(dim1)->errorMinus());
      } 
  } else {
    
    //some error message here!
  }
  return xyPointSet;
}

void LCFIAIDAPlotProcessor::CreateVertexChargeLeakagePlot(AIDA::IDataPointSet* pBJetVtxChargeDPS, AIDA::IDataPointSet* pCJetVtxChargeDPS)
{
  for (int j = 0 ; j < N_JETANGLE_BINS ; j++) {
    
    double c_numerator   = _cJet_truePlus2_recoPlus_angle[j]+_cJet_truePlus_recoPlus_angle[j]
      +_cJet_trueMinus2_recoMinus_angle[j]+_cJet_trueMinus_recoMinus_angle[j];
    double c_domininator = _cJet_truePlus2_angle[j]         +_cJet_truePlus_angle[j]         
      +_cJet_trueMinus2_angle[j]          +_cJet_trueMinus_angle[j];
    
    double b_numerator   = _bJet_truePlus2_recoPlus_angle[j]+_bJet_truePlus_recoPlus_angle[j]
      +_bJet_trueMinus2_recoMinus_angle[j]+_bJet_trueMinus_recoMinus_angle[j];
    double b_domininator = _bJet_truePlus2_angle[j]         +_bJet_truePlus_angle[j]         
      +_bJet_trueMinus2_angle[j]          +_bJet_trueMinus_angle[j];
    
    double b_leakage = 1. - b_numerator/b_domininator;
    double c_leakage = 1. - c_numerator/c_domininator;
    
    double b_leakage_error = sqrt( b_leakage * (1. -  b_leakage) / b_domininator);
    double c_leakage_error = sqrt( c_leakage * (1. -  c_leakage) / c_domininator);
    
    
    pBJetVtxChargeDPS->addPoint();
    pBJetVtxChargeDPS->point(j)->coordinate(0)->setValue(_pBJetLeakageRate->axis().binLowerEdge(j)+_pBJetLeakageRate->axis().binWidth(j)/2.);
    pBJetVtxChargeDPS->point(j)->coordinate(1)->setValue(b_leakage);
    pBJetVtxChargeDPS->point(j)->coordinate(1)->setErrorPlus(b_leakage_error);
    pBJetVtxChargeDPS->point(j)->coordinate(1)->setErrorMinus(b_leakage_error);
    
    
    pCJetVtxChargeDPS->addPoint();
    pCJetVtxChargeDPS->point(j)->coordinate(0)->setValue(_pCJetLeakageRate->axis().binLowerEdge(j)+_pCJetLeakageRate->axis().binWidth(j)/2.);
    pCJetVtxChargeDPS->point(j)->coordinate(1)->setValue(c_leakage);
    pCJetVtxChargeDPS->point(j)->coordinate(1)->setErrorPlus(c_leakage_error);
    pCJetVtxChargeDPS->point(j)->coordinate(1)->setErrorMinus(c_leakage_error);
    
    _pCJetLeakageRate->fill(_pCJetLeakageRate->axis().binLowerEdge(j)+_pCJetLeakageRate->axis().binWidth(j)/2.,c_leakage); 
    _pBJetLeakageRate->fill(_pBJetLeakageRate->axis().binLowerEdge(j)+_pBJetLeakageRate->axis().binWidth(j)/2.,b_leakage);
  }
}
#endif



void LCFIAIDAPlotProcessor::CreateVertexChargeLeakagePlot()
{
  
  for (int j = 0 ; j < N_JETANGLE_BINS ; j++) {
    
    double c_numerator   = _cJet_truePlus2_recoPlus_angle[j]+_cJet_truePlus_recoPlus_angle[j]
      +_cJet_trueMinus2_recoMinus_angle[j]+_cJet_trueMinus_recoMinus_angle[j];
    double c_domininator = _cJet_truePlus2_angle[j]         +_cJet_truePlus_angle[j]         
      +_cJet_trueMinus2_angle[j]          +_cJet_trueMinus_angle[j];
    
    double b_numerator   = _bJet_truePlus2_recoPlus_angle[j]+_bJet_truePlus_recoPlus_angle[j]
      +_bJet_trueMinus2_recoMinus_angle[j]+_bJet_trueMinus_recoMinus_angle[j];
    double b_domininator = _bJet_truePlus2_angle[j]         +_bJet_truePlus_angle[j]         
      +_bJet_trueMinus2_angle[j]          +_bJet_trueMinus_angle[j];
    
    double b_leakage = 1. - b_numerator/b_domininator;
    double c_leakage = 1. - c_numerator/c_domininator;
    
    //double b_leakage_error = sqrt( b_leakage * (1. -  b_leakage) / b_domininator);
    //double c_leakage_error = sqrt( c_leakage * (1. -  c_leakage) / c_domininator);
    
    _pCJetLeakageRate->fill(_pCJetLeakageRate->axis().binLowerEdge(j)+_pCJetLeakageRate->axis().binWidth(j)/2.,c_leakage); 
    _pBJetLeakageRate->fill(_pBJetLeakageRate->axis().binLowerEdge(j)+_pBJetLeakageRate->axis().binWidth(j)/2.,b_leakage);
    
  }
}


void LCFIAIDAPlotProcessor::FillVertexPlots(LCEvent* pEvent, unsigned int jetNumber)
{
  int jetType=FindJetType( pEvent, jetNumber );
   
  if( jetType==0 ) return;
  
  std::vector<double> mcDecayLengthVector;
  std::vector<double> bMCDecayLengthVector;
  std::vector<double> cMCDecayLengthVector;
  double b_mc_decay_length, c_mc_decay_length;

  FindTrueJetDecayLength(pEvent, jetNumber, mcDecayLengthVector, bMCDecayLengthVector, cMCDecayLengthVector);
  FindTrueJetDecayLength2(pEvent, jetNumber, b_mc_decay_length, c_mc_decay_length);


  for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection) {

    unsigned int NumVertices = FindNumVertex(pEvent, jetNumber, iTagCollection);

    if( jetType==B_JET )  {
      _pBDecayLengthAll[iTagCollection]->fill(b_mc_decay_length);
      if (NumVertices > 1)  _pBDecayLengthTwoVertices[iTagCollection]->fill(b_mc_decay_length);
    } else if( jetType==C_JET ) {
      _pCDecayLengthAll[iTagCollection]->fill(c_mc_decay_length);
      if (NumVertices > 1)  _pCDecayLengthTwoVertices[iTagCollection]->fill(c_mc_decay_length);
    }

  }


  TypesafeCollection<lcio::ReconstructedParticle> zvresjetCollection( pEvent, _ZVRESSelectedJetsCollection );
  TypesafeCollection<lcio::ReconstructedParticle> zvresdcCollection( pEvent,  _ZVRESDecayChainCollection );
 
  if (!zvresjetCollection.is_valid() || !zvresdcCollection.is_valid() ) {

    std::cerr << "Input Collections " << _ZVRESSelectedJetsCollection << " and/or " << _ZVRESDecayChainCollection << " not valid.  Will not fill vertex plots." <<std::endl;
    std::cerr << "Set MakeAdditionalPlots to false to remove this warning" << std::endl;
  }
  else 
    {
      // zvresjetCollection and zvresdcCollection is valid

      


   for (int ijet = 0 ; ijet < zvresdcCollection.getNumberOfElements(); ijet++) {

        
    int jetType = FindJetType( pEvent,  ijet );

    lcio::ReconstructedParticle*  pJet=zvresdcCollection.getElementAt(ijet);

    std::vector<Vertex*> myVertexVector;
    std::pair<Vertex*,Vertex*> myVertexPair;
    std::vector< std::pair<Vertex*,Vertex*> >  myVertexPairVector;

    std::vector<lcio::ReconstructedParticle*> decayChainPartColl = pJet->getParticles();
    for (unsigned int ipart = 0 ; ipart < decayChainPartColl.size(); ipart++) {
 
      lcio::ReconstructedParticle*  decayChainPart = dynamic_cast< lcio::ReconstructedParticle*> (decayChainPartColl[ipart]);

      //want to make an exclusive list of vertices
     
      bool alreadyContainsStart = false;
      bool alreadyContainsEnd = false;

      for (std::vector<Vertex*>::const_iterator iVert=myVertexVector.begin(); iVert!=myVertexVector.end(); iVert++) {
	if  (decayChainPart->getStartVertex() == *iVert) alreadyContainsStart = true;
	if  (decayChainPart->getEndVertex() == *iVert) alreadyContainsEnd = true;	
      }
     
      if (!alreadyContainsStart && decayChainPart->getStartVertex()!=0) myVertexVector.push_back(decayChainPart->getStartVertex());
      if (!alreadyContainsEnd && decayChainPart->getEndVertex()!=0) myVertexVector.push_back(decayChainPart->getEndVertex());


      myVertexPair.first = decayChainPart->getStartVertex();
      myVertexPair.second = decayChainPart->getEndVertex();  
      bool pairAlreadyContained = false;

      for (std::vector< std::pair<Vertex*,Vertex*> >::const_iterator iter = myVertexPairVector.begin(); iter != myVertexPairVector.end(); iter++) {
      	if ((*iter).first == myVertexPair.first && (*iter).second == myVertexPair.second) pairAlreadyContained = true;
      }
      
      if (myVertexPair.first == 0 || myVertexPair.second  == 0)  pairAlreadyContained = true;
      
      if (!pairAlreadyContained) myVertexPairVector.push_back(myVertexPair);
    }
    
    //want to find distance between Primary Vertex and next vertex (aka seconday vertex); then the secondary vertex and the tertiary vertex
    
    float reconstructedSecondaryDecayLength  = -1.;
    float reconstructedSecTerDecayLength = - 1.;
    std::vector< Vertex* > secondaryVertexVector;
    std::vector< float > secondaryDecayLengthVector;
    
    for (std::vector< std::pair<Vertex*,Vertex*> >::const_iterator iter = myVertexPairVector.begin(); 
    	 iter != myVertexPairVector.end(); iter++) {

      const Vertex* startV = (*iter).first;
    
      if (startV->isPrimary()) {
	
	const Vertex* endV = (*iter).second;

	const float* startPos = startV->getPosition();
	const float* endPos = endV->getPosition();
	
	secondaryVertexVector.push_back((*iter).second);

	reconstructedSecondaryDecayLength = CalculateDistance(startPos,endPos);
	
	_reconstructedSecondaryDecayLength->fill(reconstructedSecondaryDecayLength);
	
	
	if( jetType==B_JET ) {
	  _recoDecayLengthBJet->fill(reconstructedSecondaryDecayLength);
	  
	} else if( jetType==C_JET ) {
	  _recoDecayLengthCJet->fill(reconstructedSecondaryDecayLength);
	  
	} else {
	  _recoDecayLengthLightJet->fill(reconstructedSecondaryDecayLength);
	  
	}
      }

    }

    
    for (std::vector< std::pair<Vertex*,Vertex*> >::const_iterator iter = myVertexPairVector.begin(); 
    	 iter != myVertexPairVector.end(); iter++) {
      
      const Vertex* startV = (*iter).first;
      
      bool isSecondary=false;
      for (std::vector< Vertex* >::const_iterator iter2 = secondaryVertexVector.begin(); iter2 !=secondaryVertexVector.end(); iter2++)
	{
	  if (startV == (*iter2)) isSecondary=true;
	}
      
      if (isSecondary) { 
	const Vertex* endV = (*iter).second;
        
	const float* startPos = startV->getPosition();
	const float* endPos = endV->getPosition();
	reconstructedSecTerDecayLength = CalculateDistance(startPos,endPos);
	secondaryDecayLengthVector.push_back(reconstructedSecTerDecayLength);

	_reconstructedSecTerDecayLength->fill(reconstructedSecTerDecayLength);

	if (jetType==B_JET) _recoDecayLengthBCJet->fill(reconstructedSecTerDecayLength);
	
      }
    }

    
    int nRecoVertices = myVertexVector.size();
    
    if( jetType==B_JET ) {
      _nVerticesBJet->fill(nRecoVertices);
      _decayLengthBJet2D->fill(b_mc_decay_length,reconstructedSecondaryDecayLength);
      _decayLengthBJetCloud2D->fill(b_mc_decay_length,reconstructedSecondaryDecayLength);
           
      _decayLengthCJet2D->fill(c_mc_decay_length,reconstructedSecTerDecayLength);
      _decayLengthCJetCloud2D->fill(c_mc_decay_length,reconstructedSecTerDecayLength);

      _decayLengthBJetTrue->fill(b_mc_decay_length);
      _decayLengthBCJetTrue->fill(c_mc_decay_length);

    } else if( jetType==C_JET ) {
      _nVerticesCJet->fill(nRecoVertices);
      _decayLengthBJetTrue->fill(c_mc_decay_length);
      
    } else {
      _nVerticesLightJet->fill(nRecoVertices);
    }
    
    
  }
    }

   // Plots of the vertex quantities, if there is a vertex collection 
   TypesafeCollection<lcio::Vertex> vertexCol(pEvent, _VertexColName);
   
   if (vertexCol.is_valid()) {
     
     float primaryVertexPostion[] = {0.,0.,0.};
     
     for (int ii=0 ; ii < vertexCol.getNumberOfElements() ; ii++){
       
       Vertex* myVertex =  dynamic_cast<Vertex*>(vertexCol.getElementAt(ii));
       
       const float* vertexPosition = myVertex->getPosition();
       if (myVertex->isPrimary()) {
	 primaryVertexPostion[0] = vertexPosition[0];
	 primaryVertexPostion[1] = vertexPosition[1];
	 primaryVertexPostion[2] = vertexPosition[2];
       }
     }
     
     for (int ii=0 ; ii < vertexCol.getNumberOfElements() ; ii++){
       
       Vertex* myVertex =  dynamic_cast<Vertex*>(vertexCol.getElementAt(ii));
       const float* vertexPosition = myVertex->getPosition();
       
       double px =  double(vertexPosition[0]);
       double py =  double(vertexPosition[1]);
       double pz =  double(vertexPosition[2]);
       double ex = sqrt((double)myVertex->getCovMatrix()[0]);	  
       double ey = sqrt((double)myVertex->getCovMatrix()[2]);
       double ez = sqrt((double)myVertex->getCovMatrix()[5]);     
       
       
       if (! myVertex->isPrimary() ) {
	 
	 double distanceIP = sqrt(pow(px-primaryVertexPostion[0],2)+pow(py-primaryVertexPostion[1],2)+pow(pz-primaryVertexPostion[1],2));
	 
	 _pVertexDistanceFromIP->fill(distanceIP);
	 _pVertexPositionX->fill(px);
	 _pVertexPositionY->fill(py);
	 _pVertexPositionZ->fill(pz);
	 
       } else { //it is the primary vertex
	 
	 _pPrimaryVertexPositionX->fill(px);
	 _pPrimaryVertexPositionY->fill(py);
	 _pPrimaryVertexPositionZ->fill(pz);
	 
	 _pPrimaryVertexPullX->fill(px/ex);
	 _pPrimaryVertexPullY->fill(py/ey);
	 _pPrimaryVertexPullZ->fill(pz/ez);
       }
     }
   } else {
     
       std::cerr << "Input Collections " << _VertexColName << " is not valid.  Not filling all vertex plots." <<std::endl;
       std::cerr << "Set MakeAdditionalPlots to false to remove this warning" << std::endl;

   }
   
   
   
}

void LCFIAIDAPlotProcessor::FillZVRESTable(LCEvent* pEvent)
{
  TypesafeCollection<lcio::ReconstructedParticle> ftjetCollection( pEvent, _JetCollectionName );
  TypesafeCollection<lcio::ReconstructedParticle> zvresjetCollection( pEvent, _ZVRESSelectedJetsCollection );
  TypesafeCollection<lcio::ReconstructedParticle> zvresdcrptrackCollection( pEvent,  _ZVRESDecayChainRPTracksCollection);
  TypesafeCollection<lcio::ReconstructedParticle> zvresdcCollection( pEvent,  _ZVRESDecayChainCollection );
  
  LCCollection* relationTrackCollection = pEvent->getCollection(_TrueTracksToMCPCollection);
  
  UTIL::LCRelationNavigator* MCRelationNavigator  =  (relationTrackCollection!=NULL)  ? new LCRelationNavigator(relationTrackCollection) : NULL;
  
  
  if (MCRelationNavigator==NULL)
    std::cerr << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: " <<  pEvent->getEventNumber()<< " LCCollection " << _TrueTracksToMCPCollection << " is not valid." << std::endl;
  
  if (!ftjetCollection.is_valid())
    std::cerr << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: " <<  pEvent->getEventNumber()<< " ReconstructedParticle " << _JetCollectionName << " is not valid." << std::endl;
  if (!zvresjetCollection.is_valid()) 
    std::cerr << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: " <<  pEvent->getEventNumber()<< " ReconstructedParticle " << _ZVRESSelectedJetsCollection << " is not valid." << std::endl;
  if (!zvresdcrptrackCollection.is_valid()) 
    std::cerr << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: " <<  pEvent->getEventNumber()<< " ReconstructedParticle " << _ZVRESDecayChainRPTracksCollection << " is not valid." << std::endl;
  if (!zvresdcCollection.is_valid())
    std::cerr << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: " <<  pEvent->getEventNumber()<< " ReconstructedParticle " << _ZVRESDecayChainCollection << " is not valid." << std::endl;
 

  if (!zvresdcCollection.is_valid() || !ftjetCollection.is_valid() || !zvresjetCollection.is_valid() || !zvresdcrptrackCollection.is_valid()) {
    
    std::cerr << "Vertex-Track Association cannot be analysed.  Set MakeAdditionalPlots to false to remove this warning" << std::endl;

  } else {


  //these are to store the primary, secondary and tertiary vertices in the jet
  std::vector<Vertex*> primaryVertices; 
  std::vector<Vertex*> secondaryVertices;
  std::vector<Vertex*> tertiaryVertices;
  std::vector<Vertex*> allVertices;
  //to store pairs of vertices with a track in between
  std::vector< std::pair<Vertex*,Vertex*> > myVertexPairVector;
  primaryVertices.clear();
  secondaryVertices.clear();
  tertiaryVertices.clear();
  allVertices.clear();
  myVertexPairVector.clear();


  //loop over all the jets to find the verticies.  Only ZVRESDecayChains know about their vertices.
  
  for (int iDecayChain = 0 ; iDecayChain < zvresdcCollection.getNumberOfElements(); iDecayChain++) {
    
    lcio::ReconstructedParticle*  pDecayChain=zvresdcCollection.getElementAt(iDecayChain);
  
    //to store pairs of vertices with a track in between
    std::pair<Vertex*,Vertex*> myVertexPair;
    std::vector< std::pair<Vertex*,Vertex*> > myVertexPairVector;
    
    myVertexPairVector.clear();
    
    //these are all the particles within the jet - these guys know about their vertices
    std::vector<lcio::ReconstructedParticle*> decayChainPartColl = pDecayChain->getParticles();

    for (unsigned int iPart = 0 ; iPart < decayChainPartColl.size(); iPart++) {
 
      lcio::ReconstructedParticle*  decayChainPart = dynamic_cast< lcio::ReconstructedParticle*> (decayChainPartColl[iPart]);

      myVertexPair.first = decayChainPart->getStartVertex();
      myVertexPair.second = decayChainPart->getEndVertex();  
    
      bool pairAlreadyContained = false;

      for (std::vector< std::pair<Vertex*,Vertex*> >::const_iterator iter = myVertexPairVector.begin(); iter != myVertexPairVector.end(); iter++) {
      	if ((*iter).first == myVertexPair.first && (*iter).second == myVertexPair.second) pairAlreadyContained = true;
      }
      
      if (!pairAlreadyContained) myVertexPairVector.push_back(myVertexPair);
    }

    
    
    //first find the primary, and secondary vertices, and count up the total number of vertices
    for (std::vector< std::pair<Vertex*,Vertex*> >::const_iterator iter = myVertexPairVector.begin(); iter != myVertexPairVector.end(); iter++) {
      
      //look to see if we have a primary vertex.  It will always be in the first element, due to the way we set up the pair
      if ((*iter).first && (*iter).first->isPrimary()) 
	{
	  bool primaryAlreadyContained = false;
	  for (std::vector<Vertex*>::const_iterator piter=primaryVertices.begin(); piter != primaryVertices.end(); piter++) 
	    {
	      if (*piter == (*iter).first) primaryAlreadyContained = true; 
	    }
	  if (!primaryAlreadyContained) primaryVertices.push_back((*iter).first);
	  
      
	  // if the first vertex in the pair is a primary, the second must be a secondary
	  if ((*iter).second != 0) 
	    {
	      bool secondaryAlreadyContained = false;
	      for (std::vector<Vertex*>::const_iterator siter=secondaryVertices.begin(); siter != secondaryVertices.end(); siter++) 
		{
		  if (*siter == (*iter).second) secondaryAlreadyContained = true; 
		}	  
	      if (!secondaryAlreadyContained) secondaryVertices.push_back((*iter).second);
	    }
	}
      

      // make an exclusive list of vertices
      bool firstAlreadyContained=false, secondAlreadyContained=false;
      for (std::vector<Vertex*>::const_iterator aiter = allVertices.begin(); aiter != allVertices.end(); aiter++) 
	{
	  if ( (*iter).first  && *aiter == (*iter).first)   firstAlreadyContained=true;
	  if ( (*iter).second && *aiter == (*iter).second ) secondAlreadyContained=true;
	  
	}
      if ( !firstAlreadyContained && (*iter).first)  allVertices.push_back((*iter).first);
      if (!secondAlreadyContained && (*iter).second) allVertices.push_back((*iter).second);
    }
    
    //next find the tertiary vertices   
    for (std::vector< std::pair<Vertex*,Vertex*> >::const_iterator iter = myVertexPairVector.begin(); iter != myVertexPairVector.end(); iter++) 
      {
	for (std::vector<Vertex*>::const_iterator siter = secondaryVertices.begin(); siter != secondaryVertices.end(); siter++) 
	  {
	    if ((*iter).first && (*iter).first == (*siter) && (*iter).second) 
	      {
		bool tertiaryAlreadyContained = false;
		for (std::vector<Vertex*>::const_iterator titer=tertiaryVertices.begin(); titer != tertiaryVertices.end(); titer++) 
		  {
		    if (*titer == (*iter).second) tertiaryAlreadyContained=false;
		  }
		if (!tertiaryAlreadyContained) tertiaryVertices.push_back((*iter).second);
	      }
	  }
	
      }
    
    
  }
  

  //ok now loop over all the ZVRESSelectedJets
  for (int iJet = 0 ; iJet < zvresjetCollection.getNumberOfElements(); iJet++) {
    lcio::ReconstructedParticle*  pJet=zvresjetCollection.getElementAt(iJet);

    //get the particle within the jets
    std::vector<lcio::ReconstructedParticle*> jetParticles = pJet->getParticles();
    
    //get the true jet flavour
    int jetType = FindJetType( pEvent,  iJet );

    if (jetType == B_JET || jetType == C_JET) {  
      //loop over the particles within the jet
      for (unsigned int iPart = 0 ; iPart < jetParticles.size(); iPart++) 
	{
	  lcio::ReconstructedParticle*  myJetParticle = dynamic_cast< lcio::ReconstructedParticle*> (jetParticles[iPart]);
	  
	  //do we find a match in the decay chain?
	  bool decayChainMatch = false;
	  //do we find a match in the decay chain RPs?
	  bool decayChainRPMatch = false;
	  
	  //look for particle in ZVRESDecayChain collection
	  lcio::ReconstructedParticle* myDecayChain = 0;
	  //look for particle in ZVRESDecayChainsRPTrack collection
	  lcio::ReconstructedParticle* myDecayChainRP = 0;
	  	  
	  //count the number of vertices in the DecayChain
	  std::vector<Vertex*> myVertexVector;
	  myVertexVector.clear();
	  
	  for (int iDCJet = 0 ; iDCJet < zvresdcrptrackCollection.getNumberOfElements(); iDCJet++) {

	    lcio::ReconstructedParticle* myDCJetRP = zvresdcrptrackCollection.getElementAt(iDCJet);

	    lcio::ReconstructedParticle* myDCJetRPTrack = myDCJetRP->getParticles()[0];
	    
	    if (myDCJetRPTrack==myJetParticle) {
	      decayChainRPMatch = true;
	      myDecayChainRP = myDCJetRP;
	    }
	  }

	  //look for a matching decay chain - if we have found a matching RP - use that information
	  //otherwise look directly in ZVRESDecayChain and assume they match directly to the ZVRESSelectedJet
	  
	  if (decayChainRPMatch) {
	    
	    //need to find the decay chain the RP belongs to, in order to find out the number of vertices
	    for (int iDecayChain = 0 ; iDecayChain < zvresdcCollection.getNumberOfElements(); iDecayChain++) {
	      lcio::ReconstructedParticle*  pDecayChain=zvresdcCollection.getElementAt(iDecayChain);
	      std::vector<lcio::ReconstructedParticle*> decayChainParticles = pDecayChain->getParticles();
	      for (unsigned int iDCPart = 0 ; iDCPart < decayChainParticles.size(); iDCPart++) {
		lcio::ReconstructedParticle*  myDecayChainParticle = dynamic_cast< lcio::ReconstructedParticle*> (decayChainParticles[iDCPart]);
		if (myDecayChainParticle == myDecayChainRP) {
		  decayChainMatch = true;
		  myDecayChain = pDecayChain;
		}
	      }
	    }
	  } else {
	  
	    //we are looking at the iJet-th ZVRESSelectedJet
	    //assume that the iJet-th ZVRESDecayChain refer to the same jet
	    
	    lcio::ReconstructedParticle*  pDecayChain=zvresdcCollection.getElementAt(iJet);
	    myDecayChain = pDecayChain;
	    if (pDecayChain)  decayChainMatch = true;
	  }

	  if (decayChainMatch) {
	    
	      std::vector<lcio::ReconstructedParticle*> myDecayChainParticles = myDecayChain->getParticles();
	      for (unsigned int iDCPart = 0 ; iDCPart < myDecayChainParticles.size(); iDCPart++) {
		
		lcio::ReconstructedParticle*  myDecayChainParticle = dynamic_cast< lcio::ReconstructedParticle*> (myDecayChainParticles[iDCPart]);
		
		bool startAlreadyContained = false;
		bool endAlreadyContained = false;
		
		for (std::vector<Vertex*>::const_iterator iter=myVertexVector.begin(); iter!=myVertexVector.end(); iter++)
		  {
		    if (*iter == myDecayChainParticle->getStartVertex()) startAlreadyContained = true;
		  }
		if (!startAlreadyContained && myDecayChainParticle->getStartVertex()) myVertexVector.push_back(myDecayChainParticle->getStartVertex());
		
		
		for (std::vector<Vertex*>::const_iterator iter=myVertexVector.begin(); iter!=myVertexVector.end(); iter++)
		  {
		    if (*iter == myDecayChainParticle->getEndVertex()) endAlreadyContained = true;
		  }
		if (!endAlreadyContained && myDecayChainParticle->getEndVertex())  myVertexVector.push_back(myDecayChainParticle->getEndVertex());
		
	      }
	    }
	    
	  int numVertex=myVertexVector.size();
	  
	  bool trackFromPrimaryLight = false,  trackFromSecondaryC = false, trackFromSecondaryB = false, trackFromTertiaryC = false;
	  bool trackHasNoAssociatedMCP = false;
	  
	  //find the origin of these tracks through the MCRelationNavigator.
	  std::vector<Track*> myTracks;
	  
	  if  (decayChainRPMatch) myTracks = myDecayChainRP->getParticles()[0]->getTracks();
	  else {
	    myTracks = myJetParticle->getTracks();
	   	    
	  }
	  
	  if (myTracks.size()==0) 
	    {
	      
	      //look at the matching MC version info
	      
	      trackHasNoAssociatedMCP = true;
	    
	    }
	  
	  else if (myTracks.size()>0) 
	    {
	      
	      lcio::Track* Track = myTracks[0];
	      
	      std::vector<lcio::LCObject*> RelatedMCParticles = MCRelationNavigator->getRelatedToObjects(Track);
	      
	      //VJM 14/05/08
	      //If more than one related MC particle, seems (on a brief visual inspection) that they just all comes from the same original hadron
	      // - just from different parts of the decay chain.  
	      //So at the moment I'll just use the first MCP in the list.  It will (should) give the same answer as any other. 
	      
	      lcio::MCParticle* MCP = dynamic_cast<lcio::MCParticle*>(RelatedMCParticles[0]);
	      std::vector<lcio::MCParticle*> Parents = dynamic_cast<lcio::MCParticle*>(RelatedMCParticles[0])->getParents();
	      
	      if (Parents.empty())
		std::cerr << " In " << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: " <<  pEvent->getEventNumber()<< " : MCP has no parents" << std::endl;
	      //ERROR NO PARENT;
	      else
		{
		  int truePDGCode=MCP->getPDG();
		  int truePDGFlavour = GetPDGFlavour(truePDGCode);
		  std::vector<int> pdgCodeVector;
		  std::vector<int> pdgFlavourVector;
		  std::vector<float> pdgDecayLengths;
		  
		  pdgCodeVector.push_back(truePDGCode);
		  pdgFlavourVector.push_back(truePDGFlavour);
		  
		  std::vector<lcio::MCParticle*> ParentVec0 = MCP->getParents();
		  
		  //If there are two parents, then they point to the same original flavour
		  //if (ParentVec0.size()>1) {
		  //  std::cerr << " In " << __FILE__ << "(" << __LINE__ << "), run number: " << pEvent->getRunNumber() << ", event number: ";
		  //  std::cerr <<  pEvent->getEventNumber()<< " ;  Size of ParentVec0: " << ParentVec0.size() << std::endl;
		  //}
		  if (ParentVec0.size()>0)  {
		    lcio::MCParticle* Parent0 = ParentVec0[0];
		    
		    if (Parent0) {
		      
		      pdgCodeVector.push_back(Parent0->getPDG());
		      pdgFlavourVector.push_back(GetPDGFlavour(Parent0->getPDG()));
		      std::vector<lcio::MCParticle*> ParentVec1 = Parent0->getParents();
		      pdgDecayLengths.push_back(CalculateDistance(Parent0->getVertex(),Parent0->getEndpoint()));
		      
		      if (ParentVec1.size()>0) {
			lcio::MCParticle* Parent1 = ParentVec1[0];
			
			if (Parent1)  {
			  pdgCodeVector.push_back(Parent1->getPDG());
			  pdgFlavourVector.push_back(GetPDGFlavour(Parent1->getPDG()));
			  std::vector<lcio::MCParticle*> ParentVec2 = Parent1->getParents();
			  pdgDecayLengths.push_back(CalculateDistance(Parent1->getVertex(),Parent1->getEndpoint()));
			  
			  if (ParentVec2.size()>0) {
			    lcio::MCParticle* Parent2 = ParentVec2[0];
			    
			    if (Parent2)  {
			      pdgCodeVector.push_back(Parent2->getPDG());
			      pdgFlavourVector.push_back(GetPDGFlavour(Parent2->getPDG()));
			      std::vector<lcio::MCParticle*> ParentVec3 = Parent2->getParents();			  
			      pdgDecayLengths.push_back(CalculateDistance(Parent2->getVertex(),Parent2->getEndpoint()));
			      
			      if (ParentVec3.size()>0) {
				lcio::MCParticle* Parent3 = ParentVec3[0];
				
				if (Parent3)  {
				  pdgCodeVector.push_back(Parent3->getPDG());
				  pdgFlavourVector.push_back(GetPDGFlavour(Parent3->getPDG()));
				  std::vector<lcio::MCParticle*> ParentVec4 = Parent3->getParents();
				  pdgDecayLengths.push_back(CalculateDistance(Parent3->getVertex(),Parent3->getEndpoint()));
				  
				  if (ParentVec4.size()>0) {
				    lcio::MCParticle* Parent4 = ParentVec4[0];
				    
				    if (Parent4) {
				      pdgCodeVector.push_back(Parent4->getPDG());
				      pdgFlavourVector.push_back(GetPDGFlavour(Parent4->getPDG()));
				      std::vector<lcio::MCParticle*> ParentVec5 = Parent4->getParents();
				      pdgDecayLengths.push_back(CalculateDistance(Parent4->getVertex(),Parent4->getEndpoint()));
				      
				      if (ParentVec5.size()>0) {
					lcio::MCParticle* Parent5 = ParentVec5[0];
					
					if (Parent5) {
					  pdgCodeVector.push_back(Parent5->getPDG());
					  pdgFlavourVector.push_back(GetPDGFlavour(Parent5->getPDG()));
					  std::vector<lcio::MCParticle*> ParentVec6 = Parent5->getParents();
					  pdgDecayLengths.push_back(CalculateDistance(Parent5->getVertex(),Parent5->getEndpoint()));
					  
					  if (ParentVec6.size()>0) {
					    lcio::MCParticle* Parent6 = ParentVec6[0];
					    
					    if (Parent6) {
					      pdgCodeVector.push_back(Parent6->getPDG());
					      pdgFlavourVector.push_back(GetPDGFlavour(Parent6->getPDG()));
					      pdgDecayLengths.push_back(CalculateDistance(Parent6->getVertex(),Parent6->getEndpoint()));
					      
					    }   //eh hope that's enough...
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		  
		  
		  //assign a flavour to the jet.
		  bool isLight = false;
		  bool isC = false;
		  bool isB = false;
		  
		  
		  for (std::vector<int>::const_iterator iter = pdgFlavourVector.begin();
		       iter != pdgFlavourVector.end(); iter++) {
		    if ((*iter) == 1) isLight = true;
		    if ((*iter) == C_JET) isC = true;
		    if ((*iter) == B_JET) isB = true;
		  }
		  
		  
		  //if it's all "1"s --> then primary
		  //if there's a "4" and no "5" --> then D (secondary)
		  //if there's a "4" AND a "5"  --> then B and D (secondary and tertiary)
		  
		  if (isLight && !isC && !isB) trackFromPrimaryLight = true;
		  if (isC && !isB) trackFromSecondaryC = true;
		  if (!isC && isB) trackFromSecondaryB = true;
		  if (isC && isB)  trackFromTertiaryC = true;
		  
		}
	    }
	  
	  //find the vertex depth - only possible if we've found the RP corresponding to the track
	  
	  bool isPrimary=false, isSecondary=false, isTertiary=false;
	  
	  if (decayChainRPMatch) {
	    
	    lcio::Vertex* partVertex = myDecayChainRP->getStartVertex();
	    
	    for (std::vector<Vertex*>::const_iterator piter = primaryVertices.begin(); piter != primaryVertices.end(); piter++) 
	      {
		if ((*piter) == partVertex) isPrimary=true;
	      }
	    for (std::vector<Vertex*>::const_iterator siter = secondaryVertices.begin(); siter != secondaryVertices.end(); siter++) 
	      {
		if ((*siter) == partVertex) isSecondary=true;
	      }
	    for (std::vector<Vertex*>::const_iterator titer = tertiaryVertices.begin(); titer != tertiaryVertices.end(); titer++) 
	      {
		if ((*titer) == partVertex) isTertiary=true;
	      }
	    
	  }

	  if (numVertex == 2) 
	    {
	      if (trackHasNoAssociatedMCP) 
		{
		  if (isPrimary)   _nb_twoVertex_Primary_noMCP++;
		  else if (isSecondary) _nb_twoVertex_Secondary_noMCP++;
		  else   _nb_twoVertex_Isolated_noMCP++;                 
		}
	      
	      if (trackFromSecondaryB && jetType == B_JET ) 
		{
		  if (isPrimary)         _nb_twoVertex_bTrack_Primary++;
		  else if (isSecondary)  _nb_twoVertex_bTrack_Secondary++;
		  //else if (isTertiary)   _nb_twoVertex_bTrack_Tertiary++;
		  else                   _nb_twoVertex_bTrack_Isolated++;
		} 
	      else if ((trackFromSecondaryC || trackFromTertiaryC) && jetType == B_JET ) 
		{
		  if (isPrimary)         _nb_twoVertex_cTrack_Primary++; 
		  else if (isSecondary)  _nb_twoVertex_cTrack_Secondary++;
		  //else if (isTertiary)   _nb_twoVertex_cTrack_Tertiary++;
		  else                   _nb_twoVertex_cTrack_Isolated++;
		} 
	      else if (trackFromPrimaryLight && jetType == B_JET) 
		{
		  if (isPrimary)         _nb_twoVertex_lTrack_Primary++; 
		  else if (isSecondary)  _nb_twoVertex_lTrack_Secondary++;
		  //else if (isTertiary)   _nb_twoVertex_lTrack_Tertiary++;
		  else                   _nb_twoVertex_lTrack_Isolated++;
		}
	      
	      if (trackHasNoAssociatedMCP && jetType == C_JET) {
		
		if (isPrimary)   _nc_twoVertex_Primary_noMCP++;
		else if (isSecondary) _nc_twoVertex_Secondary_noMCP++;
		else   _nc_twoVertex_Isolated_noMCP++;                 
	      }
	      
	      if (trackFromSecondaryB && jetType == C_JET ) 
		{
		  if (isPrimary)         _nc_twoVertex_bTrack_Primary++;
		  else if (isSecondary)  _nc_twoVertex_bTrack_Secondary++;
		  //else if (isTertiary)   _nc_twoVertex_bTrack_Tertiary++;
		  else                   _nc_twoVertex_bTrack_Isolated++;
		} 
	      else if ((trackFromSecondaryC || trackFromTertiaryC) && jetType == C_JET)
		{
		  if (isPrimary)         _nc_twoVertex_cTrack_Primary++; 
		  else if (isSecondary)  _nc_twoVertex_cTrack_Secondary++;
		  //else if (isTertiary)   _nc_twoVertex_cTrack_Tertiary++;
		  else                   _nc_twoVertex_cTrack_Isolated++;
		} 
	      else if (trackFromPrimaryLight && jetType == C_JET) 
		{
		  if (isPrimary)         _nc_twoVertex_lTrack_Primary++; 
		  else if (isSecondary)  _nc_twoVertex_lTrack_Secondary++;
		  //else if (isTertiary)   _nc_twoVertex_lTrack_Tertiary++;
		  else                   _nc_twoVertex_lTrack_Isolated++;
		}
	    } 
	  else if (numVertex == 3 ) 
	    {
	      if (trackHasNoAssociatedMCP && jetType == B_JET) 
		{
		  if (isPrimary)   _nb_twoVertex_Primary_noMCP++;
		  else if (isSecondary) _nb_twoVertex_Secondary_noMCP++;
		  else if (isTertiary) _nb_twoVertex_Tertiary_noMCP++;
		  else   _nb_twoVertex_Isolated_noMCP++;                 
		}
	      
	      if (trackFromSecondaryB && jetType == B_JET) 
		{
		  if (isPrimary)         _nb_threeVertex_bTrack_Primary++; 
		  else if (isSecondary)  _nb_threeVertex_bTrack_Secondary++;
		  else if (isTertiary)   _nb_threeVertex_bTrack_Tertiary++;
		  else                   _nb_threeVertex_bTrack_Isolated++;
		} 
	      else if ((trackFromSecondaryC || trackFromTertiaryC) && jetType == B_JET)
		{
		  if (isPrimary)        _nb_threeVertex_cTrack_Primary++; 
		  else if (isSecondary) _nb_threeVertex_cTrack_Secondary++;
		  else if (isTertiary)  _nb_threeVertex_cTrack_Tertiary++;
		  else                  _nb_threeVertex_cTrack_Isolated++;
		} 
	      else if (trackFromPrimaryLight && jetType == B_JET)
		{
		  if (isPrimary)         _nb_threeVertex_lTrack_Primary++; 
		  else if (isSecondary)  _nb_threeVertex_lTrack_Secondary++;
		  else if (isTertiary)   _nb_threeVertex_lTrack_Tertiary++;
		  else                   _nb_threeVertex_lTrack_Isolated++;
		}
	      if (trackHasNoAssociatedMCP && jetType == C_JET) 
		{
		  if (isPrimary)   _nc_twoVertex_Primary_noMCP++;
		  else if (isSecondary) _nc_twoVertex_Secondary_noMCP++;
		  else if (isTertiary) _nc_twoVertex_Tertiary_noMCP++;
		  else   _nc_twoVertex_Isolated_noMCP++;                 
		}
	      else if (trackFromSecondaryB && jetType == C_JET)
		{
		  if (isPrimary)         _nc_threeVertex_bTrack_Primary++; 
		  else if (isSecondary)  _nc_threeVertex_bTrack_Secondary++;
		  else if (isTertiary)   _nc_threeVertex_bTrack_Tertiary++;
		  else                   _nc_threeVertex_bTrack_Isolated++;
		} 
	      else if ((trackFromSecondaryC || trackFromTertiaryC) && jetType == C_JET)
		{
		  if (isPrimary)        _nc_threeVertex_cTrack_Primary++; 
		  else if (isSecondary) _nc_threeVertex_cTrack_Secondary++;
		  else if (isTertiary)  _nc_threeVertex_cTrack_Tertiary++;
		  else                  _nc_threeVertex_cTrack_Isolated++;
		} 
	      else if (trackFromPrimaryLight && jetType == C_JET) 
		{
		  if (isPrimary)         _nc_threeVertex_lTrack_Primary++; 
		  else if (isSecondary)  _nc_threeVertex_lTrack_Secondary++;
		  else if (isTertiary)   _nc_threeVertex_lTrack_Tertiary++;
		  else                   _nc_threeVertex_lTrack_Isolated++;
		}
	    }  
	}
    } //jet type is B_JET or C_JET
  }
  }//collections are valid
}

float LCFIAIDAPlotProcessor::CalculateDistance(const float* pos1, const float* pos2){
  return sqrt(pow((pos1[0]-pos2[0]),2)+pow((pos1[1]-pos2[1]),2)+pow((pos1[2]-pos2[2]),2));
}

double LCFIAIDAPlotProcessor::CalculateDistance(const double* pos1, const double* pos2){
  return sqrt(pow((pos1[0]-pos2[0]),2)+pow((pos1[1]-pos2[1]),2)+pow((pos1[2]-pos2[2]),2));
}

void LCFIAIDAPlotProcessor::PrintZVRESTable()
{

  //if there is a _TrackVertexOutputFile string defined use that as the output stream, if not use std::cout
  std::filebuf* fb = new std::filebuf;  
  
  std::ostream outputFile( (!_TrackVertexOutputFile.empty()) ?                                  
			   fb->open(_TrackVertexOutputFile.c_str(),
				    std::ios_base::out|std::ios_base::trunc):  
			   std::cout.rdbuf());

  
  //no longer required
  //float twoVertex_bTrack = _nb_twoVertex_bTrack_Primary+_nb_twoVertex_bTrack_Secondary+_nb_twoVertex_bTrack_Isolated;
  //float twoVertex_cTrack = _nb_twoVertex_cTrack_Primary+_nb_twoVertex_cTrack_Secondary+_nb_twoVertex_cTrack_Isolated;
  //float twoVertex_lTrack = _nb_twoVertex_lTrack_Primary+_nb_twoVertex_lTrack_Secondary+_nb_twoVertex_lTrack_Isolated;
  
  //float threeVertex_bTrack = _nb_threeVertex_bTrack_Primary+_nb_threeVertex_bTrack_Secondary+_nb_threeVertex_bTrack_Tertiary+_nb_threeVertex_bTrack_Isolated;
  //float threeVertex_cTrack = _nb_threeVertex_cTrack_Primary+_nb_threeVertex_cTrack_Secondary+_nb_threeVertex_cTrack_Tertiary+_nb_threeVertex_cTrack_Isolated;
  //float threeVertex_lTrack = _nb_threeVertex_lTrack_Primary+_nb_threeVertex_lTrack_Secondary+_nb_threeVertex_lTrack_Tertiary+_nb_threeVertex_lTrack_Isolated;
  

  float  b_twoVertex_Primary = _nb_twoVertex_bTrack_Primary+_nb_twoVertex_cTrack_Primary+_nb_twoVertex_lTrack_Primary;
  float  b_twoVertex_Secondary = _nb_twoVertex_bTrack_Secondary+_nb_twoVertex_cTrack_Secondary+_nb_twoVertex_lTrack_Secondary;
  float  b_twoVertex_Isolated = _nb_twoVertex_bTrack_Isolated+_nb_twoVertex_cTrack_Isolated+_nb_twoVertex_lTrack_Isolated;  
  float  b_twoVertex = b_twoVertex_Primary+b_twoVertex_Secondary+b_twoVertex_Isolated;

  float  b_threeVertex_Primary = _nb_threeVertex_bTrack_Primary+_nb_threeVertex_cTrack_Primary+_nb_threeVertex_lTrack_Primary;
  float  b_threeVertex_Secondary = _nb_threeVertex_bTrack_Secondary+_nb_threeVertex_cTrack_Secondary+_nb_threeVertex_lTrack_Secondary;
  float  b_threeVertex_Tertiary = _nb_threeVertex_bTrack_Tertiary+_nb_threeVertex_cTrack_Tertiary+_nb_threeVertex_lTrack_Tertiary;
  float  b_threeVertex_Isolated = _nb_threeVertex_bTrack_Isolated+_nb_threeVertex_cTrack_Isolated+_nb_threeVertex_lTrack_Isolated; 
  float  b_threeVertex = b_threeVertex_Primary+b_threeVertex_Secondary+b_threeVertex_Tertiary+b_threeVertex_Isolated;
   


  float frac_b_twoVertex_bTrack_Primary =     100.*float(_nb_twoVertex_bTrack_Primary) / float(b_twoVertex_Primary); 
  float frac_b_twoVertex_bTrack_Secondary =   100.*float(_nb_twoVertex_bTrack_Secondary) / float(b_twoVertex_Secondary);    
  //float frac_b_twoVertex_bTrack_Tertiary =    100.*float(_nb_twoVertex_bTrack_Tertiary) / float(b_twoVertex_Tertiary);	   
  float frac_b_twoVertex_bTrack_Isolated =    100.*float(_nb_twoVertex_bTrack_Isolated) / float(b_twoVertex_Isolated);	   
  					  
  float frac_b_twoVertex_cTrack_Primary =     100.*float(_nb_twoVertex_cTrack_Primary) / float(b_twoVertex_Primary); 	   
  float frac_b_twoVertex_cTrack_Secondary =   100.*float(_nb_twoVertex_cTrack_Secondary) / float(b_twoVertex_Secondary);   
  //float frac_b_twoVertex_cTrack_Tertiary =    100.*float(_nb_twoVertex_cTrack_Tertiary) / float(b_twoVertex_Tertiary);    
  float frac_b_twoVertex_cTrack_Isolated =    100.*float(_nb_twoVertex_cTrack_Isolated) / float(b_twoVertex_Isolated);    
  					  
  float frac_b_twoVertex_lTrack_Primary =     100.*float(_nb_twoVertex_lTrack_Primary) / float(b_twoVertex_Primary); 	   
  float frac_b_twoVertex_lTrack_Secondary =   100.*float(_nb_twoVertex_lTrack_Secondary) / float(b_twoVertex_Secondary);   
  //float frac_b_twoVertex_lTrack_Tertiary =    100.*float(_nb_twoVertex_lTrack_Tertiary) / float(b_twoVertex_Tertiary);	   
  float frac_b_twoVertex_lTrack_Isolated =    100.*float(_nb_twoVertex_lTrack_Isolated) / float(b_twoVertex_Isolated);	   
  					                                 
  float frac_b_threeVertex_bTrack_Primary =     100.*float(_nb_threeVertex_bTrack_Primary) / float(b_threeVertex_Primary);   
  float frac_b_threeVertex_bTrack_Secondary =   100.*float(_nb_threeVertex_bTrack_Secondary) / float(b_threeVertex_Secondary); 
  float frac_b_threeVertex_bTrack_Tertiary =    100.*float(_nb_threeVertex_bTrack_Tertiary) / float(b_threeVertex_Tertiary);  
  float frac_b_threeVertex_bTrack_Isolated =    100.*float(_nb_threeVertex_bTrack_Isolated) / float(b_threeVertex_Isolated);  
  					   
  float frac_b_threeVertex_cTrack_Primary =     100.*float(_nb_threeVertex_cTrack_Primary) / float(b_threeVertex_Primary);   
  float frac_b_threeVertex_cTrack_Secondary =   100.*float(_nb_threeVertex_cTrack_Secondary) / float(b_threeVertex_Secondary); 
  float frac_b_threeVertex_cTrack_Tertiary =    100.*float(_nb_threeVertex_cTrack_Tertiary) / float(b_threeVertex_Tertiary);  
  float frac_b_threeVertex_cTrack_Isolated =    100.*float(_nb_threeVertex_cTrack_Isolated) / float(b_threeVertex_Isolated);  
  					   
  float frac_b_threeVertex_lTrack_Primary =     100.*float(_nb_threeVertex_lTrack_Primary) / float(b_threeVertex_Primary);   
  float frac_b_threeVertex_lTrack_Secondary =   100.*float(_nb_threeVertex_lTrack_Secondary) / float(b_threeVertex_Secondary); 
  float frac_b_threeVertex_lTrack_Tertiary =    100.*float(_nb_threeVertex_lTrack_Tertiary) / float(b_threeVertex_Tertiary);  
  float frac_b_threeVertex_lTrack_Isolated =    100.*float(_nb_threeVertex_lTrack_Isolated) / float(b_threeVertex_Isolated);  

  float frac_b_twoVertex_Primary = 100.*float(b_twoVertex_Primary) / float(b_twoVertex);
  float frac_b_twoVertex_Secondary = 100.*float(b_twoVertex_Secondary) / float(b_twoVertex);
  float frac_b_twoVertex_Isolated = 100.*float(b_twoVertex_Isolated) / float(b_twoVertex);
 
  float frac_b_threeVertex_Primary = 100.*float(b_threeVertex_Primary) / float(b_threeVertex);
  float frac_b_threeVertex_Secondary = 100.*float(b_threeVertex_Secondary) / float(b_threeVertex);
  float frac_b_threeVertex_Tertiary = 100.*float(b_threeVertex_Tertiary) / float(b_threeVertex);
  float frac_b_threeVertex_Isolated = 100.*float(b_threeVertex_Isolated) / float(b_threeVertex);
 
  

//check//  //just for checks...
//check//  frac_b_twoVertex_bTrack_Primary =     float(_nb_twoVertex_bTrack_Primary);
//check//  frac_b_twoVertex_bTrack_Secondary =   float(_nb_twoVertex_bTrack_Secondary);
//check//  //frac_b_twoVertex_bTrack_Tertiary =    float(_nb_twoVertex_bTrack_Tertiary);
//check//  frac_b_twoVertex_bTrack_Isolated =    float(_nb_twoVertex_bTrack_Isolated);
//check//  				  
//check//  frac_b_twoVertex_cTrack_Primary =     float(_nb_twoVertex_cTrack_Primary);
//check//  frac_b_twoVertex_cTrack_Secondary =   float(_nb_twoVertex_cTrack_Secondary);
//check//  //frac_b_twoVertex_cTrack_Tertiary =    float(_nb_twoVertex_cTrack_Tertiary);
//check//  frac_b_twoVertex_cTrack_Isolated =    float(_nb_twoVertex_cTrack_Isolated);
//check//  
//check//  frac_b_twoVertex_lTrack_Primary =     float(_nb_twoVertex_lTrack_Primary);
//check//  frac_b_twoVertex_lTrack_Secondary =   float(_nb_twoVertex_lTrack_Secondary);
//check//  //frac_b_twoVertex_lTrack_Tertiary =    float(_nb_twoVertex_lTrack_Tertiary);
//check//  frac_b_twoVertex_lTrack_Isolated =    float(_nb_twoVertex_lTrack_Isolated);
//check//  
//check//  frac_b_threeVertex_bTrack_Primary =     float(_nb_threeVertex_bTrack_Primary);
//check//  frac_b_threeVertex_bTrack_Secondary =   float(_nb_threeVertex_bTrack_Secondary);
//check//  frac_b_threeVertex_bTrack_Tertiary =    float(_nb_threeVertex_bTrack_Tertiary);
//check//  frac_b_threeVertex_bTrack_Isolated =    float(_nb_threeVertex_bTrack_Isolated);
//check//  
//check//  frac_b_threeVertex_cTrack_Primary =     float(_nb_threeVertex_cTrack_Primary);
//check//  frac_b_threeVertex_cTrack_Secondary =   float(_nb_threeVertex_cTrack_Secondary);
//check//  frac_b_threeVertex_cTrack_Tertiary =    float(_nb_threeVertex_cTrack_Tertiary);
//check//  frac_b_threeVertex_cTrack_Isolated =    float(_nb_threeVertex_cTrack_Isolated);
//check//  
//check//  frac_b_threeVertex_lTrack_Primary =     float(_nb_threeVertex_lTrack_Primary);
//check//  frac_b_threeVertex_lTrack_Secondary =   float(_nb_threeVertex_lTrack_Secondary);
//check//  frac_b_threeVertex_lTrack_Tertiary =    float(_nb_threeVertex_lTrack_Tertiary);
//check//  frac_b_threeVertex_lTrack_Isolated =    float(_nb_threeVertex_lTrack_Isolated);
//check//  //end of checks



  float  c_twoVertex_Primary = _nc_twoVertex_bTrack_Primary+_nc_twoVertex_cTrack_Primary+_nc_twoVertex_lTrack_Primary;
  float  c_twoVertex_Secondary = _nc_twoVertex_bTrack_Secondary+_nc_twoVertex_cTrack_Secondary+_nc_twoVertex_lTrack_Secondary;
  float  c_twoVertex_Isolated = _nc_twoVertex_bTrack_Isolated+_nc_twoVertex_cTrack_Isolated+_nc_twoVertex_lTrack_Isolated;  
  float  c_twoVertex = c_twoVertex_Primary+c_twoVertex_Secondary+c_twoVertex_Isolated;

  float  c_threeVertex_Primary = _nc_threeVertex_bTrack_Primary+_nc_threeVertex_cTrack_Primary+_nc_threeVertex_lTrack_Primary;
  float  c_threeVertex_Secondary = _nc_threeVertex_bTrack_Secondary+_nc_threeVertex_cTrack_Secondary+_nc_threeVertex_lTrack_Secondary;
  float  c_threeVertex_Tertiary = _nc_threeVertex_bTrack_Tertiary+_nc_threeVertex_cTrack_Tertiary+_nc_threeVertex_lTrack_Tertiary;
  float  c_threeVertex_Isolated = _nc_threeVertex_bTrack_Isolated+_nc_threeVertex_cTrack_Isolated+_nc_threeVertex_lTrack_Isolated; 
  float  c_threeVertex = c_threeVertex_Primary+c_threeVertex_Secondary+c_threeVertex_Tertiary+c_threeVertex_Isolated;
   


  //  float frac_c_twoVertex_bTrack_Primary =     100.*float(_nc_twoVertex_bTrack_Primary) / float(c_twoVertex_Primary); 
  //float frac_c_twoVertex_bTrack_Secondary =   100.*float(_nc_twoVertex_bTrack_Secondary) / float(c_twoVertex_Secondary);    
  //float frac_c_twoVertex_bTrack_Tertiary =    100.*float(_nc_twoVertex_bTrack_Tertiary) / float(c_twoVertex_Tertiary);	   
  //float frac_c_twoVertex_bTrack_Isolated =    100.*float(_nc_twoVertex_bTrack_Isolated) / float(c_twoVertex_Isolated);	   
  					  
  float frac_c_twoVertex_cTrack_Primary =     100.*float(_nc_twoVertex_cTrack_Primary) / float(c_twoVertex_Primary); 	   
  float frac_c_twoVertex_cTrack_Secondary =   100.*float(_nc_twoVertex_cTrack_Secondary) / float(c_twoVertex_Secondary);   
  //float frac_c_twoVertex_cTrack_Tertiary =    100.*float(_nc_twoVertex_cTrack_Tertiary) / float(c_twoVertex_Tertiary);    
  float frac_c_twoVertex_cTrack_Isolated =    100.*float(_nc_twoVertex_cTrack_Isolated) / float(c_twoVertex_Isolated);    
  					  
  float frac_c_twoVertex_lTrack_Primary =     100.*float(_nc_twoVertex_lTrack_Primary) / float(c_twoVertex_Primary); 	   
  float frac_c_twoVertex_lTrack_Secondary =   100.*float(_nc_twoVertex_lTrack_Secondary) / float(c_twoVertex_Secondary);   
  //float frac_c_twoVertex_lTrack_Tertiary =    100.*float(_nc_twoVertex_lTrack_Tertiary) / float(c_twoVertex_Tertiary);	   
  float frac_c_twoVertex_lTrack_Isolated =    100.*float(_nc_twoVertex_lTrack_Isolated) / float(c_twoVertex_Isolated);	   
  					                                 
  //float frac_c_threeVertex_bTrack_Primary =     100.*float(_nc_threeVertex_bTrack_Primary) / float(c_threeVertex_Primary);   
  //float frac_c_threeVertex_bTrack_Secondary =   100.*float(_nc_threeVertex_bTrack_Secondary) / float(c_threeVertex_Secondary); 
  //float frac_c_threeVertex_bTrack_Tertiary =    100.*float(_nc_threeVertex_bTrack_Tertiary) / float(c_threeVertex_Tertiary);  
  //float frac_c_threeVertex_bTrack_Isolated =    100.*float(_nc_threeVertex_bTrack_Isolated) / float(c_threeVertex_Isolated);  
  					   
  float frac_c_threeVertex_cTrack_Primary =     100.*float(_nc_threeVertex_cTrack_Primary) / float(c_threeVertex_Primary);   
  float frac_c_threeVertex_cTrack_Secondary =   100.*float(_nc_threeVertex_cTrack_Secondary) / float(c_threeVertex_Secondary); 
  float frac_c_threeVertex_cTrack_Tertiary =    100.*float(_nc_threeVertex_cTrack_Tertiary) / float(c_threeVertex_Tertiary);  
  float frac_c_threeVertex_cTrack_Isolated =    100.*float(_nc_threeVertex_cTrack_Isolated) / float(c_threeVertex_Isolated);  
  					   
  float frac_c_threeVertex_lTrack_Primary =     100.*float(_nc_threeVertex_lTrack_Primary) / float(c_threeVertex_Primary);   
  float frac_c_threeVertex_lTrack_Secondary =   100.*float(_nc_threeVertex_lTrack_Secondary) / float(c_threeVertex_Secondary); 
  float frac_c_threeVertex_lTrack_Tertiary =    100.*float(_nc_threeVertex_lTrack_Tertiary) / float(c_threeVertex_Tertiary);  
  float frac_c_threeVertex_lTrack_Isolated =    100.*float(_nc_threeVertex_lTrack_Isolated) / float(c_threeVertex_Isolated);  

  float frac_c_twoVertex_Primary = 100.*float(c_twoVertex_Primary) / float(c_twoVertex);
  float frac_c_twoVertex_Secondary = 100.*float(c_twoVertex_Secondary) / float(c_twoVertex);
  float frac_c_twoVertex_Isolated = 100.*float(c_twoVertex_Isolated) / float(c_twoVertex);
 
  float frac_c_threeVertex_Primary = 100.*float(c_threeVertex_Primary) / float(c_threeVertex);
  float frac_c_threeVertex_Secondary = 100.*float(c_threeVertex_Secondary) / float(c_threeVertex);
  float frac_c_threeVertex_Tertiary = 100.*float(c_threeVertex_Tertiary) / float(c_threeVertex);
  float frac_c_threeVertex_Isolated = 100.*float(c_threeVertex_Isolated) / float(c_threeVertex);
 
  

//check//  //just for checks...
//check//  frac_c_twoVertex_bTrack_Primary =     float(_nc_twoVertex_bTrack_Primary);
//check//  frac_c_twoVertex_bTrack_Secondary =   float(_nc_twoVertex_bTrack_Secondary);
//check//  //frac_c_twoVertex_bTrack_Tertiary =    float(_nc_twoVertex_bTrack_Tertiary);
//check//  frac_c_twoVertex_bTrack_Isolated =    float(_nc_twoVertex_bTrack_Isolated);
//check//  				  
//check//  frac_c_twoVertex_cTrack_Primary =     float(_nc_twoVertex_cTrack_Primary);
//check//  frac_c_twoVertex_cTrack_Secondary =   float(_nc_twoVertex_cTrack_Secondary);
//check//  //frac_c_twoVertex_cTrack_Tertiary =    float(_nc_twoVertex_cTrack_Tertiary);
//check//  frac_c_twoVertex_cTrack_Isolated =    float(_nc_twoVertex_cTrack_Isolated);
//check//  
//check//  frac_c_twoVertex_lTrack_Primary =     float(_nc_twoVertex_lTrack_Primary);
//check//  frac_c_twoVertex_lTrack_Secondary =   float(_nc_twoVertex_lTrack_Secondary);
//check//  //frac_c_twoVertex_lTrack_Tertiary =    float(_nc_twoVertex_lTrack_Tertiary);
//check//  frac_c_twoVertex_lTrack_Isolated =    float(_nc_twoVertex_lTrack_Isolated);
//check//  
//check//  frac_c_threeVertex_bTrack_Primary =     float(_nc_threeVertex_bTrack_Primary);
//check//  frac_c_threeVertex_bTrack_Secondary =   float(_nc_threeVertex_bTrack_Secondary);
//check//  frac_c_threeVertex_bTrack_Tertiary =    float(_nc_threeVertex_bTrack_Tertiary);
//check//  frac_c_threeVertex_bTrack_Isolated =    float(_nc_threeVertex_bTrack_Isolated);
//check//  
//check//  frac_c_threeVertex_cTrack_Primary =     float(_nc_threeVertex_cTrack_Primary);
//check//  frac_c_threeVertex_cTrack_Secondary =   float(_nc_threeVertex_cTrack_Secondary);
//check//  frac_c_threeVertex_cTrack_Tertiary =    float(_nc_threeVertex_cTrack_Tertiary);
//check//  frac_c_threeVertex_cTrack_Isolated =    float(_nc_threeVertex_cTrack_Isolated);
//check//  
//check//  frac_c_threeVertex_lTrack_Primary =     float(_nc_threeVertex_lTrack_Primary);
//check//  frac_c_threeVertex_lTrack_Secondary =   float(_nc_threeVertex_lTrack_Secondary);
//check//  frac_c_threeVertex_lTrack_Tertiary =    float(_nc_threeVertex_lTrack_Tertiary);
//check//  frac_c_threeVertex_lTrack_Isolated =    float(_nc_threeVertex_lTrack_Isolated);
//check//  //end of checks



  outputFile << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "                                   B-JETS                                     " << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "               Purity of reconstructed track-vertex association (%)           " << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "    Monte Carlo               Reconstructed track-vertex association         " << std::endl;
  outputFile << "    track origin        Two-vertex case               Three-vertex case       " << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "                     Pri.    Sec.    Iso.       Pri.    Sec.    Ter.    Iso. " << std::endl;   
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "    Primary       ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_lTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_lTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_lTrack_Isolated << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_lTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_lTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_lTrack_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_lTrack_Isolated << " " ;
  outputFile << std::endl;
  
  outputFile << "    B decay       ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_bTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_bTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_bTrack_Isolated << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_bTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_bTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_bTrack_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_bTrack_Isolated << " " ;
  outputFile << std::endl;
  
  outputFile << "    D decay       ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_cTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_cTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_cTrack_Isolated << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_cTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_cTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_cTrack_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_cTrack_Isolated << " " ;
  outputFile << std::endl;
  outputFile << "   ----------------------------------------------------------------------------" << std::endl;
  outputFile << "    All above     ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_Primary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_Secondary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_twoVertex_Isolated  << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_Primary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_Secondary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_b_threeVertex_Isolated  << " " ;
  outputFile << std::endl;
  outputFile << "   ----------------------------------------------------------------------------" << std::endl;
  outputFile << std::endl;
  outputFile << std::endl;
 

  outputFile << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "                                   C-JETS                                     " << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "               Purity of reconstructed track-vertex association (%)           " << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "    Monte Carlo               Reconstructed track-vertex association         " << std::endl;
  outputFile << "    track origin        Two-vertex case               Three-vertex case       " << std::endl;
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "                     Pri.    Sec.    Iso.       Pri.    Sec.    Ter.    Iso. " << std::endl;   
  outputFile << "   ---------------------------------------------------------------------------" << std::endl;
  outputFile << "    Primary       ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_lTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_lTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_lTrack_Isolated << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_lTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_lTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_lTrack_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_lTrack_Isolated << " " ;
  outputFile << std::endl;
  
//check//  outputFile << "    B decay       ";
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_twoVertex_bTrack_Primary << " " ;
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_twoVertex_bTrack_Secondary << " " ;
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_twoVertex_bTrack_Isolated << "   " ;
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_threeVertex_bTrack_Primary << " " ;
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_threeVertex_bTrack_Secondary << " " ;
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_threeVertex_bTrack_Tertiary << " " ;
//check//  outputFile.width(7);
//check//  outputFile.precision(3);
//check//  outputFile << frac_c_threeVertex_bTrack_Isolated << " " ;
//check//  outputFile << std::endl;

  outputFile << "    D decay       ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_cTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_cTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_cTrack_Isolated << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_cTrack_Primary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_cTrack_Secondary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_cTrack_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_cTrack_Isolated << " " ;
  outputFile << std::endl;
  outputFile << "   ----------------------------------------------------------------------------" << std::endl;
  outputFile << "    All above     ";
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_Primary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_Secondary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_twoVertex_Isolated  << "   " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_Primary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_Secondary  << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_Tertiary << " " ;
  outputFile.width(7);
  outputFile.precision(3);
  outputFile << frac_c_threeVertex_Isolated  << " " ;
  outputFile << std::endl;
  outputFile << "   ----------------------------------------------------------------------------" << std::endl;
  outputFile << std::endl;
  outputFile << std::endl;
 


}

void LCFIAIDAPlotProcessor::PrintNNOutput(){
  
  //if there is a _PurityEfficiencyOutputFile string defined use that as the output stream, if not use std::cout
  std::filebuf* fb = new std::filebuf;  
  
  std::ostream outputFile( (!_PurityEfficiencyOutputFile.empty()) ?                                  
		       fb->open(_PurityEfficiencyOutputFile.c_str(),
				std::ios_base::out|std::ios_base::trunc):  
			std::cout.rdbuf());

  if (outputFile.rdbuf() != fb)
      {
	delete fb;
	std::cerr << " In " << __FILE__ << "(" << __LINE__ << "): Unable to open output file " <<  _PurityEfficiencyOutputFile << "!  Redirecting output to standard out." << std::endl;
	outputFile << std::endl;
      }

 for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection)
    {
      outputFile << "\n\nRESULTS for " << iTagCollection << "-th Flavour Tag Collection named " << _FlavourTagCollectionNames[iTagCollection] << std::endl;
      outputFile << "  any number of ZVTOP vertices found   " ;
      outputFile << "   N(b) = " ;
      outputFile.width(10);
      outputFile << _pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->sumAllBinHeights();
      outputFile << "   N(c) = " ;
      outputFile.width(10);
      outputFile << _pCJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->sumAllBinHeights();
      outputFile << "   N(light) = ";
      outputFile.width(10);
      outputFile << _pLightJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->sumAllBinHeights();
      outputFile << std::endl << std::endl;

       for (unsigned int iVertexCat=1;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	
	 std::string nvname = _VertexCatNames[iVertexCat];

	 

	 outputFile << "   ";
	 outputFile.width(2);
	 outputFile << iVertexCat;
	 outputFile << (
			(iVertexCat == 1) ? "           ZVTOP vertex   found   " : "           ZVTOP vertices found   ") ;
	 
	 outputFile << "   N(b) = " ;
	 outputFile.width(10);
	 outputFile << _pBJetBTagIntegral[iTagCollection][nvname]->sumAllBinHeights();
	 outputFile << "   N(c) = " ;
	 outputFile.width(10);
	 outputFile << _pCJetBTagIntegral[iTagCollection][nvname]->sumAllBinHeights();
	 outputFile << "   N(light) = ";
	 outputFile.width(10);
	 outputFile << _pLightJetBTagIntegral[iTagCollection][nvname]->sumAllBinHeights();
	 outputFile << std::endl << std::endl;
       }

      
  
      
      outputFile.precision(5);  
      outputFile.setf(std::ios::fixed,std::ios::floatfield);  

	  
      for (unsigned int iVertexCat=1;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
	
	std::string nvname = _VertexCatNames[iVertexCat];
	
	outputFile << "numbers of jets in cuts for ";
	outputFile.width(2);
	outputFile << iVertexCat;
	outputFile << " ZVTOP vertices found" << std::endl;
	outputFile << "cut    b-tag b    b-tag other    c-tag c   c-tag other    c-tagbb c   c-tagbb other" << std::endl;
	
	int numberOfBins=_pBJetBTagIntegral[iTagCollection][nvname]->axis().bins();
	
	for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ ) {
	  
	  outputFile.width(5);
	  outputFile.precision(3);
	  outputFile << _pBJetBTagIntegral[iTagCollection][nvname]->axis().binUpperEdge(binNumber) << "   ";
	  outputFile.width(10);
	  outputFile << int(_pBJetBTagIntegral[iTagCollection][nvname]->binHeight(binNumber))  << "   ";
	  outputFile.width(10);
	  outputFile << int(_pCJetBTagIntegral[iTagCollection][nvname]->binHeight(binNumber)+_pLightJetBTagIntegral[iTagCollection][nvname]->binHeight(binNumber))<< "   ";
	  outputFile.width(10);
	  outputFile << int(_pCJetCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))  << "   ";
	  outputFile.width(10);
	  outputFile << int(_pBJetCTagIntegral[iTagCollection][nvname]->binHeight(binNumber)+_pLightJetCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))<< "   ";
	  outputFile.width(10);
	  outputFile << int(_pCJetBCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))  << "   ";
	  outputFile.width(10);
	  outputFile << int(_pBJetBCTagIntegral[iTagCollection][nvname]->binHeight(binNumber)+_pLightJetBCTagIntegral[iTagCollection][nvname]->binHeight(binNumber))<< "   ";
	  outputFile << std::endl;

	}
	outputFile << std::endl;
      }
      
       outputFile << "numbers of jets in cuts summed (for any number of vertices)" << std::endl;
       outputFile << "cut    b-tag b    b-tag other    c-tag c   c-tag other    c-tagbb c   c-tagbb other" << std::endl;
  	
       int numberOfBins=_pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->axis().bins();
	
       for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ ) {
       
	 outputFile.width(5);
	 outputFile.precision(3);
	 outputFile << _pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->axis().binUpperEdge(binNumber) << "   ";
	 outputFile.width(10);
	 outputFile << int(_pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))  << "   ";
	 outputFile.width(10);
	 outputFile << int(_pCJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber)+_pLightJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))<< "   ";
	 outputFile.width(10);
	 outputFile << int(_pCJetCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))  << "   ";
	 outputFile.width(10);
	 outputFile << int(_pBJetCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber)+_pLightJetCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))<< "   ";
	 outputFile.width(10);
	 outputFile << int(_pCJetBCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))  << "   ";
	 outputFile.width(10);
	 outputFile << int(_pBJetBCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber)+_pLightJetBCTagIntegral[iTagCollection][_VertexCatNames[0]]->binHeight(binNumber))<< "   ";
	 outputFile << std::endl; 
       }      
       outputFile << std::endl;
       outputFile << std::endl;

    }

 for (unsigned int iTagCollection=0; iTagCollection < _FlavourTagCollectionNames.size(); ++iTagCollection)
   {
     
     int numberOfBins=_pBJetBTagIntegral[iTagCollection][_VertexCatNames[0]]->axis().bins();
     
     outputFile << "\n\nRESULTS for " << iTagCollection << "-th Flavour Tag Collection named " << _FlavourTagCollectionNames[iTagCollection] << std::endl;
     
     for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
       
       std::string nvname = _VertexCatNames[iVertexCat];
       
       std::vector<double> bJetBTagEfficiency;
       std::vector<double> bJetBTagEfficiencyError;  
       std::vector<double> cJetCTagEfficiency;
       std::vector<double> cJetCTagEfficiencyError;
       std::vector<double> cJetBCTagEfficiency;
       std::vector<double> cJetBCTagEfficiencyError;
       
       CalculateTagEfficiency( _pBJetBTag[iTagCollection][nvname], bJetBTagEfficiency, bJetBTagEfficiencyError);
       CalculateTagEfficiency( _pCJetCTag[iTagCollection][nvname], cJetCTagEfficiency, cJetCTagEfficiencyError);
       CalculateTagEfficiency( _pCJetBCTag[iTagCollection][nvname], cJetBCTagEfficiency, cJetBCTagEfficiencyError);
       
       std::vector<double> bJetBTagPurity;
       std::vector<double> bJetBTagPurityError;  
       std::vector<double> cJetCTagPurity;
       std::vector<double> cJetCTagPurityError;
       std::vector<double> cJetBCTagPurity;
       std::vector<double> cJetBCTagPurityError;
       
       CalculateTagPurity( _pBJetBTag[iTagCollection][nvname], _pBTagBackgroundValues[iTagCollection][nvname], bJetBTagPurity, bJetBTagPurityError);
       CalculateTagPurity( _pCJetCTag[iTagCollection][nvname], _pCTagBackgroundValues[iTagCollection][nvname], cJetCTagPurity, cJetCTagPurityError);
       CalculateTagPurity( _pCJetBCTag[iTagCollection][nvname], _pBCTagBackgroundValues[iTagCollection][nvname], cJetBCTagPurity, cJetBCTagPurityError);
       
       outputFile << "\n Purity-Efficiencies Values for          ";
       if (iVertexCat==0)  {
	 outputFile << " any number of ZVTOP vertices found " << std::endl;
       } else {
	 outputFile << iVertexCat;
	 outputFile << (
			(iVertexCat == 1) ? " ZVTOP vertex   found   " : "  ZVTOP vertices found   ") ;
       }
       outputFile << "\n-----------------------------------------------------------------------------------------------" << std::endl;
       outputFile << "cut value    efficiency(b)       purity(b)           efficiency(c)     purity(c)          efficiency(bc)     purity(bc) " << std::endl;
   
       for (int binNumber = 0; binNumber < numberOfBins ; binNumber++ ) {
	 outputFile.width(5);
	 outputFile.precision(3);
	 outputFile << _pBJetBTagIntegral[iTagCollection][nvname]->axis().binUpperEdge(binNumber) << "      ";
	 
	 outputFile << bJetBTagEfficiency[binNumber];
	 outputFile << " +/- ";
	 outputFile << bJetBTagEfficiencyError[binNumber];
	 outputFile << "   ";
	 outputFile << bJetBTagPurity[binNumber];
	 outputFile << " +/- ";
	 outputFile << bJetBTagPurityError[binNumber];
	 outputFile << "     ";
	 
	 outputFile << cJetCTagEfficiency[binNumber];
	 outputFile << " +/- ";
	 outputFile << cJetCTagEfficiencyError[binNumber];
	 outputFile << "   ";
	 outputFile << cJetCTagPurity[binNumber];
	 outputFile << " +/- ";
	 outputFile << cJetCTagPurityError[binNumber];
	 outputFile << "     ";
	 
	 outputFile << cJetBCTagEfficiency[binNumber];
	 outputFile << " +/- ";
	 outputFile << cJetBCTagEfficiencyError[binNumber];
	 outputFile << "   ";
	 outputFile << cJetBCTagPurity[binNumber];
	 outputFile << " +/- ";
	 outputFile << cJetBCTagPurityError[binNumber];
	 outputFile << "     ";
	 
	 outputFile << std::endl; 
       }
    }
     
     outputFile << std::endl;
     outputFile << std::endl;
     
   }
 
}

void LCFIAIDAPlotProcessor::InternalVectorInitialisation()
{
  //sets the size of some vectors and initialises some counters to zero

  _cJet_truePlus2_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus2_angle.resize(N_JETANGLE_BINS);
  
  _cJet_truePlus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_truePlus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus2_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_truePlus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_truePlus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueNeut_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_trueNeut_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueNeut_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_trueMinus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _cJet_trueMinus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _cJet_trueMinus2_recoMinus_angle.resize(N_JETANGLE_BINS);
  
  _bJet_truePlus2_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus_angle.resize(N_JETANGLE_BINS);	
  _bJet_trueNeut_angle.resize(N_JETANGLE_BINS);	
  _bJet_trueMinus_angle.resize(N_JETANGLE_BINS);	
  _bJet_trueMinus2_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_truePlus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus2_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_truePlus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_truePlus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_trueNeut_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_trueNeut_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_trueNeut_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_trueMinus_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus_recoMinus_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus2_recoPlus_angle.resize(N_JETANGLE_BINS); 
  _bJet_trueMinus2_recoNeut_angle.resize(N_JETANGLE_BINS);
  _bJet_trueMinus2_recoMinus_angle.resize(N_JETANGLE_BINS);
 

  for (int j = 0 ; j < N_JETANGLE_BINS ; j++) {
    _cJet_truePlus2_recoPlus_angle[j]=0;  
    _cJet_truePlus2_recoNeut_angle[j]=0;  
    _cJet_truePlus2_recoMinus_angle[j]=0; 
    _cJet_truePlus_recoPlus_angle[j]=0;   
    _cJet_truePlus_recoNeut_angle[j]=0;   
    _cJet_truePlus_recoMinus_angle[j]=0;  
    _cJet_trueNeut_recoPlus_angle[j]=0;   
    _cJet_trueNeut_recoNeut_angle[j]=0;   
    _cJet_trueNeut_recoMinus_angle[j]=0;  
    _cJet_trueMinus_recoPlus_angle[j]=0;  
    _cJet_trueMinus_recoNeut_angle[j]=0;  
    _cJet_trueMinus_recoMinus_angle[j]=0; 
    _cJet_trueMinus2_recoPlus_angle[j]=0; 
    _cJet_trueMinus2_recoNeut_angle[j]=0; 
    _cJet_trueMinus2_recoMinus_angle[j]=0;
    
    _bJet_truePlus2_angle[j]=0;	       
    _bJet_truePlus_angle[j]=0;	       
    _bJet_trueNeut_angle[j]=0;	       
    _bJet_trueMinus_angle[j]=0;	       
    _bJet_trueMinus2_angle[j]=0;	       
    _bJet_truePlus2_recoPlus_angle[j]=0;  
    _bJet_truePlus2_recoNeut_angle[j]=0;  
    _bJet_truePlus2_recoMinus_angle[j]=0; 
    _bJet_truePlus_recoPlus_angle[j]=0;   
    _bJet_truePlus_recoNeut_angle[j]=0;   
    _bJet_truePlus_recoMinus_angle[j]=0;  
    _bJet_trueNeut_recoPlus_angle[j]=0;   
    _bJet_trueNeut_recoNeut_angle[j]=0;   
    _bJet_trueNeut_recoMinus_angle[j]=0;  
    _bJet_trueMinus_recoPlus_angle[j]=0;  
    _bJet_trueMinus_recoNeut_angle[j]=0;  
    _bJet_trueMinus_recoMinus_angle[j]=0; 
    _bJet_trueMinus2_recoPlus_angle[j]=0; 
    _bJet_trueMinus2_recoNeut_angle[j]=0; 
    _bJet_trueMinus2_recoMinus_angle[j]=0;
  }

  
  _cJet_truePlus2=0;
  _cJet_truePlus=0;
  _cJet_trueNeut=0;
  _cJet_trueMinus=0;
  _cJet_trueMinus2=0;
  _cJet_truePlus2_recoPlus=0; 
  _cJet_truePlus2_recoNeut=0;
  _cJet_truePlus2_recoMinus=0;
  _cJet_truePlus_recoPlus=0; 
  _cJet_truePlus_recoNeut=0;
  _cJet_truePlus_recoMinus=0;
  _cJet_trueNeut_recoPlus=0; 
  _cJet_trueNeut_recoNeut=0;
  _cJet_trueNeut_recoMinus=0;
  _cJet_trueMinus_recoPlus=0; 
  _cJet_trueMinus_recoNeut=0;
  _cJet_trueMinus_recoMinus=0;
  _cJet_trueMinus2_recoPlus=0; 
  _cJet_trueMinus2_recoNeut=0;
  _cJet_trueMinus2_recoMinus=0;
  
  _bJet_truePlus2=0;
  _bJet_truePlus=0;	
  _bJet_trueNeut=0;	
  _bJet_trueMinus=0;	
  _bJet_trueMinus2=0;
  _bJet_truePlus2_recoPlus=0; 
  _bJet_truePlus2_recoNeut=0;
  _bJet_truePlus2_recoMinus=0;
  _bJet_truePlus_recoPlus=0; 
  _bJet_truePlus_recoNeut=0;
  _bJet_truePlus_recoMinus=0;
  _bJet_trueNeut_recoPlus=0; 
  _bJet_trueNeut_recoNeut=0;
  _bJet_trueNeut_recoMinus=0; 
  _bJet_trueMinus_recoPlus=0;  
  _bJet_trueMinus_recoNeut=0;  
  _bJet_trueMinus_recoMinus=0;
  _bJet_trueMinus2_recoPlus=0; 
  _bJet_trueMinus2_recoNeut=0;
  _bJet_trueMinus2_recoMinus=0;
        
  _nb_twoVertex_bTrack_Primary=0;	  
  _nb_twoVertex_bTrack_Secondary=0;  
  _nb_twoVertex_bTrack_Tertiary=0;
  _nb_twoVertex_bTrack_Isolated=0;
  _nb_twoVertex_cTrack_Primary=0; 	  
  _nb_twoVertex_cTrack_Secondary=0;	  
  _nb_twoVertex_cTrack_Tertiary=0; 	  
  _nb_twoVertex_cTrack_Isolated=0; 
  _nb_twoVertex_lTrack_Primary=0; 	  
  _nb_twoVertex_lTrack_Secondary=0;	  
  _nb_twoVertex_lTrack_Tertiary=0;	  
  _nb_twoVertex_lTrack_Isolated=0;
  _nb_threeVertex_bTrack_Primary=0;   
  _nb_threeVertex_bTrack_Secondary=0;  
  _nb_threeVertex_bTrack_Tertiary=0;    
  _nb_threeVertex_bTrack_Isolated=0;  
  _nb_threeVertex_cTrack_Primary=0;  
  _nb_threeVertex_cTrack_Secondary=0; 
  _nb_threeVertex_cTrack_Tertiary=0; 
  _nb_threeVertex_cTrack_Isolated=0; 
  _nb_threeVertex_lTrack_Primary=0;   
  _nb_threeVertex_lTrack_Secondary=0;  
  _nb_threeVertex_lTrack_Tertiary=0;   
  _nb_threeVertex_lTrack_Isolated=0;
 
  _nb_threeVertex_Primary_noMCP=0;
  _nb_threeVertex_Secondary_noMCP=0;
  _nb_threeVertex_Tertiary_noMCP=0;
  _nb_threeVertex_Isolated_noMCP=0;
  _nb_twoVertex_Primary_noMCP=0;
  _nb_twoVertex_Secondary_noMCP=0;
  _nb_twoVertex_Tertiary_noMCP=0;
  _nb_twoVertex_Isolated_noMCP=0;
  
  _nc_twoVertex_bTrack_Primary=0;	  
  _nc_twoVertex_bTrack_Secondary=0;  
  _nc_twoVertex_bTrack_Tertiary=0;
  _nc_twoVertex_bTrack_Isolated=0;
  _nc_twoVertex_cTrack_Primary=0; 	  
  _nc_twoVertex_cTrack_Secondary=0;	  
  _nc_twoVertex_cTrack_Tertiary=0; 	  
  _nc_twoVertex_cTrack_Isolated=0; 
  _nc_twoVertex_lTrack_Primary=0; 	  
  _nc_twoVertex_lTrack_Secondary=0;	  
  _nc_twoVertex_lTrack_Tertiary=0;	  
  _nc_twoVertex_lTrack_Isolated=0;
  _nc_threeVertex_bTrack_Primary=0;   
  _nc_threeVertex_bTrack_Secondary=0;  
  _nc_threeVertex_bTrack_Tertiary=0;    
  _nc_threeVertex_bTrack_Isolated=0;  
  _nc_threeVertex_cTrack_Primary=0;  
  _nc_threeVertex_cTrack_Secondary=0; 
  _nc_threeVertex_cTrack_Tertiary=0; 
  _nc_threeVertex_cTrack_Isolated=0; 
  _nc_threeVertex_lTrack_Primary=0;   
  _nc_threeVertex_lTrack_Secondary=0;  
  _nc_threeVertex_lTrack_Tertiary=0;   
  _nc_threeVertex_lTrack_Isolated=0;
 
  _nc_threeVertex_Primary_noMCP=0;
  _nc_threeVertex_Secondary_noMCP=0;
  _nc_threeVertex_Tertiary_noMCP=0;
  _nc_threeVertex_Isolated_noMCP=0;
  _nc_twoVertex_Primary_noMCP=0;
  _nc_twoVertex_Secondary_noMCP=0;
  _nc_twoVertex_Tertiary_noMCP=0;
  _nc_twoVertex_Isolated_noMCP=0;
  


}


void LCFIAIDAPlotProcessor::FindTrueJetDecayLength( LCEvent* pEvent, unsigned int jetNumber, std::vector<double>&  decaylengthvector, 
						    std::vector<double>&  bjetdecaylengthvector, std::vector<double>&  cjetdecaylengthvector)
{

  //looks in the MC to find the longest decay length for the b and c jets 
  //returns -99 for light jets

  int jettype = FindJetType( pEvent,  jetNumber );

  if (abs(jettype)!= C_JET && abs(jettype)!=B_JET) 
    {
      return;
    } else {
  
    int pdgcode =  FindJetPDGCode( pEvent,  jetNumber );
    
    TypesafeCollection<lcio::MCParticle> mcParticleCollection( pEvent, _MCParticleColName );  
    
    if (!mcParticleCollection.is_valid()) 
      { 
	std::cerr << "In " << __FILE__ << "(" << __LINE__ << "): for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << "mcParticleCollection is not valid" << std::endl;
	return;
      } 
    else 
      {
	for (int imc = 0 ; imc < mcParticleCollection.getNumberOfElements(); imc++ ) {
	  
	  lcio::MCParticle* myMCParticle = mcParticleCollection.getElementAt( imc );
	  
	  if (myMCParticle->getPDG() == pdgcode) {
	    
	    //v./ std::cout << "1: " << pdgcode << " " << GetPDGFlavour(pdgcode) << std::endl;

	    //get decay length
	    const double* vertex = myMCParticle->getVertex();
	    const double* endpoint = myMCParticle->getEndpoint();	 	  
	    double decay_length = CalculateDistance(vertex, endpoint);
	    
	    //if (decay_length>0) {

	      decaylengthvector.push_back(decay_length);
	      if (GetPDGFlavour(pdgcode) == C_JET ) cjetdecaylengthvector.push_back(decay_length);
	      if (GetPDGFlavour(pdgcode) == B_JET ) bjetdecaylengthvector.push_back(decay_length);
	      //}
	    //get the daughters
	    const MCParticleVec& daughters = myMCParticle->getDaughters();
	    
	    for ( MCParticleVec::const_iterator idaughter = daughters.begin(); idaughter!=daughters.end(); idaughter++)
	      {
		int code = (*idaughter)->getPDG();
		
		//std::cout << "PDG code of daughter 1: " << code << std::endl;
		//check if they are still bs or ds
		int flavour = GetPDGFlavour(code);
		//std::cout << "Daughter 1 flavour is: " << flavour << std::endl;
		if (abs(flavour) == C_JET || abs(flavour) == B_JET) {
		  
		  //v./  std::cout << "2: " << code << " " << flavour << std::endl;


		  //get decay length
		  const double* vertex = (*idaughter)->getVertex();
		  const double* endpoint = (*idaughter)->getEndpoint();	 	  
		  double decay_length = CalculateDistance(vertex, endpoint);
		
		  //if (decay_length>0) {
		    decaylengthvector.push_back(decay_length);  
		    if (abs(flavour) == C_JET ) cjetdecaylengthvector.push_back(decay_length);
		    if (abs(flavour) == B_JET ) bjetdecaylengthvector.push_back(decay_length);

		    
		    // }
		  
		  
		  // std::cout << "Decay length 2 is:     " << decay_length << std::endl;
		  
		  //look for daughters again - just in case...
		  
		  const MCParticleVec& daughters2 = (*idaughter)->getDaughters();
		  
		  for ( MCParticleVec::const_iterator idaughter2 = daughters2.begin(); idaughter2!=daughters2.end(); idaughter2++)
		    {
		      int code2 = (*idaughter2)->getPDG();
		      
		      //check if they are still bs or ds
		      int flavour2 = GetPDGFlavour(code2);
		      
		      if (abs(flavour2) == C_JET || abs (flavour2) == B_JET) {
			
			//v./ std::cout << "3: " << code2 << " " << flavour2 << std::endl;

			//get decay length
			const double* vertex = (*idaughter2)->getVertex();
			const double* endpoint = (*idaughter2)->getEndpoint();	 	  
			double decay_length = CalculateDistance(vertex, endpoint);

			//if (decay_length>0) {
			  decaylengthvector.push_back(decay_length);
			  if (abs(flavour2) == C_JET ) cjetdecaylengthvector.push_back(decay_length);
			  if (abs(flavour2) == B_JET ) bjetdecaylengthvector.push_back(decay_length);
			  //}

			const MCParticleVec& daughters3 = (*idaughter2)->getDaughters();

			for ( MCParticleVec::const_iterator idaughter3 = daughters3.begin(); idaughter3!=daughters3.end(); idaughter3++)
			  {
			    int code3 = (*idaughter3)->getPDG();
		      
			    //check if they are still bs or ds
			    int flavour3 = GetPDGFlavour(code3);
			    
			    if (abs(flavour3)==C_JET || abs(flavour3)==B_JET) {
			   
			      //v./ std::cout << "4: " << code3 << " " << flavour3 << std::endl;
   
			      //get decay length
			      const double* vertex = (*idaughter3)->getVertex();
			      const double* endpoint = (*idaughter3)->getEndpoint();	 	  
			      double decay_length = CalculateDistance(vertex, endpoint);
			
			      //if (decay_length>0) {
				decaylengthvector.push_back(decay_length);
				if (abs(flavour3) == C_JET ) cjetdecaylengthvector.push_back(decay_length);
				if (abs(flavour3) == B_JET ) bjetdecaylengthvector.push_back(decay_length);
				//}

			      const MCParticleVec& daughters4 = (*idaughter3)->getDaughters();

			      for ( MCParticleVec::const_iterator idaughter4 = daughters4.begin(); idaughter4!=daughters4.end(); idaughter4++)
				{
				  int code4 = (*idaughter4)->getPDG();
				  
				  //check if they are still bs or ds
				  int flavour4 = GetPDGFlavour(code4);
				  
				  if (abs(flavour4)==C_JET || abs(flavour4)==B_JET) {
				    
				    //v./ std::cout << "5: " << code4 << " " << flavour4 << "   " << std::endl;

				    //get decay length
				    const double* vertex = (*idaughter4)->getVertex();
				    const double* endpoint = (*idaughter4)->getEndpoint();	 	  
				    double decay_length = CalculateDistance(vertex, endpoint);
				    //if (decay_length>0) {
				      decaylengthvector.push_back(decay_length);
				      if (abs(flavour4) == C_JET ) cjetdecaylengthvector.push_back(decay_length);
				      if (abs(flavour4) == B_JET ) bjetdecaylengthvector.push_back(decay_length);
				      //}
				  }
				}
			      
			    }
			  }
			

		    
		      }//b or c
		      
		    }//loop over daughter's daugthers
		  
		}//daughter is b or c
	      }//loop over daughters
	  }
	}
      }
  }
  
  return;
}

void LCFIAIDAPlotProcessor::FindTrueJetDecayLength2( LCEvent* pEvent, unsigned int jetNumber, double& bjetdecaylength, double& cjetdecaylength)
{
  double b_decay_length0 = -1.;
  double c_decay_length0 = -1.;
  
  //looks in the MC to find the longest decay length for the b and c jets 
 
  int jettype = FindJetType( pEvent,  jetNumber );

  if (abs(jettype)!=B_JET  && abs(jettype)!=C_JET )
    {
      return;
    } else {
  
    int pdgcode =  FindJetPDGCode( pEvent,  jetNumber );
    
    TypesafeCollection<lcio::MCParticle> mcParticleCollection( pEvent, _MCParticleColName );  
    
    if (!mcParticleCollection.is_valid()) 
      { 
	std::cerr << __FILE__ << "(" << __LINE__ << "): MCParticle collection " <<  _MCParticleColName << " is not valid"	 
		  << " for event " << pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << std::endl;
	return;
      } 
    else 
      {
	for (int imc = 0 ; imc < mcParticleCollection.getNumberOfElements(); imc++ ) {
	  
	  lcio::MCParticle* myMCParticle = mcParticleCollection.getElementAt( imc );
	  if (myMCParticle->getPDG() == pdgcode) {
	    
	    
	    std::vector<double> BJetDecayLengthVector;
	    std::vector<double> CJetDecayLengthVector;
	    std::vector<int> BJetDecayVector;
	    std::vector<int> CJetDecayVector;
	    //std::cout << "1: " << pdgcode << " " << GetPDGFlavour(pdgcode) << std::endl;

	    int code0 = pdgcode;
	    int flavour0 = GetPDGFlavour(code0);
	    

	    //get start vertex
	    const double* vertex0 = myMCParticle->getVertex();
	    
	    //get the daughters
	    const MCParticleVec& daughters = myMCParticle->getDaughters();
	    
	    for ( MCParticleVec::const_iterator idaughter = daughters.begin(); idaughter!=daughters.end(); idaughter++)
	      {
		int code1 = (*idaughter)->getPDG();		
		int flavour1 = GetPDGFlavour(code1);
		
		if (abs(flavour1) != B_JET && abs(flavour0) == B_JET) {
		
		  const double* vertex1 = (*idaughter)->getVertex();
		  double decay_length = CalculateDistance(vertex0, vertex1);
		  BJetDecayLengthVector.push_back(decay_length);
		}

		if (abs(flavour1) != C_JET && abs(flavour0) == C_JET) {
		  const double* vertex1 = (*idaughter)->getVertex();
		  double decay_length = CalculateDistance(vertex0, vertex1);
		  CJetDecayLengthVector.push_back(decay_length);
		}

		//look for daughters again - just in case...
		
		const MCParticleVec& daughters2 = (*idaughter)->getDaughters();
		
		for ( MCParticleVec::const_iterator idaughter2 = daughters2.begin(); idaughter2!=daughters2.end(); idaughter2++)
		  {
		    int code2 = (*idaughter2)->getPDG();
		    int flavour2 = GetPDGFlavour(code2);
		    
		    if (abs(flavour2) != B_JET && abs(flavour1) == B_JET) {
		      const double* vertex2 = (*idaughter2)->getVertex();
		      double decay_length = CalculateDistance(vertex0, vertex2);
		      BJetDecayLengthVector.push_back(decay_length);
		    }
		    
		    if (abs(flavour2) != C_JET && abs(flavour1) == C_JET) {
		      const double* vertex2 = (*idaughter2)->getVertex();
		      double decay_length = CalculateDistance(vertex0, vertex2);
		      CJetDecayLengthVector.push_back(decay_length);
		    } 
		    
		    const MCParticleVec& daughters3 = (*idaughter2)->getDaughters();
		    
		    for ( MCParticleVec::const_iterator idaughter3 = daughters3.begin(); idaughter3!=daughters3.end(); idaughter3++)
		      {
			int code3 = (*idaughter3)->getPDG();
			int flavour3 = GetPDGFlavour(code3);
			
			if (abs(flavour3) != B_JET && abs(flavour2) == B_JET) {
			  const double* vertex3 = (*idaughter3)->getVertex();
			  double decay_length = CalculateDistance(vertex0, vertex3);
			  BJetDecayLengthVector.push_back(decay_length);
			}

			if (abs(flavour3) != C_JET && abs(flavour2) == C_JET) {
			  const double* vertex3 = (*idaughter3)->getVertex();
			  double decay_length = CalculateDistance(vertex0, vertex3);
			  CJetDecayLengthVector.push_back(decay_length);
			} 

			const MCParticleVec& daughters4 = (*idaughter3)->getDaughters();
			
			for ( MCParticleVec::const_iterator idaughter4 = daughters4.begin(); idaughter4!=daughters4.end(); idaughter4++)
			  {
			    int code4 = (*idaughter4)->getPDG();
			    int flavour4 = GetPDGFlavour(code4);
			    
			    if (abs(flavour4) != B_JET && abs(flavour3) == B_JET) {
			      const double* vertex4 = (*idaughter4)->getVertex();
			      double decay_length = CalculateDistance(vertex0, vertex4);
			      BJetDecayLengthVector.push_back(decay_length);
			    }
			     if (abs(flavour4) != C_JET && abs(flavour3) == C_JET) {
			       const double* vertex4 = (*idaughter4)->getVertex();
			       double decay_length = CalculateDistance(vertex0, vertex4);
			       CJetDecayLengthVector.push_back(decay_length);
			     }
			     
			     
			     const MCParticleVec& daughters5 = (*idaughter4)->getDaughters();
			    
			     for ( MCParticleVec::const_iterator idaughter5 = daughters5.begin(); idaughter5!=daughters5.end(); idaughter5++)
			       {
				 int code5 = (*idaughter5)->getPDG();
				 int flavour5 = GetPDGFlavour(code5);
				 
				 if (abs(flavour5) != B_JET && abs(flavour4) == B_JET) {
				   const double* vertex5 = (*idaughter5)->getVertex();
				   double decay_length = CalculateDistance(vertex0, vertex5);
				   BJetDecayLengthVector.push_back(decay_length);
				 }
				 if (abs(flavour5) != C_JET && abs(flavour4) == C_JET) {
				   const double* vertex5 = (*idaughter5)->getVertex();
				   double decay_length = CalculateDistance(vertex0, vertex5);
				   CJetDecayLengthVector.push_back(decay_length);
				 }
				 
				 //std::cout << code0 << " " << code1 << " " << code2 << " " << code3 << " " << code4 << " " << code5 << std::endl;

			       }
			  }
		      }
		  }
	      }
	    //most of the time we will only have one entry per vector
	    for (std::vector<double>::const_iterator iter = BJetDecayLengthVector.begin(); iter != BJetDecayLengthVector.end() ; iter++)
	      {
		double decay_length = *iter;
		if (decay_length > b_decay_length0 ) b_decay_length0 = decay_length;
	      }
	    
	     for (std::vector<double>::const_iterator iter = CJetDecayLengthVector.begin(); iter != CJetDecayLengthVector.end() ; iter++)
	      {
		double decay_length = *iter;
		if (decay_length > c_decay_length0 ) c_decay_length0 = decay_length;
	      }
	     
	  }
	}
      }
  }
 
  bjetdecaylength = b_decay_length0;
  cjetdecaylength = c_decay_length0;
  
 return;
}
  
  
int LCFIAIDAPlotProcessor::GetPDGFlavour(int code){
  //this little addition will take care of wierder mesons 

  code = abs(code);


  //only really want the weak decaying particles

  if (code==411) return C_JET;
  if (code==421) return C_JET;
  if (code==431) return C_JET;
  
  if (code==511) return B_JET;
  if (code==521) return B_JET;
  if (code==531) return B_JET;
  if (code==541) return B_JET;
  
  if (code==4122) return C_JET;
  if (code==4222) return C_JET;
  if (code==4112) return C_JET;
  if (code==4212) return C_JET;
  if (code==4232) return C_JET;
  if (code==4132) return C_JET;
  
    
  if (code==5122) return B_JET;
  if (code==5222) return B_JET;
  if (code==5112) return B_JET;
  if (code==5212) return B_JET;
  if (code==5132) return B_JET;
  if (code==5232) return B_JET;
  if (code==5332) return B_JET;
  
  return 1;
  

  if( code  > 10000   )
    {
      code  = code%1000;
    }
  
  if( code  > 1000   )
    {
      
      if( code/ 1000 == C_JET   )
	{
	  return C_JET;
	}
      
      if( code/ 1000 == B_JET   )
	{
	  return B_JET;
	}		
    }
  
  if( code  > 100   )
    {
      if( code / 100  == C_JET   )
	{
	  return C_JET;
	  
	}    
      
      if( code/100 == B_JET   )
	{
	  return B_JET;
	}
    }
  
  return 1;
}



#endif // endif of the check to see if MARLIN_USE_AIDA has been defined
