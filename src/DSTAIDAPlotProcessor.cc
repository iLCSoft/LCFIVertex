// If AIDA is being used then this will have been defined
#ifndef MARLIN_USE_AIDA

#warning "-------------------------------------------------------------------------------"
#warning "- DSTAIDAPlotProcessor requires MARLIN_USE_AIDA to be defined. Did you enable -"
#warning "- AIDA when compiling Marlin? DSTAIDAPlotProcessor will not be compiled.      -"
#warning "-------------------------------------------------------------------------------"

// Can't do anything else.
#else

#include "DSTAIDAPlotProcessor.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>


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

// #ifdef USING_RAIDA
// #pragma message "USING_RAIDA defined"
// #else
// #define USING_JAIDA
// #pragma message "USING_JAIDA defined"
// #endif


#ifdef USING_JAIDA//Data point sets aren't implemented in RAIDA - which is a shame as they have functionality not given by histograms
//such as the facility to set the error
#include <AIDA/IDataPointSet.h>
#include <AIDA/IDataPointSetFactory.h>
#include <AIDA/IDataPoint.h>
#endif



#include "EVENT/LCCollection.h"
#include "UTIL/PIDHandler.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "EVENT/LCParameters.h"
#include "EVENT/LCFloatVec.h"
#include "EVENT/LCIntVec.h"

#include "util/inc/vector3.h"
#include "util/inc/util.h"

using std::includes;
using std::map;
using std::string;
using std::endl;
using std::cout;
using std::vector;
using std::stringstream;
using std::abs;
using std::pair;
using std::ofstream;
using std::cerr;
using std::set;
using namespace lcio;

#include "TypesafeCollection.h"

//Needs to be instantiated for Marlin to know about it (I think)
DSTAIDAPlotProcessor aDSTAIDAPlotProcessor;

DSTAIDAPlotProcessor::DSTAIDAPlotProcessor() : marlin::Processor("DSTAIDAPlotProcessor")
{
	_description = "Plots various outputs from the flavour tag" ;

	// register steering parameters: name, description, class-variable, default value
	//The name of the collection of ReconstructedParticles that is the jet
	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" , 
				"Name of the collection of ReconstructedParticles that is the jet"  ,
				_JetCollectionName ,
				string("FTSelectedJets") ) ;	
	

}

DSTAIDAPlotProcessor::~DSTAIDAPlotProcessor()
{
}

void DSTAIDAPlotProcessor::init()
{
	printParameters();
	cout << _description << endl
		<< "-------------------------------------------------" << endl
		<< endl;
		
	_nRun=0;


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
      
       
      CreateTagPlots();

      CreateFlavourTagTuple();
      
    } else {

    std::cerr  << "### " << __FILE__ << "(" << __LINE__ << "): Unable to get the histogram factory! No histograms will be made."<< std::endl;
  }
  


}

void DSTAIDAPlotProcessor::processRunHeader( LCRunHeader* pRun )
{

  
  CreateFlavourTagInputPlots(pRun);
	
}

void DSTAIDAPlotProcessor::processEvent( LCEvent* pEvent )
{
  //check that what we need is here  
  //get the jet collection
  LCCollection* jetCollection;
  try{
    jetCollection=pEvent->getCollection( _JetCollectionName );
  }
  catch(DataNotAvailableException &e){
    std::cout << "Collection " <<_JetCollectionName << " is unavailable, can't plot anything " << std::endl;
    return;
  }

  PIDHandler pidh( jetCollection ) ;
  //get the algorithm id for the MC info
  try{
    pidh.getAlgorithmID("MCTruth");
  }
  catch(UnknownAlgorithm  &e)
    {
      std::cout << "MCTruth Algorithm is unavailable, can't plot anything " <<  std::endl;
      return; 
    }
  //get the algorithm id for the MC info
  try{
    pidh.getAlgorithmID("LCFIFlavourTag");
  }
  catch(UnknownAlgorithm  &e){
    std::cout << "LCFIFlavourTag Algorithm is unavailable, can't plot anything " << std::endl;
    return;
  }  
  
  
  //apply any cuts on the event here
  if( _passesEventCuts(pEvent) )
    {
      
      ReconstructedParticle* pJet;
      //loop over the jets
      for( int jetNumber=0; jetNumber<jetCollection->getNumberOfElements(); ++jetNumber )
	{
	  pJet=dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt(jetNumber) );
	  
	  //only do anything if the jet passes the jet cuts
	  if( _passesJetCuts(pJet) )
	    {
	      FillTagPlots( pEvent, jetNumber );
	      FillInputsPlots( pEvent, jetNumber );
	    }
	} 
    }



}




void DSTAIDAPlotProcessor::end()
{

  CalculateIntegralAndBackgroundPlots();
  CalculateEfficiencyPurityPlots();

}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently applies no cuts at all*/
bool DSTAIDAPlotProcessor::_passesEventCuts( LCEvent* )
{

	return true;
}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently selects jets for which the jet polar angle theta is -0.95< cos(theta) <0.95. 
*/
bool DSTAIDAPlotProcessor::_passesJetCuts( ReconstructedParticle* pJet )
{
	//
	// This cut added on the suggestion of Sonja Hillert 12/Jan/07.
	//
	// Selects jets for which the cosine of the jet polar
	// angle theta for all jets is not too large.
	//
	// Make something that's easy to search for to track down erroneous cuts:
	// GREPTAG_CUT : Jet cut on abs(cos(theta)) of jet axis
	//
	
	double CThJ_lower=0;	//lower cut on cos(theta) of the jet axis
	double CThJ_upper=0.95;	//upper cut (Sho Ryu Ken!)

	vertex_lcfi::util::Vector3 zAxis( 0, 0, 1 ); //work out theta from dot product with jet axis

	const double* mom=pJet->getMomentum();
	vertex_lcfi::util::Vector3 jetMomentum( mom[0], mom[1], mom[2] );

	jetMomentum.makeUnit();
	double cosTheta=jetMomentum.dot( zAxis );
	if( fabs(cosTheta)<=CThJ_lower || fabs(cosTheta)>=CThJ_upper ) return false;


	// If control gets to this point then the jet has passed
	return true;
}


void DSTAIDAPlotProcessor::_displayCollectionNames( LCEvent* pEvent )
{
	const vector<string>* pCollectionNames=pEvent->getCollectionNames();
	
	cout << "The available collections are: (name - type)" << endl;
	for( vector<string>::const_iterator i=pCollectionNames->begin(); i<pCollectionNames->end(); ++i )
	{
		LCCollection* pCollection=pEvent->getCollection( (*i) );
		const string typeName=pCollection->getTypeName();
		cout << "  " << (*i) << " - " << typeName << endl;
	}
	cout << endl;
}


void DSTAIDAPlotProcessor::CreateTagPlots()
{
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

 
 
  for (unsigned int iVertexCat=0;  iVertexCat <  N_VERTEX_CATEGORIES+1; ++iVertexCat ){
    _pLightJetBTag[_VertexCatNames[iVertexCat]]=0; 
    _pLightJetCTag[_VertexCatNames[iVertexCat]]=0;
    _pLightJetBCTag[_VertexCatNames[iVertexCat]]=0;
    
    _pBJetBTag[_VertexCatNames[iVertexCat]]=0;	 
    _pBJetCTag[_VertexCatNames[iVertexCat]]=0;
    _pBJetBCTag[_VertexCatNames[iVertexCat]]=0;
    
    _pCJetBTag[_VertexCatNames[iVertexCat]]=0;	 
    _pCJetCTag[_VertexCatNames[iVertexCat]]=0;  
    _pCJetBCTag[_VertexCatNames[iVertexCat]]=0;  
    
    _pBTagBackgroundValues[_VertexCatNames[iVertexCat]]=0; 
    _pCTagBackgroundValues[_VertexCatNames[iVertexCat]]=0; 
    _pBCTagBackgroundValues[_VertexCatNames[iVertexCat]]=0; 
     
  }

      
  

  for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ){
    
    std::string nvname = _VertexCatNames[iVertexCat];
    
    
    pTree->cd( "/" + name() + "/" + "/");
	    
    if (!pTree->cd( _NumVertexCatDir[iVertexCat]+"/")) {
      pTree->mkdir( _NumVertexCatDir[iVertexCat]+"/");
      pTree->cd( _NumVertexCatDir[iVertexCat]+"/");
    }	      
	    
    _pLightJetBTag[nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by B-tag NN value. ("+ nvname +")",_numberOfPoints , 0, 1.0 );
    _pLightJetCTag[nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
    _pLightJetBCTag[nvname] = pHistogramFactory->createHistogram1D( "Numbers of light jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
    _pBJetBTag[nvname]     = pHistogramFactory->createHistogram1D( "Numbers of B jets by B-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
    _pBJetCTag[nvname]     = pHistogramFactory->createHistogram1D( "Numbers of B jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
    _pBJetBCTag[nvname]    = pHistogramFactory->createHistogram1D( "Numbers of B jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
    _pCJetBTag[nvname]     = pHistogramFactory->createHistogram1D( "Numbers of C jets by B-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 ); 
    _pCJetCTag[nvname]     = pHistogramFactory->createHistogram1D( "Numbers of C jets by C-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );
    _pCJetBCTag[nvname]    = pHistogramFactory->createHistogram1D( "Numbers of C jets by BC-tag NN value. ("+ nvname +")", _numberOfPoints, 0, 1.0 );    
    
  }
	
}

void DSTAIDAPlotProcessor::CreateFlavourTagTuple()
{
#ifdef USING_JAIDA
  //AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

  //something dosen't work for me with the tuples in RAIDA
  AIDA::ITupleFactory* pTupleFactory=marlin::AIDAProcessor::tupleFactory( this );



    pTree->cd( "/" + name());
    
    //make the ntuple
    //this breaks the paradigm of reading these in from the flavour tag collections themselves
    std::string columnNames="int TrueJetFlavour=-1,  int NumberOfVertices=-1, int NumberOfTracksInVertices=-1,  float DecayLength = -999., float DecayLengthSignificance= -999., float JointProbRPhi= -999., float JointProbZ= -999., float PTCorrectedMass= -999., float RawMomentum= -999., float SecondaryVertexProbability= -999. ";

   	
    if (!pTree->cd( "/" + name() + "/" )) {
      pTree->cd( "/" + name() + "/"); 
      if (!pTree->cd("/")) {
	pTree->mkdir(  +"/");
	pTree->cd(  "/");
      }
    }
    
    if (!pTree->cd( "TupleDir/")) {
      pTree->mkdir( "TupleDir/");
      pTree->cd( "TupleDir/");
    }
    
    _pMyTuple=pTupleFactory->create( "FlavourTagInputsTuple","FlavourTagInputsTuple", columnNames);
 

#else
  //just in case this has side-effects..., previously the output was assigned to
  //a variable but not used if JAIDA was not used. Now we call this just in case
  //this function has side effects even though we do not use the return value
  marlin::AIDAProcessor::tree( this );
#endif	
 
 
    
}


void DSTAIDAPlotProcessor::CreateFlavourTagInputPlots(LCRunHeader* )
{

  
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  std::vector<std::string> VarNames;

  VarNames.push_back("NumVertices");
  VarNames.push_back("JointProbRPhi"); 
  VarNames.push_back("JointProbZ");
  VarNames.push_back("NumTracksInVertices");
  VarNames.push_back("DecayLength");
  VarNames.push_back("DecayLengthSignificance"); 
  VarNames.push_back("RawMomentum");
  VarNames.push_back("PTCorrectedMass");
  VarNames.push_back("SecondaryVertexProbability");

  for (size_t i = 0;i < VarNames.size();++i)
    {
      
      // If there is no histogram for this name then create one
      if( _inputsHistogramsBJets[VarNames[i]]==0 )
	{
	  pTree->cd( "/" + name() + "/");
	  if( !pTree->cd("NNInputs" ) )
	    {
	      pTree->mkdir( "NNInputs" ) ; 
	      pTree->cd("NNInputs" ) ;
	    }
	  
	  int numberOfPoints=_numberOfPoints/4;
	  double lowerBin=-1;
	  double higerBin=1;
	  
	  
	  if( VarNames[i]=="NumVertices" )
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
	  
	  
		
	  else if (VarNames[i]=="DecayLengthSignificance") 
	    {
	      numberOfPoints=100;
	      lowerBin=0.;
	      higerBin=100.;
	    }
	  else if (VarNames[i]=="DecayLength")
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
	  else if (VarNames[i]=="RawMomentum" ) 
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
		    
	  pTree->cd( "/" + name() + "/" + "NNInputs");
	  if( !pTree->cd( "bJets" ) )
	    {
	      pTree->mkdir( "bJets" ) ; 
	      pTree->cd(    "bJets" ) ;
	    }
		    
	  _inputsHistogramsBJets[VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
	  
	  pTree->cd( "/" + name() + "/" + "NNInputs");
	  if( !pTree->cd( "cJets" ) )
	    {
	      pTree->mkdir( "cJets" ) ; 
	      pTree->cd(    "cJets" ) ;
	    }
	  
	  _inputsHistogramsCJets[VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		    
	  pTree->cd( "/" + name() + "/" + "NNInputs");
	  if( !pTree->cd( "udsJets/" ) )
	    {
	      pTree->mkdir( "udsJets/" ) ; 
	      pTree->cd(    "udsJets/" ) ;
	    }
	  
	  _inputsHistogramsUDSJets[VarNames[i]]=pHistogramFactory->createHistogram1D( VarNames[i], numberOfPoints, lowerBin, higerBin );
		  
	}//end of histogram creation
    }
	    

}

 



void DSTAIDAPlotProcessor::FillTagPlots( LCEvent* pEvent, unsigned int jetNumber)
{
 //  int jetType=FindTrueJetType( pEvent, jetNumber );
//   if( jetType==0 ) return;
  
 //get the jet collection
  LCCollection* jetCollection=pEvent->getCollection( _JetCollectionName );
  //  TypesafeCollection<lcio::ReconstructedParticle> *jetCollection( pEvent, _JetCollectionName );

  //get the jet
  ReconstructedParticle* pJet=dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt(jetNumber) );

  //get the PIDHandler fot this jet collection
  PIDHandler pidh( jetCollection ) ;

  
  //get the algorithm id for the MC info
       
  int mcid = pidh.getAlgorithmID("MCTruth");

  //get the Particle id object containing the MC info
  const ParticleID& mcpid = pidh.getParticleID(pJet,mcid);
  
  //get the paramters for the MC info
  FloatVec mcparams = mcpid.getParameters();
  
  //get the true jet flavour
  int jetType = int(mcparams[pidh.getParameterIndex(mcid,"TrueJetFlavour")]);

  //get the algorithm id for the FT info
  int alid = pidh.getAlgorithmID("LCFIFlavourTag");
  
  //get the Particle id object containing the FT info
  const ParticleID& pid = pidh.getParticleID(pJet,alid);

  //get the paramters for the FT info
  FloatVec params = pid.getParameters();
  
  //get the actual flavour tags
  double bTag= params[pidh.getParameterIndex(alid,"BTag")];
  double cTag= params[pidh.getParameterIndex(alid,"CTag")];
  double cTagBBack= params[pidh.getParameterIndex(alid,"BCTag")];
  
  unsigned int NumVertices = int(params[pidh.getParameterIndex(alid,"NumVertices")]);
  //int CQVtx =  FindCQVtx(pEvent, jetNumber);
  //int BQVtx =  FindBQVtx(pEvent, jetNumber);
  //int trueJetCharge = int(FindTrueJetHadronCharge(pEvent,jetNumber));
  
  std::string nvname = _VertexCatNames[ (NumVertices>=N_VERTEX_CATEGORIES) ? (N_VERTEX_CATEGORIES) : (NumVertices)];


	      
  if( jetType==B_JET )  {
    
    if( bTag<=1 && bTag>=0 )
      {
	_pBJetBTag[nvname]->fill( bTag );
      } 
    else 
      {
	_pBJetBTag[nvname]->fill( -0.005 );
      }
    
    if( cTag<=1 && cTag>=0 )
      {
	_pBJetCTag[nvname]->fill( cTag );
      }
    else 
      {
	_pBJetCTag[nvname]->fill( -0.005 );
      }
    if( cTagBBack<=1 && cTagBBack>=0 ) 
      {
	_pBJetBCTag[nvname]->fill( cTagBBack );
      }
    else 
      {
	_pBJetBCTag[nvname]->fill( -0.005 );
      }
    
  } else if( jetType==C_JET ) {
    
    if( bTag<=1 && bTag>=0 )
      {
	_pCJetBTag[nvname]->fill( bTag );
      } 
    else 
      {
	_pCJetBTag[nvname]->fill( -0.005 );
      }
    
    if( cTag<=1 && cTag>=0 )
      {
	_pCJetCTag[nvname]->fill( cTag );
      }
    else 
      {
	_pCJetCTag[nvname]->fill( -0.005 );
      }
    
    if( cTagBBack<=1 && cTagBBack>=0 ) 
      {
	_pCJetBCTag[nvname]->fill( cTagBBack );
      }
    else 
      {
	_pCJetBCTag[nvname]->fill( -0.005 );
      }
  } else {
    if( bTag<=1 && bTag>=0 )
      {
	_pLightJetBTag[nvname]->fill( bTag );
      } 
    else 
      {
	_pLightJetBTag[nvname]->fill( -0.005 );
      }
    
    if( cTag<=1 && cTag>=0 )
      {
	_pLightJetCTag[nvname]->fill( cTag );
      }
    else 
      {
	_pLightJetCTag[nvname]->fill( -0.005 );
      } 
    
    if( cTagBBack<=1 && cTagBBack>=0 ) 
      {
	_pLightJetBCTag[nvname]->fill( cTagBBack );
      }
    else 
      {
	_pLightJetBCTag[nvname]->fill( -0.005 );
      }
  }
	
}

void DSTAIDAPlotProcessor::FillInputsPlots( LCEvent* pEvent, unsigned int jetNumber )
{  


  
//jet the jet collection
  LCCollection* jetCollection=pEvent->getCollection( _JetCollectionName );
  //  TypesafeCollection<lcio::ReconstructedParticle> *jetCollection( pEvent, _JetCollectionName );

  //get the jet
  ReconstructedParticle* pJet=dynamic_cast<ReconstructedParticle*>( jetCollection->getElementAt(jetNumber) );

  //get the PIDHandler fot this jet collection
  PIDHandler pidh( jetCollection ) ;

  
  //get the algorithm id for the MC info
  int mcid = pidh.getAlgorithmID("MCTruth");

  //get the Particle id object containing the MC info
  const ParticleID& mcpid = pidh.getParticleID(pJet,mcid);
  
  //get the paramters for the MC info
  FloatVec mcparams = mcpid.getParameters();
  
  //get the true jet flavour
  int jetType = int(mcparams[pidh.getParameterIndex(mcid,"TrueJetFlavour")]);
  

  int alid = pidh.getAlgorithmID("LCFIFlavourTag");
  
  //get the Particle id object containing the FT info
  const ParticleID& pid = pidh.getParticleID(pJet,alid);

  //get the paramters for the FT info
  FloatVec params = pid.getParameters();
  
  //get the actual flavour tags

  
  int NumVertices = int(params[pidh.getParameterIndex(alid,"NumVertices")]);

  
 
  float JointProbRPhi=-999.;
  float JointProbZ=-999.;
  int NumTracksInVertices=-999;
  float DecayLength=-999.;
  float DecayLengthSignificance=-999.;
  float RawMomentum=-999.;
  float PTCorrectedMass=-999.;
  float SecondaryVertexProbability=-999.;

  if(NumVertices > 1)
    {
      JointProbRPhi=	    params[pidh.getParameterIndex(alid,"JointProbRPhi")];
      JointProbZ=	       params[pidh.getParameterIndex(alid,"JointProbZ")];
      NumTracksInVertices= int(params[pidh.getParameterIndex(alid,"NumTracksInVertices")]);
      DecayLength=params[pidh.getParameterIndex(alid,"DecayLength")];
      DecayLengthSignificance=  params[pidh.getParameterIndex(alid,"DecayLengthSignificance")];
      RawMomentum=	       params[pidh.getParameterIndex(alid,"RawMomentum")];
      PTCorrectedMass=    params[pidh.getParameterIndex(alid,"PTCorrectedMass")];
      SecondaryVertexProbability=params[pidh.getParameterIndex(alid,"SecondaryVertexProbability")];
    }
		
  if (_pMyTuple) {
    
    _pMyTuple->fill( 0, jetType );
    _pMyTuple->fill( 1, NumVertices );
    if(NumTracksInVertices > 1)
      {
	_pMyTuple->fill( 2, NumTracksInVertices );
	_pMyTuple->fill( 3, DecayLength);
	_pMyTuple->fill( 4, DecayLengthSignificance);
	_pMyTuple->fill( 5, JointProbRPhi);
	_pMyTuple->fill( 6, JointProbZ);
	_pMyTuple->fill( 7, PTCorrectedMass);
	_pMyTuple->fill( 8, RawMomentum);
	_pMyTuple->fill( 9, SecondaryVertexProbability);
      }		  
    _pMyTuple->addRow();
    
  }//endif _pMyTuple


 if( jetType==B_JET ) _inputsHistogramsBJets["NumVertices"]->fill(NumVertices);
  else if( jetType==C_JET ) _inputsHistogramsCJets["NumVertices"]->fill(NumVertices);
  else _inputsHistogramsUDSJets["NumVertices"]->fill(NumVertices);
 

  if(NumVertices > 1)
    {
      if( jetType==B_JET ) _inputsHistogramsBJets["NumTracksInVertices"]->fill(NumTracksInVertices);
      else if( jetType==C_JET ) _inputsHistogramsCJets["NumTracksInVertices"]->fill(NumTracksInVertices);
      else _inputsHistogramsUDSJets["NumTracksInVertices"]->fill(NumTracksInVertices);
      
      if( jetType==B_JET ) _inputsHistogramsBJets["DecayLength"]->fill(DecayLength);
      else if( jetType==C_JET ) _inputsHistogramsCJets["DecayLength"]->fill(DecayLength);
      else _inputsHistogramsUDSJets["DecayLength"]->fill(DecayLength);

      if( jetType==B_JET ) _inputsHistogramsBJets["DecayLengthSignificance"]->fill(DecayLengthSignificance);
      else if( jetType==C_JET ) _inputsHistogramsCJets["DecayLengthSignificance"]->fill(DecayLengthSignificance);
      else _inputsHistogramsUDSJets["DecayLengthSignificance"]->fill(DecayLengthSignificance);

      if( jetType==B_JET ) _inputsHistogramsBJets["JointProbRPhi"]->fill(JointProbRPhi);
      else if( jetType==C_JET ) _inputsHistogramsCJets["JointProbRPhi"]->fill(JointProbRPhi);
      else _inputsHistogramsUDSJets["JointProbRPhi"]->fill(JointProbRPhi);

      if( jetType==B_JET ) _inputsHistogramsBJets["JointProbZ"]->fill(JointProbZ);
      else if( jetType==C_JET ) _inputsHistogramsCJets["JointProbZ"]->fill(JointProbZ);
      else _inputsHistogramsUDSJets["JointProbZ"]->fill(JointProbZ);
 
      if( jetType==B_JET ) _inputsHistogramsBJets["PTCorrectedMass"]->fill(PTCorrectedMass);
      else if( jetType==C_JET ) _inputsHistogramsCJets["PTCorrectedMass"]->fill(PTCorrectedMass);
      else _inputsHistogramsUDSJets["PTCorrectedMass"]->fill(PTCorrectedMass);
      
      if( jetType==B_JET ) _inputsHistogramsBJets["RawMomentum"]->fill(RawMomentum);
      else if( jetType==C_JET ) _inputsHistogramsCJets["RawMomentum"]->fill(RawMomentum);
      else _inputsHistogramsUDSJets["RawMomentum"]->fill(RawMomentum);

      if( jetType==B_JET ) _inputsHistogramsBJets["SecondaryVertexProbability"]->fill(SecondaryVertexProbability);
      else if( jetType==C_JET ) _inputsHistogramsCJets["SecondaryVertexProbability"]->fill(SecondaryVertexProbability);
      else _inputsHistogramsUDSJets["SecondaryVertexProbability"]->fill(SecondaryVertexProbability);


    }

		
  
	      
}



void DSTAIDAPlotProcessor::CalculateIntegralAndBackgroundPlots() {
  
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );
  
      
  for (unsigned int iVertexCat=1;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
    //sum over the different vertex catagories, this information goes into the "AnyNumberOfVertices" directory
    
    _pBJetBTag[_VertexCatNames[0]] -> 	add(*_pBJetBTag[_VertexCatNames[iVertexCat]]);
    _pBJetCTag[_VertexCatNames[0]] -> 	add(*_pBJetCTag[_VertexCatNames[iVertexCat]]);
    _pBJetBCTag[_VertexCatNames[0]] -> 	add(*_pBJetBCTag[_VertexCatNames[iVertexCat]]);
    _pCJetBTag[_VertexCatNames[0]] -> 	add(*_pCJetBTag[_VertexCatNames[iVertexCat]]);
    _pCJetCTag[_VertexCatNames[0]] -> 	add(*_pCJetCTag[_VertexCatNames[iVertexCat]]);
    _pCJetBCTag[_VertexCatNames[0]] -> 	add(*_pCJetBCTag[_VertexCatNames[iVertexCat]]);
    _pLightJetBTag[_VertexCatNames[0]] -> 	add(*_pLightJetBTag[_VertexCatNames[iVertexCat]]);
    _pLightJetCTag[_VertexCatNames[0]] -> 	add(*_pLightJetCTag[_VertexCatNames[iVertexCat]]);
	_pLightJetBCTag[_VertexCatNames[0]] -> add(*_pLightJetBCTag[_VertexCatNames[iVertexCat]]);
  }

  for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {
    //add up all the background values
    
    pTree->cd( "/" + name() + "/" +_NumVertexCatDir[iVertexCat]);
    
    _pBTagBackgroundValues[_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-B jets by B-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetBTag[_VertexCatNames[iVertexCat]],*_pCJetBTag[_VertexCatNames[iVertexCat]]);
    _pCTagBackgroundValues[_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-C jets by C-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetCTag[_VertexCatNames[iVertexCat]],*_pBJetCTag[_VertexCatNames[iVertexCat]]); 
	_pBCTagBackgroundValues[_VertexCatNames[iVertexCat]] = pHistogramFactory->add("Numbers of non-C jets by BC-tag NN value.  ("+ _VertexCatNames[iVertexCat]+")",*_pLightJetBCTag[_VertexCatNames[iVertexCat]],*_pBJetBCTag[_VertexCatNames[iVertexCat]]);
  }
  
  
}

void DSTAIDAPlotProcessor::CalculateEfficiencyPurityPlots()
{
  AIDA::IHistogramFactory* pHistogramFactory=marlin::AIDAProcessor::histogramFactory( this );
  AIDA::ITree* pTree=marlin::AIDAProcessor::tree( this );

#ifdef USING_JAIDA
  AIDA::IDataPointSetFactory* pDataPointSetFactory=marlin::AIDAProcessor::dataPointSetFactory(this);
#endif


  for (unsigned int iVertexCat=0;  iVertexCat <=  N_VERTEX_CATEGORIES; ++iVertexCat ) {

    pTree->cd( "/" + name() + "/" +_NumVertexCatDir[iVertexCat] );
    
    std::string nvname = _VertexCatNames[iVertexCat];
	
#ifdef USING_JAIDA
    AIDA::IDataPointSet* _pBJetBTagEfficiency = CreateEfficiencyPlot( _pBJetBTag[nvname] , pDataPointSetFactory->create("B-Tag efficiency  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pCJetCTagEfficiency = CreateEfficiencyPlot( _pCJetCTag[nvname] , pDataPointSetFactory->create("C-Tag efficiency  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pCJetBCTagEfficiency = CreateEfficiencyPlot( _pCJetBCTag[nvname] , pDataPointSetFactory->create("BC-Tag efficiency  ("+ nvname +")",2));
#endif

	
    _pBJetBTagIntegral[nvname] =   	
      CreateIntegralHistogram( _pBJetBTag[nvname], 
			       pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pBJetBTag[nvname]->axis().bins(),_pBJetBTag[nvname]->axis().lowerEdge(),_pBJetBTag[nvname]->axis().upperEdge()));
    
    _pCJetBTagIntegral[nvname] =     
      CreateIntegralHistogram( _pCJetBTag[nvname], 
			       pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pCJetBTag[nvname]->axis().bins(),_pCJetBTag[nvname]->axis().lowerEdge(),_pCJetBTag[nvname]->axis().upperEdge()));
    
    _pLightJetBTagIntegral[nvname] = 
      CreateIntegralHistogram( _pLightJetBTag[nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing B-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pLightJetBTag[nvname]->axis().bins(),_pLightJetBTag[nvname]->axis().lowerEdge(),_pLightJetBTag[nvname]->axis().upperEdge()));
    
    _pBJetCTagIntegral[nvname] = 
      CreateIntegralHistogram( _pBJetCTag[nvname], 
			       pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pBJetCTag[nvname]->axis().bins(),_pBJetCTag[nvname]->axis().lowerEdge(),_pBJetCTag[nvname]->axis().upperEdge()));
    
    _pCJetCTagIntegral[nvname] =     
      CreateIntegralHistogram( _pCJetCTag[nvname], 
			       pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pCJetCTag[nvname]->axis().bins(),_pCJetCTag[nvname]->axis().lowerEdge(),_pCJetCTag[nvname]->axis().upperEdge()));
    
    _pLightJetCTagIntegral[nvname] = 
      CreateIntegralHistogram( _pLightJetCTag[nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing C-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pLightJetCTag[nvname]->axis().bins(),_pLightJetCTag[nvname]->axis().lowerEdge(),_pLightJetCTag[nvname]->axis().upperEdge()));
    
    _pBJetBCTagIntegral[nvname] =    
      CreateIntegralHistogram( _pBJetBCTag[nvname], 
				   pHistogramFactory->createHistogram1D("B-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
									_pBJetBCTag[nvname]->axis().bins(),_pBJetBCTag[nvname]->axis().lowerEdge(),_pBJetBCTag[nvname]->axis().upperEdge()));
    
    _pCJetBCTagIntegral[nvname] =   
      CreateIntegralHistogram( _pCJetBCTag[nvname], 
			       pHistogramFactory->createHistogram1D("C-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +") (DON'T TRUST ERRORS!)",
								    _pCJetBCTag[nvname]->axis().bins(),_pCJetBCTag[nvname]->axis().lowerEdge(),_pCJetBCTag[nvname]->axis().upperEdge()));
    
    _pLightJetBCTagIntegral[nvname] = 
      CreateIntegralHistogram( _pLightJetBCTag[nvname], 
			       pHistogramFactory->createHistogram1D("Light-Jets: Numbers of events passing BC-Tag NN Cut  ("+ nvname +")",
								    _pLightJetBCTag[nvname]->axis().bins(),_pLightJetBCTag[nvname]->axis().lowerEdge(),_pLightJetBCTag[nvname]->axis().upperEdge()));
    
    
	
#ifdef USING_JAIDA
    AIDA::IDataPointSet* _pBJetBTagPurity =  CreatePurityPlot( _pBJetBTag[nvname],  _pBTagBackgroundValues[nvname] , pDataPointSetFactory->create("B-Jet purity for B-Tag  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pCJetCTagPurity =  CreatePurityPlot( _pCJetCTag[nvname],  _pCTagBackgroundValues[nvname] , pDataPointSetFactory->create("C-Jet purity for C-Tag  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pCJetBCTagPurity = CreatePurityPlot( _pCJetBCTag[nvname], _pBJetBCTag[nvname], pDataPointSetFactory->create("C-Jet purity for BC-Tag  ("+ nvname +")",2));      
    
    AIDA::IDataPointSet* _pCJetBTagLeakage =      CreateLeakageRatePlot( _pCJetBTag[nvname],      pDataPointSetFactory->create("C-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pLightJetBTagLeakage =  CreateLeakageRatePlot( _pLightJetBTag[nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pBJetCTagLeakage =      CreateLeakageRatePlot( _pBJetCTag[nvname],      pDataPointSetFactory->create("B-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pLightJetCTagLeakage =  CreateLeakageRatePlot( _pLightJetCTag[nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pBJetBCTagLeakage =     CreateLeakageRatePlot( _pBJetBCTag[nvname],     pDataPointSetFactory->create("B-Jets: Leakage Rate into BC-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pLightJetBCTagLeakage = CreateLeakageRatePlot( _pLightJetBCTag[nvname], pDataPointSetFactory->create("Light-Jets: Leakage Rate into BC-Tag Sample  ("+ nvname +")",2));     
    AIDA::IDataPointSet* _pNonBJetBTagLeakage =   CreateLeakageRatePlot( _pBTagBackgroundValues[nvname],      pDataPointSetFactory->create("C-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pNonCJetCTagLeakage =   CreateLeakageRatePlot( _pCTagBackgroundValues[nvname],  pDataPointSetFactory->create("Light-Jets: Leakage Rate into B-Tag Sample  ("+ nvname +")",2));
    AIDA::IDataPointSet* _pNonCJetBCTagLeakage =  CreateLeakageRatePlot( _pBCTagBackgroundValues[nvname],      pDataPointSetFactory->create("B-Jets: Leakage Rate into C-Tag Sample  ("+ nvname +")",2));
    
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
  
  

}


AIDA::IHistogram1D* DSTAIDAPlotProcessor::CreateIntegralHistogram(const AIDA::IHistogram1D* pNN, AIDA::IHistogram1D* pIntegral)
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



#ifdef USING_JAIDA
AIDA::IDataPointSet* DSTAIDAPlotProcessor::CreateEfficiencyPlot(const AIDA::IHistogram1D* pSignal, AIDA::IDataPointSet* pDataPointSet)
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

AIDA::IDataPointSet* DSTAIDAPlotProcessor::CreatePurityPlot(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet)  
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
AIDA::IDataPointSet* DSTAIDAPlotProcessor::CreateLeakageRatePlot(const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet)
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


AIDA::IDataPointSet* DSTAIDAPlotProcessor::CreateXYPlot(const AIDA::IDataPointSet* pDataPointSet0, const AIDA::IDataPointSet* pDataPointSet1, AIDA::IDataPointSet* xyPointSet, const int dim0, const int dim1 )
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

#endif


#endif // end of "#ifndef MARLIN_USE_AIDA"

