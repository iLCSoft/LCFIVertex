
#include "DSTPlotProcessor.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <set>
#ifdef USEROOT
/////////////////////
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
/////////////////////
#endif // of #ifdef USEROOT

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

//Needs to be instantiated for Marlin to know about it (I think)
DSTPlotProcessor aDSTPlotProcessor;

DSTPlotProcessor::DSTPlotProcessor() : marlin::Processor("DSTPlotProcessor")
{
	_description = "Plots various outputs from the flavour tag" ;

	// register steering parameters: name, description, class-variable, default value
	//The name of the collection of ReconstructedParticles that is the jet
	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"JetCollectionName" , 
				"Name of the collection of ReconstructedParticles that is the jet"  ,
				_JetCollectionName ,
				string("FTSelectedJets") ) ;

	//The output filename
	registerProcessorParameter( "OutputFilename" , 
				"Filename for the output"  ,
				_OutputFilename ,
				string("DSTPlotProcessorOutput") ) ;

	registerProcessorParameter( "CheckDSTParameters" , 
				"Checks the DST parameters against the full parameters - gives a lot of print out"  ,
				_checkDST ,
				int(0) ) ;

	

	//Add these collections id you want to check the DST collections against the full collections
	registerOptionalParameter("FlavourTagCollection" , 
				 "Name of the LCFloatVec Collection containing the flavour tags, only need this if you want to check DST parameters"  ,
				 _FlavourTagCollectionName,
				 std::string("FlavourTag") ) ;
	
	registerOptionalParameter("FlavourTagInputsCollection" , 
				 "Name of the LCFloatVec Collection that contains the flavour tag inputs, only need this if you want to check DST parameters"  ,
				 _FlavourTagInputsCollectionName,
				 std::string("FlavourTagInputs") ) ;
	

	
	registerOptionalParameter("TrueJetFlavourCollection" , 
				 "Name of the LCIntVec collection containing the true flavour of the jets, only need this if you want to check DST parameters"  ,
				 _TrueJetFlavourColName ,
				 std::string("TrueJetFlavour") ) ;
	

}

DSTPlotProcessor::~DSTPlotProcessor()
{
}

void DSTPlotProcessor::init()
{
	printParameters();
	cout << _description << endl
		<< "-------------------------------------------------" << endl
		<< endl;
		
	_nRun=0;
	_jetEMax=0;

}

void DSTPlotProcessor::processRunHeader( LCRunHeader* pRun )
{
	//

  //map the parameter names to the indexes, for the flavout tags, NN input vales, and mc truth info, if the DST parameters are to be checked
  if(_checkDST==1)
    {
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
    }

	
}

void DSTPlotProcessor::processEvent( LCEvent* pEvent )
{
  
  if(_checkDST==1)
    _checkDSTParameters( pEvent);

	//Get the collection of jets. Can't do anything if the collection isn't there
	//so don't bother catching the exception and terminate.
	LCCollection* pJetCollection=pEvent->getCollection( _JetCollectionName );
	
	//make sure the collection is of the right type
	if( pJetCollection->getTypeName()!=LCIO::RECONSTRUCTEDPARTICLE )
	{
		stringstream message;
		message << endl
			<< "########################################################################################\n"
			<< "# DSTPlotProcessor -                                                                      #\n"
			<< "#   The jet collection requested (\"" << _JetCollectionName << "\") is not of the type \"" << LCIO::RECONSTRUCTEDPARTICLE << "\"  #\n"
			<< "########################################################################################" << endl;
		throw EventException( message.str() );
	}
	
	//apply any cuts on the event here
	if( _passesEventCuts(pEvent) )
	{
		ReconstructedParticle* pJet;
		//loop over the jets
		for( int a=0; a<pJetCollection->getNumberOfElements(); ++a )
		{
			//Dynamic casts are not the best programming practice in the world, but I can't see another way of doing this
			//in the LCIO framework.  This cast should be safe though because we've already tested the type.
			pJet=dynamic_cast<ReconstructedParticle*>( pJetCollection->getElementAt(a) );
			
			if( _passesJetCuts(pJet) )
			{
				_fillPlots(pEvent, a );
				_jetEnergy.add(pJet->getEnergy());
				if (pJet->getEnergy() > _jetEMax) _jetEMax = pJet->getEnergy();
			}
		}
	}
}


void DSTPlotProcessor::end()
{
	_outputDataToFile( _OutputFilename );

	cout << "The largest jet energy was " << _jetEMax << endl;
}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently applies no cuts at all*/
bool DSTPlotProcessor::_passesEventCuts( LCEvent* pEvent )
{
	//
	// No event cuts at present
	//

	return true;
}

// IMPORTANT - If you change the cuts make sure you change the line below to show the changes in the docs
/*! Currently selects jets for which the jet polar angle theta is -0.95< cos(theta) <0.95. 
*/
bool DSTPlotProcessor::_passesJetCuts( ReconstructedParticle* pJet )
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

void DSTPlotProcessor::_fillPlots( LCEvent* pEvent, unsigned int jet)
{

  //jet the jet collection
  LCCollection* pJetCollection=pEvent->getCollection( _JetCollectionName );

  //get the jet
  ReconstructedParticle* pJet=dynamic_cast<ReconstructedParticle*>( pJetCollection->getElementAt(jet) );
  
  //get the PIDHandler fot this jet collection
  PIDHandler pidh( pJetCollection ) ;

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
  
  //get the algorithm id for the MC info
  int mcid = pidh.getAlgorithmID("MCTruth");

  //get the Particle id object containing the MC info
  const ParticleID& mcpid = pidh.getParticleID(pJet,mcid);
  
  //get the paramters for the MC info
  FloatVec mcparams = mcpid.getParameters();
  
  //get the true jet flavour
  float jetType = mcparams[pidh.getParameterIndex(mcid,"TrueJetFlavour")];

       
  //cout<<"JetType = "<<jetType<<endl;
  
  if( jetType==B_JET )
    {
      if( bTag<=1 && bTag>=0 ) _BTagEfficiencyPurity.add_signal( bTag );
      if( cTag<=1 && cTag>=0 ) _CTagEfficiencyPurity.add_background( cTag );
      if( cTagBBack<=1 && cTagBBack>=0 ) _BCTagEfficiencyPurity.add_background( cTagBBack );
    }
  else if( jetType==C_JET )
    {
      if( bTag<=1 && bTag>=0 ) _BTagEfficiencyPurity.add_background( bTag );
      if( cTag<=1 && cTag>=0 ) _CTagEfficiencyPurity.add_signal( cTag );
      if( cTagBBack<=1 && cTagBBack>=0 ) _BCTagEfficiencyPurity.add_signal( cTagBBack );
    }
  else
    {
      if( bTag<=1 && bTag>=0 ) _BTagEfficiencyPurity.add_background( bTag );
      if( cTag<=1 && cTag>=0 ) _CTagEfficiencyPurity.add_background( cTag );
      //don't add background for the bcnet result because we're only considering the b background
    }
	
}

void DSTPlotProcessor::_outputDataToFile( string filename )
{
#ifdef USEROOT
	stringstream filenameStream;
	filenameStream << filename << ".root";

	TFile rootFile( filenameStream.str().c_str(), "RECREATE", "Various plots from the plot processor" );

	//have to put all the results into normal c arrays for root to deal with
	//there may be an easier way of doing this conversion but I don't know what it is
	
	
	int numberOfPoints=100;
	
	//get the results in STL form
	vector< pair<double,double> > BTagefficiencyPurityPairs=_BTagEfficiencyPurity.eff_pur( numberOfPoints );
	vector< pair<double,double> > CTagefficiencyPurityPairs=_CTagEfficiencyPurity.eff_pur( numberOfPoints );
	vector< pair<double,double> > BCTagefficiencyPurityPairs=_BCTagEfficiencyPurity.eff_pur( numberOfPoints );
	
	//have to put all the results into normal c arrays for root to deal with
	//there may be an easier way of doing this conversion but I don't know what it is
	double BNetEfficiency[100];
	double BNetPurity[100];
	double CNetEfficiency[100];
	double CNetPurity[100];
	double BCNetEfficiency[100];
	double BCNetPurity[100];
		
	for( int a=0; a<numberOfPoints; ++a )
	  {
	    BNetEfficiency[a]=BTagefficiencyPurityPairs[a].first;
	    BNetPurity[a]=BTagefficiencyPurityPairs[a].second;
	    CNetEfficiency[a]=CTagefficiencyPurityPairs[a].first;
	    CNetPurity[a]=CTagefficiencyPurityPairs[a].second;
	    BCNetEfficiency[a]=BCTagefficiencyPurityPairs[a].first;
	    BCNetPurity[a]=BCTagefficiencyPurityPairs[a].second;
	  }
	
	TGraph BNetGraph( 99, BNetEfficiency, BNetPurity );
	BNetGraph.SetName( (" B Tag") );
	TGraph CNetGraph( 99, CNetEfficiency, CNetPurity );
	CNetGraph.SetName( ( " C Tag") );
	TGraph BCNetGraph( 99, BCNetEfficiency, BCNetPurity );
	BCNetGraph.SetName( ( " BC Tag") );
	
	TMultiGraph TagGraphs;
	TagGraphs.SetName( "TagGraphs" );
	TagGraphs.Add( &BNetGraph );
	TagGraphs.Add( &CNetGraph );
	TagGraphs.Add( &BCNetGraph );
	
	TAxis* pAxis;
	TagGraphs.SetTitle( (string("Efficiency-purity for ") ).c_str() );
		
	//why the hell doesn't this work?
	pAxis=TagGraphs.GetXaxis();
	if(pAxis) pAxis->SetTitle( "Efficiency" );
	
	pAxis=TagGraphs.GetYaxis();
	if(pAxis) pAxis->SetTitle( "Purity" );
	
	rootFile.Add( &TagGraphs );

	rootFile.Write();

	
	// now add in the jet energies
	TH1F jetEnergyHistogram( "jetEnergyHistogram", "Jet energies", 200, 0, 120 );
	const vector<double> jetEnergyData=_jetEnergy.sorted_data();
	for( vector<double>::const_iterator i=jetEnergyData.begin(); i<jetEnergyData.end(); ++i ) jetEnergyHistogram.Fill( (*i) );
	
	rootFile.Write();
	rootFile.Close();
#endif //of

}

void DSTPlotProcessor::_displayCollectionNames( LCEvent* pEvent )
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


void DSTPlotProcessor::_checkDSTParameters( LCEvent* pEvent)
{
  try{
    LCCollection* JetCollection=pEvent->getCollection( _JetCollectionName );
    
    try{
     LCCollection* trueJetFlavourCollection = pEvent->getCollection(_TrueJetFlavourColName);   
     
     try{
       LCCollection* flavourTagInputsCollection = pEvent->getCollection( _FlavourTagInputsCollectionName);
       
       try{
	 LCCollection* TagCollection=pEvent->getCollection( _FlavourTagCollectionName );


	 
	 for( int jet=0; jet<JetCollection->getNumberOfElements(); jet++ )
	   {


	     //get all the info from the full collections
	     ReconstructedParticle*  pJet  = dynamic_cast<ReconstructedParticle*>( JetCollection->getElementAt(jet) );
	     
	     //get the flavour tags
	     LCFloatVec* flavourTags = dynamic_cast<LCFloatVec*>(TagCollection->getElementAt( jet ));
	     float fullBTag = (*flavourTags)[_IndexOfForEachTag["BTag"]];
	     float fullCTag = (*flavourTags)[_IndexOfForEachTag["CTag"]];
	     float fullBCTag =(*flavourTags)[_IndexOfForEachTag["BCTag"]]; 
	     
	     LCFloatVec* pJetFlavour = dynamic_cast<LCFloatVec*>(trueJetFlavourCollection->getElementAt( jet ));
	     float fullTruePDGCode = (*pJetFlavour)[_FlavourIndex["TruePDGCode"]];
	     float fullTruePartonCharge =(*pJetFlavour)[_FlavourIndex["TruePartonCharge"]];
	     float fullTrueHadronCharge= (*pJetFlavour)[_FlavourIndex["TrueHadronCharge"]];
	     float fullTrueJetFlavour = (*pJetFlavour)[_FlavourIndex["TrueJetFlavour"]];

	     LCFloatVec* tagInputs = dynamic_cast<LCFloatVec*>(flavourTagInputsCollection->getElementAt( jet ));
	     int  fullNumVertices = int((*tagInputs)[_InputsIndex["NumVertices"]]);
	     
	     float fullJointProbRPhi=999;
	     float fullJointProbZ=999;
	     float fullNumTracksInVertices=999;
	     float fullDecayLength=999;
	     float fullDecayLengthSignificance=999;
	     float fullRawMomentum=999;
	     float fullPTCorrectedMass=999;
	     float fullSecondaryVertexProbability=999;

	     if(fullNumVertices > 1)
	       {
		 fullJointProbRPhi=(*tagInputs)[_InputsIndex["JointProbRPhi"]];
		 fullJointProbZ=(*tagInputs)[_InputsIndex["JointProbZ"]];
		 fullNumTracksInVertices=(*tagInputs)[_InputsIndex["NumTracksInVertices"]];
		 fullDecayLength=(*tagInputs)[_InputsIndex["DecayLength"]];
		 fullDecayLengthSignificance=(*tagInputs)[_InputsIndex["DecayLengthSignificance"]];
		 fullRawMomentum=(*tagInputs)[_InputsIndex["RawMomentum"]];
		 fullPTCorrectedMass=(*tagInputs)[_InputsIndex["PTCorrectedMass"]];
		 fullSecondaryVertexProbability=(*tagInputs)[_InputsIndex["SecondaryVertexProbability"]];
	       }


	     //get all the info from the DST collections
	     PIDHandler pidh( JetCollection ) ;
	     int alid = pidh.getAlgorithmID("LCFIFlavourTag");
	     const ParticleID& pid = pidh.getParticleID(pJet,alid);
	     FloatVec params = pid.getParameters();
	     
	     float dstBTag = params[pidh.getParameterIndex(alid,"BTag")];
	     float dstCTag = params[pidh.getParameterIndex(alid,"CTag")];
	     float dstBCTag =params[pidh.getParameterIndex(alid,"BCTag")]; 

	     int mcid = pidh.getAlgorithmID("MCTruth");
	     const ParticleID& mcpid = pidh.getParticleID(pJet,mcid);
	     FloatVec mcparams = mcpid.getParameters();
	     
	     float dstTruePDGCode = mcparams[pidh.getParameterIndex(mcid,"TruePDGCode")];
	     float dstTruePartonCharge =mcparams[pidh.getParameterIndex(mcid,"TruePartonCharge")];
	     float dstTrueHadronChange= mcparams[pidh.getParameterIndex(mcid,"TrueHadronCharge")];
	     float dstTrueJetFlavour = mcparams[pidh.getParameterIndex(mcid,"TrueJetFlavour")];

	     int  dstNumVertices = int(params[pidh.getParameterIndex(alid,"NumVertices")]);
	     
	     float dstJointProbRPhi=999;
	     float dstJointProbZ=999;
	     float dstNumTracksInVertices=999;
	     float dstDecayLength=999;
	     float dstDecayLengthSignificance=999;
	     float dstRawMomentum=999;
	     float dstPTCorrectedMass=999;
	     float dstSecondaryVertexProbability=999;

	     if(dstNumVertices > 1)
	       {
		 dstJointProbRPhi=params[pidh.getParameterIndex(alid,"JointProbRPhi")];
		 dstJointProbZ=params[pidh.getParameterIndex(alid,"JointProbZ")];
		 dstNumTracksInVertices=params[pidh.getParameterIndex(alid,"NumTracksInVertices")];
		 dstDecayLength=params[pidh.getParameterIndex(alid,"DecayLength")];
		 dstDecayLengthSignificance=params[pidh.getParameterIndex(alid,"DecayLengthSignificance")];
		 dstRawMomentum=params[pidh.getParameterIndex(alid,"RawMomentum")];
		 dstPTCorrectedMass=params[pidh.getParameterIndex(alid,"PTCorrectedMass")];
		 dstSecondaryVertexProbability=params[pidh.getParameterIndex(alid,"SecondaryVertexProbability")];
	       }
	     

	     //compare the values
	     cout<< dstBTag <<" "<<fullBTag<<endl;;
	     cout<< dstCTag <<" "<<fullCTag<<endl;;
	     cout<< dstBCTag <<" "<<fullBCTag<<endl;; 
	     
	     
	     cout<< dstTruePDGCode <<" "<<fullTruePDGCode<<endl;;
	     cout<< dstTruePartonCharge <<" "<<fullTruePartonCharge<<endl;;
	     cout<< dstTrueHadronChange<<" "<<fullTrueHadronCharge<<endl;;
	     cout<< dstTrueJetFlavour <<" "<<fullTrueJetFlavour<<endl;;
	     cout<< dstNumVertices <<" "<<fullNumVertices<<endl;;

	     if(dstNumVertices > 1)
	       {
		 cout<< dstJointProbRPhi <<" "<<fullJointProbRPhi<<endl;;
		 cout<< dstJointProbZ <<" "<<fullJointProbZ<<endl;;
		 cout<< dstNumTracksInVertices <<" "<<fullNumTracksInVertices<<endl;;
		 cout<< dstDecayLength <<" "<<fullDecayLength<<endl;;
		 cout<< dstDecayLengthSignificance <<" "<<fullDecayLengthSignificance<<endl;;
		 cout<< dstRawMomentum <<" "<<fullRawMomentum<<endl;;
		 cout<< dstPTCorrectedMass <<" "<<fullPTCorrectedMass<<endl;;
		 cout<< dstSecondaryVertexProbability <<" "<<fullSecondaryVertexProbability<<endl;;
	       }

	     
	     float ep = 0.00001;
	     if(fabs( dstBTag -fullBTag) > ep)
	       cout<<"B Parameters not equal"<<endl;
	     if(fabs( dstCTag -fullCTag) > ep)
	       cout<<"C Parameters not equal"<<endl;
	     if(fabs( dstBCTag -fullBCTag) > ep) 
	       cout<<"BC Parameters not equal"<<endl;
	       
	       
	     if(fabs( dstTruePDGCode -fullTruePDGCode) > ep)
	       cout<<"TPDG Parameters not equal"<<endl;
	     if(fabs( dstTruePartonCharge -fullTruePartonCharge) > ep)
	       cout<<"TP Parameters not equal"<<endl;
	     if(fabs( dstTrueHadronChange-fullTrueHadronCharge) > ep)
	       cout<<"TH Parameters not equal"<<endl;
	     if(fabs( dstTrueJetFlavour -fullTrueJetFlavour) > ep)
	       cout<<"TJ Parameters not equal"<<endl;
	     if(fabs( dstNumVertices -fullNumVertices) > ep)
	       cout<<"NV Parameters not equal"<<endl;

	     if(dstNumVertices > 1)
	       {
		 if(fabs( dstJointProbRPhi -fullJointProbRPhi) > ep)
		   cout<<"JPRP Parameters not equal"<<endl;
		 if(fabs( dstJointProbZ -fullJointProbZ) > ep)
		   cout<<"JPZ Parameters not equal"<<endl;
		 if(fabs( dstNumTracksInVertices -fullNumTracksInVertices) > ep)
		   cout<<"NT Parameters not equal"<<endl;
		 if(fabs( dstDecayLength -fullDecayLength) > ep)
		   cout<<"DL Parameters not equal"<<endl;
		 if(fabs( dstDecayLengthSignificance -fullDecayLengthSignificance) > ep)
		   cout<<"DLS Parameters not equal"<< dstDecayLengthSignificance -fullDecayLengthSignificance<<endl;
		 if(fabs( dstRawMomentum -fullRawMomentum) > ep)
		   cout<<"raw p Parameters not equal"<<endl;
		 if(fabs( dstPTCorrectedMass -fullPTCorrectedMass) > ep)
		   cout<<"pt cor Parameters not equal"<<endl;
		 if(fabs( dstSecondaryVertexProbability -fullSecondaryVertexProbability) > ep)
		   cout<<"sv Parameters not equal "<<endl;
	       }



	   }

       }

	 catch(DataNotAvailableException &e){
	     std::cout << "Collection " <<_JetCollectionName << " is unavailable " << std::endl;
	 }
       }
       catch(DataNotAvailableException &e){
	   std::cout << "Collection " <<_TrueJetFlavourColName << " is unavailable " << std::endl;
       }
     }
     catch(DataNotAvailableException &e){
	 std::cout << "Collection " << _FlavourTagInputsCollectionName << " is unavailable " << std::endl;
     }
   }
   catch(DataNotAvailableException &e){
       std::cout << "Collection " <<_FlavourTagCollectionName << " is unavailable " << std::endl;
   }
}



