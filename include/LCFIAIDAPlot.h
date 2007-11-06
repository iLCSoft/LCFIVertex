//===========================================================================================================
// LCFIAIDAPlot Class - make plots of the LCFI flavour tag and vertex charge code
// 
// Please note that as of the time of writing (November 2007)  LCFIAIDAPlot will only compile with AIDAJNI 
// and not RAIDA - this is as not all of the methods required are defined in RAIDA v01-03
//
// To compile with cmake you may need to add  -DBUILD_WITH="ROOT;AIDAJNI" -DAIDAJNI_HOME=${AIDAJNI_HOME} to your 
//
//




#ifndef MARLIN_USE_AIDA
// This check is in the C++ file, but do it again just in case this is
// included elsewhere for whatever reason.
#warning "LCFIAIDAPlot.h has been included but MARLIN_USE_AIDA is not defined. All code in LCFIAIDAPlot.h will be ignored"

#else // "else" from "ifndef MARLIN_USE_AIDA"


#ifndef LCFIAIDAPlot_h 
#define LCFIAIDAPlot_h 

//Marlin and LCIO includes 
#include "marlin/Processor.h" 
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCIntVec.h"
#include "EVENT/LCFloatVec.h"

#include <iostream>
#include <fstream>

//AIDA includes...
#include <AIDA/IHistogram1D.h>
#include <AIDA/IDataPointSet.h>
#include <AIDA/ITuple.h>


template<class T>
class TypesafeCollection
{
public:
	TypesafeCollection( lcio::LCEvent* pEvent, std::string collectionName )
	{
		// Try and get the collection
		try
		{
			_pCollection=pEvent->getCollection( collectionName );
		}
		catch(...)
		{
			_error.str(""); // clear whatever was there beforehand
			_error << __FILE__ << "(" << __LINE__ << "): Unable to get the collection \"" << collectionName << "\" for event "
					<< pEvent->getEventNumber() << " in run " << pEvent->getRunNumber() << ".";

			_pCollection=0;
			return;
		}
		
		if( !_checkCollectionType() )  //The collection is not of the expected type
		{
			_error.str("");  //clear whatever was there beforehand
			_error << __FILE__ << "(" << __LINE__ << "): The jet collection \"" << collectionName << "\" for event " << pEvent->getEventNumber()
					<< " in run " << pEvent->getRunNumber() << " is not the correct type.";

			_pCollection=0;
		}
	}

	virtual ~TypesafeCollection(){}

	bool is_valid()
	{
		if( _pCollection==0 ) return false;
		else return true;
	}

	std::string last_error()
	{
		return _error.str();
	}

	int getNumberOfElements()
	{
		if( _pCollection==0 ) return 0;
		else return _pCollection->getNumberOfElements();
	}

	T* getElementAt( int element )
	{
		if( _pCollection==0 ) return 0;
		//TODO Modify the last error message to say way these are returning 0
		if( element<0 || element>=_pCollection->getNumberOfElements() ) return 0;

		return dynamic_cast<T*>( _pCollection->getElementAt( element ) );
	}
private:
	lcio::LCCollection* _pCollection;
	std::stringstream _error; ///< Holds a string explaining the last error encountered.
	
//	template<class T>
	bool _checkCollectionType();
};

template<>
bool TypesafeCollection<lcio::ReconstructedParticle>::_checkCollectionType()
{
	return ( _pCollection->getTypeName()==LCIO::RECONSTRUCTEDPARTICLE ); //return true if they're equal, false otherwise
}

template<>
bool TypesafeCollection<lcio::LCFloatVec>::_checkCollectionType()
{
	return ( _pCollection->getTypeName()==LCIO::LCFLOATVEC ); //return true if they're equal, false otherwise
}

template<>
bool TypesafeCollection<lcio::LCIntVec>::_checkCollectionType()
{
	return ( _pCollection->getTypeName()==LCIO::LCINTVEC ); //return true if they're equal, false otherwise
}

class LCFIAIDAPlot : public marlin::Processor
{ 
public: 
	//this bit has to be done for all Marlin processors 
	virtual Processor* newProcessor() { return new LCFIAIDAPlot; } 
 
	LCFIAIDAPlot(); 
	virtual ~LCFIAIDAPlot(); 
 
	virtual void init(); 
 
	virtual void processRunHeader( LCRunHeader* pRun ); 
 
	virtual void processEvent( LCEvent* pEvent ); 
 
	virtual void check( LCEvent* pEvent ); 
 
	virtual void end(); 
protected: 
	std::vector<std::string> _FlavourTagCollectionNames;
	std::vector<std::string> _FlavourTagInputsCollectionNames;
	std::string _TrueJetFlavourColName;
	std::string _TrueJetHadronChargeColName;
	std::string _TrueJetPDGCodeColName;
	std::string _TrueJetPartonChargeColName;
	std::string _JetCollectionName;
	std::string _VertexColName;
	std::string _CVertexChargeCollection;
	std::string _BVertexChargeCollection;

	double _CosThetaJetMax;
	double _CosThetaJetMin;
	double _PTJetMin;
	double _PTJetMax;
	double _BTagNNCut;
	double _CTagNNCut;

	bool _PrintNeuralNetOutput;
	bool _MakeTuple;

	int _iVertexChargeTagCollection;
	int _myVertexChargeTagCollection;
	std::string _NeuralNetOutputFile;

	std::vector<std::string> _VertexCatNames;
	std::vector<std::string>  _NumVertexCatDir;
	std::vector<std::string> _ZoomedVarNames;
	std::string _MCParticleColName;

	AIDA::IHistogram2D* _pBJetCharge2D;
	AIDA::IHistogram2D* _pCJetCharge2D;

	AIDA::IHistogram1D* _pBJetLeakageRate;
	AIDA::IHistogram1D* _pCJetLeakageRate;
	AIDA::IHistogram1D* _pBJetVertexCharge;
	AIDA::IHistogram1D* _pCJetVertexCharge;

	std::stringstream _myStringStream; ///< Just a string stream for knocking up error/warning messages. Only implemented as a member because it's used a lot.
	std::vector< std::map<std::string,unsigned int> > _IndexOfForEachTag;
	std::vector< std::map<std::string,unsigned int> > _InputsIndex;
	std::vector< std::map<std::string,unsigned int> > _ZoomedInputsIndex;

	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _inputsHistogramsBJets;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _inputsHistogramsCJets;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _inputsHistogramsUDSJets;

	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _zoomedInputsHistogramsBJets;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _zoomedInputsHistogramsCJets;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _zoomedInputsHistogramsUDSJets;

	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetCTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetCTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetCTag;     
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBCTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBCTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBCTag;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBTagBackgroundValues;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCTagBackgroundValues;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBCTagBackgroundValues;

	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBTagIntegral;    
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetCTagIntegral;  
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetCTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetCTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pBJetBCTagIntegral; 
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pCJetBCTagIntegral;
	std::vector< std::map<std::string,AIDA::IHistogram1D*> > _pLightJetBCTagIntegral; 
					
					
	AIDA::ITuple* _pMyTuple;

	int _numberOfPoints;

	int _lastRunHeaderProcessed;
	int _suppressOutputForRun;

	bool PassesEventCuts( LCEvent* pEvent );
	bool PassesJetCuts( ReconstructedParticle* pJet );
	void FillInputsPlots( LCEvent* pEvent, unsigned int jetNumber );
	void FillTagPlots( LCEvent* pEvent, unsigned int jetNumber );
	
	static const int B_JET=5;
	static const int C_JET=4;
	static const int N_JETANGLE_BINS=10; // this is something that could potentially become a user defined parameter

	//number of different vertex categories we want to look at: 1 vertex, 2 verticies, >=3 verticies
	static const unsigned int N_VERTEX_CATEGORIES=3;  
	
	float CalculateDistance(const float* pos1, const float* pos2);
	int FindJetType( LCEvent* pEvent, unsigned int jetNumber );
	float FindJetHadronCharge(LCEvent* pEvent, unsigned int jetNumber);
	int FindJetPDGCode( LCEvent* pEvent, unsigned int jetNumber );
	float FindJetPartonCharge(LCEvent* pEvent, unsigned int jetNumber);
	
	int FindNumVertex( LCEvent* pEvent, unsigned int jetNumber, unsigned int iInputsCollection);
	int FindCQVtx( LCEvent* pEvent, unsigned int jetNumber);
	int FindBQVtx( LCEvent* pEvent, unsigned int jetNumber);
	void PrintNNOutput();
	void InternalVectorInitialisation();

	AIDA::IDataPointSet* CreateEfficiencyPlot(const AIDA::IHistogram1D* pSignal, AIDA::IDataPointSet* pDataPointSet);
	AIDA::IDataPointSet* CreateIntegralPlot(const AIDA::IHistogram1D* pNN, AIDA::IDataPointSet* pIntegral);
	AIDA::IDataPointSet* CreatePurityPlot(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet);
	AIDA::IDataPointSet* CreateLeakageRatePlot(const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet);
	
	AIDA::IDataPointSet* CreateXYPlot(const AIDA::IDataPointSet* pDataPointSet0, const AIDA::IDataPointSet* pDataPointSet1, AIDA::IDataPointSet* xyPointSet, const int dim0=0, const int dim1=0 );
	AIDA::IHistogram1D* CreateIntegralHistogram(const AIDA::IHistogram1D* pNN, AIDA::IHistogram1D* pIntegral);

	void CreateVertexChargeLeakagePlot(AIDA::IDataPointSet* pBJetVtxChargeDPS, AIDA::IDataPointSet* pCJetVtxChargeDPS);


	AIDA::IHistogram1D* _pVertexDistanceFromIP;
	AIDA::IHistogram1D* _pVertexPositionX;
	AIDA::IHistogram1D* _pVertexPositionY;
	AIDA::IHistogram1D* _pVertexPositionZ;
	
	AIDA::IHistogram1D* _pPrimaryVertexPullX;
	AIDA::IHistogram1D* _pPrimaryVertexPullY;
	AIDA::IHistogram1D* _pPrimaryVertexPullZ;
	AIDA::IHistogram1D* _pPrimaryVertexPositionX;
	AIDA::IHistogram1D* _pPrimaryVertexPositionY;
	AIDA::IHistogram1D* _pPrimaryVertexPositionZ;
	

	int _cJet_truePlus2;
	int _cJet_truePlus;
	int _cJet_trueNeut;
	int _cJet_trueMinus;
	int _cJet_trueMinus2;
	int _cJet_truePlus2_recoPlus; 
	int _cJet_truePlus2_recoNeut;
	int _cJet_truePlus2_recoMinus;
	int _cJet_truePlus_recoPlus; 
	int _cJet_truePlus_recoNeut;
	int _cJet_truePlus_recoMinus;
	int _cJet_trueNeut_recoPlus; 
	int _cJet_trueNeut_recoNeut;
	int _cJet_trueNeut_recoMinus;
	int _cJet_trueMinus_recoPlus; 
	int _cJet_trueMinus_recoNeut;
	int _cJet_trueMinus_recoMinus;
	int _cJet_trueMinus2_recoPlus; 
	int _cJet_trueMinus2_recoNeut;
	int _cJet_trueMinus2_recoMinus;
	int _bJet_truePlus2;
	int _bJet_truePlus;	
	int _bJet_trueNeut;	
	int _bJet_trueMinus;	
	int _bJet_trueMinus2;
	int _bJet_truePlus2_recoPlus; 
	int _bJet_truePlus2_recoNeut;
	int _bJet_truePlus2_recoMinus;
	int _bJet_truePlus_recoPlus; 
	int _bJet_truePlus_recoNeut;
	int _bJet_truePlus_recoMinus;
	int _bJet_trueNeut_recoPlus; 
	int _bJet_trueNeut_recoNeut;
	int _bJet_trueNeut_recoMinus;
	int _bJet_trueMinus_recoPlus; 
	int _bJet_trueMinus_recoNeut;
	int _bJet_trueMinus_recoMinus;
	int _bJet_trueMinus2_recoPlus; 
	int _bJet_trueMinus2_recoNeut;
	int _bJet_trueMinus2_recoMinus;


	std::vector< unsigned int>  _cJet_truePlus2_angle;
	std::vector< unsigned int>  _cJet_truePlus_angle;
	std::vector< unsigned int>  _cJet_trueNeut_angle;
	std::vector< unsigned int>  _cJet_trueMinus_angle;
	std::vector< unsigned int>  _cJet_trueMinus2_angle;
		     
	std::vector< unsigned int>  _cJet_truePlus2_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_truePlus2_recoNeut_angle;
	std::vector< unsigned int>  _cJet_truePlus2_recoMinus_angle;
	std::vector< unsigned int>  _cJet_truePlus_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_truePlus_recoNeut_angle;
	std::vector< unsigned int>  _cJet_truePlus_recoMinus_angle;
	std::vector< unsigned int>  _cJet_trueNeut_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_trueNeut_recoNeut_angle;
	std::vector< unsigned int>  _cJet_trueNeut_recoMinus_angle;
	std::vector< unsigned int>  _cJet_trueMinus_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_trueMinus_recoNeut_angle;
	std::vector< unsigned int>  _cJet_trueMinus_recoMinus_angle;
	std::vector< unsigned int>  _cJet_trueMinus2_recoPlus_angle; 
	std::vector< unsigned int>  _cJet_trueMinus2_recoNeut_angle;
	std::vector< unsigned int>  _cJet_trueMinus2_recoMinus_angle;
		     
	std::vector< unsigned int>  _bJet_truePlus2_angle;
	std::vector< unsigned int>  _bJet_truePlus_angle;	
	std::vector< unsigned int>  _bJet_trueNeut_angle;	
	std::vector< unsigned int>  _bJet_trueMinus_angle;	
	std::vector< unsigned int>  _bJet_trueMinus2_angle;
	std::vector< unsigned int>  _bJet_truePlus2_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_truePlus2_recoNeut_angle;
	std::vector< unsigned int>  _bJet_truePlus2_recoMinus_angle;
	std::vector< unsigned int>  _bJet_truePlus_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_truePlus_recoNeut_angle;
	std::vector< unsigned int>  _bJet_truePlus_recoMinus_angle;
	std::vector< unsigned int>  _bJet_trueNeut_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_trueNeut_recoNeut_angle;
	std::vector< unsigned int>  _bJet_trueNeut_recoMinus_angle;
	std::vector< unsigned int>  _bJet_trueMinus_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_trueMinus_recoNeut_angle;
	std::vector< unsigned int>  _bJet_trueMinus_recoMinus_angle;
	std::vector< unsigned int>  _bJet_trueMinus2_recoPlus_angle; 
	std::vector< unsigned int>  _bJet_trueMinus2_recoNeut_angle;
	std::vector< unsigned int>  _bJet_trueMinus2_recoMinus_angle;
	 
}; 



 
#endif // endif for "ifndef LCFIAIDAPlot_h"

#endif // endif for "ifndef MARLIN_USE_AIDA" ... "else"

