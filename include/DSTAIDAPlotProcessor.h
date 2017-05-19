#ifndef DSTAIDAPlotProcessor_h
#define DSTAIDAPlotProcessor_h

#include <vector>
#include <algorithm>

#include "../vertex_lcfi/util/inc/util.h"

//Marlin and LCIO includes
#include "marlin/Processor.h"
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"


//AIDA includes...
#include <AIDA/IHistogram1D.h>
#include <AIDA/IDataPointSet.h>
#include <AIDA/ITuple.h>


class DSTAIDAPlotProcessor : public marlin::Processor
{
public:
	//The usual Marlin processor methods
	virtual Processor* newProcessor() { return new DSTAIDAPlotProcessor; }
	DSTAIDAPlotProcessor();
	DSTAIDAPlotProcessor(const DSTAIDAPlotProcessor&) = delete;
	DSTAIDAPlotProcessor& operator=(const DSTAIDAPlotProcessor&) = delete;
	virtual ~DSTAIDAPlotProcessor();
	virtual void init();
	virtual void processRunHeader( LCRunHeader* pRun );
	virtual void processEvent( LCEvent* pEvent );
	//don't need this
	//virtual void check( LCEvent* pEvent );
	virtual void end();
protected:
	std::string _JetCollectionName{};	/**< @internal The name of the collection of ReconstructedParticles that is the jet (comes from the steering file).*/


	int _nRun=-1; /**< @internal The current run number.*/



	//!Histograms of the neural net B-tag outputs for true light-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pLightJetBTag{};
	//!Histograms of the neural net C-tag outputs for true light-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pLightJetCTag{};
	//!Histograms of the neural net B-tag outputs for true B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pBJetBTag{};
	//!Histograms of the neural net C-tag outputs for true B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pBJetCTag{};
	//!Histograms of the neural net B-tag outputs for true C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pCJetBTag{};
	//!Histograms of the neural net C-tag outputs for true C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pCJetCTag{}; 
	//!Histograms of the neural net BC-tag outputs for true B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)   
	std::map<std::string,AIDA::IHistogram1D*> _pBJetBCTag{};
	//!Histograms of the neural net BC-tag outputs for true C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pCJetBCTag{};
	//!Histograms of the neural net BC-tag outputs for true light-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pLightJetBCTag{};
	//!Histograms of the neural net B-tag outputs for non B-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pBTagBackgroundValues{};
	//!Histograms of the neural net C-tag outputs for non C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pCTagBackgroundValues{};
	//!Histograms of the neural net BC-tag outputs for non C-jets - seperately for different number of vertices in the jets, 1, 2, >=3, any (sum of previous)
	std::map<std::string,AIDA::IHistogram1D*> _pBCTagBackgroundValues{};
	
	//!Histograms of the neural net tags - number of events that pass a given cut: jet NN value > given NN value for the three tags - B-tag, C-tag, BC-tag
	//!  - separately for true B jets, true C jets & true light jets and different number of vertices in the jets, 1, 2 or >=3 & any (sum of previous three)
	//! See comments above
	std::map<std::string,AIDA::IHistogram1D*> _pBJetBTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pCJetBTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pLightJetBTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pBJetCTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pCJetCTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pLightJetCTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pBJetBCTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pCJetBCTagIntegral{};
	std::map<std::string,AIDA::IHistogram1D*> _pLightJetBCTagIntegral{};


	std::vector<std::string> _VertexCatNames{};
	std::vector<std::string>  _NumVertexCatDir{};
	
	//!Histograms of the neural net inputs for true B-jets
	std::map<std::string,AIDA::IHistogram1D*>  _inputsHistogramsBJets{};
	//!Histograms of the neural net inputs for true C-jets
	std::map<std::string,AIDA::IHistogram1D*>  _inputsHistogramsCJets{};
	//!Histograms of the neural net inputs for light B-jets
	std::map<std::string,AIDA::IHistogram1D*>  _inputsHistogramsUDSJets{};



	//!number of different vertex categories we want to look at: 1 vertex, 2 vertices, >=3 vertices
	static const unsigned int N_VERTEX_CATEGORIES=3;  

	int _numberOfPoints=0;

	//!Tuple of the input variables - only filled for one input collection - selected with UseFlavourTagCollectionForVertexCharge
	AIDA::ITuple* _pMyTuple=nullptr;

	//useful constants
	static const int C_JET=4;/**< @internal Useful constant for the jet flavour*/
	static const int B_JET=5;/**< @internal Useful constant for the jet flavour*/

	void _displayCollectionNames( lcio::LCEvent* pEvent );/**< @internal Just prints out the available collections in the LCIO file to standard output.*/
	bool _passesEventCuts( lcio::LCEvent* pEvent );	///< @internal A function that contains all the event cuts - returns true if the event passes all of the cuts, false otherwise.
	bool _passesJetCuts( lcio::ReconstructedParticle* pJet ); ///< @internal A function that contains all the jet cuts - returns true if the event passes all of the cuts, false otherwise.


       	void FillTagPlots( LCEvent* pEvent, unsigned int jetNumber );

	void CreateTagPlots();


	void CreateFlavourTagInputPlots(LCRunHeader* pRun );
	void CreateFlavourTagTuple();
	void CalculateIntegralAndBackgroundPlots();
	void CalculateEfficiencyPurityPlots();

	void FillInputsPlots( LCEvent* pEvent, unsigned int jetNumber );


	AIDA::IDataPointSet* CreateEfficiencyPlot(const AIDA::IHistogram1D* pSignal, AIDA::IDataPointSet* pDataPointSet);
	AIDA::IHistogram1D* CreateIntegralHistogram(const AIDA::IHistogram1D* pNN, AIDA::IHistogram1D* pIntegral);
	AIDA::IDataPointSet* CreatePurityPlot(const AIDA::IHistogram1D* pSignal, const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet);
	AIDA::IDataPointSet* CreateLeakageRatePlot(const AIDA::IHistogram1D* pBackground, AIDA::IDataPointSet* pDataPointSet);
	AIDA::IDataPointSet* CreateXYPlot(const AIDA::IDataPointSet* pDataPointSet0, const AIDA::IDataPointSet* pDataPointSet1, AIDA::IDataPointSet* xyPointSet, const int dim0=0, const int dim1=0 );





};

#endif //ifndef DSTAIDAPlotProcessor_h
