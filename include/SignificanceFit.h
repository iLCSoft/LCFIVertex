#ifndef SignificanceFitProcessor_h
#define SignificanceFitProcessor_h

// include file for the maps ana processor
#include "marlin/Processor.h"
#include "lcio.h"

#include "EVENT/LCIO.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/CalorimeterHit.h"
#include "IMPL/LCFlagImpl.h"


// standard c++ includes
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>

//AIDA includes...
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IDataPointSet.h>
#include <AIDA/ITuple.h>

/** Calculates the resolution parameters needed for the joint probability tag. 
 *
 * The aim of the processor is to calculate the resolution parameters used for the joint probability flavour tag. This is done by considering the tracks with negative impact parameters. The curve is modeled as a gaussian + two exponentials.
   
 * The processor works only when using JAIDA! This processor effectively does nothing if the  USING_RAIDA flag is defined when the code has been compiled!

 * <H4>Input</H4>
 * - A collection of ReconstructedParticles that represents the jets in the event (obtained from a jet
 * finder, say SatoruJetFinderProcessor).
 * - All the other optional parameters define the beginning and the end of the local fitting of the gaussian and exponential curves. Additionally there are also parameters defining the initial input values used in the fitting.
 *
 * <H4>Output</H4>
 * The processor writes to screen the values that need to be imported ( in the same order ) into the flavour tag inputs processor. Additionally the processor displays on screen the distribution and the fitted functions. 
 *
 *  @author Erik Devetak(erik.devetak1@physics.ox.ac.uk)
*/



// namespaces using it for lcio and marlin
using namespace lcio ;
using namespace marlin ;

class SignificanceFit : public Processor 
{	
	public:
		virtual Processor*  newProcessor() { return new SignificanceFit ; }
		

		SignificanceFit();
//initialises everything		
		virtual void init() ;

		// Called for every run.
		virtual void processRunHeader( LCRunHeader* run ) ;

		// Called for every event - the working horse.		
		virtual void processEvent( LCEvent * evt ) ;


		virtual void check( LCEvent * evt ) ;


		//Called after data processing for clean up.

  		virtual void end() ;


	protected:

 /** Input collection name and type.
   */
		AIDA::IHistogram1D* rphisto;
		AIDA::IHistogram1D* zhisto;
		AIDA::IHistogram1D* rphistogauss;
		AIDA::IHistogram1D* zhistogauss;
		AIDA::IHistogram1D* rphistoexp1;
		AIDA::IHistogram1D* zhistoexp1;
		AIDA::IHistogram1D* rphistoexp2;
		AIDA::IHistogram1D* zhistoexp2;
		AIDA::IAnalysisFactory* _analysisFactory; 

		int _nRun ;
		int _nEvt ;
		std::string _JetRPColName;
		std::string _IPVertexCollectionName;
		std::string _DecayChainRPColName;
		double _Cutoffgauss;
		double _Cutoffexp1;
		double _Cutoffexp2;
		double _GaussAmp;
		double _GaussSigma;
		double _ExpAmp;
		double _ExpLambda;
	private:		
// use this for LCIO naming convention		

		
};
		
#endif
