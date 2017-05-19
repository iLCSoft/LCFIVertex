// First of all, make sure that AIDA was enabled by setting the environment
// // variable MARLIN_USE_AIDA when Marlin was compiled. The makefile will then
// // have done the setup and defined this macro.
//
#ifndef MARLIN_USE_AIDA

#warning "--------------------------------------------------------------------------------"
#warning "- SignificanceFit requires MARLIN_USE_AIDA to be defined. Did you enable -"
#warning "- AIDA when compiling Marlin? SignificanceFit will not be compiled.      -"
#warning "--------------------------------------------------------------------------------"

// Can't do anything else.
#else

//change USING_RAIDA to USING_JAIDA if you are interested in using this processor! In RAIDA the processor is only a place holder! This is really due to the limited functionality of RAIDA!!!
#define USING_RAIDA

// #ifdef USING_RAIDA
// #pragma message "USING_RAIDA defined"
// #else
// #define USING_JAIDA
// #pragma message "USING_JAIDA defined"
// #endif



#include"SignificanceFit.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <string>


// LCIO includes
#include <EVENT/LCCollection.h>
#include "EVENT/LCFloatVec.h"

#include <EVENT/Vertex.h>

//LCFI includes
#include <inc/jet.h>
#include <inc/track.h>
#include <inc/event.h>
#include <util/inc/memorymanager.h>
#include <inc/lciointerface.h>
#include <util/inc/vector3.h>
#include <util/inc/projection.h>


// AIDA includes...


#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IAxis.h>
#include <AIDA/IFitFactory.h>
#include <AIDA/IFitter.h>
#include <AIDA/IPlotter.h>
#include <AIDA/IPlotterFactory.h>
#include <AIDA/IPlotterRegion.h>
#include <AIDA/IFitResult.h>
#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IFunctionFactory.h>
#include <AIDA/IFunction.h>
#include <AIDA/IFitParameterSettings.h>


#include <UTIL/LCRelationNavigator.h>


#define PI 3.14159265

using namespace lcio ;
using namespace marlin ;
using namespace vertex_lcfi;

SignificanceFit aSignificanceFit;


SignificanceFit::SignificanceFit() : Processor("SignificanceFit")
{
  
  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			   "JetRPCollection" , 
			   "Name of the ReconstructedParticle collection that represents jets",
			   _JetRPColName ,
			   std::string("Jets") ) ;
  registerInputCollection( lcio::LCIO::VERTEX,
			   "IPVertexCollection" , 
			   "Name of the Vertex collection that contains the primary vertex (Optional)"  ,
			   _IPVertexCollectionName ,
			   std::string("IPVertex") ) ;
   registerOptionalParameter( "Cutoffgauss" ,
   			      "Maximum value for the local fitting of the gaussian, starting point 0"  ,
			     _Cutoffgauss,
			      double(5));
 registerOptionalParameter( "Cutoffexp1" ,
   			      "Maximum value for the local fitting of the first exponential, starting point Cutoffgauss"  ,
			    _Cutoffexp1,
			      double(40));

  registerOptionalParameter( "Cutoffexp2" ,
			     "Maximum value for the local fitting of the second exponential, starting point Cutoffexp1",
			     _Cutoffexp2,
			      double(200));
  registerOptionalParameter( "GaussAmlitudeInit" ,
			     "Initialization value for Gauss amplitude in the Gaussian local fit",
			     _GaussAmp,
			      double(5000));
  registerOptionalParameter( "GaussSigmaInit" ,
			     "Initialization value for Gauss sigma in the Gaussian local fit",
			     _GaussSigma,
			      double(2));

  registerOptionalParameter( "ExpAmlitudeInit" ,
			     "Initialization value for amplitude of both exponentials in the exponentials local fit",
			     _ExpAmp,
			      double(100));
  registerOptionalParameter( "ExpLambdaInit" ,
			     "Initialization value for decay constant of both exponentials in the exponentials local fit",
			     _ExpLambda,
			      double(0));
}


void SignificanceFit::init() 
{

#ifdef USING_JAIDA
  
  AIDA::IHistogramFactory* histofact = AIDAProcessor::histogramFactory(this);
  
  rphisto =histofact->createHistogram1D("-d_0 over sigma(d_0)",400,0,_Cutoffexp2);
  zhisto = histofact->createHistogram1D("-z_0 over sigma(z_0)",400,0,_Cutoffexp2);


  // A work around doing local fitting given that AIDAJNI does not yet support this. AIDA 3.3 tough does.  

  rphistogauss =histofact->createHistogram1D("-d_0 over sigma(d_0)",10,0,_Cutoffgauss);
  zhistogauss =histofact->createHistogram1D("-z_0 over sigma(z_0)",10,0,_Cutoffgauss);

  rphistoexp1 =histofact->createHistogram1D("-d_0 over sigma(d_0)",10,_Cutoffgauss,_Cutoffexp1);
  zhistoexp1 =histofact->createHistogram1D("-z_0 over sigma(z_0)",10,_Cutoffgauss,_Cutoffexp1);

  rphistoexp2 =histofact->createHistogram1D("-d_0 over sigma(d_0)",10,_Cutoffexp1,_Cutoffexp2);
  zhistoexp2 =histofact->createHistogram1D("-z_0 over sigma(z_0)",10,_Cutoffexp1,_Cutoffexp2);

  printParameters() ;
#endif	  
	
	_nRun = 0 ;
	_nEvt = 0 ;

}


void SignificanceFit::processEvent( LCEvent* )
{
#ifdef USING_JAIDA

	LCCollection* JetRPCol = evt->getCollection( _JetRPColName );
	LCCollection* VertexCol;
	VertexCol = evt->getCollection( _IPVertexCollectionName );

	Vector3 IPPos;
	Matrix3x3 IPErr;

		int nVerts = VertexCol->getNumberOfElements()  ;
		bool done = 0;

		for(int i=0; i< nVerts ; i++)
		  {
		    lcio::Vertex* iVertex = dynamic_cast<lcio::Vertex*>(VertexCol->getElementAt(i));
		    if (iVertex->isPrimary())
		      {
			IPPos.x() = iVertex->getPosition()[0];
			IPPos.y() = iVertex->getPosition()[1];
			IPPos.z() = iVertex->getPosition()[2];
			IPErr(0,0) = iVertex->getCovMatrix()[0];
			IPErr(1,0) = iVertex->getCovMatrix()[1];
			IPErr(1,1) = iVertex->getCovMatrix()[2];
			IPErr(2,0) = iVertex->getCovMatrix()[3];
			IPErr(2,1) = iVertex->getCovMatrix()[4];
			IPErr(2,2) = iVertex->getCovMatrix()[5];
			done = 1;
		      }
		    if (done) break;
		  }
		if (!done) done = 1;//TODO Throw something
		
		
		Event* MyEvent = new Event(IPPos,IPErr);

	//Event* MyEvent = new Event(Vector3(),SymMatrix3x3());

	MemoryManager<Event>::Event()->registerObject(MyEvent);
	std::vector< vertex_lcfi::Jet*> ThisVector;
	
	int nRCP = JetRPCol->getNumberOfElements()  ;
	
	if(nRCP ==0 ) std::cerr<<"Warning: SignificanceFit:336 : NO jets present "<<std::endl; 
	
	for(int i=0; i< nRCP ; i++)
	{
	  ReconstructedParticle* JetRP = dynamic_cast<ReconstructedParticle*>(JetRPCol->getElementAt(i)); 
	  vertex_lcfi::Jet* ThisJet = jetFromLCIORP(MyEvent,JetRP);
	  ThisVector.push_back(ThisJet);
	}
       	
	if(ThisVector.size()>0)
	  {
	    for(std::vector<Jet*>::iterator iJet2=ThisVector.begin();iJet2 != ThisVector.end();++iJet2)
	      {
		
	      for (std::vector<vertex_lcfi::Track*>::const_iterator iTrack = (*iJet2)->tracks().begin(); iTrack != (*iJet2)->tracks().end() ;++iTrack)   
		{

		  double SignificanceRPhi = (*iTrack)->signedSignificance(RPhi, *iJet2);
		  double SignificanceZ = (*iTrack)->signedSignificance(Z, *iJet2);
		  if(SignificanceRPhi<0)
		    {
		  				  
		      rphistogauss->fill(fabs(SignificanceRPhi));
		      rphistoexp1->fill(fabs(SignificanceRPhi));
		      rphistoexp2->fill(fabs(SignificanceRPhi));
		      rphisto->fill(fabs(SignificanceRPhi));
		    
		    }
		  if(SignificanceZ<0)
		    {
		      zhistogauss->fill(fabs(SignificanceZ));
		      zhistoexp1->fill(fabs(SignificanceZ));
		      zhistoexp2->fill(fabs(SignificanceZ));
		      zhisto->fill(fabs(SignificanceZ));
			      
		    } 
		}

	      }

	  }
	MetaMemoryManager::Event()->delAllObjects();
#endif	  

	_nEvt++;

}


void SignificanceFit::processRunHeader( LCRunHeader* )
{
  _nRun++ ;
}

void SignificanceFit::check( LCEvent* )
{



}



void SignificanceFit::end(){ 
 
#ifdef USING_JAIDA

  _analysisFactory = AIDAProcessor::GetIAnalysisFactory(this);
  AIDA::IFitFactory* fitFactory =_analysisFactory->createFitFactory();
  AIDA::IPlotter* plotter = _analysisFactory->createPlotterFactory()->create("Plot");
  AIDA::IPlotter* plotterZ = _analysisFactory->createPlotterFactory()->create("Plot");
  AIDA::IFitter* fitter      = fitFactory->createFitter("Chi2");
  AIDA::IFunctionFactory* funFactory =_analysisFactory->createFunctionFactory(*AIDAProcessor::tree(this));
 
  std::vector<double> initialParsgauss;
  std::vector<double> initialPars;
  std::vector<double> globalPars;
  std::vector<double> outPars;

  std::vector<double> zinitialParsgauss;
  std::vector<double> zinitialPars;
  std::vector<double> zglobalPars;
  std::vector<double> zoutPars;
  
  AIDA::IFunction* functiongauss = funFactory->createFunctionByName("g","g");

  AIDA::IFunction* functionexp1 =  funFactory->createFunctionByName("e","e");
  AIDA::IFunction* functionexp2 =  funFactory->createFunctionByName("e","e");

  AIDA::IFunction* functionfinal = funFactory->createFunctionByName("g","g+e+e");

  
  initialPars.push_back(_ExpAmp);
  initialPars.push_back(_ExpLambda);


  initialParsgauss.push_back(_GaussAmp);  
  initialParsgauss.push_back(0);
  initialParsgauss.push_back(_GaussSigma);

  functiongauss->setParameters( initialParsgauss );

  functionexp1->setParameters( initialPars );
  functionexp2->setParameters( initialPars );


  fitter->fitParameterSettings("mean").setBounds(-0.0000000001,0.0000000001);

  AIDA::IFitResult* resultgauss = fitter->fit(*rphistogauss, *functiongauss);
  AIDA::IFitResult* resultexp1 = fitter->fit(*rphistoexp1, *functionexp1);
  AIDA::IFitResult* resultexp2 = fitter->fit(*rphistoexp2, *functionexp2);

  std::vector<std::string> fParNamesgauss = resultgauss->fittedParameterNames();
  std::vector<std::string> fParNamesexp1 = resultexp1->fittedParameterNames();
  std::vector<std::string> fParNamesexp2 = resultexp2->fittedParameterNames();
  std::vector<double> fParsgauss     = resultgauss->fittedParameters();
  std::vector<double> fParsexp1     = resultexp1->fittedParameters();
  std::vector<double> fParsexp2     = resultexp2->fittedParameters();

  std::cout<<" RPHI - Plane    "<<std::endl;
  std::cout<<std::endl;
 
  std::cout<<" Local Fit of Gaussian   "<<std::endl;
  for(unsigned int i=0; i<fParsgauss.size();i++)
    {
      std::cout<<" Parameter   "<< fParNamesgauss[i]<<"    "<<fParsgauss[i]<<std::endl; 
      globalPars.push_back(fParsgauss[i]);
    }
  std::cout<<"Chi^2/ndf:   "<<resultgauss->quality()<<std::endl;

  std::cout<<" Local Fit of first Exponential   "<<std::endl;
  for(unsigned int i=0; i<fParsexp1.size();i++)
    {
      std::cout<<" Parameter   "<< fParNamesexp1[i]<<"    "<<fParsexp1[i]<<std::endl; 
      globalPars.push_back(fParsexp1[i]);
    }
  std::cout<<"Chi^2/ndf:   "<<resultexp1->quality()<<std::endl;

  std::cout<<" Local Fit of second Exponential   "<<std::endl;
  for(unsigned int i=0; i<fParsexp2.size();i++)
    {
      std::cout<<" Parameter   "<< fParNamesexp1[i]<<"    "<<fParsexp2[i]<<std::endl; 
      globalPars.push_back(fParsexp2[i]);	
    }

  std::cout<<"Chi^2/ndf:   "<<resultexp2->quality()<<std::endl;


  functionfinal->setParameters(globalPars);
  AIDA::IFitResult* resultfinal = fitter->fit(*rphisto, *functionfinal);

  std::vector<std::string> fParNamesglobal = resultfinal->fittedParameterNames();
  std::vector<double> fParsglobal     = resultfinal->fittedParameters();

  std::cout<<" Global Fit  "<<std::endl;
  for(unsigned int i=0; i<fParsglobal.size();i++)
    {
      std::cout<<" Parameter   "<< fParNamesglobal[i]<<"    "<<fParsglobal[i]<<std::endl; 
    }
  std::cout<<"Chi^2/ndf:   "<<resultfinal->quality()<<std::endl;

  plotter->createRegions(2,2,0);
  plotter->region(0)->plot(*rphistogauss);
  plotter->region(1)->plot(*rphistoexp1);
  plotter->region(2)->plot(*rphistoexp2);
  plotter->region(3)->plot(*rphisto);
  plotter->region(0)->plot(resultgauss->fittedFunction());
  plotter->region(1)->plot(resultexp1->fittedFunction());
  plotter->region(2)->plot(resultexp2->fittedFunction());
  plotter->region(3)->plot(resultfinal->fittedFunction());

  outPars.push_back(fParsglobal[2]);
  outPars.push_back(-sqrt(2/(PI))*fParsglobal[3]/(fParsglobal[0]*fParsglobal[2]*fParsglobal[4]));
  outPars.push_back(-(fParsglobal[4]));
  outPars.push_back(-sqrt(2/(PI))*fParsglobal[5]/(fParsglobal[0]*fParsglobal[2]*fParsglobal[6]));
  outPars.push_back(-(fParsglobal[6]));
 
  plotter->show();

  
  zinitialPars.push_back(100);
  zinitialPars.push_back(0);
  
  zinitialParsgauss.push_back(5000);  
  zinitialParsgauss.push_back(0);
  zinitialParsgauss.push_back(2);

  functiongauss->setParameters( zinitialParsgauss );

  functionexp1->setParameters( zinitialPars );
  functionexp2->setParameters( zinitialPars );

 
  fitter->fitParameterSettings("mean").setBounds(-0.0000000001,0.0000000001);

  AIDA::IFitResult* zresultgauss = fitter->fit(*zhistogauss, *functiongauss);
  AIDA::IFitResult* zresultexp1 = fitter->fit(*zhistoexp1, *functionexp1);
  AIDA::IFitResult* zresultexp2 = fitter->fit(*zhistoexp2, *functionexp2);


  std::vector<std::string> zfParNamesgauss = zresultgauss->fittedParameterNames();
  std::vector<std::string> zfParNamesexp1 = zresultexp1->fittedParameterNames();
  std::vector<std::string> zfParNamesexp2 = zresultexp2->fittedParameterNames();
  std::vector<double> zfParsgauss     = zresultgauss->fittedParameters();
  std::vector<double> zfParsexp1     = zresultexp1->fittedParameters();
  std::vector<double> zfParsexp2     = zresultexp2->fittedParameters();
  
  std::cout<<std::endl;  
  std::cout<<std::endl;  
  std::cout<<std::endl;

  std::cout<<" Z - axis    "<<std::endl;
  std::cout<<std::endl;

  std::cout<<" Local Fit of Gaussian   "<<std::endl;
  for(unsigned int i=0; i<fParsgauss.size();i++)
    {
      std::cout<<" Parameter   "<< zfParNamesgauss[i]<<"    "<<zfParsgauss[i]<<std::endl; 
      zglobalPars.push_back(zfParsgauss[i]);
    }
  std::cout<<"Chi^2/ndf:   "<<zresultgauss->quality()<<std::endl;

  std::cout<<" Local Fit of first Exponential   "<<std::endl;
  for(unsigned int i=0; i<zfParsexp1.size();i++)
    {
      std::cout<<" Parameter   "<< zfParNamesexp1[i]<<"    "<<zfParsexp1[i]<<std::endl; 
      zglobalPars.push_back(zfParsexp1[i]);
    }
  std::cout<<"Chi^2/ndf:   "<<zresultexp1->quality()<<std::endl;

  std::cout<<" Local Fit of second Exponential   "<<std::endl;
  for(unsigned int i=0; i<zfParsexp2.size();i++)
    {
      std::cout<<" Parameter   "<< zfParNamesexp1[i]<<"    "<<zfParsexp2[i]<<std::endl; 
      globalPars.push_back(zfParsexp2[i]);	
    }

  std::cout<<"Chi^2/ndf:   "<<zresultexp2->quality()<<std::endl;


  functionfinal->setParameters(zglobalPars);
  AIDA::IFitResult* zresultfinal = fitter->fit(*zhisto, *functionfinal);

  std::vector<std::string> zfParNamesglobal = zresultfinal->fittedParameterNames();
  std::vector<double> zfParsglobal     = zresultfinal->fittedParameters();

  std::cout<<" Global Fit  "<<std::endl;
  for(unsigned int i=0; i<zfParsglobal.size();i++)
    {
      std::cout<<" Parameter   "<< zfParNamesglobal[i]<<"    "<<zfParsglobal[i]<<std::endl; 
    }
  std::cout<<"Chi^2/ndf:   "<<zresultfinal->quality()<<std::endl;

  plotterZ->createRegions(2,2,0);
  plotterZ->region(0)->plot(*zhistogauss);
  plotterZ->region(1)->plot(*zhistoexp1);
  plotterZ->region(2)->plot(*zhistoexp2);
  plotterZ->region(3)->plot(*zhisto);
  plotterZ->region(0)->plot(zresultgauss->fittedFunction());
  plotterZ->region(1)->plot(zresultexp1->fittedFunction());
  plotterZ->region(2)->plot(zresultexp2->fittedFunction());
  plotterZ->region(3)->plot(zresultfinal->fittedFunction());

  plotterZ->show();

  zoutPars.push_back(zfParsglobal[2]);
  zoutPars.push_back(-sqrt(2/(PI))*zfParsglobal[3]/(zfParsglobal[0]*zfParsglobal[2]*zfParsglobal[4]));
  zoutPars.push_back(-1*zfParsglobal[4]);
  zoutPars.push_back(-sqrt(2/(PI))*zfParsglobal[5]/(zfParsglobal[0]*zfParsglobal[2]*zfParsglobal[6]));
  zoutPars.push_back(-1*zfParsglobal[6]);


  std::cout<<"PARAMETERS FOR RPHI Joint Probability"<<std::endl;
  for(unsigned int i=0;i<5;i++)
    {
      std::cout<<outPars[i]<<std::endl;
    }
  std::cout<<"PARAMETERS FOR Z Joint Probability"<<std::endl;
  for(unsigned int i=0;i<5;i++)
    {
      std::cout<<zoutPars[i]<<std::endl;
    }

  MetaMemoryManager::Run()->delAllObjects();

#endif	  
  std::cout << "SignificanceFitProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}



#endif // end if ifdef MARLIN_USE_AIDA


