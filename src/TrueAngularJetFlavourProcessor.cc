#include "TrueAngularJetFlavourProcessor.h"
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
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/LCCollectionVec.h>

#include <vector>
#include <string>
#include <map>

using namespace marlin ;
using namespace lcio;
using std::vector;
using std::string;
using std::map;
using EVENT::Track;

TrueAngularJetFlavourProcessor aTrueAngularJetFlavourProcessor ;

TrueAngularJetFlavourProcessor::TrueAngularJetFlavourProcessor() : Processor("TrueAngularJetFlavourProcessor") {
  
 // modify processor description
  _description = "TrueAngularJetFlavourProcessor - Determines the true flavour of a jet from the MC Paticles associated to the Jets RP" ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( lcio::LCIO::RECONSTRUCTEDPARTICLE,
			      "JetRPCollection" , 
			      "Name of the ReconstructedParticle collection that represents jets"  ,
			      _JetRPColName ,
			      std::string("Jets") ) ;
  registerInputCollection( lcio::LCIO::MCPARTICLE,
  			      "MCParticleCollection" , 
			      "Name of the collection that holds all MC particles. "  ,
			      _MCParticleColName ,
			      std::string("MCParticle") ) ;
  registerOutputCollection( lcio::LCIO::LCINTVEC,
  			      "TrueJetFlavourCollection" , 
			      "Name of the output collection of LCIntVec (same order as jets)"  ,
			      _TrueJetFlavourColName ,
			      std::string("TrueJetFlavour") ) ;
  registerOptionalParameter( "MaximumAngle" ,
   			      "Maximum value allowed between MCParticle and jet momentum expressed in degrees"  ,
			     _MaximumAngle,
			      double(180));
  }

void TrueAngularJetFlavourProcessor::init() 
{ 
	// usually a good idea to
	printParameters() ;
	
	_nRun = 0 ;
	_nEvt = 0 ;
}

void TrueAngularJetFlavourProcessor::processRunHeader( LCRunHeader* run) { 
	_nRun++ ;
} 

void TrueAngularJetFlavourProcessor::processEvent( LCEvent * evt ) { 

	LCCollection* JetRPCol = evt->getCollection( _JetRPColName );
       	LCCollection* MCParticleCol = evt->getCollection( _MCParticleColName );

	std::vector<ReconstructedParticle*> MyJets;
	std::vector<int> MCflavour;
	
	LCCollectionVec* OutCollection = new LCCollectionVec("LCIntVec");
	evt->addCollection(OutCollection,_TrueJetFlavourColName);	

       	int nRCP = JetRPCol->getNumberOfElements()  ;

	if(nRCP ==0 ) std::cerr<<"Warning: TrueAngularFetFlavourProcessor.cc:88 : NO jets present "<<std::endl; 

	for(int i=0; i< nRCP ; i++)
	{
	  MyJets.push_back(dynamic_cast<ReconstructedParticle*>(JetRPCol->getElementAt(i))); 
	}
	
	vector<int> jetflavourdata(MyJets.size(),1);
	std::vector<std::vector<int> > Ntrkjet;
	bool originalhadron = 0;
	int code = 0;

	std::vector<MCParticle*> HeavyMCs;
	std::vector<MCParticle*> JetMCs;
	std::vector<double>  TrueFlavour;
	
	std::vector<int> typepart;


	if( MCParticleCol->getNumberOfElements() == 0 ) std::cerr<<"Warning: TrueAngularFetFlavourProcessor.cc:107 : NO MC data presentjets present "<<std::endl; 


	for(int nMCpart = 0; nMCpart< MCParticleCol->getNumberOfElements();nMCpart++ )
	  {
	    
	    MCParticle* MCused = dynamic_cast<MCParticle*> (MCParticleCol->getElementAt(nMCpart));
	    
	    //note this whole process telies on 551/100 = 5 in integer cast!!! (which is true here)
	      
	    code = abs( MCused->getPDG() );
	    
	    //this little addition will take care of wierder mesons 
	    if( code  > 10000   )
	      {
		code  = code%1000;
	      }
	    
	    if( code  > 1000   )
	      {
		
		if( code/ 1000 == 4   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(4);
		  }
		
		if( code/ 1000 == 5   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(5);
		  }		
	      }
	    
	    if( code  > 100   )
	      {
		if( code / 100  == 4   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(4);
		  }
		
		if( code/100 == 5   )
		  {
		    HeavyMCs.push_back(MCused);   
		    MCflavour.push_back(5);
		  }
	      }
	    
	  }


	for (std::vector<MCParticle*>::const_iterator iParticle = HeavyMCs.begin(); iParticle != HeavyMCs.end() ;++iParticle)	
	  {
	    originalhadron =0;
	    MCParticle* MCused = *iParticle;
	    while(!originalhadron)
	      {
		// at hadron level particles can have only 1 parent only at parton level 
		// there might be more than 1 parent. if we are at parton level we are too far anyway.
		// refer to LCIO manual for details.
		if( MCused->getParents().size() == 1) 
		  {

		    //		    std::cout<<" ParentalPDG    "<<MCused->getParents()[0]->getPDG()<<std::endl;
		    if (abs( MCused->getParents()[0]->getPDG()) >100  || (abs(MCused->getParents()[0]->getPDG())>10 && abs(MCused->getParents()[0]->getPDG())<81))
		      {
			MCused = MCused->getParents()[0];
		      }
		    else
		      {

			//check if already in the vector

			//unique does not work
			std::vector<MCParticle*>::const_iterator Doubles = find( JetMCs.begin() ,JetMCs.end() , MCused);
			if(Doubles == JetMCs.end()|| JetMCs.size() ==0 )
			  {
			    //			    std::cout<<"number    "<< MCused->getPDG() <<std::endl;
			    JetMCs.push_back(MCused);
			  }  
			originalhadron = 1;
		      }
		  }
		else
		  {
		    //check if already in the vector
			//unique does not work
		    std::vector<MCParticle*>::const_iterator Doubles = find(JetMCs.begin(),JetMCs.end(), MCused);
		    if(Doubles == JetMCs.end()|| JetMCs.size() ==0)
		      {
			//			std::cout<<"NUMBER    "<<MCused->getPDG()<<std::endl;
			JetMCs.push_back(MCused);
		      }
			originalhadron = 1;
		  }
	      }
	  }

	if(MyJets.size()>0)
	  {

	    
	    for(std::vector<MCParticle*>::const_iterator iParticle3 = JetMCs.begin(); iParticle3 != JetMCs.end() ;++iParticle3)
	      {
		double dist = -999999999; 
		int  MCcounter =0;
		int  JetCounter =0;
		int JetValue =-1;
		double angle = 0;

		vector<double> normMomMC;
		const double* MomentumMC = (*iParticle3)->getMomentum();
		double length = sqrt(MomentumMC[0]*MomentumMC[0]+MomentumMC[1]*MomentumMC[1]+MomentumMC[2]*MomentumMC[2]);
		normMomMC.push_back(MomentumMC[0]/length);
		normMomMC.push_back(MomentumMC[1]/length);
		normMomMC.push_back(MomentumMC[2]/length);
		
		//		std::cout<<"numberUSED    "<< (*iParticle3)->getPDG() <<std::endl;

		for(std::vector<ReconstructedParticle*>::const_iterator iParticle2 = MyJets.begin(); iParticle2 != MyJets.end() ;++iParticle2)
		  {
		    vector<double> normMom;
		    const double* Momentum = (*iParticle2)->getMomentum();
		    double length = sqrt(Momentum[0]*Momentum[0]+Momentum[1]*Momentum[1]+Momentum[2]*Momentum[2]);
		    normMom.push_back(Momentum[0]/length);
		    normMom.push_back(Momentum[1]/length);
		    normMom.push_back(Momentum[2]/length);
		    
		    double disttemp= sqrt(((normMomMC[0]-normMom[0])*(normMomMC[0]-normMom[0]))+
					  ((normMomMC[1]-normMom[1])*(normMomMC[1]-normMom[1]))+
					  ((normMomMC[2]-normMom[2])*(normMomMC[2]-normMom[2])));

		    //		    std::cout<<"Disttemp    "<<JetCounter<< "      "<<disttemp<<std::endl;

		    if (disttemp < fabs(dist))
		      {
			angle = 180*asin( disttemp/2 )*2/3.14159265;
	
			if(angle<_MaximumAngle)
			  {
			    dist = disttemp;
			    JetValue = JetCounter;
			  }
			else
			  {
      			    std::cout<< "Heavy MC Particle not assigned due to angle cut   "<<std::endl;
			  }
		      }
		    
		    JetCounter++;
		  };

		if(dist>0 && JetValue >=0)
		  {
		    if(jetflavourdata[JetValue]>3)
		      {
			if( jetflavourdata[JetValue] < MCflavour[MCcounter])
			  {
			    jetflavourdata[JetValue] =  MCflavour[MCcounter];
			  }
		      }
		    else
		      {
			jetflavourdata[JetValue] =  MCflavour[MCcounter];
		      }
		  }

		MCcounter++;

	      }
	  }
	
	std::cout << "TJ:";
	for(unsigned int iii=0; iii<jetflavourdata.size(); iii++)
	  {
       	    std::cout <<jetflavourdata[iii]<< " ";
	    LCIntVec *OutVec = new LCIntVec();    
	    OutVec->push_back(jetflavourdata[iii]);
	    OutCollection->addElement(OutVec);
	  };
	std::cout << ",";

	//Create the collection to store the result

	_nEvt ++ ;
}



void TrueAngularJetFlavourProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrueAngularJetFlavourProcessor::end(){ 
	
	std::cout << "TrueAngularJetFlavourProcessor::end()  " << name() 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

}
