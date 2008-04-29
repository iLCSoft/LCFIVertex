#include "KnowYourInputs.h"

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <set>
#include <cmath>

// Marlin includes
#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>

// MarlinUtil includes
#include "HelixClass.h"

// LCIO includes
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/LCObject.h>
#include <EVENT/SimTrackerHit.h>

// GEAR includes
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gear/VXDParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/BField.h>


using namespace lcio;
using namespace marlin;
using namespace std;


KnowYourInputs aKnowYourInputs;


KnowYourInputs::KnowYourInputs() : Processor("KnowYourInputs") {
  
  _description = "KnowYourInputs processor evaluates hits+clusters+tracks" ;
   
  registerOptionalParameter("doTrack",
			    "Enable Track diagnostics",
			    _trackEnable, bool(1));

  registerOptionalParameter("doReconstructedParticle",
			    "Enable ReconstructedParticle diagnostics",
			    _recParEnable, bool(1));

  registerOptionalParameter("doSimTrackerHit",
			    "Enable SimTrackerHit diagnostics",
			    _simTrkHitEnable, bool(1));
  
  registerOptionalParameter("doTrackerHit",
			    "Enable TrackerHit diagnostics",
			    _trkHitEnable, bool(1));

  registerOptionalParameter("doCalorimeterHit",
			    "Enable CalorimeterHit diagnostics",
			    _calHitEnable, bool(1));
  
  registerOptionalParameter("doSimCalorimeterHit",
			    "Enable SimCalorimeterHit diagnostics",
			    _simCalHitEnable, bool(1));
}


void KnowYourInputs::init() { 

  printParameters() ;

  // GEAR parameters used within this code
  _bField = Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();
  message<marlin::MESSAGE>(log() << "using GEAR B field " << _bField << " T");

  _knownCollections.clear();
  _histograms.clear();

}


void KnowYourInputs::processEvent( LCEvent * evt ) { 


  message<marlin::MESSAGE>(log() << "detector " << evt->getDetectorName() 
			   << ", run " << evt->getRunNumber()
			   << ", event " << evt->getEventNumber());

  // get all collection names in this event and make sure we did not lose any
  // collections since the last event. (we'll deal with new collections later)
  size_t last_num_of_collections=_knownCollections.size();
  const StringVec* colNames = evt->getCollectionNames();
  for (size_t i=0; i<_knownCollections.size(); i++) {
    bool found=false;
    for (size_t k=0; k<colNames->size(); k++) {
      found |= (_knownCollections[i]==(*colNames)[k]);
    }
    if (!found && _knownCollectionTypes[i]!=LCIO::SIMCALORIMETERHIT) {
      // a collection that was encountered in a previous event is now gone.
      // this usually indicates serious trouble, except when it happens for
      // SimCalorimeterHit collections where it seems to be normal
      message<marlin::ERROR>(log() << "collection " << _knownCollections[i]
			     << " not present in this event!");
    }
  }

  // loop over all collections and make standard plots
  for (size_t i=0; i<colNames->size(); i++) {

    // read collection
    string collName=(*colNames)[i];
    LCCollection *collection=evt->getCollection(collName);
    string collType=collection->getTypeName();

    // many plots are a bit specific to individual subdetectors
    // (e.g. as far as resolutions go). therefore, try to identify
    // known subdetectors from the collection name where applicable
    string id = collName;
    for (size_t n=0; n<id.size(); n++) id[n]=toupper(id[n]);
    string subdet="N/A";
    if (id.find("VXD")!=string::npos) {
      subdet="VXD";
    } else if (id.find("VTX")!=string::npos) {
      subdet="VXD";
    } else if (id.find("TPC")!=string::npos) {
      subdet="TPC";
    } else if (id.find("SIT")!=string::npos) {
      subdet="SIT";
    } else if (id.find("FTD")!=string::npos) {
      subdet="FTD";
    } else if (id.find("ETD")!=string::npos) {
      subdet="ETD";
    } else if (id.find("ECAL")!=string::npos) {
      subdet="ECAL";
    } else if (id.find("HCAL")!=string::npos) {
      subdet="HCAL";
    } else if (id.find("MUON")!=string::npos) {
      subdet="MUON";
    }

    // check whether this collection is already known.
    // changes in the list of collections usually indicate serious trouble.
    bool found=false;
    for (size_t k=0; k<_knownCollections.size(); k++) {
      found |= (_knownCollections[k]==collName);
    }
    if (!found) {
      if (last_num_of_collections==0 || collType==LCIO::SIMCALORIMETERHIT) {
	// two cases are normal:
	// a) there are no previous collections recorded at all
	//    -> we are looking at the very first event.
	// b) SimCalorimeterHit collections drop in and out anyway
	message<marlin::MESSAGE>(log() << "collection found: " << collName
				 << " (" << collType << ")");
      } else {
	// if this happens in any other situation, it is a problem!
	message<marlin::ERROR>(log() << "new " << collType << " collection"
			       << " appearing in this event: " << collName);
      }
      // now declare this a known collection
      _knownCollections.push_back(collName);
      _knownCollectionTypes.push_back(collType);
    }

    // now branch out to specialist routines
    if ((collection->getTypeName()==LCIO::TRACK)
	&& _trackEnable ) {
      trackerPlots(evt,collName,subdet);
    } else if ((collection->getTypeName()==LCIO::RECONSTRUCTEDPARTICLE)
	       && _recParEnable) {
      trackerPlots(evt,collName,subdet);
    } else if ((collection->getTypeName()==LCIO::SIMTRACKERHIT)
	       && _simTrkHitEnable) {
      simTrackerHitPlots(evt,collName,subdet);
    } else if ((collection->getTypeName()==LCIO::TRACKERHIT)
	       && _trkHitEnable) {
      trackerHitPlots(evt,collName,subdet);
    } else if ((collection->getTypeName()== LCIO::CALORIMETERHIT)
	       && _calHitEnable) {
    } else if ((collection->getTypeName()== LCIO::SIMCALORIMETERHIT)
	       && _simCalHitEnable) {
      calPlots(evt, collName,subdet);
    }
  }
  setReturnValue(true);
}


void KnowYourInputs::trackerHitPlots( const LCEvent *evt,
				      const string id, const string subdet ) {

  // create or select histogram directory
  if (!_histograms[id])
    _histograms[id]= new HistMap(this,"TrackerHit/"+id);
  HistMap* histos= _histograms[id];

  // put histogram name together
  stringstream histname,histtitle;
  histname << "resolution_" << id;
  histtitle << " distance (mm) between rawhits and trackerhits, " << id;

  // guess range of expected deviations from subdetector name
  double histrange;
  if (subdet=="VXD") histrange=0.05;
  else if (subdet=="TPC") histrange=2;
  else if (subdet=="SIT") histrange=0.1;
  else if (subdet=="FTD") histrange=0.1;
  else {
    histrange=0.2;
    message<marlin::WARNING>(log() << "cannot identify subdetector from"
			     << " TrackerHit collection name " << id);
  }

  // loop over all hits
  LCCollection *coll=evt->getCollection(id);
  for (int ihit=0; ihit<coll->getNumberOfElements(); ihit++) {
    EVENT::TrackerHit* trkHit=dynamic_cast<EVENT::TrackerHit*>
      ( coll->getElementAt(ihit) );

    // check distance to matching SimTrackerHits (i.e. tracker resolution)
    LCObjectVec simHits= trkHit->getRawHits();
    for(unsigned int simElement=0;simElement<simHits.size();simElement++){
      SimTrackerHit* simHit=dynamic_cast<EVENT::SimTrackerHit*>
	( simHits[simElement]);
      
      const double* recPos=trkHit->getPosition();
      const double* simPos=simHit->getPosition();
      double dist=sqrt(pow(recPos[2]-simPos[2],2)
		       +pow(recPos[1]-simPos[1],2)
		       +pow(recPos[0]-simPos[0],2));
      histos->fill(histname.str(),dist,1,histtitle.str(),100,0,histrange);
    }
  }
}


void KnowYourInputs::simTrackerHitPlots( const LCEvent *evt,
					 const string id, const string subdet) {

  // create or select histogram directory
  if (!_histograms[id])
    _histograms[id]= new HistMap(this,"SimTrackerHit/"+id);
  HistMap* histos= _histograms[id];
  
  // sort all hits with same cell ID (=same layer) into vectors
  std::map<int,std::vector<EVENT::SimTrackerHit*> > cellMap;
  LCCollection *coll=evt->getCollection(id);
  for (int ihit=0; ihit<coll->getNumberOfElements(); ihit++) {
    EVENT::SimTrackerHit* trkHit=dynamic_cast<EVENT::SimTrackerHit*>
      ( coll->getElementAt(ihit) );
    
    int cellID=trkHit->getCellID();
    // cell ID 0 usually means no cell ID available, so skip these hits
    if (cellID==0) continue;

    cellMap[cellID].push_back(trkHit);    
  }

  // now loop over all cell IDs we found to find closest hits within each layer
  for(std::map<int,std::vector<EVENT::SimTrackerHit*> >::iterator
	iter=cellMap.begin(); iter!=cellMap.end();iter++){
    int cellID=(*iter).first;

    // define histogram properties for hits with this cell ID
    stringstream histname,histtitle;
    histname << "mindist_" << id << "_layer" << abs(cellID);
    histtitle << " min distance (mm) between hits from different particles, "
	      << id << ", layer " << abs(cellID);
    if (subdet=="TPC" || subdet=="SIT" || subdet=="VXD") {
      double r=sqrt(pow(cellMap[cellID][0]->getPosition()[0],2)
		    +pow(cellMap[cellID][0]->getPosition()[1],2));
      histtitle << ", radius " << r;
    }
    double histrange;
    if (subdet=="VXD") histrange=0.2;
    else if (subdet=="TPC") histrange=100;
    else if (subdet=="FTD") histrange=1;
    else if (subdet=="SIT") histrange=0.5;
    else if (subdet=="ETD") histrange=5;
    else {
      histrange=100;
      message<marlin::WARNING>(log() << "cannot identify subdetector from"
			       << " SimTrackerHit collection name " << id);
    }

    // now loop over all combinations of hits with same cell ID
    for(unsigned int ihit=0; ihit<cellMap[cellID].size();ihit++){
      double mindist=999999.;
      for(unsigned int ihit2=ihit+1;ihit2<cellMap[cellID].size();ihit2++){
	// check hits from different particle only
	if (cellMap[cellID][ihit]->getMCParticle()
	    !=cellMap[cellID][ihit2]->getMCParticle()){
	  const double* pos1 = cellMap[cellID][ihit]->getPosition();
	  const double* pos2 = cellMap[cellID][ihit2]->getPosition();
	  double dist=sqrt(pow(pos2[2]-pos1[2],2)
			   +pow(pos2[1]-pos1[1],2)
			   +pow(pos2[0]-pos1[0],2));
	  if (dist<mindist) mindist=dist;
	}
      }
      histos->fill(histname.str(),mindist,1,histtitle.str(),100,0,histrange);
    }
  }
}


void KnowYourInputs::calPlots( const LCEvent *evt,
			       const string collectionName,
			       const string subdet ) {

  if (!_histograms[collectionName])
    _histograms[collectionName]= new HistMap(this,"CalorimeterHit/"+collectionName);
  HistMap* histos= _histograms[collectionName];

  
  LCCollection *coll=evt->getCollection(collectionName);
  LCRelationNavigator relColl(evt->getCollection("RelationCaloHit"));
  
  string id=collectionName;
    
  double rMinBarrel = Global::GEAR->getEcalBarrelParameters().getExtent()[0];
  double rMaxBarrel = Global::GEAR->getHcalBarrelParameters().getExtent()[1];
  double zMaxBarrel = max(Global::GEAR->getHcalBarrelParameters().getExtent()[3],Global::GEAR->getEcalBarrelParameters().getExtent()[3]);
  double zMinBarrel = -zMaxBarrel;
  
  double rMinEndcap = min(Global::GEAR->getHcalEndcapParameters().getExtent()[0],Global::GEAR->getEcalEndcapParameters().getExtent()[0]);
  double rMaxEndcap = max(Global::GEAR->getHcalEndcapParameters().getExtent()[1],Global::GEAR->getEcalEndcapParameters().getExtent()[1]);
  double zMaxEndcap = Global::GEAR->getHcalEndcapParameters().getExtent()[3];
  double zMinEndcap = - zMaxEndcap;
      
  for (int ihit=0; ihit<coll->getNumberOfElements(); ihit++) {
     
    double energy=0;
    double positionX=0;
    double positionY=0;
    double positionZ=0;
    double positionR=0;

    if(coll->getTypeName()== LCIO::CALORIMETERHIT){
      EVENT::CalorimeterHit *calHit=dynamic_cast<EVENT::CalorimeterHit*>
	( coll->getElementAt(ihit) );

      energy = calHit->getEnergy();
      positionX = calHit->getPosition()[0];
      positionY = calHit->getPosition()[1];
      positionZ = calHit->getPosition()[2];
      positionR = sqrt(pow(positionX,2) + pow(positionY,2));
      
      
      LCObjectVec relSimObject = (relColl.getRelatedToObjects(calHit));
      
      for(unsigned int vecElem=0; vecElem<relSimObject.size();vecElem++){	
	SimCalorimeterHit *simCalRelationHit = dynamic_cast<EVENT::SimCalorimeterHit*>
	  (relSimObject[vecElem]);
       
	double positionX1=simCalRelationHit->getPosition()[0];
	double positionY1=simCalRelationHit->getPosition()[1];
	double positionZ1=simCalRelationHit->getPosition()[2];

	double dist=sqrt(pow((positionX-positionX1),2)+ pow((positionY-positionY1),2)+ pow((positionZ-positionZ1),2));
	double distZ=positionZ-positionZ1;

	histos->fill("Dist_"+id,dist,1,"3D distance between the CaloHit and SimCaloHit ("+id+")",100,0,20);
	histos->fill("Dist2D_"+id,positionZ1,distZ,1,"Distance between the CaloHit and SimCaloHit wrt Z ("+id+")",200,zMinEndcap-200,zMaxEndcap+200,200,-20,50);
	histos->fill("Dist_1D"+id,distZ,1,"Distance between Calohit and SimCaloHit ("+id+")",200,-10,20);
      }
    }
    
    else if(coll->getTypeName()== LCIO::SIMCALORIMETERHIT){
      EVENT::SimCalorimeterHit *simCalHit=dynamic_cast<EVENT::SimCalorimeterHit*>
	( coll->getElementAt(ihit) );

      energy = simCalHit->getEnergy();
      positionX = simCalHit->getPosition()[0];
      positionY = simCalHit->getPosition()[1];
      positionZ = simCalHit->getPosition()[2];
      positionR = sqrt(pow(positionX,2) + pow(positionY,2));
    }
   
      histos->fill("Energy_"+id,energy,1,"Energy of hit ("+id+")",1000,0,5);
      histos->fill("PositionXvsY_"+id,positionX,positionY,1,"Position of hit in X vs Y coordinates ("+id+")",500,-3000,3000,500,-3000,3000);
      histos->fill("PositionZvsR_Barrel_"+id,positionZ,positionR,1,"Position of hit in Z vs R coordinates in Barrel ("+id+")",1000,zMinBarrel-200,zMaxBarrel+200,1000,rMinBarrel-200,rMaxBarrel+200);
      histos->fill("PositionZvsR_Endcap_"+id,positionZ,positionR,1,"Position of hit in Z vs R coordinates in Endcap ("+id+")",1000,zMinEndcap-200,zMaxEndcap+200,1000,rMinEndcap-200,rMaxEndcap+200);
      
  }
  histos->fill("numHits_"+id,coll->getNumberOfElements(),1,"Number of hits ("+id+")",100,100,1200);
   
}


void KnowYourInputs::trackerPlots( const LCEvent *evt,
				   const string collectionName,
				   const string subdet) {


  if (!_histograms[collectionName])
    _histograms[collectionName]= new HistMap(this,"Track/"+collectionName);
  HistMap* histos= _histograms[collectionName];

  LCCollection *coll=evt->getCollection(collectionName);

  // find matching collection with MCParticle references
  LCCollection *refColl=0;
  const StringVec* colNames = evt->getCollectionNames();
  for (size_t i=0; i<colNames->size(); i++) {
    if ((*colNames)[i]==collectionName+"MCP") {
      refColl=evt->getCollection(collectionName+"MCP");
    }
  }

  if (!refColl) {
    message<marlin::WARNING>("collection "+collectionName+"MCP not found.");
  } else {

  }

  string id=collectionName;
  string type="";
  
  HelixClass * helix = new HelixClass();

  for (int itrk=0; itrk<coll->getNumberOfElements(); itrk++) {
    
    double pt=0;
    double pz=0;
    double charge=0;
    double phi=0;
    double costheta=0;
    double halfpi = acos(0.0);
    double Pi = 2*halfpi;

    if(coll->getTypeName()== LCIO::TRACK){
      EVENT::Track *trk=dynamic_cast<EVENT::Track*>
	( coll->getElementAt(itrk) );

     type="Track";

      phi=trk->getPhi();
      while (phi>2*Pi) phi-=2*Pi;
      while (phi<0) phi+=2*Pi;
      double d0=trk->getD0();
      double z0=trk->getZ0();
      double omega=trk->getOmega();
      double tanlambda=trk->getTanLambda();
      double theta=halfpi-atan(tanlambda);

      costheta=cos(theta);

      // convert curvature to momentum
      helix->Initialize_Canonical(phi,d0,z0,omega,tanlambda,_bField);
      pt=sqrt(helix->getMomentum()[0]*helix->getMomentum()[0]
	      +helix->getMomentum()[1]*helix->getMomentum()[1]);
      pz=helix->getMomentum()[2];
      double psquare=helix->getMomentum()[0]*helix->getMomentum()[0]
	+helix->getMomentum()[1]*helix->getMomentum()[1]
	+helix->getMomentum()[2]*helix->getMomentum()[2];
      charge=helix->getCharge();


      //-----------------
      // basic quantities
      //-----------------
      
      histos->fill("d0Track_"+id,d0,1,"d_{0} ("+id+")", 100,-1,1);
      histos->fill("z0Track_"+id,z0,1,"z_{0} ("+id+")", 100,-100,100);
      
      histos->fill("curvatureTrack_"+id,omega,1,"curvature ("+id+")", 100,-1,1);
      histos->fill("ndfTrack_"+id,trk->getNdf(),1,
		"track fit degrees of freedom ("+id+")", 400,0,400);
      histos->fill("chi2Track_"+id,trk->getChi2()/trk->getNdf(),1,
		"#Chi^{2}/ndf ("+id+")", 100,0,50);
      histos->fill("innerRadiusTrack_"+id,trk->getRadiusOfInnermostHit(),1,
		"Radius of innermost hit ("+id+")",100,0,300);
      
            
      //-------------
      // correlations
      //-------------
      // compare cost(theta) from momentum vector with tan(lambda)
      if (psquare>0) {
	histos->fill("pcost_vs_cost_"+id,pz/sqrt(psquare),cos(theta),1,
		  "momentum cos(theta) vs cos(theta) from tan(#Lambda), "+id,
		  100,-1.01,1.01,100,-1.01,1.01);
      }
      
      //--------------
      // special cases
      //--------------
      if (fabs(pz)>1000) {
	histos->fill("phi_for_bad_pz_"+id,trk->getPhi(),1,"phi for bad pz ("+id+")",
		  100,-Pi,3.5*Pi);
      }
      
      //--------------------
      // error matrix checks
      //--------------------
      histos->fill("errd0_"+id,sqrt(trk->getCovMatrix()[0]),1,
		"#sigma(d_{0}), "+id,100,0,0.1);
      histos->fill("errz0_"+id,sqrt(trk->getCovMatrix()[10]),1,
		"#sigma(z_{0}), "+id,100,0,0.002);
      double xscat=pow(sin(theta),3)/psquare;
      histos->fill("errd0_vs_mat_"+id,xscat,trk->getCovMatrix()[0],1,
		"d_{0} error vs (1/p sin^{3}#Theta)^{2}, "+id,
		100,0,25,100,0,0.001);
      histos->fill("errz0_vs_mat_"+id,xscat,sqrt(trk->getCovMatrix()[10]),1,
		"z_{0} error vs (1/p sin^{3}#Theta)^{2}, "+id,
		100,0,25,100,0,0.001);
      
      if (psquare>100 and fabs(cos(theta))<0.5) {
	vector<int> hitnums = trk->getSubdetectorHitNumbers();
	if (hitnums.size()>=4)
	  histos->fill("numTPChits_"+id,hitnums[3],1,"number of TPC hits, tracks above 10GeV and within |cos(#theta)|<0.5",40,0,400);
      }
    }
    else if(coll->getTypeName()== LCIO::RECONSTRUCTEDPARTICLE){
      EVENT::ReconstructedParticle *recPar=dynamic_cast<EVENT::ReconstructedParticle*>
	( coll->getElementAt(itrk) );
      
      type="RecParticles";

      double momentumX = recPar->getMomentum()[0];
      double momentumY = recPar->getMomentum()[1];
      double momentumZ = recPar->getMomentum()[2];
      double energy = recPar->getEnergy();

      phi=atan2(momentumY,momentumX);
      while (phi>2*Pi) phi-=2*Pi;
      while (phi<0) phi+=2*Pi;

      pt=sqrt(pow(momentumX,2)+pow(momentumY,2));
      pz=momentumZ;
      double psquare=(pow(momentumX,2)+pow(momentumY,2)+pow(momentumZ,2));

      costheta=pz/sqrt(psquare);
      charge=recPar->getCharge();

      histos->fill("energyRec_"+id,energy,1,"Energy of the reconstructed particle("+id+")", 200,0,100);

      cout << collectionName << "RECOPARTICLE:"
	   << " id=" << recPar->getType()
           << " numtrack=" << recPar->getTracks().size()
	   << " numpart=" << recPar->getParticles().size()
	   << " numclust=" << recPar->getClusters().size() << endl;
     }

    //-----------------
    // basic quantities
    //-----------------
   
    histos->fill("pt"+type+"_"+id,pt,1,"p_{T} ("+id+")", 100,0,100);
    histos->fill("pz"+type+"_"+id,pz,1,"p_{z} ("+id+")", 100,-100,100);
    histos->fill("charge_"+id,charge,1,"Charge of the particle ("+id+")", 50,-2,2);
    histos->fill("phi"+type+"_"+id,phi,1,"phi of "+type+" ("+id+")",
	      50,-Pi,3.5*Pi);
    histos->fill("phiRestr"+type+"_"+id,phi,1,"phi of "+type+", forced to 0..2#pi ("
	      +id+")",50,0,2*Pi);
    histos->fill("costheta_"+id,costheta,1,
	      "cos(#Theta) ("+id+")",102,-1.01,1.01);
  }

  histos->fill("num"+type+"_"+id,coll->getNumberOfElements(),1,
	    "Number of "+type+" ("+id+")",100,0,100);
}



void KnowYourInputs::end(){ 
}
