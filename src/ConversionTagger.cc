#include "ConversionTagger.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include <UTIL/LCTOOLS.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include <marlin/Global.h>
#include <marlin/VerbosityLevels.h>
#include <MarlinUtil.h>
#include <HelixClass.h>

#include <gear/BField.h>

#define M_ELECTRON  0.00051
#define M_PIPLUS    0.13957
#define M_KZERO     0.49767
#define M_PROTON    0.93827
#define M_LAMBDA    1.11568

using namespace lcio;
using namespace marlin;
using namespace std;


ConversionTagger aConversionTagger ;


ConversionTagger::ConversionTagger() : Processor("ConversionTagger") {
  
  // modify processor description
  _description = "ConversionTagger processor does conversion and V0 tagging" ;
  
  _twopi=2*acos(-1.0);

  // input collections to run the conversion/V0 tagging on.
  // default is empty list, which means the code will run over
  // every single collection of type TRACK or RECONSTRUCTEDPARTICLE
  std::vector<std::string> inputCollections;
  inputCollections.clear();
  registerInputCollections( LCIO::RECONSTRUCTEDPARTICLE, 
			    "InputCollections" , 
			    "Input Collection Names (TRACK or RECONSTRUCTEDPARTICLE)" ,
			    _InputCollections ,
			    inputCollections);

  registerOptionalParameter("MassRangePhoton",
			    "upper mass limit for photon conversion candidates"
			    " (GeV)",
			    _massRangePhoton,0.01);
  registerOptionalParameter("MassRangeKaon",
			    "Kshort candidates: max deviation of candidate mass"
			    " from PDG Kshort mass (GeV)",
			    _massRangeKaon,0.02);
  registerOptionalParameter("MassRangeLambda",
			    "Lambda candidates: max deviation of candidate mass"
			    " from PDG Lambda mass",
			    _massRangeLambda,0.02);
  registerOptionalParameter("DistCut",
			    "max distance of closest approach between two "
			    " helices (mm)",
			    _distCut,1.0);
  registerOptionalParameter("MinDistFromIP",
			    "min distance of V0 candidates from IP (mm)",
			    _minDistFromIP,1.0);

  vector<int> tag_all;
  tag_all.push_back(22);
  tag_all.push_back(310);
  tag_all.push_back(3122);
  registerOptionalParameter("PdgToTag",
			    "list of particles types to look for (described"
			    " in terms of positive PDG value)",
			    _PdgToTag,tag_all);

  registerOptionalParameter("CheatMode",
			    "whether or not to use MC info for conv/V0 finding",
			    _cheatMode,false);

  registerOptionalParameter("CheatEvenMore",
			    "even tag MC conv/V0 when only one track reconstructed?",
			    _cheatEvenMore,false);

}


void ConversionTagger::init() { 

  printParameters() ;

  _BField = Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();

  for (size_t i=0; i<_PdgToTag.size(); i++) _TagPDG[_PdgToTag[i]]=true;

  histos = new HistMap(this);
}


void ConversionTagger::processEvent( LCEvent * evt ) { 

  const StringVec* colNames = evt->getCollectionNames();

  if (_InputCollections.size()==0) {
    // run over all track and reconstructed particle collections in the event
    for (size_t i=0; i<colNames->size(); i++) {
      LCCollection *collection=evt->getCollection((*colNames)[i]);
      if ((collection->getTypeName()==LCIO::TRACK) ||
	  (collection->getTypeName()==LCIO::RECONSTRUCTEDPARTICLE)) {
	tagger(evt,(*colNames)[i]);
      }
    }
  } else {
    for (size_t i=0; i<_InputCollections.size(); i++) {
      tagger(evt,_InputCollections[i]);
    }
  }

  setReturnValue(true);

}


MCParticle* ConversionTagger::isFromV0(ReconstructedParticle* rp,
				       vector<LCRelationNavigator*>relCols) {

  vector<MCParticle*> mcp;
  vector<double> mcp_weight;

  // loop over all associated tracks
  for (size_t itrk=0; itrk<rp->getTracks().size(); itrk++) {
    Track* trk=rp->getTracks()[itrk];

    // does this track have any MC correspondence in the LCRelations?
    for (size_t irel=0; irel<relCols.size(); irel++) {
      for (size_t i=0; i<relCols[irel]->getRelatedToObjects(trk).size(); i++) {
	mcp.push_back(dynamic_cast<EVENT::MCParticle*>
		      (relCols[irel]->getRelatedToObjects(trk)[i]));
	mcp_weight.push_back(relCols[irel]->getRelatedToWeights(trk)[i]);
      }
      for (size_t i=0; i<relCols[irel]->getRelatedFromObjects(trk).size(); i++) {
	mcp.push_back(dynamic_cast<EVENT::MCParticle*>
		      (relCols[irel]->getRelatedFromObjects(trk)[i]));
	mcp_weight.push_back(relCols[irel]->getRelatedToWeights(trk)[i]);
      }
    }
  }

  // find all different conv/V0 that contributed to this particle (if any)
  map<MCParticle*,double> mcmap;
  for (size_t i=0; i<mcp.size(); i++) {
    if (mcp[i]->getParents().size()==0) continue;
    MCParticle* parent = mcp[i]->getParents()[0];
    int parpdg = abs(parent->getPDG());
    while (parpdg!=22 && parpdg!=310 && parpdg!=3122
	   && parent->getParents().size()>0) {
      parent=parent->getParents()[0];
      parpdg=abs(parent->getPDG());
    }
    if (parpdg==22 || parpdg==310 || parpdg==3122) mcmap[parent]+=mcp_weight[i];
  }

  // check which conv/V0 had the largest contribution to this particle
  double maxweight=0;
  MCParticle* maxparent=NULL;
  for (map<MCParticle*,double>::iterator it=mcmap.begin();
       it!=mcmap.end(); it++) {
    if (it->second>maxweight) {
      maxweight=it->second;
      maxparent=it->first;
    }
  }
  if (mcmap.size()>1) streamlog_out(ERROR) << "more than one conv/V0 parent"
					   << " found. this means we might not"
					   << " recognize pairs properly"
					   << " in cheat mode!" << endl;
  
  // in case the weight corresponds to the fraction of hits on the track
  // coming from a given MCParticle, we require a certain minimum contribution
  // of 80% in order to declare this a conversion/V0 track.
  // in practice it seems though as if all weights are either 0 or 1.
  if (maxweight>0.8) {
    return maxparent;
  } else {
    return NULL;
  }
}


void ConversionTagger::tagger( LCEvent *evt,
			       const string collectionName) {

  // get input collection and convert to ReconstructedParticle type
  // if necessary
  LCCollection *coll=evt->getCollection(collectionName);
  bool delete_coll=false;
  if (coll->getTypeName()==LCIO::TRACK) {
    LCCollectionVec* newcoll = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    delete_coll=true;
    for (int i=0; i<coll->getNumberOfElements(); i++) {
      Track *trk=dynamic_cast<Track*>(coll->getElementAt(i));
      ReconstructedParticle *newpart = CreateRecoPart(trk);
      newcoll->addElement(newpart);
    }
    coll = newcoll;
  }

  // for cheater mode, get LCRelation of MC truth
  vector<LCRelationNavigator*> relCols;
  if (_cheatMode) {
    relCols.clear();
    const StringVec* colNames = evt->getCollectionNames();
    // find all LCRelation collections in the event
    for (size_t i=0; i<colNames->size(); i++) {
      LCCollection *collection=evt->getCollection((*colNames)[i]);
      if (collection->getTypeName()==LCIO::LCRELATION) {
	relCols.push_back(new LCRelationNavigator(collection));
      }
    }
  }

  // create output collection
  LCCollectionVec* recocoll = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  // bookkeeping of reconstructed particles tagged as conversions
  vector<bool> tagged;
  for (int i=0; i<coll->getNumberOfElements(); i++) tagged.push_back(false);

  for (int irp1=0; irp1<coll->getNumberOfElements(); irp1++) {
    ReconstructedParticle *rp1=dynamic_cast<ReconstructedParticle*>
      (coll->getElementAt(irp1));

    MCParticle* parentV0=NULL;
    if (_cheatMode) parentV0=isFromV0(rp1,relCols);
    if (_cheatEvenMore) tagged[irp1]=(parentV0!=NULL);

    for (int irp2=irp1+1; irp2<coll->getNumberOfElements(); irp2++) {
      ReconstructedParticle *rp2=dynamic_cast<ReconstructedParticle*>
	(coll->getElementAt(irp2));

      // require exactly one track each
      if (rp1->getTracks().size()!=1) continue;
      if (rp2->getTracks().size()!=1) continue;

      // get all information required for finding potential vertex
      int cand_type=0;

      if (_cheatMode) {
	// cheat mode: take all true conversions+V0 according to MC tree.
	// this is quite CPU inefficient because we start from reconstructed
	// particles rather than only checking each true conv/V0 whether it
	// has reconstructed tracks. however, starting from all combinations
	// of reconstructed particles allows this piece of code
	// to fit in better with the realistic reconstruction code.

	MCParticle* parent2V0=isFromV0(rp2,relCols);
	if (parentV0 && parent2V0) {
	  if (parentV0==parent2V0) {
	    cand_type=abs(parentV0->getPDG());
	    streamlog_out(DEBUG) << "cheat mode found cand_type="
				 << cand_type << endl;
	  } else {
	    streamlog_out(DEBUG) << "two particles from different V0 parents!"
				 << parentV0->getPDG() << " and "
				 << parent2V0->getPDG() << endl;
	    continue;
	  }
	} else {
	  continue;
	}

      } else {
        // realistic conversion/V0 reconstruction

	// require opposite charges
	if (rp1->getCharge()*rp2->getCharge()>=0) continue;

	// get helix representation of these ReconstructedParticles
	Track* trk1=rp1->getTracks()[0];
	Track* trk2=rp2->getTracks()[0];
	HelixClass helix1,helix2;
	helix1.Initialize_Canonical(trk1->getPhi(),trk1->getD0(),
				    trk1->getZ0(),trk1->getOmega(),
				    trk1->getTanLambda(),_BField);
	helix2.Initialize_Canonical(trk2->getPhi(),trk2->getD0(),
				    trk2->getZ0(),trk2->getOmega(),
				    trk2->getTanLambda(),_BField);
	const double* mom1=rp1->getMomentum();
	const double* mom2=rp2->getMomentum();

	// find distance of closest approach of helices, plus vertex properties
	float vpos[3],vmom[3];
	double dist=helix1.getDistanceToHelix(&helix2,vpos,vmom);
	histos->fill("helix_dist",dist,1,"distance between helices",100,0,20);
	if (dist>_distCut) continue;

	// for some reason the extrapolation sometimes ends up with NaN momentum
	if ( ( !(vmom[0]<0) && !(vmom[0]>=0) ) ||
	     ( !(vmom[1]<0) && !(vmom[1]>=0) ) ||
	     ( !(vmom[2]<0) && !(vmom[2]>=0) ) ) {
	  streamlog_out(ERROR) << "extrapolated momenta are NaN. "
			       << " dist=" << dist << endl;
	  continue;
	}

	// remove candidates that are very close to the IP
	double vertex_radius=sqrt(vpos[0]*vpos[0]
				  +vpos[1]*vpos[1]
				  +vpos[2]*vpos[2]);
	if (vertex_radius<_minDistFromIP) continue;

	// invariant mass of the combination: either around 0 or K0 mass?
	double conv_mass = diParticleMass(mom1,mom2,M_ELECTRON,M_ELECTRON,vmom);
	double K0_mass = diParticleMass(mom1,mom2,M_PIPLUS,M_PIPLUS,vmom);
	double Lambda_mass;
	if (mom1[0]*mom1[0]+mom1[1]*mom1[1]+mom1[2]*mom1[2]<
	    mom2[0]*mom2[0]+mom2[1]*mom2[1]+mom2[2]*mom2[2]) {
	  // lambda pion has smaller momentum than lambda proton.
	  // thus in this case we assume mom1 to be the pion
	  Lambda_mass = diParticleMass(mom1,mom2,M_PIPLUS,M_PROTON,vmom);
	} else {
	  Lambda_mass = diParticleMass(mom2,mom1,M_PIPLUS,M_PROTON,vmom);
	}
	histos->fill("conv_mass",conv_mass,1,"conv_mass",100,0,0.3);
	histos->fill("K0_mass",K0_mass,1,"K0_mass",100,0,1);
	histos->fill("Lambda_mass",Lambda_mass,1,"Lambda_mass",100,1,1.3);

	// check whether our candidate is either close to photon mass
	// or K0 mass or Lambda0 mass
	if (conv_mass>_massRangePhoton && fabs(K0_mass-M_KZERO)>_massRangeKaon
	    && fabs(Lambda_mass-M_LAMBDA)>_massRangeLambda) continue;


	// vertex probability (cut on distance of closest approach first?)
	
	// can we improve mass resolution by track refit with vertex constraint?

	// we should make sure that no track is used in more than one
	// candidate. reason: at least one of the partners will be a false
	// conversion tag then!

	histos->fill("radius",vertex_radius,1,
		     "3d radius of conv/V0 vertex candidates",100,0,2000);

	// determine type of candidate
	if (conv_mass<=_massRangePhoton) {
	  cand_type=22;
	} else if (fabs(K0_mass-M_KZERO)<=_massRangeKaon) {
	  cand_type=310;
	} else {
	  cand_type=3122;
	}
      }

      if (!_TagPDG[cand_type]) continue;

      // whatever is left here will be stored as conversion candidate
      ReconstructedParticleImpl* recopart = new ReconstructedParticleImpl();
      recopart->setType(cand_type);
      recopart->addTrack(rp1->getTracks()[0]);
      recopart->addTrack(rp2->getTracks()[0]);
      recocoll->addElement(recopart);
      tagged[irp1]=true;
      tagged[irp2]=true;
    }
  }

  // another pass: look for composite particles that might correspond
  // to a previously identified conversion or V0 candidate
  for (int irp=0; irp<coll->getNumberOfElements(); irp++) {
    ReconstructedParticleImpl *rp=dynamic_cast<ReconstructedParticleImpl*>
      (coll->getElementAt(irp));

    // we cannot check charge here for the time being because I found
    // PandoraPFA is creating Kshorts with charge +-1...
    //if (rp->getCharge()!=0) continue;
    if (rp->getTracks().size()<2) continue;
    int pdg=abs(rp->getType());
    if (pdg!=22 && pdg!=130 && pdg!=310 && pdg!=3122) continue;
    if (!_TagPDG[pdg]) continue;

    tagged[irp]=true;
    ReconstructedParticleImpl* recopart = new ReconstructedParticleImpl(*rp);
    recocoll->addElement(recopart);
  }

  evt->addCollection( recocoll , collectionName+"Conv") ;

  // create new ReconstructedParticle collection without tagged particles
  LCCollectionVec* remainder = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  remainder->setSubset();
  for (int irp=0; irp<coll->getNumberOfElements(); irp++) {
    if (!tagged[irp]) {
      ReconstructedParticleImpl* newpart=new ReconstructedParticleImpl(*dynamic_cast<ReconstructedParticleImpl*>(coll->getElementAt(irp)));
      remainder->addElement(newpart);
    }
  }
  streamlog_out(DEBUG0) << "ConversionTagger: adding V0-purged collection "
			<< "V0Veto"+collectionName << " to the event" << endl;
  evt->addCollection(remainder, "V0Veto"+collectionName) ;

  streamlog_out(DEBUG0) << "original number of entries: "
			<< coll->getNumberOfElements()
			<< ", new number of entries: "
			<< remainder->getNumberOfElements() << endl;

  // delete our temporary collection (if any).
  // destructor of LCCollectionVec automatically deletes all elements
  if (delete_coll) delete coll;

  // clean up other temporary stuff
  if (_cheatMode) {
    for (size_t i=0; i<relCols.size(); i++) delete relCols[i];
  }
}


void ConversionTagger::end(){ 
}


ReconstructedParticle* ConversionTagger::CreateRecoPart(Track* trk){

  double d0 = trk->getD0();
  double z0 = trk->getZ0();
  double omega = trk->getOmega();
  double phi0  = trk->getPhi();
  double tanlambda = trk->getTanLambda();
  HelixClass* helix = new HelixClass();
  helix->Initialize_Canonical(phi0, d0, z0, omega, tanlambda, _BField);
  double* p = new double[3];
  for (int k=0; k < 3; ++k) p[k] = helix->getMomentum()[k];
  delete helix;
  helix = 0;


  ReconstructedParticleImpl* recopart = new ReconstructedParticleImpl();
  recopart->setType(211);
  recopart->addTrack(trk);

  recopart->setCharge(trk->getOmega()/fabs(trk->getOmega()));
  recopart->setMass(0.13957);
  double momentum[3] = {0.0, 0.0, 0.0};
  double energy2=recopart->getMass()*recopart->getMass();
  for (int i = 0; i < 3; ++i) {
    momentum[i] = p[i];
    energy2+=momentum[i]*momentum[i];
  }
  recopart->setMomentum(momentum);
  recopart->setEnergy(sqrt(energy2));

  // to do
  /*
 	setCovMatrix (const float *cov)
 	setCovMatrix (const EVENT::FloatVec &)
 	setReferencePoint (const float *reference)
  */

  delete [] p;
  return recopart;
}


double ConversionTagger::diParticleMass(const double* mom1,
					const double* mom2,
					const double mass1, const double mass2,
					const float* vertmom) {

  double e1=sqrt(mass1*mass1+mom1[0]*mom1[0]+mom1[1]*mom1[1]+mom1[2]*mom1[2]);
  double e2=sqrt(mass2*mass2+mom2[0]*mom2[0]+mom2[1]*mom2[1]+mom2[2]*mom2[2]);
  double sqmass=(e1+e2)*(e1+e2);
  for (int i=0; i<3; i++) sqmass-=vertmom[i]*vertmom[i];
  //streamlog_out( DEBUG0 ) << "diParticleMass: part1=" << mom1[0] << ", "
  //			  << mom1[1] << ", " << mom1[2] << "; energy "
  //			  << e1 << endl;
  //streamlog_out( DEBUG0 ) << "                part2=" << mom2[0] << ", "
  //			  << mom2[1] << ", " << mom2[2] << "; energy "
  //			  << e2 << endl;
  //streamlog_out( DEBUG0 ) << "                mass=" << sqrt(sqmass) << endl;
  return sqrt(sqmass);
}
