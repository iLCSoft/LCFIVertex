#include "V0Performance.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

#include <UTIL/LCTOOLS.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/LCCollectionVec.h>

#include <marlin/VerbosityLevels.h>

// fix for mac osx (K. Fujii)
#ifdef __APPLE__
#define abs(a) (a >= 0 ? a : -a)
#endif

using namespace lcio;
using namespace marlin;
using namespace std;

//===============================
// what are we trying to do here?
//===============================
//
// - find actual conversions and V0s in Mokka output
// - check which of their tracks are reconstructed
// - check which of their tracks are identified as particle flow objects
// - check which conversions+V0 are identified as such


V0Performance aV0Performance ;


V0Performance::V0Performance() : Processor("V0Performance") {
  
  // modify processor description
  _description = "V0Performance processor does performance analysis of conversion and V0 tagging" ;
  
  // minimum number of SimTrackerHits that at least two conversion/V0 decay
  // particles have to have in order for the V0/conversion to be taken as
  // inside the tracking system
  registerOptionalParameter("minHits",
			    "minimum hit numbers for conv/V0 decay particles",
			    _minHits,-999);
  vector<string> allColls;
  allColls.clear();
  registerOptionalParameter("collectionsToPrint",
			    "names of LCCollections to print statistics for",
			    _CollectionsToPrint,allColls);

  V0Name[V0Gamma]="photon conversion";
  V0Name[V0K0]="Kshort";
  V0Name[V0Lambda]="Lambda";
}


void V0Performance::init() { 

  printParameters() ;

  // bookkeeping for efficiency/purity determination at track level
  for (int i=0; i<V0LastType; i++) {
    mc_num[i]=0;
    mc_tracks[i]=0;
  }
  num_tracks_total.clear();

  histos = new HistMap(this);
}


void V0Performance::processRunHeader( LCRunHeader* run) { 

} 

void V0Performance::processEvent( LCEvent * evt ) { 

  V0Candidates.clear();
  numHits.clear();

  const StringVec* colNames = evt->getCollectionNames();

  streamlog_out(DEBUG0) << "detector " << evt->getDetectorName() 
			 << ", run " << evt->getRunNumber()
			 << ", event " << evt->getEventNumber()
			 << endl;

  _maxCollNameLength=15;
  vector<string> mccols,trackcols,relationcols,recocols;
  for (size_t i=0; i<colNames->size(); i++) {
    string colname=(*colNames)[i];
    // check whether this collection is among the ones to investigate
    bool useThis=true;
    if (_CollectionsToPrint.size()>0) {
      useThis=false;
      for (vector<string>::iterator it=_CollectionsToPrint.begin();
	   it!=_CollectionsToPrint.end(); it++) {
	useThis|=((*it)==colname);
      }
    }
    LCCollection *collection=evt->getCollection(colname);
    if (collection->getTypeName()==LCIO::MCPARTICLE) {
      mccols.push_back(colname);
    } else if (collection->getTypeName()==LCIO::TRACK && useThis) {
      trackcols.push_back(colname);
      if (colname.length()>_maxCollNameLength)
	_maxCollNameLength=colname.length();
    } else if (collection->getTypeName()==LCIO::LCRELATION) {
      relationcols.push_back(colname);
      LCTOOLS::printParameters(evt->getParameters());
      LCRelationNavigator nav(collection);
      streamlog_out(DEBUG0) << "relation collection " << nav.getFromType() 
	   << "   " << nav.getToType() << " named " << colname << endl;
    } else if (collection->getTypeName()==LCIO::RECONSTRUCTEDPARTICLE && useThis) {
      recocols.push_back(colname);
      if (colname.length()>_maxCollNameLength)
	_maxCollNameLength=colname.length();
    }
  }
  streamlog_out(DEBUG0) << "collections: " << mccols.size() << " MCPARTICLE" << endl;
  for (size_t i=0; i<mccols.size(); i++) {
    streamlog_out(DEBUG0) << "                " << mccols[i] << ":  "
	 << evt->getCollection(mccols[i])->getNumberOfElements()
	 << " entries" << endl;
  }
  streamlog_out(DEBUG0) << "             " << trackcols.size() << " TRACK" << endl;
  for (size_t i=0; i<trackcols.size(); i++) {
    streamlog_out(DEBUG0) << "                " << trackcols[i] << ":  "
	 << evt->getCollection(trackcols[i])->getNumberOfElements()
	 << " entries" << endl;
  }
  streamlog_out(DEBUG0) << "             " << recocols.size() << " RECOPART" << endl;
  for (size_t i=0; i<recocols.size(); i++) {
    streamlog_out(DEBUG0) << "                " << recocols[i] << ":  "
	 << evt->getCollection(recocols[i])->getNumberOfElements()
	 << " entries" << endl;
  }
  streamlog_out(DEBUG0) << "             " << relationcols.size() << " LCRELATION" << endl;
  for (size_t i=0; i<relationcols.size(); i++) {
    streamlog_out(DEBUG0) << "                " << relationcols[i] << ":  "
	 << evt->getCollection(relationcols[i])->getNumberOfElements()
	 << " entries" << endl;
  }


  // provide access to all LCRelation collections
  lcRel.clear();
  for (size_t i=0; i<relationcols.size(); i++) {
    lcRel[relationcols[i]]= new LCRelationNavigator(evt->getCollection(relationcols[i]));
  }

  // get number of tracker hits for each MCParticle
  // (we will need this to define which V0/conversions happened well within
  // the tracking detectors)
  hitAnalysis(evt);


  // find MCParticle collection to find true conversions+V0s
  for (size_t i=0; i<mccols.size(); i++) {
    treeAnalysis(evt,mccols[i]);
  }

  // then, get tracking output to see what was reconstructed
  for (size_t i=0; i<trackcols.size(); i++) {
    trackAnalysis(evt,trackcols[i]);
  }

  // finally, get particle flow output to see what was identified as V0
  for (size_t i=0; i<recocols.size(); i++) {
    recoAnalysis(evt,recocols[i]);
  }


  for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
    dump_V0Candidate(V0Candidates[iv0]);
    delete V0Candidates[iv0];
  }

  for (size_t i=0; i<relationcols.size(); i++) {
    delete lcRel[relationcols[i]];
  }

  setReturnValue(true);

}


void V0Performance::treeAnalysis( const LCEvent *evt, const string collectionName ) {

  LCCollection *mcColl=evt->getCollection(collectionName);
  for (int imc=0; imc<mcColl->getNumberOfElements(); imc++) {
    EVENT::MCParticle *mcp=dynamic_cast<EVENT::MCParticle*>
      ( mcColl->getElementAt(imc) );

    int pdg = mcp->getPDG();
    const MCParticleVec &daughters = mcp->getDaughters();
    int numDaughters=daughters.size();
    double radius=0,z=0;
    if (numDaughters>=2) {
      const double* vtx = daughters[0]->getVertex();
      radius = sqrt(vtx[0]*vtx[0]+vtx[1]*vtx[1]);
      z = vtx[2];
    }
    if (abs(pdg)==22 && numDaughters==2 && mcp->isDecayedInTracker()
	&& (daughters[0]->getPDG()*daughters[1]->getPDG()==-11*11)
	&& numHits[daughters[0]]>=_minHits && numHits[daughters[1]]>=_minHits) {
      // this is a photon conversion in the tracking system
      V0Candidate_type *newconv = new V0Candidate_type();
      newconv->V0Type=V0Gamma;
      newconv->radius=radius;
      newconv->z=z;
      newconv->mother = mcp;
      newconv->daughters = daughters;
      newconv->tracks.clear();
      newconv->recopart.clear();
      newconv->numTrackerHits=0;
      for (size_t i=0; i<daughters.size(); i++)
	newconv->numTrackerHits+=numHits[daughters[i]];
      V0Candidates.push_back(newconv);
      ++mc_num[V0Gamma];
      mc_tracks[V0Gamma]+=2;
    } else if (mcp->isDecayedInTracker()
	       && numDaughters>=2 && (abs(pdg)==310 || abs(pdg)==3122)) {

      // long lived neutral particles
      int thisV0Type = V0K0;
      if (abs(pdg)!=310) thisV0Type=V0Lambda;
	
      // first check how many charged daughters there are
      int nCharged=0;
      double totalCharge=0;
      bool has_nucleon=false;
      for (int idaughter=0; idaughter<numDaughters; idaughter++) {
	if (daughters[idaughter]->getCharge()==-1000) {
	  streamlog_out(ERROR) << "Mokka failed to find charge for PDG "
			       << daughters[idaughter]->getPDG()
			       << "; assuming this to be neutral" << endl;
	} else if (daughters[idaughter]->getCharge()!=0
		   && numHits[daughters[idaughter]]>=_minHits) {
	  ++nCharged;
	  totalCharge+=daughters[idaughter]->getCharge();
	  has_nucleon|=(daughters[idaughter]->getPDG()==2212);
	}
      }
      if (totalCharge!=0) {
	if (has_nucleon && totalCharge!=0) {
	  // we are probably looking at a nuclear interaction.
	  // ignore this for now.
	} else {
	  streamlog_out(ERROR) << "charges of decay products do not sum up!"
			       << endl;
	  streamlog_out(ERROR) << "BAD THING in tracker: pdg=" << pdg
			       << ", numdaughters=" << numDaughters
			       << ", radius=" << radius
			       << ", z=" << z << ", total charge="
			       << totalCharge;
	  for (int i=0; i<numDaughters; i++) {
	    streamlog_out(ERROR) << ", PDG=" << daughters[i]->getPDG()
				   << "(" << daughters[i]->getCharge() << ")";
	  }
	  streamlog_out(ERROR) << endl;
	}
      }
      if (nCharged>=2 && totalCharge==0) {
	if (nCharged!=numDaughters) {
	  streamlog_out(WARNING) << V0Name[thisV0Type] << " has " << nCharged
				 << " charged and " << numDaughters-nCharged
				 << " neutral daughters in MC."
				 << " mass reconstruction will fail." << endl;
	}
	if (thisV0Type==V0Lambda) {
	  double pimom=0.0;
	  double protmom=0.0;
	  for (int i=0; i<numDaughters; i++) {
	    double mom2=
	      daughters[i]->getMomentum()[0]*daughters[i]->getMomentum()[0]
	      +daughters[i]->getMomentum()[1]*daughters[i]->getMomentum()[1]
	      +daughters[i]->getMomentum()[2]*daughters[i]->getMomentum()[2];
	    if (abs(daughters[i]->getPDG())==211) pimom=sqrt(mom2);
	    if (abs(daughters[i]->getPDG())==2212) protmom=sqrt(mom2);
	  }
	  histos->fill("lambda_daughter_mom_pi",pimom,
		       1,"momentum of pions from Lambda decay",100,0,20);
	  histos->fill("lambda_daughter_mom_proton",protmom,
		       1,"momentum of protons from Lambda decay",100,0,20);
	  histos->fill("lambda_daughters_mom",pimom,protmom,1,
		       "lambda pion vs proton momentum",100,0,20,100,0,20);
	}

	V0Candidate_type *newconv = new V0Candidate_type();
	newconv->V0Type=thisV0Type;
	newconv->radius=radius;
	newconv->z=z;
	newconv->mother = mcp;
	newconv->daughters = daughters;
	newconv->tracks.clear();
	newconv->recopart.clear();
	newconv->numTrackerHits=0;
	for (size_t i=0; i<daughters.size(); i++)
	  newconv->numTrackerHits+=numHits[daughters[i]];
	V0Candidates.push_back(newconv);
	++mc_num[thisV0Type];
	mc_tracks[thisV0Type]+=nCharged;
      }
    }
  }
  streamlog_out(DEBUG0) << "number of conv/V0 in this event: "
			 << V0Candidates.size() << endl;
}


void V0Performance::hitAnalysis( const LCEvent *evt ) {

  const StringVec* colNames = evt->getCollectionNames();
  for (size_t k=0; k<colNames->size(); k++) {
    LCCollection *collection=evt->getCollection((*colNames)[k]);
    if (collection->getTypeName()!=LCIO::SIMTRACKERHIT) continue;

    for (int i=0; i<collection->getNumberOfElements(); i++) {
      EVENT::SimTrackerHit* hit = dynamic_cast<EVENT::SimTrackerHit*>
	( collection->getElementAt(i) );
      ++numHits[hit->getMCParticle()];
    }
  }
  
}



void V0Performance::dump_V0Candidate(V0Candidate_type* cand) {

  // general information
  streamlog_out(MESSAGE) << V0Name[cand->V0Type] << " found: " << endl;
  streamlog_out(MESSAGE) << "radius=" << cand->radius << " mm, z=" << cand->z << " mm" << endl;
  streamlog_out(MESSAGE) << "daughter momenta:" << endl;
  for (size_t i=0; i<cand->daughters.size(); i++) {
    streamlog_out(MESSAGE)
      << "   momentum="
      << cand->daughters[i]->getMomentum()[0] << ", "
      << cand->daughters[i]->getMomentum()[1] << ", "
      << cand->daughters[i]->getMomentum()[2] << endl;
  }

  // MC decay products
  streamlog_out(MESSAGE) << "daughter particles:"
			 << setw(_maxCollNameLength-12) << " ";
  for (size_t i=0; i<cand->daughters.size(); i++) {
    if (cand->daughters[i]->getCharge()!=0
	&& cand->daughters[i]->getCharge()!=-1000) {
      streamlog_out(MESSAGE) << setw(14) << cand->daughters[i]->getPDG();
    } else {
      stringstream pid;
      pid << "(" << cand->daughters[i]->getPDG() << ")";
      streamlog_out(MESSAGE) << setw(14) << pid.str();
    }
  }
  streamlog_out(MESSAGE) << endl;

  // reconstructed tracks
  for (map<string,vector<Track*> >::iterator it=cand->tracks.begin();
       it!=cand->tracks.end(); it++) {
    streamlog_out(MESSAGE) << "tracks in " << setw(_maxCollNameLength)
			   << it->first << ":";
    for (size_t itrk=0; itrk<it->second.size(); itrk++) {
      if (it->second[itrk]) {
	streamlog_out(MESSAGE) << setw(7) << "yes";
      } else {
	streamlog_out(MESSAGE) << setw(7) << "no";
      }
      streamlog_out(MESSAGE) << " (" << setw(4)
	   << int(cand->track_weights[it->first][itrk]*100+0.5)/100.
	   << ")";
    }
    streamlog_out(MESSAGE) << endl;
  }

  // particle flow objects
  for (map<string,vector<ReconstructedParticle*> >::iterator
	 it=cand->recopart.begin(); it!=cand->recopart.end(); it++) {
    streamlog_out(MESSAGE) << "PFOs   in " << setw(_maxCollNameLength)
			   << it->first << ":";
    for (size_t itrk=0; itrk<it->second.size(); itrk++) {
      if (it->second[itrk]) {
	if (it->second[itrk]->getTracks().size()>1) {
	  streamlog_out(MESSAGE) << setw(7) << "comp";
	  streamlog_out(MESSAGE) << " (" << setw(4)
				 << it->second[itrk]->getType()
				 << ")";
	} else {
	  streamlog_out(MESSAGE) << setw(7) << "yes";
	  streamlog_out(MESSAGE) << " (" << setw(4)
				 << it->second[itrk]->getType()
				 << ")";
	}
      } else {
	streamlog_out(MESSAGE) << setw(14) << "no ( N/A)";
      }
    }
    streamlog_out(MESSAGE) << endl;
  }
  streamlog_out(MESSAGE) << endl;

  // histograms for this candidate
  if (cand->V0Type==V0Gamma) {
    histos->fill("conv_radius_mc",cand->radius,1,
		 "radius of MC conversion",100,0,2000);
    histos->fill("conv_z_mc",cand->z,1,"z of MC conversion",100,-2500,2500);
  } else {
    histos->fill("v0_radius_mc",cand->radius,1,"radius of MC V0",100,0,2000);
    histos->fill("v0_z_mc",cand->z,1,"z of MC V0",100,-2500,2500);
  }
  histos->fill("mass_mc",cand->mother->getMass(),1,"true mass of V0/conversion",
	       100,0,5);
}


void V0Performance::trackAnalysis( const LCEvent *evt,
				      const string collectionName) {

  LCCollection *trkColl=evt->getCollection(collectionName);
  num_tracks_total[collectionName]+=trkColl->getNumberOfElements();

  for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
    // initialize track vectors for this LCCollection if necessary
    if (!V0Candidates[iv0]->tracks[collectionName].size()) {
      for (size_t idght=0;
	   idght<V0Candidates[iv0]->daughters.size(); idght++) {
	V0Candidates[iv0]->tracks[collectionName].push_back(NULL);
	V0Candidates[iv0]->track_weights[collectionName].push_back(0);
      }
    }
  }

  for (int itrk=0; itrk<trkColl->getNumberOfElements(); itrk++) {
    EVENT::Track *trk=dynamic_cast<EVENT::Track*>(trkColl->getElementAt(itrk));

    // trace this track back to corresponding MCParticle.
    // since we do not make any assumptions about which LCRelation collection
    // connects where and which way, we just look everywhere.
    vector<MCParticle*> mcp;
    vector<float> mcp_weight;
    for (map<string,LCRelationNavigator*>::iterator it=lcRel.begin();
	 it!=lcRel.end(); it++) {
      LCRelationNavigator *relate = it->second;
      for (size_t i=0; i<relate->getRelatedToObjects(trk).size(); i++) {
	mcp.push_back(dynamic_cast<EVENT::MCParticle*>(relate->getRelatedToObjects(trk)[i]));
	mcp_weight.push_back(relate->getRelatedToWeights(trk)[i]);
      }
      for (size_t i=0; i<relate->getRelatedFromObjects(trk).size(); i++) {
	mcp.push_back(dynamic_cast<EVENT::MCParticle*>(relate->getRelatedFromObjects(trk)[i]));
	mcp_weight.push_back(relate->getRelatedFromWeights(trk)[i]);
      }
    }

    // now check whether we found any track associated with any of our
    // conversion or V0 candidates.
    // now, it might have been easier to follow the LCRelations from the
    // V0 candidates and just get matching tracks for each daughter. However,
    // how would we then find out which collection a track belongs to?
    for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
      for (size_t idght=0; idght<V0Candidates[iv0]->daughters.size(); idght++) {
	for (size_t imc=0; imc<mcp.size(); imc++) {
	  if (V0Candidates[iv0]->daughters[idght]==mcp[imc]) {
	    // we do have a matching track here. if there is more than one,
	    // store the one with highest weight in LCRelation to MCParticle
	    if (mcp_weight[imc]
		>=V0Candidates[iv0]->track_weights[collectionName][idght]) {
	      V0Candidates[iv0]->tracks[collectionName][idght]=trk;
	      V0Candidates[iv0]->track_weights[collectionName][idght]
		=mcp_weight[imc];
	    }
	  }
	}
      }
    }

  }

  // book-keeping: count how many times we did (not) do well
  for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
    size_t numReconstructed=0;
    for (size_t idght=0; idght<V0Candidates[iv0]->tracks[collectionName].size();
	 idght++) {
      if (V0Candidates[iv0]->tracks[collectionName][idght]) ++numReconstructed;
    }
    if (numReconstructed==0) {
      ++foundNOtracks[V0Candidates[iv0]->V0Type][collectionName];
    } else if (numReconstructed==1) {
      ++foundONEtrack[V0Candidates[iv0]->V0Type][collectionName];
    } else if (numReconstructed>=2) {
      ++foundTWOtracks[V0Candidates[iv0]->V0Type][collectionName];
    }
    if (numReconstructed==V0Candidates[iv0]->daughters.size()) {
      ++foundALLtracks[V0Candidates[iv0]->V0Type][collectionName];
    }
  }
}


void V0Performance::recoAnalysis( const LCEvent *evt, const string collectionName ) {

  LCCollection *recoColl=evt->getCollection(collectionName);

  for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
    // initialize recopart vectors for this LCCollection if necessary
    if (!V0Candidates[iv0]->recopart[collectionName].size()) {
      for (size_t idght=0;
	   idght<V0Candidates[iv0]->daughters.size(); idght++) {
	V0Candidates[iv0]->recopart[collectionName].push_back(NULL);
      }
    }
  }

  vector<vector<int> > matches;

  for (int ipart=0; ipart<recoColl->getNumberOfElements(); ipart++) {
    EVENT::ReconstructedParticle *part
      =dynamic_cast<EVENT::ReconstructedParticle*>
      (recoColl->getElementAt(ipart));

    if (!part) {
      // this is weird
      streamlog_out(ERROR) << "collection " << collectionName
			   << " contains NULL object at position "
			   << ipart << " of " << recoColl->getNumberOfElements()
			   << ". skipping." << endl;
      vector<int> dummy;
      matches.push_back(dummy);
      continue;
    }

    // is there more than one track associated with this object?
    if (part->getTracks().size()>0) {
      num_tracks_total[collectionName]+=part->getTracks().size();
    }
    if (part->getTracks().size()>1) {
      ++num_composites_total[collectionName];
      if (part->getType()==22) {
	++num_composites_by_type[V0Gamma][collectionName];
      } else if (part->getType()==130 || part->getType()==310) {
	++num_composites_by_type[V0K0][collectionName];
      } else if (part->getType()==3122) {
	++num_composites_by_type[V0Lambda][collectionName];
      }
      //streamlog_out(DEBUG0) << "collection " << collectionName
      //		     << " contains object of type " << part->getType()
      //		     << " with " << part->getTracks().size()
      //		     << " tracks" << endl;
    }


    // for all tracks associated with this particle, check whether they
    // point back to any track used in our V0 candidates
    vector<int> matches_this_recopart;
    for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
      int matches_this_v0=0;
      for (size_t itrk=0; itrk<part->getTracks().size(); itrk++) {
	for (map<string,vector<Track*> >::iterator it
	       =V0Candidates[iv0]->tracks.begin();
	     it!=V0Candidates[iv0]->tracks.end(); it++) {
	  for (size_t idght=0; idght<it->second.size(); idght++) {
	    if (it->second[idght]==part->getTracks()[itrk]) {
	      // we do have a matching track here.
	      ++matches_this_v0;
	      if (V0Candidates[iv0]->recopart[collectionName][idght]) {
                // we have more than one match. since we store only the
		// first match, some information about whether or not
		// a composite object was compiled properly will be lost!
		streamlog_out(ERROR) << "more than one matching recoparticle!" << endl;
	      } else {
		V0Candidates[iv0]->recopart[collectionName][idght]=part;
	      }
	    }
	  }
	}
      }
      matches_this_recopart.push_back(matches_this_v0);
    }
    matches.push_back(matches_this_recopart);
  }

  // book-keeping: count how many times we did (not) do well
  for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {

    // check how many daughters appeared in any recoparticle at all
    size_t numReconstructed=0;
    for (size_t idght=0; idght<V0Candidates[iv0]->recopart[collectionName].size();
	 idght++) {
      if (V0Candidates[iv0]->recopart[collectionName][idght])
	++numReconstructed;
    }

    // check whether there were composite particles with >1 daughter particle
    int max_tracks_in_single_recopart=0;
    for (int ipart=0; ipart<recoColl->getNumberOfElements(); ipart++) {
      if (matches[ipart][iv0]>max_tracks_in_single_recopart)
	max_tracks_in_single_recopart=matches[ipart][iv0];
    }

    if (numReconstructed==0) {
      ++foundNOtracks[V0Candidates[iv0]->V0Type][collectionName];
    } else if (numReconstructed==1) {
      ++foundONEtrack[V0Candidates[iv0]->V0Type][collectionName];
    } else if (numReconstructed>=2) {
      ++foundTWOtracks[V0Candidates[iv0]->V0Type][collectionName];
      if (max_tracks_in_single_recopart>=2) {
	++foundComposite[V0Candidates[iv0]->V0Type][collectionName];
      }
    }
    if (numReconstructed==V0Candidates[iv0]->daughters.size()) {
      ++foundALLtracks[V0Candidates[iv0]->V0Type][collectionName];
    }
  }

}


void V0Performance::check( LCEvent * evt ) { 

}


void V0Performance::end(){ 

  streamlog_out(MESSAGE) << "PERFORMANCE SUMMARY:" << endl;
  double ntrue_tracks=0;
  for (int i=0; i<V0LastType; i++) {
    streamlog_out(MESSAGE) << " actual " << V0Name[i] << " in MCParticles: "
			 << mc_num[i] << " with "
			 << mc_tracks[i] << " tracks" << endl;
    ntrue_tracks+=mc_tracks[i];
  }
  for (map<string,int>::iterator it=num_tracks_total.begin();
       it!=num_tracks_total.end(); it++) {
    for (int i=0; i<V0LastType; i++) {
      if (mc_num[i]==0) continue;
      streamlog_out(MESSAGE) << " number of MC " <<  V0Name[i]
			     << " with no reconstructed tracks in "
			     << it->first << ": " << foundNOtracks[i][it->first]
			     << " of " << mc_num[i] << " ("
			     << int(0.5+1000.*foundNOtracks[i][it->first]
				 /mc_num[i])/10. << "%)" << endl;
      streamlog_out(MESSAGE) << " number of MC " <<  V0Name[i]
			     << " with one reconstructed track  in "
			     << it->first << ": " << foundONEtrack[i][it->first]
			     << " of " << mc_num[i] << " ("
			     << int(0.5+1000.*foundONEtrack[i][it->first]
				 /mc_num[i])/10. << "%)" << endl;
      streamlog_out(MESSAGE) << " number of MC " <<  V0Name[i]
			     << " with >=2 reconstructed tracks in "
			     << it->first << ": " << foundTWOtracks[i][it->first]
			     << " of " << mc_num[i] << " ("
			     << int(0.5+1000.*foundTWOtracks[i][it->first]
				 /mc_num[i])/10. << "%)" << endl;
      streamlog_out(MESSAGE) << " number of MC " <<  V0Name[i]
			     << " with all daughters reconstructed in "
			     << it->first << ": " << foundALLtracks[i][it->first]
			     << " of " << mc_num[i] << " ("
			     << int(0.5+1000.*foundALLtracks[i][it->first]
				 /mc_num[i])/10. << "%)" << endl;
      streamlog_out(MESSAGE) << " number of MC " <<  V0Name[i]
			     << " reconstructed as composite in "
			     << it->first << ": " << foundComposite[i][it->first]
			     << " of " << mc_num[i] << " ("
			     << int(0.5+1000.*foundComposite[i][it->first]
				 /mc_num[i])/10. << "%)" << endl;
      streamlog_out(MESSAGE) << " total number of " << V0Name[i]
			     << " candidates with >=2 tracks in "
			     << it->first << ": "
			     << num_composites_by_type[i][it->first] << endl;
    }
    streamlog_out(MESSAGE) << " total number of composite objects"
			   << " with >=2 tracks in "
			   << it->first << ": "
			   << num_composites_total[it->first] << endl;
  }
}
