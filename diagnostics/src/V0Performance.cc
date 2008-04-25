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
  
}


void V0Performance::init() { 

  // minimum number of SimTrackerHits that at least two conversion/V0 decay
  // particles have to have in order for the V0/conversion to be taken as
  // inside the tracking system
  minHits=-999;

  printParameters() ;

  // bookkeeping for efficiency/purity determination at track level
  total_num_v0=0;
  ntrue_tracks=0;
  ncand_tracks.clear();
  ngood_tracks.clear();
  ngood_tracks_in_composites.clear();
  ntracks_in_composites.clear();

  histos = new HistMap(this);
}


void V0Performance::processRunHeader( LCRunHeader* run) { 

} 

void V0Performance::processEvent( LCEvent * evt ) { 

  V0Candidates.clear();

  const StringVec* colNames = evt->getCollectionNames();

  streamlog_out(MESSAGE) << "detector " << evt->getDetectorName() 
			 << ", run " << evt->getRunNumber()
			 << ", event " << evt->getEventNumber()
			 << endl;

  vector<string> mccols,trackcols,relationcols,recocols;
  for (size_t i=0; i<colNames->size(); i++) {
    string colname=(*colNames)[i];
    LCCollection *collection=evt->getCollection(colname);
    if (collection->getTypeName()==LCIO::MCPARTICLE) {
      mccols.push_back(colname);
    } else if (collection->getTypeName()==LCIO::TRACK) {
      trackcols.push_back(colname);
    } else if (collection->getTypeName()==LCIO::LCRELATION) {
      relationcols.push_back(colname);
      LCTOOLS::printParameters(evt->getParameters());
      LCRelationNavigator nav(collection);
      streamlog_out(MESSAGE) << "relation collection " << nav.getFromType() 
	   << "   " << nav.getToType() << " named " << colname << endl;
    } else if (collection->getTypeName()==LCIO::RECONSTRUCTEDPARTICLE) {
      recocols.push_back(colname);
    }
  }
  streamlog_out(MESSAGE) << "collections: " << mccols.size() << " MCPARTICLE" << endl;
  for (size_t i=0; i<mccols.size(); i++) {
    streamlog_out(MESSAGE) << "                " << mccols[i] << ":  "
	 << evt->getCollection(mccols[i])->getNumberOfElements()
	 << " entries" << endl;
  }
  streamlog_out(MESSAGE) << "             " << trackcols.size() << " TRACK" << endl;
  for (size_t i=0; i<trackcols.size(); i++) {
    streamlog_out(MESSAGE) << "                " << trackcols[i] << ":  "
	 << evt->getCollection(trackcols[i])->getNumberOfElements()
	 << " entries" << endl;
  }
  streamlog_out(MESSAGE) << "             " << recocols.size() << " RECOPART" << endl;
  for (size_t i=0; i<recocols.size(); i++) {
    streamlog_out(MESSAGE) << "                " << recocols[i] << ":  "
	 << evt->getCollection(recocols[i])->getNumberOfElements()
	 << " entries" << endl;
  }
  streamlog_out(MESSAGE) << "             " << relationcols.size() << " LCRELATION" << endl;
  for (size_t i=0; i<relationcols.size(); i++) {
    streamlog_out(MESSAGE) << "                " << relationcols[i] << ":  "
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
  for (size_t i=0; i<colNames->size(); i++) {
    LCCollection *collection=evt->getCollection((*colNames)[i]);
    if (collection->getTypeName()==LCIO::MCPARTICLE) {
      treeAnalysis(evt,(*colNames)[i]);
    }
  }

  // then, get tracking output to see what was reconstructed
  for (size_t i=0; i<colNames->size(); i++) {
    LCCollection *collection=evt->getCollection((*colNames)[i]);
    if (collection->getTypeName()==LCIO::TRACK) {
      trackAnalysis(evt,(*colNames)[i]);
    }
  }

  // relation between tracks and MC particles
  for (size_t i=0; i<colNames->size(); i++) {
    LCCollection *collection=evt->getCollection((*colNames)[i]);
    if (collection->getTypeName()==LCIO::LCRELATION) {
      
    }
  }

  // finally, get particle flow output to see what was identified as V0
  for (size_t i=0; i<colNames->size(); i++) {
    LCCollection *collection=evt->getCollection((*colNames)[i]);
    if (collection->getTypeName()==LCIO::RECONSTRUCTEDPARTICLE) {
      recoAnalysis(evt,(*colNames)[i]);
    }
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
    bool not_from_IP=false;
    if (numDaughters>=2) {
      const double* vtx = daughters[0]->getVertex();
      radius = sqrt(vtx[0]*vtx[0]+vtx[1]*vtx[1]);
      z = vtx[2];
      // now we assume that the IP is always at 0,0,0. for future event
      // generation with IP smearing at generator level this might not
      // work.
      not_from_IP=(radius>0.001) || (fabs(z)>0.001); // is 1 micron ok?
    }
    if (abs(pdg)==22 && numDaughters==2 && mcp->isDecayedInTracker()
	&& (daughters[0]->getPDG()*daughters[1]->getPDG()==-11*11)
	&& numHits[daughters[0]]>=minHits && numHits[daughters[1]]>=minHits) {
      // this is a photon conversion in the tracking system
      V0Candidate_type *newconv = new V0Candidate_type();
      newconv->is_conversion=true;
      newconv->radius=radius;
      newconv->z=z;
      newconv->mother = mcp;
      newconv->daughters = daughters;
      newconv->tracks.clear();
      newconv->recopart.clear();
      V0Candidates.push_back(newconv);
      ntrue_tracks+=2;
      cout << "CONVERSION: number of hits from daughter particles: "
	   << numHits[daughters[0]] << ", " << numHits[daughters[1]] << endl;
    } else if (not_from_IP && mcp->isDecayedInTracker()
	       && numDaughters>=2 && (abs(pdg)==130 || abs(pdg)==3122)) {
      // long lived neutral particles
      // first check how many charged daughters there are
      int nCharged=0;
      double totalCharge=0;
      bool has_proton=false;
      for (int idaughter=0; idaughter<numDaughters; idaughter++) {
	if (daughters[idaughter]->getCharge()==-1000) {
	  streamlog_out(ERROR) << "Mokka failed to find charge for PDG "
			       << daughters[idaughter]->getPDG()
			       << "; assuming this to be neutral" << endl;
	} else if (daughters[idaughter]->getCharge()!=0
		   && numHits[daughters[idaughter]]>=minHits) {
	  ++nCharged;
	  totalCharge+=daughters[idaughter]->getCharge();
	  has_proton|=(daughters[idaughter]->getPDG()==2212);
	}
      }
      if (totalCharge!=0) {
	if (has_proton && totalCharge==1) {
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
      if (nCharged>=2) {
	V0Candidate_type *newconv = new V0Candidate_type();
	newconv->is_conversion=false;
	newconv->radius=radius;
	newconv->z=z;
	newconv->mother = mcp;
	newconv->daughters = daughters;
	newconv->tracks.clear();
	newconv->recopart.clear();
	V0Candidates.push_back(newconv);
	ntrue_tracks+=nCharged;
      }
    }
  }
  streamlog_out(MESSAGE) << "number of conv/V0 in this event: "
			 << V0Candidates.size() << endl;
  total_num_v0+=V0Candidates.size();
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
  if (cand->is_conversion) {
    streamlog_out(MESSAGE) << "photon conversion found: " << endl;
  } else {
    streamlog_out(MESSAGE) << "V0 found (PDG " << cand->mother->getPDG() << "):" << endl;
  }
  streamlog_out(MESSAGE) << "radius=" << cand->radius << " mm, z=" << cand->z << " mm" << endl;
  streamlog_out(MESSAGE) << "daughter momenta:" << endl;
  for (size_t i=0; i<cand->daughters.size(); i++) {
    streamlog_out(MESSAGE)
      << "   mass=" << cand->daughters[i]->getMass() << ", momentum="
      << cand->daughters[i]->getMomentum()[0] << ", "
      << cand->daughters[i]->getMomentum()[1] << ", "
      << cand->daughters[i]->getMomentum()[2] << endl;
  }

  // MC decay products
  streamlog_out(MESSAGE) << "daughter particles:";
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
    streamlog_out(MESSAGE) << "tracks in " << setw(12) << it->first << ":";
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
    streamlog_out(MESSAGE) << "PFOs   in " << setw(12) << it->first << ":";
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
    streamlog_out(MESSAGE) << endl << endl;
  }

  // histograms for this candidate
  if (cand->is_conversion) {
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

  ncand_tracks[collectionName]+=trkColl->getNumberOfElements();
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
      // initialize track vectors for this LCCollection if necessary
      if (!V0Candidates[iv0]->tracks[collectionName].size()) {
	for (size_t idght=0;
	     idght<V0Candidates[iv0]->daughters.size(); idght++) {
	  V0Candidates[iv0]->tracks[collectionName].push_back(NULL);
	  V0Candidates[iv0]->track_weights[collectionName].push_back(0);
	}
      }
      for (size_t idght=0; idght<V0Candidates[iv0]->daughters.size(); idght++) {
	for (size_t imc=0; imc<mcp.size(); imc++) {
	  if (V0Candidates[iv0]->daughters[idght]==mcp[imc]) {
	    // if there was no previous track match, count this one as found now
	    if (V0Candidates[iv0]->track_weights[collectionName][idght]==0) {
	      ++ngood_tracks[collectionName];
	    }
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
}


void V0Performance::recoAnalysis( const LCEvent *evt, const string collectionName ) {

  LCCollection *recoColl=evt->getCollection(collectionName);

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
      continue;
    }

    // is there more than one track associated with this object?
    if (part->getTracks().size()>0) {
      ncand_tracks[collectionName]+=part->getTracks().size();
    }
    if (part->getTracks().size()>1) {
      ntracks_in_composites[collectionName]+=part->getTracks().size();
      streamlog_out(MESSAGE) << "collection " << collectionName
			     << " contains object of type " << part->getType()
			     << " with " << part->getTracks().size()
			     << " tracks" << endl;
    }


    // for all tracks associated with this particle, check whether they
    // point back to any track used in our V0 candidates
    for (size_t itrk=0; itrk<part->getTracks().size(); itrk++) {
      for (size_t iv0=0; iv0<V0Candidates.size(); iv0++) {
	// initialize recopart vectors for this LCCollection if necessary
	if (!V0Candidates[iv0]->recopart[collectionName].size()) {
	  for (size_t idght=0;
	       idght<V0Candidates[iv0]->daughters.size(); idght++) {
	    V0Candidates[iv0]->recopart[collectionName].push_back(NULL);
	  }
	}
	for (map<string,vector<Track*> >::iterator it
	       =V0Candidates[iv0]->tracks.begin();
	     it!=V0Candidates[iv0]->tracks.end(); it++) {
	  for (size_t idght=0; idght<it->second.size(); idght++) {
	    if (it->second[idght]==part->getTracks()[itrk]) {
	      // we do have a matching track here.
	      if (V0Candidates[iv0]->recopart[collectionName][idght]) {
		streamlog_out(ERROR) << "more than one matching recoparticle!" << endl;
	      } else {
		V0Candidates[iv0]->recopart[collectionName][idght]=part;
		++ngood_tracks[collectionName];
		if (part->getTracks().size()>1) {
		  ++ngood_tracks_in_composites[collectionName];
		}
	      }
	    }
	  }
	}
      }
    }
  }
}


void V0Performance::check( LCEvent * evt ) { 

}


void V0Performance::end(){ 

  streamlog_out(MESSAGE) << "PERFORMANCE SUMMARY:" << endl;
  streamlog_out(MESSAGE) << " actual conv/V0 in MCParticles: " << total_num_v0
			 << " with " << ntrue_tracks << " tracks" << endl;
  for (map<string,int>::iterator it=ncand_tracks.begin();
       it!=ncand_tracks.end(); it++) {
    streamlog_out(MESSAGE) << " " << it->first << " collection contains "
			   << it->second << " tracks, including "
			   << ngood_tracks[it->first]
			   << " ones actually coming from V0/conversion"
			   << endl;
    if (ntrue_tracks>0) {
      streamlog_out(MESSAGE) << "    ==> efficiency "
			     << ngood_tracks[it->first]/double(ntrue_tracks)
			     << endl;
    }
    if (it->second>0) {
      streamlog_out(MESSAGE) << "    ==> purity     "
			     << ngood_tracks[it->first]/double(it->second)
			     << endl;
    }      

    streamlog_out(MESSAGE) << " " << it->first << " collection contains "
			   << ngood_tracks_in_composites[it->first]
			   << " correct tracks in composite objects"
			   << " and " << ntracks_in_composites[it->first]
			   << " tracks in composite objects with >1 track"
			   << endl;
    if (ntrue_tracks>0) {
      streamlog_out(MESSAGE) << "    ==> composite efficiency "
			     << ngood_tracks_in_composites[it->first]
	/double(ntrue_tracks) << endl;
    }
    if (ntracks_in_composites[it->first]>0) {
      streamlog_out(MESSAGE) << "    ==> composite purity     "
			     << ngood_tracks_in_composites[it->first]
	/double(ntracks_in_composites[it->first])
			     << endl;
    }      
  }
}

