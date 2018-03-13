/*
 * TruthEval.C
 *
 *  Created on: Oct 29, 2017
 *      Author: yuhw@nmsu.edu
 */


#include "TruthEval.h"
#include "TruthTrack.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>


#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <limits>
#include <boost/lexical_cast.hpp>
#include <map>

#define LogInfo(exp)		    std::cout<<"INFO: "   <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogDebug(exp)		    std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		    std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	    std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

TruthEval::TruthEval(const std::string& name, const std::string &out) :
SubsysReco(name),
_event(0),
_g4truth_container(nullptr),
_out_name(out)
{
	ResetEvalVars();
	InitEvalTree();
}

int TruthEval::Init(PHCompositeNode* topNode) {
	return Fun4AllReturnCodes::EVENT_OK;
}

int TruthEval::InitRun(PHCompositeNode* topNode) {

	int ret = GetNodes(topNode);
	if(ret != Fun4AllReturnCodes::EVENT_OK) return ret;

	return Fun4AllReturnCodes::EVENT_OK;
}

int TruthEval::process_event(PHCompositeNode* topNode) {

	if(Verbosity() >= Fun4AllBase::VERBOSITY_SOME)
		std::cout << "Entering TruthEval::process_event: " << _event << std::endl;

	ResetEvalVars();

	multimap<int, float> _m_track_edp;
	for(auto iter=_g4hit_coil->getHits().first; iter!=_g4hit_coil->getHits().second; ++iter) {
		PHG4Hit* hit = iter->second;
		int track_id = hit->get_trkid();
		float edep = hit->get_edep();
		auto iterm = _m_track_edp.find(track_id);
		if( iterm != _m_track_edp.end()) {
			iterm->second += edep;
		} else {
			_m_track_edp.insert(multimap<int, float>::value_type(track_id, edep));
		}
	}

	int iarr = 0;
	auto range = _g4truth_container->GetParticleRange();
	for(auto iter = range.first; iter!= range.second; ++iter) {
		PHG4Particle *particle =  iter->second;

		int track_id = particle->get_track_id()>0;
		if(track_id) continue; ///input proton

		TruthTrack track;
		track.pid = particle->get_pid();
		track.px = particle->get_px();
		track.py = particle->get_py();
		track.pz = particle->get_pz();
		track.e = particle->get_e();
		int vtx_id = particle->get_vtx_id();

		PHG4VtxPoint *vtx = _g4truth_container->GetVtx(vtx_id);

		track.vx = vtx->get_x();
		track.vy = vtx->get_y();
		track.vz = vtx->get_z();
		track.t = vtx->get_t();

		if(sqrt(track.vx*track.vx+track.vy*track.vy) < 1 and abs(track.vz) < 3.95) track.det_id = 0;
		else track.det_id = 1;

		track.total_edep_in_coil = 0;
		auto iterm = _m_track_edp.find(track_id);
		if(iterm != _m_track_edp.end()) {
			track.total_edep_in_coil = iterm->second;
		}

		new ((*_tca_truthtracks)[iarr++]) TruthTrack(track);
	}

	_tout->Fill();

	if(Verbosity() >= Fun4AllBase::VERBOSITY_SOME)
		std::cout << "Leaving TruthEval::process_event: " << _event << std::endl;
	++_event;

	return Fun4AllReturnCodes::EVENT_OK;
}

int TruthEval::End(PHCompositeNode* topNode) {
	if(Verbosity() >= Fun4AllBase::VERBOSITY_SOME)
		std::cout << "TruthEval::End" << std::endl;

	PHTFileServer::get().cd(_out_name.c_str());
	_tout->Write();

	return Fun4AllReturnCodes::EVENT_OK;
}

int TruthEval::InitEvalTree() {
	PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

	if (!_tca_truthtracks)
		_tca_truthtracks = new TClonesArray("TruthTrack");

	_tout = new TTree("T", "TruthEval");
	_tout->Branch("TruthParticle", _tca_truthtracks);
	return 0;
}

int TruthEval::ResetEvalVars() {
	_tca_truthtracks->Clear();
	return 0;
}

int TruthEval::GetNodes(PHCompositeNode* topNode) {

	_g4truth_container = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
	if (!_g4truth_container) {
		cerr << PHWHERE << " ERROR: Can't find node G4TruthInfo" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}

	_g4hit_coil = findNode::getClass<PHG4HitContainer>(topNode,"G4HIT_Coil");
	if (!_g4hit_coil) {
		cerr << PHWHERE << " ERROR: Can't find node G4HIT_Coil" << endl;
		return Fun4AllReturnCodes::ABORTEVENT;
	}


	return Fun4AllReturnCodes::EVENT_OK;
}







