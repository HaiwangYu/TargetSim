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
	InitEvalTree();
	ResetEvalVars();
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

	map<int, float> _m_track_edep;
	map<int, float> _m_track_path;

	for(auto iter=_g4hit_coil->getHits().first; iter!=_g4hit_coil->getHits().second; ++iter) {
		PHG4Hit* hit = iter->second;
		int track_id = hit->get_trkid();

		float edep = hit->get_edep();

		_total_edep += edep;
		
		float x0 = hit->get_x(0);
		float x1 = hit->get_x(1);
		float y0 = hit->get_y(0);
		float y1 = hit->get_y(1);
		float z0 = hit->get_z(0);
		float z1 = hit->get_z(1);

		float path = sqrt(
				pow(x1-x0,2)+
				pow(y1-y0,2)+
				pow(z1-z0,2)
				);

		auto iter_edep = _m_track_edep.find(track_id);
		if( iter_edep != _m_track_edep.end()) {
			iter_edep->second += edep;
		} else {
			_m_track_edep.insert(multimap<int, float>::value_type(track_id, edep));
		}

		auto iter_path = _m_track_path.find(track_id);
		if( iter_path != _m_track_path.end()) {
			iter_path->second += path;
		} else {
			_m_track_path.insert(multimap<int, float>::value_type(track_id, path));
		}

		if(verbosity > 2) {
			cout << __FILE__ << ": " << __LINE__ << endl;
			hit->identify();
			cout
			<< " edep: " << edep
			<< " path: " << path
			<< endl;

			cout
			<< "_m_track_path: "
			<< "{ " << track_id
			<< " -> " << _m_track_edep[track_id] << "}"
			<< endl;

			cout
			<< "_m_track_path: "
			<< "{ " <<track_id
			<< " -> " << _m_track_path[track_id] << "}"
			<< endl;
		}
	}

	int iarr = 0;
	auto range = _g4truth_container->GetParticleRange();
	for(auto iter = range.first; iter!= range.second; ++iter) {
		PHG4Particle *particle =  iter->second;

		int track_id = particle->get_track_id();
		if(track_id>0) continue; ///input proton

		TruthTrack track;
		track.parentid = particle->get_parent_id();
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


		float target_r = 1;
		float target_z = 3.95;

		float coil_in_r = 6;
		float coil_ot_r = 22.225;
	 	float coil_min_y = 2;
		float coil_max_y = 24.7;	

		if(
				sqrt(track.vx*track.vx+track.vy*track.vy) < target_r
				and abs(track.vz) < target_z
				) track.det_id = 0;
		else if(
				sqrt(track.vx*track.vx+track.vz*track.vz) > coil_in_r and
				sqrt(track.vx*track.vx+track.vz*track.vz) < coil_ot_r and
				abs(track.vy) > coil_min_y and
				abs(track.vy) < coil_max_y
				) track.det_id = 1;
		else track.det_id = 9999;

		track.edep_coil = 0;
		auto iter_edep = _m_track_edep.find(track_id);
		if(iter_edep != _m_track_edep.end()) {
			track.edep_coil = iter_edep->second;
		}

		track.path_coil = 0;
		auto iter_path = _m_track_path.find(track_id);
		if(iter_path != _m_track_path.end()) {
			track.path_coil = iter_path->second;
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

	_tout1->Fill();
	_tout1->Write();

	return Fun4AllReturnCodes::EVENT_OK;
}

int TruthEval::InitEvalTree() {
	PHTFileServer::get().open(_out_name.c_str(), "RECREATE");

	if (!_tca_truthtracks)
		_tca_truthtracks = new TClonesArray("TruthTrack");

	_tout = new TTree("T", "TruthEval");
	_tout->Branch("TruthParticle", _tca_truthtracks);

	_total_edep = 0;
	_tout1 = new TTree("T1", "RunLevel");
	_tout1->Branch("total_edep", &_total_edep, "total_edep/F");

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







