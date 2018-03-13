/*
 * TruthEval.C
 *
 *  Created on: Oct 29, 2017
 *      Author: yuhw@nmsu.edu
 */


#include "TruthEval.h"
#include "TruthTrack.h"

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
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

#define LogInfo(exp)		    std::cout<<"INFO: "   <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogDebug(exp)		    std::cout<<"DEBUG: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogError(exp)		    std::cout<<"ERROR: "  <<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl
#define LogWarning(exp)	    std::cout<<"WARNING: "<<__FILE__<<": "<<__LINE__<<": "<< exp << std::endl

using namespace std;

TruthEval::TruthEval(const std::string& name) :
SubsysReco(name),
_event(0),
_g4truth_container(nullptr),
_out_name("eval.root")
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
	return Fun4AllReturnCodes::EVENT_OK;
}







