/*
 * TruthEval.h
 *
 *  Created on: Oct 29, 2017
 *      Author: yuhw@nmsu.edu
 */

#ifndef _H_TruthEval_H_
#define _H_TruthEval_H_

// Fun4All includes
#include <fun4all/SubsysReco.h>

// STL includes
#include <vector>
#include <string>
#include <iostream>
#include <list>
#include <map>

class PHG4TruthInfoContainer;
class PHG4HitContainer;

class TFile;
class TTree;

class TClonesArray;

class TruthEval: public SubsysReco {

public:

	TruthEval(const std::string &name = "TruthEval");
	virtual ~TruthEval() {
	}

	int Init(PHCompositeNode *topNode);
	int InitRun(PHCompositeNode *topNode);
	int process_event(PHCompositeNode *topNode);
	int End(PHCompositeNode *topNode);

	int InitEvalTree();
	int ResetEvalVars();

	const std::string& get_out_name() const {
		return _out_name;
	}

	void set_out_name(const std::string& outName) {
		_out_name = outName;
	}

private:

	int GetNodes(PHCompositeNode *topNode);

	int _event;


	PHG4TruthInfoContainer *_g4truth_container;

	PHG4HitContainer *_g4hit_coil;

	std::string _out_name;

	TTree *_tout;
	TClonesArray *_tca_truthtracks;

};


#endif /* _H_TruthEval_H_ */
