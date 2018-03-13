/*
 * TruthEval.h
 *
 *  Created on: Oct 29, 2017
 *      Author: yuhw@nmsu.edu
 */

#ifndef _H_TruthTrack_H_
#define _H_TruthTrack_H_

#include <TObject.h>

class TruthTrack : public TObject{
public:
	TruthTrack () :
		pid(0),
		vx(-9999), vy(-9999), vz(-9999),
		px(-9999), py(-9999), pz(-9999),
		edep_in_coil(-9999){}
	int pid;
	double vx, vy, vz, px, py, pz;
	double edep_in_coil;

	ClassDef(TruthTrack, 1)
};

#endif
