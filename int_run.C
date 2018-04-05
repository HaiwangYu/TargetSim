void int_run() {
	gROOT->ProcessLine(".x Fun4All_G4_BNL_R1.C\(-1\)");
	gROOT->ProcessLine(".L DisplayOn.C");
	PHG4Reco* g4 = DisplayOn();
	g4->ApplyCommand("/vis/viewer/set/background white");
	g4->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 270 0");
	Fun4AllServer *se = Fun4AllServer::instance();
	se->run(1);
	se->End();
}
