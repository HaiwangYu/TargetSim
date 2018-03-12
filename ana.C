void ana(
		const char* in = "G4.root",
		const char* out = "hist.root"
		){

	TFile *outf = TFile::Open(out, "recreate");

	TFile *inf = TFile::Open(in, "read");
	TTree *T = (TTree*) inf->Get("T");


	T->SetAlias("x", "0.5*(G4HIT_Coil.x[][0]+G4HIT_Coil.x[][1])");
	T->SetAlias("y", "0.5*(G4HIT_Coil.y[][0]+G4HIT_Coil.y[][1])");
	T->SetAlias("z", "0.5*(G4HIT_Coil.z[][0]+G4HIT_Coil.z[][1])");
	T->SetAlias("e", "G4HIT_Coil.edep");


	TH2D *h2d_ezy = new TH2D("h2d_ezy","e [GeV]; z [cm]; y [cm]", 100, -25, 25, 100, -25, 25);
	T->Project("h2d_ezy", "y:z", "e");
	outf->cd();
	h2d_ezy->Write();

	double ybin[51];
	for(int i=0;i<51;++i) ybin[i] = i - 25.;

	//for(int i=27;i<28;++i) {
	for(int i=0;i<50;++i) {
		const char* hname = Form("h2d_ezx_%02d",i);
		TH2D *h2d_ezx = new TH2D(hname,
				Form("edep [GeV], %02.0f<y<%02.0f cm; z [cm]; x [cm]", ybin[i], ybin[i+1]),
				100, -25, 25, 100, -25, 25);
		T->Project(hname, "x:z" , Form("e&&(y>%f&&y<%f)",ybin[i],ybin[i+1]));
		outf->cd();
		h2d_ezx->Write();
	}

	outf->Close();
}
