
//brief DY anlalysis module based on picoDST, generate femtoDST contains NTuple and histos


void Fun4Drellyan(
		int nEvents = 0, //
		char *input_file = "368042_0.root", //
		std::string out_name = "fpdsts", 
		bool is_sim = false,
		int run_number = 395402
		)
{

	gSystem->Load("libMWGOO");
	gSystem->Load("libfvtx_subsysreco.so");
	gSystem->Load("librpc_subsysreco");

	//gSystem->Load("/phenix/spin/phnxsp01/yuhw/install/lib/libpicodst_object.so");
	gSystem->Load("libpicodst_object.so");
	gSystem->Load("/direct/phenix+spin2/darshana/TestLibrary/test/install/lib/libdnpfilter.so");

	recoConsts *rc = recoConsts::instance();
	rc->set_IntFlag("RUNNUMBER", run_number);

	Fun4AllServer *se = Fun4AllServer::instance();
	se->Verbosity(0);

	se->registerSubsystem(new FvtxReadbackDST());
	se->registerSubsystem(new RpcReadbackDST());


	dnpfilter *mdyevt = new dnpfilter(out_name);
	//mdyevt->Verbosity(10);
	mdyevt->set_is_sim(is_sim);
	se->registerSubsystem(mdyevt);


	Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin");
	in->fileopen(input_file);
	//in->AddListFile(input_file);
	se->registerInputManager(in);

	se->run(nEvents);
	se->End();

	//se->dumpHistos("se.root");

	delete se;

	cout << "Completed reconstruction." << endl;

}
