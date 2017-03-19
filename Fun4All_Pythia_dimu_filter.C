// Read in a DST and save a subset of nodes or event out
// Jin and Darshana March_2015

void
Fun4All_Pythia_dimu_filter(
    int nEvents = 1000, //
    char *input_file = "phpythia_1.root",
    const int nSkip = 0,
    char *dst_file = "phpythia.root")
{

  // load libraries
  gSystem->Load("libfvtx_subsysreco.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libmutoo_subsysreco");
  gSystem->Load("libfun4allfuncs_muons");
  gSystem->Load("libMWGOO");
  gSystem->Load("libmutrg");
  gSystem->Load("librpc_subsysreco");
  gSystem->Load("librpc_muotrackreco");
  gSystem->Load("librecal");
  //gSystem->Load("libg4analysis.so");
  gSystem->Load("libpicodst_object");


  gSystem->Load("/direct/phenix+spin2/darshana/Libraries/event_gen/src/PHPythia/install/lib/libPHPythiaEventGen.so");
  gSystem->Load("/direct/phenix+spin2/darshana/Libraries/event_gen/src/PHPythia/install/lib/libPHPythia.so");
  gSystem->Load("/direct/phenix+spin2/darshana/Libraries/event_gen/src/PHPythia/install/lib/libPHPyTrigger.so");		// 
  gSystem->Load("/direct/phenix+spin2/darshana/Libraries/event_gen/src/PHPythia/install/lib/libPHPyParticleSelect.so");
  gSystem->Load("libsimreco.so");
  //gSystem->Load("libfun4allfuncs.so");	// framework only
 /* gSystem->Load("libPHPythiaEventGen.so");
  gSystem->Load("libPHPythia.so");
  gSystem->Load("libPHPyTrigger.so");		// For PHPyTrigger derived classes
  gSystem->Load("libPHPyParticleSelect.so");	// For PHPyParticleSelect derived classes
  gSystem->Load("libsimreco.so");	*/// framework + reco modules

  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);


  // Counter
  MuonCounter * muc = new MuonCounter();
  muc->set_log_scale_dump();
  se->registerSubsystem(muc);

  //se->registerSubsystem(new MuonReadbackDST());
 // se->registerSubsystem(new FvtxReadbackDST());

  gSystem->ListLibraries();

// input manager

  Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin");
  in->fileopen(input_file);
  se->registerInputManager(in);


  PHPyDimuonTrigger *dimuonTrig = new PHPyDimuonTrigger();
  dimuonTrig->set_arm_selection(PHPyDimuonTrigger::MUON_ARM);
  dimuonTrig->set_minimum_pz(2);
  se->registerSubsystem(dimuonTrig);

// output manager
  if (dst_file)
    {
      std::cout << "registering Fun4AllDstOutputManager" << std::endl;
      Fun4AllDstOutputManager *dstManager = new Fun4AllDstOutputManager("DSTOUT", dst_file);

      se->registerOutputManager(dstManager);
    }

  ///////////////////////////////////////////
  // Analyze the Data.
  //////////////////////////////////////////


  se->skip(nSkip);
  se->run(nEvents);
  se->End();

  delete se;

  cout << "Completed reconstruction." << endl;
}
