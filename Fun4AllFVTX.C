/*!
  creates Real Case pdsts for run 16 Jin, Darshana, Xiaorong, Haiwang and Kathy
*/

void Fun4AllFVTX(
  int nEvents = 0,
  bool use_svx_cluster = false,
  const char *inputfile = "PISAEvent.root",
  const char *evtgen = "phpythia.root",
  const char *dstfile   = "dst_pisa.root",
  const char *dstreaderfile = "dst_reader.root",
  const char *inputfile2 = "/phenix/u/jinhuang/links/fvtx_tmp/miliped_work/Filter_Clock/DST.list", 
  const char *ana_file = "muon_ana_ntuples.root",
  const char *eval_file = "mutoo_eval_ntuple.root",
  const char *singlepdstout = "singlemuon_pdst.root",
  const char *dimuonpdstout = "dimuon_pdst.root",
  const char *hit_dstout = NULL, 
  const char *qafile = "qaRoot.root",
  int run_number =  367593,
)

{
  bool use_lvl2 = false;
  bool use_fvtx = true;
  bool write_ndst = false;
  bool write_pdst = true;
  bool write_qa = false;
  bool read_pythia = true;
  bool use_svx = true; 
  bool use_rpc = true;
  // bool do_eval = false;
  //bool is_pp = true;
  //bool is_sim = true;
  // flags
  bool use_muon_arm = true;
  bool use_dstreader = false;
  bool do_embedding = false; 
  bool write_align_dst = false;
  bool fill_mutr_mc_info = true;
  bool do_multiplicity = false;
  bool use_perfect_detector = false;

  const char *embed_src_topnode = "IMBED_SRC_TOP";

  // load libraries
  gSystem->Load("libPHGeant" );
  gSystem->Load("libfun4all");
  gSystem->Load("libsimreco_base");
  gSystem->Load("libmutoo_subsysreco");
  gSystem->Load("libmuon_subsysreco");
  gSystem->Load("libfun4allfuncs_muons");
  gSystem->Load("libfvtx_subsysreco");
  gSystem->Load("liblvl2");
  gSystem->Load("libPythia6.so");
  gSystem->Load("libPHPythia.so");
  gSystem->Load("libPHPythiaEventGen.so");
  gSystem->Load("libpicodst_object.so");

  cerr << "libraries loaded.\n";

//For Muon Tracker Dead Maps 
 if( use_perfect_detector )
    {

      // tell muid tube efficiency must be set to 1
      // rc->set_DoubleFlag("MUIOO_TUBE_EFF",1.0);
      // make perfect mutr detector
      TMutDatabaseCntrl::set_database_access("disable_HV",false);
      TMutDatabaseCntrl::set_database_access("attenuated_channels",false);
      TMutDatabaseCntrl::set_database_access("dead_channels",false);
       
    } else 
{

  TMutDatabaseCntrl::set_database_access("disable_HV",true);
  TMutDatabaseCntrl::set_database_access( "use_local_landau_parameters_file", true );
  TMutDatabaseCntrl::set_filename("use_local_landau_parameters_file","landau_parameters.txt" );
  TMutDatabaseCntrl::set_database_access("attenuated_channels",true);
  TMutDatabaseCntrl::set_database_access("dead_channels",true);
}

// This includes the MUID tube efficiencies
  bool use_local_muid_tube_eff = true;
  if {use_local_muid_tub_eff} {
    
    std::string muid_eff_south("twopack_eff_south.dat");
    std::string muid_eff_north("twopack_eff_north.dat");

    std::cout << "Fun4Muons_Pisa - using local two pack efficiciency files:" << muid_eff_south << ", "  << muid_eff_north <<std::endl;

    TMuiHVMask::set_mode(TMuiHVMask::FROM_FILE);
    TMuiHVMask::set_filename_south(muid_eff_south.c_str());
    TMuiHVMask::set_filename_north(muid_eff_north.c_str());
  }

  bool use_local_hv_files = false;
	if( use_local_hv_files )
	{
		std::stringstream hv_file;
		hv_file << "/phenix/spin/phnxsp01/keyaaron/Run13MuTrHV/mut.disabledAnodes.dat_run" << run_number;
		TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",true);
		TMutDatabaseCntrl::set_filename("use_local_dead_HV_file", hv_file.str().c_str());      
	} 
       else TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",false);

   mMfmMT::setMapFileFlag( mMfmMT::MAP_3D_PLUS_PLUS );
   mMfmMT::setMapFileScale( 1.0 );
   MuonUtil::set_check_mapfile_scale( false );


  // recoConsts setup

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("SVXACTIVE", use_svx);
  rc->set_IntFlag("SIMVERTEXFLAG",2);
  rc->set_FloatFlag("SIMZ0VERTEX",0);
  rc->set_FloatFlag("SIMT0VERTEX",0);
  rc->set_FloatFlag("SIMZ0VERTEXWIDTH",0.00);
  rc->set_FloatFlag("SIMT0VERTEXWIDTH",0);
  rc->set_IntFlag("SIMULATIONFLAG",2);
  rc->set_IntFlag("EMBEDFLAG",0);
  rc->set_IntFlag("PRINT_MUTOO_PARAMETERS",1);
  rc->set_IntFlag("RUNNUMBER", run_number);
  rc->set_CharFlag("FVTX_EMBED_SRC_TOPNODE",embed_src_topnode);//"TOP"

  if (use_lvl2) 
 {
  // Set Lvl2 flags
  rc->set_IntFlag("LVL2_REAL_DATA",1);
  rc->set_IntFlag("LVL2_YEAR",4);
  rc->set_IntFlag("FORCE_LVL2",1);
  rc->set_IntFlag("LVL2_USE_ASCII_DB",1);
  rc->set_CharFlag("LVL2_DB_DIR","/afs/rhic.bnl.gov/phenix/users/frawley/lvl2_db/RUN4_REAL");
  rc->set_CharFlag("Run2Lvl2AlgoName", "");
 }

  if (use_fvtx)
    {
  TFvtxGlobalParCntrl::set_bool_par("is_pp",true);
  TFvtxGlobalParCntrl::set_bool_par("is_sim",false);
  TFvtxGlobalParCntrl::set_bool_par("use_svx",use_svx);
  TFvtxGlobalParCntrl::set_bool_par("deadmap_auto_load",false);//make this true or remove it if you need fvtx deadmaps
  TFvtxGlobalParCntrl::set_bool_par("deadmap_use_calibration_database", true); 
  TFvtxGlobalParCntrl::set_bool_par("geom_use_calibration_database", false);
  TFvtxGlobalParCntrl::set_string_par("geom_root_file_path","/phenix/u/jinhuang/links/fvtx_data/geometry/simulation_base/");
  TFvtxGlobalParCntrl::set_string_par("geom_root_file_name", "fvtxgeom.root");
}


  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  ///////////////////////////////////////////
  // Subsystems
  //////////////////////////////////////////

  // run header and trigger setting
  se->registerSubsystem( new HeadSimreco() );

  // vertex simulation
  // puts the vertex from the pisa header node into vtxOut object
  gSystem->Load( "libsimreco_base.so" );
  VtxSimreco* vtx_reco = new VtxSimreco();
  vtx_reco->SmearZ( true );//false
  vtx_reco->UseXY( false );
  vtx_reco->OverwriteBBC( true );
  vtx_reco->XVertexSigma( 0.00 ); 
  vtx_reco->YVertexSigma( 0.00 );//0.01
  vtx_reco->ZVertexSigma( 2.00 ); 
  se->registerSubsystem( vtx_reco );
  
  // counter
  MuonCounter* counter = new MuonCounter();
  counter->set_event_dump(1000);
  se->registerSubsystem( counter );

if (use_lvl2) {
		gSystem->Load("liblvl2");
		SubsysReco *lvl2reco = new Lvl2Reco();
		lvl2reco->Verbosity(0);
		se->registerSubsystem(lvl2reco);

		SubsysReco *lvl2stats = new Lvl2StatsEval();
		lvl2stats->Verbosity(0);
		se->registerSubsystem(lvl2stats);

		Lvl2RunTrigSelect *lvl2runtrigselect = new Lvl2RunTrigSelect();
		lvl2runtrigselect->AddTrigger("L2MutrDimuonSouthTrigger");
		lvl2runtrigselect->AddTrigger("L2MutrDimuonNorthTrigger");
		lvl2runtrigselect->AddTrigger("L2MuidDimuonSouthTrigger");
		lvl2runtrigselect->AddTrigger("L2MuidDimuonNorthTrigger");

		// tell level2 to reject events which do not trigger on any of the previous
		lvl2runtrigselect->SetReturnCode("ABORT");
		lvl2runtrigselect->Verbosity(0);
		se->registerSubsystem(lvl2runtrigselect);
	}

  if ( use_svx ) 
    {
      SvxParManager *svxpar = new SvxParManager();
      svxpar->Verbosity(0);
      svxpar->set_BeamCenter(0.,0.);
      svxpar->set_OffsetVtxToCnt(0.0, 0.0, 0.0);
      svxpar->set_OffsetEastToWest(0.0, 0.0, 0.0);
      svxpar->set_ReadGeoParFromFile(true);//1
      svxpar->set_GeometryFileName("/direct/phenix+subsys+fvtx/jinhuang/geometry/simulation_base/svxPISA.par.ideal"); 

      bool use_local_vtx_eff = true;//ideal false
  
  if( use_local_vtx_eff ){
      const int vtx_runnumber = 421822;
      const char *pixel_refmap = "/direct/phenix+hhj2/dcm07e/run14MiniProd/fieldon/blank_pixel_refmap.txt";

      string pixel_diffmap =Form("/direct/phenix+prod01/phnxreco/millepede/deadmaps_run15pp200/zerofield/woraw/pixels/pixel_deadmap_run15pp200woraw_run%d.dat",vtx_runnumber);

      string pixel_chipmap =Form("/direct/phenix+prod01/phnxreco/millepede/deadmaps_run15pp200/zerofield/woraw/pixels/chip_deadmap_run15pp200woraw_run%d.dat",vtx_runnumber);

      string strip_deadchannel = Form("/direct/phenix+prod01/phnxreco/millepede/deadmaps_run15pp200/zerofield/woraw/strips/run%d_strip_hotdeadChannels.txt",vtx_runnumber);

      string strip_deadRCC =Form("/direct/phenix+prod01/phnxreco/millepede/deadmaps_run15pp200/zerofield/woraw/strips/run%d_strip_hotdeadReadouts.txt",vtx_runnumber);

      svxpar->set_UseRefDiffPixelMap(1);
      svxpar->set_ReadRefDiffPixelMapFromFile(1);
      svxpar->set_RefDiffPixelMapFiles(pixel_refmap, pixel_diffmap,pixel_chipmap);


      svxpar->set_ReadStripHotDeadFromFile(1);
      svxpar->set_StripHotDeadFileName(strip_deadchannel);
      svxpar->set_StripHotDeadReadoutsFileName(strip_deadRCC);
}
      se->registerSubsystem(svxpar);

      SvxSimulator *svxsim = new SvxSimulator();
      svxsim->set_StripixelNoise(0.0);
      svxsim->Verbosity(0);
      se->registerSubsystem(svxsim);
      
      SvxApplyHotDead *svxhotdead = new SvxApplyHotDead();
      svxhotdead->Verbosity(0);
      se->registerSubsystem(svxhotdead);
      
      SvxReco *svxrec = new SvxReco();
      svxrec->Verbosity(0);
      svxrec->set_ThisIsSimulation();
      svxrec->set_UseStripThresholdDatbase(true);
      svxrec->set_StripixelAdcSumThreshold(0);
      se->registerSubsystem(svxrec);

      SvxPriVertexSeedFinder *svxvtxseedfinder = new SvxPriVertexSeedFinder();
      svxvtxseedfinder->Verbosity(0);
      se->registerSubsystem(svxvtxseedfinder);
      
      SvxStandAloneReco *svxstandalone = new SvxStandAloneReco();
      svxstandalone->Verbosity(0);
      svxstandalone->setVertexRecoFlag(2);
      se->registerSubsystem(svxstandalone);

      SvxPrimVertexFinder *svxprimvtxfinder = new SvxPrimVertexFinder();
      svxprimvtxfinder->Verbosity(0);
      se->registerSubsystem(svxprimvtxfinder);
      ////////////////////////////////////////////////
    }


      MuonUnpackPisa* muon_unpack_pisa( new MuonUnpackPisa() );
      muon_unpack_pisa->Verbosity( 1 );
      muon_unpack_pisa->set_flag(MuonUnpackPRDF::SKIP_ZERO_SUPPRESSION,1);
      muon_unpack_pisa->set_flag(MuonUnpackPisa::DO_RESPONSE,1);
      se->registerSubsystem( muon_unpack_pisa );
     // se->registerSubsystem( new MuonAnaTuples("MUONANATUPLES", ana_file ));
     // se->registerSubsystem( new MuonEval("MuonEval", eval_file ));
      // mutoo reconstruction

//major work PISA need to  be unpacked

  if (use_rpc)
    {
                //RPC
      gSystem->Load( "librpc_subsysreco" );
      se->registerSubsystem( new RpcUnpackPisa());
      se->registerSubsystem( new RpcReco());
    }

  // fvtx prdf unpacker
  if (use_fvtx)
    {
      FvtxUnpackPisa *fvtx_unpack = new FvtxUnpackPisa();
      fvtx_unpack->set_do_response( true );
      fvtx_unpack->Verbosity(0);
      se->registerSubsystem( fvtx_unpack );
      

      FvtxReco* fvtxreco = new FvtxReco();
      fvtxreco->set_do_embedding(do_embedding);
      se->registerSubsystem(fvtxreco);
      
      // primary vertex from FVTX
      FvtxPrimVertex* fvtxprimvtx = new FvtxPrimVertex();
      se->registerSubsystem(fvtxprimvtx);                         // (**) in development

//to get the primary vertex from the VTX (saved in Evt_vtx*)not all events have vtx vertex. only few events has it
  FvtxPrimVertex* fvtxprimvtx1 = new FvtxPrimVertex();
  fvtxprimvtx1->set_source( FvtxPrimVertex::Segments );
  fvtxprimvtx1->set_clustering( FvtxPrimVertex::AllInOne );
  fvtxprimvtx1->set_max_vertexes(1);
  fvtxprimvtx1->set_vertex_output(0,"VTX_ONLY_FvtxPrimVertex",20);
  fvtxprimvtx1->set_arms_in_use(false,false,true,true); // North, South, East, West
  se->registerSubsystem(fvtxprimvtx1);



 /* if ( is_pp )
      {
	TMutNode<mFvtxFindTrackPar>::find_node(se->topNode(),"mFvtxFindTrackPar")->set_allowTwoHitTracks(true);
	TMutNode<mFvtxFindTrackPar>::find_node(se->topNode(),"mFvtxFindTrackPar")->set_filterTwoHitTracks(true);      
      }
*/
      se->registerSubsystem( new MuiooReco() );
      se->registerSubsystem( new MuonDev() );

      se->registerSubsystem( new FvtxRecoWithMut() );
	
/*	// primary vertex from FVTX
      FvtxPrimVertex* fvtxprimvtx = new FvtxPrimVertex();
      fvtxprimvtx->set_clustering(FvtxPrimVertex::Graph);//AllInOne, Cesars, Graph, AntiKt
      se->registerSubsystem(fvtxprimvtx);

     /*! module to fill the container for multiplicity study
*/
  if(do_multiplicity){
			mFillFvtxPrimVtxContainer* mfpvc = new mFillFvtxPrimVtxContainer(1);
			mfpvc->Verbosity(0);//0: minimum out put, 1 is normal, 2 is debug mode
			mfpvc->set_cut_tracklet_vtx_dca(999999.);//default is 1. cm
			mfpvc->set_use_tracklet_vtx_dca_cut(false); //default is true
			mfpvc->set_use_chi2_tracklet_cut(true);//default is true
			mfpvc->set_use_2_hit_tracklet(true); // default is true
			se->registerSubsystem(mfpvc);
		      }


    }  
  

   if(fill_mutr_mc_info)  
	      {
	        // muon evaluation module, used to get MC information in pico DST
	        MuonEval* mueval = new MuonEval();
		mueval->SetSignalNodeName("TOP");
	        mueval->set_flags(0); // no ntuple output needed
	        se->registerSubsystem (mueval);
	      }


  if (write_qa)
    {
      // QA histograms
      gSystem->Load("libdstqa_muons.so");
      se->registerSubsystem( new QAMut() );
      se->registerSubsystem( new QAFvtx() );
    }

    // picoDST
  if( singlepdstout && write_pdst )
    {
        // global Reco
      se->registerSubsystem( new GlobalReco() );
      
      // MWG
      gSystem->Load("libMWGOO.so");
      PHInclusiveNanoCuts *MWGcuts = new MWGInclusiveNanoCutsv2();
      se->registerSubsystem(new MWGFvtxReco(MWGcuts));
      
      mFillTriggerEmulatorContainer* triggerEmulator = new mFillTriggerEmulatorContainer( );
      se->registerSubsystem(triggerEmulator);

      // module which counts tracklets and clusters withing 8 cone ranges
      FvtxConeTracklets* fvtxcone = new FvtxConeTracklets();
      se->registerSubsystem(fvtxcone);     

      mFillSingleMuonContainer* msngl = new mFillSingleMuonContainer();
      msngl->set_vtx_vertex_name("VTX_ONLY_FvtxPrimVertex");
      msngl->set_bbcz_cut(100);
      se->registerSubsystem(msngl);

      mFillRpcSingleMuonContainer* msngl_rpc = new mFillRpcSingleMuonContainer();
      se->registerSubsystem(msngl_rpc);      
      
      mFillDiMuonContainer* mdi = new mFillDiMuonContainer(false); // do not make mixed events
      mdi->set_is_pp(true);
      mdi->set_mass_cut(0.5);
      se->registerSubsystem(mdi);

      mFillRpcDiMuonContainer* mdi_rpc = new mFillRpcDiMuonContainer();
      se->registerSubsystem(mdi_rpc);

 if(fill_mutr_mc_info)
    	{
      	  //********************* Mutr MC information Module ******************//
    	  //mFillMCSingleMuonContainer* msngl_mc = new mFillMCSingleMuonContainer();
    	 mFillMCSingleMuonFvtxContainer* msngl_mc = new mFillMCSingleMuonFvtxContainer();
    	 se->registerSubsystem(msngl_mc);
   
         mFillMCSingleMuonContainer* mcsngl = new mFillMCSingleMuonContainer();
         se->registerSubsystem(mcsngl);

         mFillMCDiMuonContainer* mcdi = new mFillMCDiMuonContainer();
         se->registerSubsystem(mcdi);

	}

      //se->registerOutputManager(outdimu);
   

    
      if ( write_ndst )
        {
          Fun4AllOutputManager *outndst = new Fun4AllDstOutputManager("Outndst", ndstfile);
          outndst->AddNode("Sync");
          outndst->AddNode("TrigLvl1");
          outndst->AddNode("VtxOut");
          outndst->AddNode("PHGlobal");
          outndst->AddNode("PHPythiaHeader");
          outndst->AddNode("PHPythia");
          outndst->AddNode("PHMuoTracksOO");
          outndst->AddNode("TFvtxCompactTrk");
          se->registerOutputManager(outndst);
        }

      
      mDSTReader * dr = new mDSTReader();
      dr -> add_node(mDSTReader::FvtxCompactTrk);
      dr -> add_node(mDSTReader::FvtxHit);
      se->registerSubsystem(dr);

      Fun4AllOutputManager *outdimu = new Fun4AllDstOutputManager("Outdimu",dimuonpdstout);
      outdimu->AddNode("Sync");
      outdimu->AddNode("PHGlobal");
      outdimu->AddNode("DiMuonContainer");
      outdimu->AddNode("MCDiMuonContainer");
      outdimu->AddNode("SingleMuonContainer");
      outdimu->AddNode("MCSingleMuonContainer");
      outdimu->AddNode("RpcSingleMuonContainer");
      outdimu->AddNode("RpcDiMuonContainer");
      outdimu->AddNode("DSTReader");
      outdimu->AddNode("PHPythiaHeader");
      outdimu->AddNode("PHPythia");
      outdimu->AddNode("VtxOut");
      outdimu->AddNode("TriggerEmulatorContainer");
      if(do_multiplicity)outdimu->AddNode("FvtxPrimVtxContainer");
      if (fill_mutr_mc_info) outdimu->AddNode("MCSingleMuonFvtxContainer"); // Mutr MC information container

      outdimu->AddEventSelector("mFillDiMuonContainer");

      se->registerOutputManager(outdimu);

}


  ///////////////////////////////////////////
  // IOManagers...
  ///////////////////////////////////////////

/*  if( dstfile ) 
    {
      Fun4AllDstOutputManager *dstManager  = new Fun4AllDstOutputManager("DSTOUT",  dstfile);
    
      dstManager->AddNode("RunHeader");
      dstManager->AddNode("EventHeader");
      dstManager->AddNode("VtxOut");
      dstManager->AddNode("BbcOut");
      dstManager->AddNode("BbcRaw");
      dstManager->AddNode("ZdcOut");
      dstManager->AddNode("ZdcRaw");

      dstManager->AddNode("TMCPrimary");
      dstManager->AddNode("PHPythiaHeader");
      dstManager->AddNode("PHPythia");

      
   if ( use_muon_arm )
     {
      // Muioo nodes
      dstManager->AddNode("TMuiHitO");
      dstManager->AddNode("TMuiMCHitO");
          
      // Mutoo nodes
      dstManager->AddNode("TMutHit");
      dstManager->AddNode("TMutMCHit");
      dstManager->AddNode("TMutMCTrk");

      //    dstManager->AddNode("PHMuoTracksOO");
     }
      
   if ( use_fvtx )
     {
       // FVTX nodes
       dstManager->AddNode("TFvtxHit");
       dstManager->AddNode("TFvtxMCHit");
       dstManager->AddNode("TFvtxPisaHit");
     }

      // From EVA node
      dstManager->AddNode("header");
      dstManager->AddNode("fkin");
      dstManager->AddNode("primary");
      dstManager->AddNode("pythia");
      
      se->registerOutputManager(dstManager);
    }
*/

  if ( use_dstreader )
    {
      mDSTReader * dr = new mDSTReader();
      dr -> add_node(mDSTReader::FvtxHit);
      se->registerSubsystem(dr);
      Fun4AllOutputManager *dr_out = new Fun4AllDstOutputManager("OutDST",dstreaderfile);
      dr_out->AddNode("Sync");
      dr_out->AddNode("DSTReader");
      se->registerOutputManager(dr_out);
    }


  ///////////////////////////////////////////
  // Analyze the Data.
  //////////////////////////////////////////


  Fun4AllPisaInputManager *inMan = new Fun4AllPisaInputManager("PisaIn");
  se->registerInputManager(inMan);
  se->fileopen(inMan->Name(),inputfile);


  Fun4AllDstInputManager *ipythia = new Fun4AllNoSyncDstInputManager("DSTin2","DST");
  se->registerInputManager(ipythia);
  se->fileopen(ipythia->Name(),evtgen);

 /* if (inputfile2)
    {

      Fun4AllNoSyncDstInputManager *indst = new Fun4AllNoSyncDstInputManager("DSTIn", "DST", embed_src_topnode);
      se->registerInputManager(indst);
      indst->AddListFile(inputfile2); // for one file

    }*/

  se->run(nEvents);
  se->End();

  if ( write_qa )
	{
	Fun4AllHistoManager *hm = se->getHistoManager("QA");
	if( hm ) hm->setOutfileName( qafile );
	else cout << "Fun4Muons_Pisa - unable to find QA histograms" << endl;
	se->dumpHistos();
	}
  cout << "Completed reconstruction." << endl;
}
