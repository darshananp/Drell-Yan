/*!
  creates Ideal Case pdsts for run 16 Jin and Darshana
*/

void Fun4FvtxIdeal(
  int nEvents = 100,
  bool use_svx_cluster = false,
  const char *inputfile = "PISAEvent.root",
  const char *dstfile   = "dst_pisa.root",
  const char *dstreaderfile = "dst_reader.root",
  const char *inputfile2 = "/phenix/u/jinhuang/links/fvtx_tmp/miliped_work/Filter_Clock/DST.list", //
  const char *ana_file = "muon_ana_ntuples.root",
  const char *singlepdstout = "singlemuon_pdst.root",
  const char *dimuonpdstout = "dimuon_pdst.root",
  const char *hit_dstout = NULL, 
  const char *qafile = "qaRoot.root",
  int run_number =  367593
)

{
  // flags
  bool use_muon_arm = true;
  bool use_rpc = true;
  bool use_fvtx = true;
  bool use_svx = true; 

  bool write_qa = false;
  bool write_ndst = false;
  bool write_pdst = true;
  bool use_dstreader = false;
  bool do_embedding = false; 

  bool use_lvl2 = false;
  bool do_eval = false;
  bool write_align_dst = false;
  bool fill_mutr_mc_info = true;

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

  TMutDatabaseCntrl::set_database_access("disable_HV",true);
  TMutDatabaseCntrl::set_database_access( "use_local_landau_parameters_file", true );
  TMutDatabaseCntrl::set_filename("use_local_landau_parameters_file","landau_parameters.txt" );

  bool use_local_muid_tube_eff = false;
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
      char* hv_file = "mut.disabledAnodes.dat_run307084";
      TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",true);
      TMutDatabaseCntrl::set_filename("use_local_dead_HV_file", hv_file );      
    } 
  else TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",false);

  TMutDatabaseCntrl::set_database_access("attenuated_channels",true);
  TMutDatabaseCntrl::set_database_access("dead_channels",true);
 
  ///////////////////////////////////////////
  // recoConsts setup
  ////////////////////////////////////////// Wgroups macro 
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("SVXACTIVE", 1);
  rc->set_IntFlag("SIMVERTEXFLAG",2);
  rc->set_FloatFlag("SIMZ0VERTEX",0);
  rc->set_FloatFlag("SIMT0VERTEX",0);
  rc->set_FloatFlag("SIMZ0VERTEXWIDTH",0.02);
  rc->set_FloatFlag("SIMT0VERTEXWIDTH",0);
  rc->set_IntFlag("SIMULATIONFLAG",2);
  rc->set_IntFlag("EMBEDFLAG",0);
  rc->set_IntFlag("PRINT_MUTOO_PARAMETERS",1);
  rc->set_IntFlag("RUNNUMBER", run_number);
  rc->set_CharFlag("FVTX_EMBED_SRC_TOPNODE",embed_src_topnode);

  TFvtxGlobalParCntrl::set_bool_par("is_pp",true);
  TFvtxGlobalParCntrl::set_bool_par("is_sim",true);
  TFvtxGlobalParCntrl::set_bool_par("use_svx",use_svx);

  TFvtxGlobalParCntrl::set_bool_par("geom_use_calibration_database", false);
  TFvtxGlobalParCntrl::set_string_par("geom_root_file_path","/direct/phenix+subsys+fvtx/jinhuang/geometry/simulation_base/");
  
  TMutExtVtx::get().set_verbosity( MUTOO::SOME );
  
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
  vtx_reco->SmearZ( true );
  vtx_reco->UseXY( false );
  vtx_reco->OverwriteBBC( true );
  vtx_reco->XVertexSigma( 0.00 ); 
  vtx_reco->YVertexSigma( 0.00 );
  vtx_reco->ZVertexSigma( 0.01 ); 
  se->registerSubsystem( vtx_reco );
  

  if ( use_svx ) 
    {
      SvxParManager *svxpar = new SvxParManager();
      svxpar->Verbosity(0);
      svxpar->set_BeamCenter(0.,0.);
      svxpar->set_OffsetVtxToCnt(0.0, 0.0, 0.0);
      svxpar->set_OffsetEastToWest(0.0, 0.0, 0.0);
      svxpar->set_ReadGeoParFromFile(false);
      svxpar->set_GeometryFileName("svxPISA.par.ideal");

/*  bool use_local_vtx_eff = false;//ideal false
  
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
}*/
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

  if ( use_muon_arm ) 
    {
      // counter
      MuonCounter* counter = new MuonCounter();
      counter->set_event_dump(1000);
      se->registerSubsystem( counter );

      MuonUnpackPisa* muon_unpack_pisa( new MuonUnpackPisa() );
      muon_unpack_pisa->Verbosity( 1 );
      muon_unpack_pisa->set_flag(MuonUnpackPRDF::SKIP_ZERO_SUPPRESSION,1);
      se->registerSubsystem( muon_unpack_pisa );      

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

    }  


  if ( use_muon_arm )
    {
      // mutoo reconstruction
      se->registerSubsystem( new MuiooReco() );
      se->registerSubsystem( new MuonDev() );
    }

  if (use_rpc)
    {
      //RPC
      gSystem->Load( "librpc_subsysreco" );
      se->registerSubsystem( new RpcUnpackPRDF());
      se->registerSubsystem( new RpcReco());
    }

  if (use_fvtx)
    {

      se->registerSubsystem( new FvtxRecoWithMut() );

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
      
      // module which counts tracklets and clusters withing 8 cone ranges
      FvtxConeTracklets* fvtxcone = new FvtxConeTracklets();
      se->registerSubsystem(fvtxcone);     
	
      
      mFillSingleMuonContainer* msngl = new mFillSingleMuonContainer();
      se->registerSubsystem(msngl);
      msngl->set_bbcz_cut(100);

      mFillRpcSingleMuonContainer* msngl_rpc = new mFillRpcSingleMuonContainer();
      se->registerSubsystem(msngl_rpc);      
      
      mFillMCSingleMuonContainer* mcsngl = new mFillMCSingleMuonContainer();
      se->registerSubsystem(mcsngl);

      mFillDiMuonContainer* mdi = new mFillDiMuonContainer(false); // do not make mixed events
      mdi->set_is_pp(true);
      se->registerSubsystem(mdi);
      mdi->set_mass_cut(0.5);

      mFillRpcDiMuonContainer* mdi_rpc = new mFillRpcDiMuonContainer();
      se->registerSubsystem(mdi_rpc);

      mFillMCDiMuonContainer* mcdi = new mFillMCDiMuonContainer();
      se->registerSubsystem(mcdi); 
    
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
      outdimu->AddEventSelector("mFillDiMuonContainer");
      se->registerOutputManager(outdimu);

}

  ///////////////////////////////////////////
  // IOManagers...
  ///////////////////////////////////////////

  // dst
  if( dstfile ) 
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

  if (inputfile2)
    {

      Fun4AllNoSyncDstInputManager *indst = new Fun4AllNoSyncDstInputManager("DSTIn", "DST", embed_src_topnode);
      se->registerInputManager(indst);
      indst->AddListFile(inputfile2); // for one file

    }

  Fun4AllPisaInputManager *inMan = new Fun4AllPisaInputManager("PisaIn");
  se->registerInputManager(inMan);
  se->fileopen(inMan->Name(),inputfile);
  se->run(nEvents);
  se->End();

  Fun4AllHistoManager *hm = se->getHistoManager("QA");
  if( hm ) hm->setOutfileName( qafile );
  else cout << "Fun4Muons_Pisa - unable to find QA histograms" << endl;

  cout << "Completed reconstruction." << endl;
}
