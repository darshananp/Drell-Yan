// $Id: Fun4FVTX_Pisa.C,v 1.1 2014/07/10 17:52:38 darshana Exp $
/*!
  creates a unreconstructed DST
*/

void Fun4FVTX_Pisa(
  int nEvents = 0,
  const char *inputfile = "PISAEvent.root",
  const char *dstfile   = "dst_pisa.root",
  const char *dstreaderfile = "dst_reader.root",
  const char *ana_file = "muon_ana_ntuples.root",
  const char *singlepdstout = "singlemuon_pdst.root",
  const char *dimuonpdstout = "dimuon_pdst.root",
  const char *qafile = "qaRoot.root",
  int run_number =  392296 // 392296 Run13 //375910 - for CuAu// 360502 - this run number will give perfect dead map
)

{
  // flags
  bool use_muon_arm = true;
  bool use_rpc = true;
  bool use_fvtx = true;
  bool use_svx = false;

  bool write_qa = false;
  bool write_ndst = false;
  bool write_pdst = true;
  bool use_dstreader = false;

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
  TMutDatabaseCntrl::set_filename("use_local_landau_parameters_file",
                                  "landau_parameters.txt" );

  bool use_local_hv_files = false;
  if( use_local_hv_files )
    {
      // generate hv file from run number
      char* hv_file = "mut.disabledAnodes.dat_run307084";
      TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",true);
      TMutDatabaseCntrl::set_filename("use_local_dead_HV_file", hv_file );      
    } 
  else TMutDatabaseCntrl::set_database_access("use_local_dead_HV_file",false);

  TMutDatabaseCntrl::set_database_access("attenuated_channels",true);
  TMutDatabaseCntrl::set_database_access("dead_channels",true);
  //TMutDatabaseCntrl::set_database_access("dead_FEMs",true);

  /*
    disable check of mapFileScale since its decided by pisa,
    not by the run number.
  */
  mMfmMT::setMapFileFlag( mMfmMT::MAP_3D_PLUS_PLUS );
  mMfmMT::setMapFileScale( 1.0 );
  MuonUtil::set_check_mapfile_scale( false );
  
  ///////////////////////////////////////////
  // recoConsts setup
  //////////////////////////////////////////
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

  // Read the FVTX dead channel map and geometry from the database:
  // https://www.phenix.bnl.gov/WWW/offline/wikioff/index.php/FVTX/Global_Parameters
  TFvtxGlobalParCntrl::set_bool_par("is_pp",true);
  TFvtxGlobalParCntrl::set_bool_par("is_sim",true);
  TFvtxGlobalParCntrl::set_bool_par("use_svx",use_svx);
  // TFvtxDatabaseCntrl::set_flag("deadmap_use_calibration_database", true); it is default

  // temporarily setup load default PISA geometry for simulation.
  // the default one in database require PISA to run with loading read misalignment
  TFvtxGlobalParCntrl::set_bool_par("geom_use_calibration_database", false);
  TFvtxGlobalParCntrl::set_string_par("geom_root_file_path",
      "/direct/phenix+subsys+fvtx/jinhuang/geometry/simulation_base/");

  //See what happens if we change the b-field scale:
  //mMfmMT::setMapFileScale( 1.00 );
  
  // mutoo vertex source configuration
  // this allows to print which vertex is used and its value
  //TMutExtVtx::get().set_vtx_source( TMutExtVtx::MC );
  //TMutExtVtx::get().set_smear_z( false );
  //TMutExtVtx::get().set_vertex_name("SVX_PRECISE"); //if using a detector to reconstruct vertex, chose which detector
  
  //TMutExtVtx::get().set_verbosity( MUTOO::SOME );
  
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
  //    se->registerSubsystem( new TrigSimreco() );

  // vertex simulation
  // puts the vertex from the pisa header node into vtxOut object
  gSystem->Load( "libsimreco_base.so" );
  VtxSimreco* vtx_reco = new VtxSimreco();
  vtx_reco->SmearZ( false );
  vtx_reco->UseXY( false );
  vtx_reco->OverwriteBBC( true );
  vtx_reco->XVertexSigma( 0.01 ); // FVTX resolution of 100 microns
  vtx_reco->YVertexSigma( 0.01 );
  vtx_reco->ZVertexSigma( 0.01 ); // typically 2(0.5)cm for pp, 0.5,0.002,0.01 for HI when using BBC,VTX,FVTX
  se->registerSubsystem( vtx_reco );
  
  // global detectors subsystem
  //se->registerSubsystem( new HeadReco() );
  //se->registerSubsystem( new TrigReco( ));
  //se->registerSubsystem( new BbcReco() );
  //se->registerSubsystem( new ZdcReco() );

  // SVX reconstruction"
  if ( use_svx ) 
    {
      SvxParManager *svxpar = new SvxParManager();
      svxpar->Verbosity(0);
      svxpar->set_BeamCenter(0.,0.);
      svxpar->set_OffsetVtxToCnt(0.0, 0.0, 0.0);
      svxpar->set_OffsetEastToWest(0.0, 0.0, 0.0);
      svxpar->set_ReadGeoParFromFile(false);
      svxpar->set_GeometryFileName("svxPISA.par.ideal");
      se->registerSubsystem(svxpar);

      SvxSimulator *svxsim = new SvxSimulator();
      svxsim->set_StripixelNoise(0.0);
      svxsim->Verbosity(0);
      se->registerSubsystem(svxsim);
      
      SvxApplyHotDead *svxhotdead = new SvxApplyHotDead();
      svxhotdead->Verbosity(0);
      se->registerSubsystem(svxhotdead);
      
      ////////////// this step is not necessary if you want to embed  //////////
      SvxReco *svxrec = new SvxReco();
      svxrec->Verbosity(0);
      // svxrec->Load_ThresholdFile("threshold.h");
      svxrec->set_ThisIsSimulation();
      svxrec->set_UseStripThresholdDatbase(true);
      //svxrec->Load_ThresholdFile("threshold_ideal.h");
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

      // muon prdf unpacker
      MuonUnpackPisa* muon_unpack_pisa( new MuonUnpackPisa() );
      muon_unpack_pisa->Verbosity( 1 );
      muon_unpack_pisa->set_flag(MuonUnpackPRDF::SKIP_ZERO_SUPPRESSION,1);
      se->registerSubsystem( muon_unpack_pisa );      

      // mutoo reconstruction moved to after FVTX reco
//      se->registerSubsystem( new MuiooReco() );
//      se->registerSubsystem( new MuonDev() );
    }

  // fvtx prdf unpacker
  if (use_fvtx)
    {
      FvtxUnpackPisa *fvtx_unpack = new FvtxUnpackPisa();
      fvtx_unpack->set_do_response( true );
      fvtx_unpack->Verbosity(0);
      se->registerSubsystem( fvtx_unpack );
      
      ///////////////  this step is not necessary if you want to embed  //////////
      FvtxReco* fvtxreco = new FvtxReco();
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
      // Perform FVTX-Mutr track matching and refit track in case fvtxreco->set_do_mutr_matching(false)
      ///////////////  this step is not necessary if you want to embed  //////////
      se->registerSubsystem( new FvtxRecoWithMut() );
      /////////////////////////

    }
  
  if (write_qa)
    {
      // QA histograms
      gSystem->Load("libdstqa_muons.so");
      //  se->registerSubsystem( new QAMui() );  // unstable
      se->registerSubsystem( new QAMut() );
      se->registerSubsystem( new QAFvtx() );
    }

    // picoDST
  if( singlepdstout && write_pdst )
    {
        //      se->registerSubsystem( new MpcReco() );
        //      gSystem->Load("librxnp_subsysreco.so");
        //      se->registerSubsystem( new RxnpReco() );
        //      se->registerSubsystem( new RpSumXYReco() ); // recalibrator and rp doesn't work together!
      
        // global Reco
      se->registerSubsystem( new GlobalReco() );
      
      // MWG
      gSystem->Load("libMWGOO.so");
      PHInclusiveNanoCuts *MWGcuts = new MWGInclusiveNanoCutsv2();
      se->registerSubsystem(new MWGFvtxReco(MWGcuts));
      
      // module which counts tracklets and clusters withing 8 cone ranges
       FvtxConeTracklets* fvtxcone = new FvtxConeTracklets();
       se->registerSubsystem(fvtxcone);     
      
      mFillMCSingleMuonContainer* mcsngl = new mFillMCSingleMuonContainer();
      se->registerSubsystem(mcsngl);	
      
      mFillMCDiMuonContainer* mcdi = new mFillMCDiMuonContainer();
      se->registerSubsystem(mcdi);	
      
      mFillSingleMuonContainer* msngl = new mFillSingleMuonContainer();
      se->registerSubsystem(msngl);
      msngl->set_bbcz_cut(100);
      
      mFillDiMuonContainer* mdi = new mFillDiMuonContainer(false); // do not make mixed events
      mdi->set_is_sim(true);
      mdi->set_is_pp(true);
      se->registerSubsystem(mdi);
      mdi->set_mass_cut(0.5);
      
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
      
      Fun4AllOutputManager *outsmu = new Fun4AllDstOutputManager("Outsmu",singlepdstout);
      outsmu->AddNode("Sync");
      outsmu->AddNode("SingleMuonContainer");
      outsmu->AddNode("PHGlobal");
      outsmu->AddNode("MCSingleMuonContainer");
      outsmu->AddEventSelector("mFillSingleMuonContainer");//I have to comment out this when running embedding with single muon generator to have two files with same event entries one with and other without embedding.
      se->registerOutputManager(outsmu);
      
      Fun4AllOutputManager *outdimu = new Fun4AllDstOutputManager("Outdimu",dimuonpdstout);
      outdimu->AddNode("Sync");
      outdimu->AddNode("DiMuonContainer");
      outdimu->AddNode("MCDiMuonContainer");
      outdimu->AddNode("SingleMuonContainer");
      outdimu->AddNode("PHGlobal");
      outdimu->AddNode("MCSingleMuonContainer");
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

      //    dstManager->AddNode("TrigLvl1");
      //    dstManager->AddNode("L2Decision");
      //    dstManager->AddNode("Lvl2OutArray");
      
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

      // PHGlobal
      //    dstManager->AddNode("PHGlobal");
      //    dstManager->AddNode("PHGlobal_MUON");
      
      se->registerOutputManager(dstManager);
    }

  if ( use_dstreader )
    {
      mDSTReader * dr = new mDSTReader();
      dr -> add_node(mDSTReader::FvtxHit);
      // you can add more nodes, available nodes here can be found at mDSTReader::node_id
      // dr -> add_node( ... );
      se->registerSubsystem(dr);
      Fun4AllOutputManager *dr_out = new Fun4AllDstOutputManager("OutDST",dstreaderfile);
      dr_out->AddNode("Sync");
      dr_out->AddNode("DSTReader");
      // add more nodes with out->AddNode();
      se->registerOutputManager(dr_out);
    }


  ///////////////////////////////////////////
  // Analyze the Data.
  //////////////////////////////////////////
  
  //pfileopen(inputfile);
  //prun(nEvents);
  Fun4AllPisaInputManager *inMan = new Fun4AllPisaInputManager("PisaIn");
  se->registerInputManager(inMan);
  se->fileopen(inMan->Name(),inputfile);
  se->run(nEvents);
  se->End();

  Fun4AllHistoManager *hm = se->getHistoManager("QA");
  if( hm ) hm->setOutfileName( qafile );
  else cout << "Fun4Muons_Pisa - unable to find QA histograms" << endl;
  
  //se->dumpHistos();

  cout << "Completed reconstruction." << endl;
}
