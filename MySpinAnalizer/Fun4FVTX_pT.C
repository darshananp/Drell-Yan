#include <exception>

using namespace std;

void
Fun4FVTX_pT(int nEvents = 100000,char *input_file = "afterqa.list", char *dst_file = "dst_out.root")
{

  // load libraries
  gSystem->Load("libmutoo_subsysreco.so");
  gSystem->Load("libfun4all.so");
  gSystem->Load("librecal.so");
  gSystem->Load("libfun4allfuncs.so");
  gSystem->Load("liblpc.so");
  gSystem->Load("libcompactCNT.so");
  gSystem->Load("libfun4allfuncs_muons");
  gSystem->Load("libMWGOO");
  gSystem->Load("libmutrg");
  gSystem->Load("librpc_subsysreco");
  gSystem->Load("libSvxDstQA.so");
  gSystem->Load("libpicodst_object");
 // gSystem->Load("libspin_analyzer");
  gSystem->Load("/gpfs/mnt/gpfs02/phenix/spin3/spin2/darshana/Analysis/Run13/New_Lib/install/lib/libspin_analyzer.so");
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //! load a QA text file marking good run and good corssings
  //saModLoadSpinInfo * mspin = new saModLoadSpinInfo("SpinInfoWithAdditionalQA");
  //mspin -> load_crossing_flag("spinQA.txt");

  saSpinAnalyzer * sa = new saSpinAnalyzer("SpinAnalyzer");

  sa->set_auto_save_hist(
      string(input_file) + string("_") + string("pt_hist.root"));
  sa->Verbosity(1);

  sa->RegisterModule(new saModuleDimuonDYDarshana());
  se->registerSubsystem(sa);

  ///////////////////////////////////////////
  // DST
  //////////////////////////////////////////
  if (dst_file)
    {
      std::cout << "registering Fun4AllDstOutputManager" << std::endl;

      Fun4AllDstOutputManager *dstManager = new Fun4AllDstOutputManager(
          "DSTOUT", string(input_file) + string("_") + string(dst_file));

      dstManager->AddNode("saEventProperty");
      dstManager->AddNode("DiMuonContainer");
      dstManager->AddNode("SimpleDimuonFlag");
      dstManager->AddNode("Sync");
      dstManager->AddNode("TrigLvl1");


      dstManager->AddEventSelector("SpinAnalyzer");
      se->registerOutputManager(dstManager);
    }
  ///////////////////////////////////////////
  // Analyze the Data.
  //////////////////////////////////////////

  gSystem->ListLibraries();

  Fun4AllDstInputManager *in = new Fun4AllDstInputManager("DSTin");
//  in->fileopen(input_file);
  in->AddListFile(input_file);
  se->registerInputManager(in);

  se->run(nEvents);
  se->End();

  delete se;

  cout << "Completed reconstruction." << endl;
}
