#include <cassert>

using namespace std;

void
Test_SpinAna_mass()
{
  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  TVirtualFitter::SetDefaultFitter("Minuit2");
  gSystem->Load("/gpfs/mnt/gpfs02/phenix/spin3/spin2/darshana/Analysis/Run13/New_Lib/install/lib/libspin_analyzer.so");
  //gSystem->Load("/direct/phenix+spin2/darshana/Analysis/Run13/Drell_Yan_ALL/install/lib/libspin_analyzer.so");

  ReadTest();
}

void
ReadTest()
{

  TFile *_file0 = TFile::Open("afterqa.list_mass_hist.root");
  //TFile *_file0 = TFile::Open("7to14.root");
  //TFile *_file0 = TFile::Open("14.root");

  saHist * Ms_os_sig_n_saHist = (saHist *) _file0->GetObjectChecked("Ms_os_sig_n_saHist", "saHist");
  assert(Ms_os_sig_n_saHist);

  saHist * Ms_os_sig_s_saHist = (saHist *) _file0->GetObjectChecked("Ms_os_sig_s_saHist", "saHist");
  assert(Ms_os_sig_s_saHist);

  saHist * Ms_ss_n_saHist = (saHist *) _file0->GetObjectChecked("Ms_ss_n_saHist", "saHist");
  assert(Ms_ss_n_saHist);

  saHist * Ms_ss_s_saHist = (saHist *) _file0->GetObjectChecked("Ms_ss_s_saHist", "saHist");
  assert(Ms_ss_s_saHist);


  saHist * RelLumi_saHist = (saHist *) _file0->GetObjectChecked("RelLumi_saHist", "saHist");
  assert(RelLumi_saHist);

  Ms_os_sig_n_saHist->AutoLoad();
  Ms_os_sig_s_saHist->AutoLoad();
  Ms_ss_n_saHist->AutoLoad();
  Ms_ss_s_saHist->AutoLoad();

  RelLumi_saHist->AutoLoad();

  ///////////////////////

  TCanvas *c1 = new TCanvas("Test_SpinAna_ReadTest", "Test_SpinAna_ReadTest",1000, 900);
  c1->Divide(2, 2);
  c1->cd(1)->SetGrid();
  c1->Update();
  Ms_os_sig_s_saHist->getHisto("A_LL")->Draw("E1X0");
  Ms_os_sig_s_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");
  Ms_os_sig_s_saHist->getHisto("A_LL")->SetMarkerStyle(21);
  Ms_os_sig_s_saHist->getHisto("A_LL")->SetMarkerColor( kRed );
  Ms_os_sig_s_saHist->getHisto("A_LL")->SetLineWidth(3);
  Ms_os_sig_s_saHist->getHisto("A_LL")->SetMarkerSize(1.2);

double ALL_os_bin1_s= (Ms_os_sig_s_saHist->getHisto("A_LL"))->GetBinContent(1);
double ERR_os_bin1_s= (Ms_os_sig_s_saHist->getHisto("A_LL"))->GetBinError(1);
double ALL_os_bin2_s= (Ms_os_sig_s_saHist->getHisto("A_LL"))->GetBinContent(2);
double ERR_os_bin2_s= (Ms_os_sig_s_saHist->getHisto("A_LL"))->GetBinError(2);


  c1->cd(2)->SetGrid();
  c1->Update();
  Ms_os_sig_n_saHist->getHisto("A_LL")->Draw("E1X0");
  Ms_os_sig_n_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");
  Ms_os_sig_n_saHist->getHisto("A_LL")->SetMarkerStyle(21);
  Ms_os_sig_n_saHist->getHisto("A_LL")->SetMarkerColor( kRed );
  Ms_os_sig_n_saHist->getHisto("A_LL")->SetLineWidth(3);
  Ms_os_sig_n_saHist->getHisto("A_LL")->SetMarkerSize(1.2);

double ALL_os_bin1_n= (Ms_os_sig_n_saHist->getHisto("A_LL"))->GetBinContent(1);
double ERR_os_bin1_n= (Ms_os_sig_n_saHist->getHisto("A_LL"))->GetBinError(1);
double ALL_os_bin2_n= (Ms_os_sig_n_saHist->getHisto("A_LL"))->GetBinContent(2);
double ERR_os_bin2_n= (Ms_os_sig_n_saHist->getHisto("A_LL"))->GetBinError(2);


  c1->cd(3)->SetGrid();
  c1->Update();
  Ms_ss_s_saHist->getHisto("A_LL")->Draw("E1X0");
  Ms_ss_s_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");
  Ms_ss_s_saHist->getHisto("A_LL")->SetMarkerStyle(21);
  Ms_ss_s_saHist->getHisto("A_LL")->SetMarkerColor( kRed );
  Ms_ss_s_saHist->getHisto("A_LL")->SetLineWidth(3);
  Ms_ss_s_saHist->getHisto("A_LL")->SetMarkerSize(1.2);

double ALL_ss_bin1_s= (Ms_ss_s_saHist->getHisto("A_LL"))->GetBinContent(1);
double ERR_ss_bin1_s= (Ms_ss_s_saHist->getHisto("A_LL"))->GetBinError(1);
double ALL_ss_bin2_s= (Ms_ss_s_saHist->getHisto("A_LL"))->GetBinContent(2);
double ERR_ss_bin2_s= (Ms_ss_s_saHist->getHisto("A_LL"))->GetBinError(2);

  c1->cd(4)->SetGrid();
  c1->Update();
  Ms_ss_n_saHist->getHisto("A_LL")->Draw("E1X0");
  Ms_ss_n_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");
  Ms_ss_n_saHist->getHisto("A_LL")->SetMarkerStyle(21);
  Ms_ss_n_saHist->getHisto("A_LL")->SetMarkerColor( kRed );
  Ms_ss_n_saHist->getHisto("A_LL")->SetLineWidth(3);
  Ms_ss_n_saHist->getHisto("A_LL")->SetMarkerSize(1.2);

double ALL_ss_bin1_n= (Ms_ss_n_saHist->getHisto("A_LL"))->GetBinContent(1);
double ERR_ss_bin1_n= (Ms_ss_n_saHist->getHisto("A_LL"))->GetBinError(1);
double ALL_ss_bin2_n= (Ms_ss_n_saHist->getHisto("A_LL"))->GetBinContent(2);
double ERR_ss_bin2_n= (Ms_ss_n_saHist->getHisto("A_LL"))->GetBinError(2);



  cout<<"ALL OS bin1_S = "<< ALL_os_bin1_s << "+/-" <<ERR_os_bin1_s<<endl;
  cout<<"ALL OS bin2_S = "<< ALL_os_bin2_s << "+/-" <<ERR_os_bin2_s<<endl;
  cout<<"ALL OS bin1_N = "<< ALL_os_bin1_n << "+/-" <<ERR_os_bin1_n<<endl;
  cout<<"ALL OS bin2_N = "<< ALL_os_bin2_n << "+/-" <<ERR_os_bin2_n<<endl;
  cout<<"ALL SS bin1_S = "<< ALL_ss_bin1_s << "+/-" <<ERR_ss_bin1_s<<endl;
  cout<<"ALL SS bin2_S = "<< ALL_ss_bin2_s << "+/-" <<ERR_ss_bin2_s<<endl;
  cout<<"ALL SS bin1_N = "<< ALL_ss_bin1_n << "+/-" <<ERR_ss_bin1_n<<endl;
  cout<<"ALL SS bin2_N = "<< ALL_ss_bin2_n << "+/-" <<ERR_ss_bin2_n<<endl;
}
