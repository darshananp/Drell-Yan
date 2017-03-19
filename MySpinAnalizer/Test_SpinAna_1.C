#include <cassert>

using namespace std;

void
Test_SpinAna_1()
{
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  TVirtualFitter::SetDefaultFitter("Minuit2");

  gSystem->Load("/direct/phenix+spin2/darshana/Analysis/Run13/Drell_Yan_ALL/install/lib/libspin_analyzer.so");

  ReadTest();
}

void
ReadTest()
{

  TFile *_file0 = TFile::Open("good_filtered.list_hist.root");

  saHist * InvMass_saHist = (saHist *) _file0->GetObjectChecked("InvMass_saHist", "saHist");
  assert(InvMass_saHist);
  saHist * RelLumi_saHist = (saHist *) _file0->GetObjectChecked("RelLumi_saHist", "saHist");
  assert(RelLumi_saHist);

  InvMass_saHist->AutoLoad();

  RelLumi_saHist->AutoLoad();

  ///////////////////////

  TCanvas *c1 = new TCanvas("Test_SpinAna_ReadTest", "Test_SpinAna_ReadTest",1000, 900);
  c1->Divide(2, 2);

  c1->cd(1)->SetGrid();
  c1->Update();

  InvMass_saHist->getHisto("YIELD")->Draw("E1X0");
  InvMass_saHist->getHisto("YIELD")->SetYTitle("Yield");
  InvMass_saHist->getHisto("YIELD")->SetMarkerStyle(21);
  InvMass_saHist->getHisto("YIELD")->SetMarkerColor( kRed );
  InvMass_saHist->getHisto("YIELD")->SetLineWidth(3);
  InvMass_saHist->getHisto("YIELD")->SetMarkerSize(1.2);
  c1->cd(2)->SetGrid();
  c1->Update();

  InvMass_saHist->getHisto("A_LL")->Draw("E1X0");
  InvMass_saHist->getHisto("A_LL")->SetYTitle("Double spin asymmetry");
  InvMass_saHist->getHisto("A_LL")->SetMarkerStyle(21);
  InvMass_saHist->getHisto("A_LL")->SetMarkerColor( kRed );
  InvMass_saHist->getHisto("A_LL")->SetLineWidth(3);
  InvMass_saHist->getHisto("A_LL")->SetMarkerSize(1.2);

  c1->cd(3)->SetGrid();
  c1->Update();

  InvMass_saHist->getHisto("A_L_Blue")->Draw("E1X0");
  InvMass_saHist->getHisto("A_L_Blue")->SetYTitle("Single spin asymmetry (Blue)");
  InvMass_saHist->getHisto("A_L_Blue")->SetMarkerStyle(21);
  InvMass_saHist->getHisto("A_L_Blue")->SetMarkerColor( kRed );
  InvMass_saHist->getHisto("A_L_Blue")->SetLineWidth(3);
  InvMass_saHist->getHisto("A_L_Blue")->SetMarkerSize(1.2);
  c1->cd(4)->SetGrid();
  c1->Update();

  InvMass_saHist->getHisto("A_L_Yellow")->Draw("E1X0");
  InvMass_saHist->getHisto("A_L_Yellow")->SetYTitle("Single spin asymmetry (Yelllow)");
  InvMass_saHist->getHisto("A_L_Yellow")->SetMarkerStyle(21);
  InvMass_saHist->getHisto("A_L_Yellow")->SetMarkerColor( kRed );
  InvMass_saHist->getHisto("A_L_Yellow")->SetLineWidth(3);
  InvMass_saHist->getHisto("A_L_Yellow")->SetMarkerSize(1.2);
  c1->Paint();

}
