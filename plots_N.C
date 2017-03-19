#include "plots_N.h"
#include <math.h>
void plots_N()

{

  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);

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
  gSystem->Load("libpicodst_object.so");

  TFile *file1 = new TFile("/gpfs/mnt/gpfs02/phenix/spin3/spin2/darshana/Analysis/Run13/Simulation_outputs/filtered_DY_13.root");
  TFile *file2 = new TFile("/gpfs/mnt/gpfs02/phenix/spin3/spin2/darshana/Analysis/Run13/Simulation_outputs/filtered_MB_data_13.root");
  TFile *file3 = new TFile("/gpfs/mnt/gpfs02/phenix/spin3/spin2/darshana/Analysis/Run13/Filtered_outputs/fvtx_Run13.root");
  TFile *file4 = new TFile("/gpfs/mnt/gpfs02/phenix/spin3/spin2/darshana/Analysis/Run13/Simulation_outputs/primedatafinal.root");

  t1 = (TTree*)file1->Get("T");
  t2 = (TTree*)file2->Get("T");
  t3 = (TTree*)file3->Get("T");
  t4 = (TTree*)file4->Get("T");

 
 TCut muonCuts("Tr0_DG0 < 20 && Tr1_DG0 < 20 && Tr0_DDG0 < 8 && Tr1_DDG0 < 8"
                "&& Tr0_nidhits > 6 && Tr1_nidhits > 6 && Tr0_ntrhits > 10"
                "&& Tr1_ntrhits > 10 && Tr0_lastgap > 3 && Tr1_lastgap > 3"
                "&& Tr0_pz*Tr1_pz > 0 && abs(Evt_vtxchi2) < 4"
		"&& same_event == 1"
                "&& abs(Tr0_pz)>2 && abs(Tr1_pz)>2"
		"&& Tr0_trchi2<10 && Tr1_trchi2<10");

 TCut fvtxCuts_mass_cut ("mass_fvtxmutr != 0");

 TCut fvtxCuts_0("abs(Tr0_dphi_fvtx) < 10 && abs(Tr0_dr_fvtx) < 10 && abs(Tr0_dtheta_fvtx) < 10");

 TCut fvtxCuts_1("abs(Tr1_dphi_fvtx) < 10 && abs(Tr1_dr_fvtx) < 10 && abs(Tr1_dtheta_fvtx) < 10");

 TCut fvtxCuts_2("abs(Tr0_dphi_fvtx) < 10 && abs(Tr0_dr_fvtx) < 10 && abs(Tr0_dtheta_fvtx) < 10"
	       "&& abs(Tr1_dphi_fvtx) < 10 && abs(Tr1_dr_fvtx) < 10 && abs(Tr1_dtheta_fvtx) < 10 && mass_fvtxmutr != 0");

 TCut fvtxCuts_3 = fvtxCuts_1 || fvtxCuts_0;

 TCut EvertexCut("abs(Evt_fvtxZ)< 10 ");

 TCut charge0_cut (" charge==0");  TCut charge1_cut (" charge!=0");
 TCut charge2_cut (" charge==2");  TCut charge3_cut (" charge==-2");

 TCut massCut("mass > 4 && mass < 8 ");
 TCut massCut2("mass > 3.75 && mass < 8 ");
 TCut massCut3("mass > 4.5 && mass < 8 ");

 TCut armCut_N("Tr0_pz>0");  TCut armCut_S("Tr0_pz<0");

 TCut Trig_N ("fLL1Mu2DNorth==1");  TCut Trig_S ("fLL1Mu2DSouth==1 ");

 TCut cut_eta("abs(rapidity) > 1.2 && abs(rapidity)<2.4");

 TCut bb_cut ("flavor_index==1");  TCut cc_cut ("flavor_index==2");

 TCut mixxed_ct("flavor_index==1 || flavor_index==2 || flavor_index==0");

 TCut RTrig ("(lvl1_trigscaled&0x00010000)>0");

  cut_dy_pm = EvertexCut  && charge0_cut  && massCut && muonCuts && cut_eta && armCut_N ;// && fvtxCuts_2 ;
  cut_bb_pm = EvertexCut  && charge0_cut  && massCut && muonCuts && cut_eta && armCut_N && bb_cut;// && fvtxCuts_2;
  cut_cc_pm = EvertexCut  && charge0_cut  && massCut && muonCuts && cut_eta && armCut_N && cc_cut;// && fvtxCuts_2;
  cut_rd_pm = EvertexCut  && charge0_cut  && massCut3 && muonCuts && cut_eta && armCut_N  && RTrig;// && fvtxCuts_2;
  cut_rd_pp = EvertexCut  && charge2_cut  && massCut && muonCuts && cut_eta && armCut_N  && RTrig;// && fvtxCuts_2;
  cut_rd_mm = EvertexCut  && charge3_cut  && massCut && muonCuts && cut_eta && armCut_N  && RTrig;// && fvtxCuts_2;

  cut_bb_pp = EvertexCut  && charge2_cut  && massCut && muonCuts && cut_eta && armCut_N && bb_cut;// && fvtxCuts_2;
  cut_bb_mm = EvertexCut  && charge3_cut  && massCut && muonCuts && cut_eta && armCut_N && bb_cut;// && fvtxCuts_2;
  cut_pr_pm = EvertexCut  && charge0_cut  && massCut2 && muonCuts && cut_eta && armCut_N;// && fvtxCuts_2 ;
//-------------------------------------   Define Template  Functions ------------------------------
  TCanvas *c_2d_tmp = new TCanvas(TString("c_2d_tmp"), "Define Templates",1600 * .7, 1600 * .7);
  c_2d_tmp->Divide(2,2);
  c_2d_tmp->cd(1);//->SetGrid();
  TH2D *h_2d_dy = new TH2D("h_2d_dy","", 16, 4, 8,50,0,50);
  h_2d_dy->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_dy->GetYaxis()->SetTitle("#of Tracklets");
  h_2d_dy->SetMarkerStyle(21);
  h_2d_dy->SetMarkerSize(.5);
  h_2d_dy->SetMarkerColor( kBlue );
  t1->Project("h_2d_dy","nTracklets:mass",cut_dy_pm);
  h_2d_dy->Sumw2();

  TF2 *f_2d_dy = new TF2("f_2d_dy",_func_2D_dy,4,8,0,50,9);
  f_2d_dy->SetParameter(0,-1.30231);
  f_2d_dy->SetParameter(1,0.054889);
  f_2d_dy->SetParameter(2,0.615495);
  f_2d_dy->SetParameter(3,-0.322594);
  f_2d_dy->SetParameter(4,-1220.4);
  f_2d_dy->SetParameter(5,5265.2);
  f_2d_dy->SetParameter(6,2644.42);
  f_2d_dy->SetParameter(7,-0.0632658);
  f_2d_dy->SetParameter(8,2.10978); 

  //h_2d_dy->Fit(f_2d_dy, "MRLS");
  TFitResultPtr r_dy = h_2d_dy->Fit(f_2d_dy, "MRLS");
  cout<<"Cov matrix determination"<<endl;
  r_dy->Print("V");  
  gMinuit->mnmatu(1);
  f_p0_dy = f_2d_dy->GetParameter(0); 
  f_p1_dy = f_2d_dy->GetParameter(1);
  f_p2_dy = f_2d_dy->GetParameter(2); 
  f_p3_dy = f_2d_dy->GetParameter(3); 
  f_p4_dy = f_2d_dy->GetParameter(4); 
  f_p5_dy = f_2d_dy->GetParameter(5); 
  f_p6_dy = f_2d_dy->GetParameter(6); 
  f_p7_dy = f_2d_dy->GetParameter(7); 
  f_p8_dy = f_2d_dy->GetParameter(8);

  h_2d_dy->SetDrawOption("lego2");
  h_2d_dy->GetXaxis()->SetTitleOffset(1.8);   h_2d_dy->GetYaxis()->SetTitleOffset(1.8);


  c_2d_tmp->cd(2);//->SetGrid();

  TH2D *h_2d_bb_pp = new TH2D("h_2d_bb_pp","", 16, 4, 8,50,0,50);
  h_2d_bb_pp->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_bb_pp->GetYaxis()->SetTitle("#of Tracklets");
  h_2d_bb_pp->SetMarkerStyle(21);
  h_2d_bb_pp->SetMarkerSize(.5);
  h_2d_bb_pp->SetMarkerColor( kBlue );
  t2->Project("h_2d_bb_pp","nTracklets:mass",cut_bb_pp);
  h_2d_bb_pp->Sumw2();

  TH2D *h_2d_bb_mm = new TH2D("h_2d_bb_mm","", 16, 4, 8,50,0,50);
  h_2d_bb_mm->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_bb_mm->GetYaxis()->SetTitle("#of Tracklets");
  h_2d_bb_mm->SetMarkerStyle(21);
  h_2d_bb_mm->SetMarkerSize(.5);
  h_2d_bb_mm->SetMarkerColor( kBlue );
  t2->Project("h_2d_bb_mm","nTracklets:mass",cut_bb_mm);
  h_2d_bb_mm->Sumw2();

  TH2D *h_2d_bb = new TH2D("h_2d_bb","", 16, 4, 8,50,0,50);
  h_2d_bb->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_bb->GetYaxis()->SetTitle("#of Tracklets");
  h_2d_bb->SetMarkerStyle(21);
  h_2d_bb->SetMarkerSize(.5);
  h_2d_bb->SetMarkerColor( kBlue );
  t2->Project("h_2d_bb","nTracklets:mass",cut_bb_pm);
  h_2d_bb->Sumw2();

  TF2 *f_2d_bb = new TF2("f_2d_bb",_func_2D_bb,4,8,0,50,6);
  f_2d_bb->SetParameter(0,-0.4341);
  f_2d_bb->SetParameter(1,-0.02556);
  f_2d_bb->SetParameter(2,9.808);
  f_2d_bb->SetParameter(3,2.656);
  f_2d_bb->SetParameter(4,-178.4);
  f_2d_bb->SetParameter(5,393.2);
  //h_2d_bb->Fit( f_2d_bb , "MRL");

  TFitResultPtr r_bb = h_2d_bb->Fit(f_2d_bb, "MRLS");

  gMinuit->mnmatu(1);
  f_p0_bb = f_2d_bb->GetParameter(0);
  f_p1_bb = f_2d_bb->GetParameter(1);
  f_p2_bb = f_2d_bb->GetParameter(2);
  f_p3_bb = f_2d_bb->GetParameter(3);
  f_p4_bb = f_2d_bb->GetParameter(4);
  f_p5_bb = f_2d_bb->GetParameter(5);

  h_2d_bb->SetDrawOption("lego2");
  h_2d_bb->GetXaxis()->SetTitleOffset(1.8);   h_2d_bb->GetYaxis()->SetTitleOffset(1.8);


  c_2d_tmp->cd(3);//->SetGrid();
  TH2D *h_2d_cc = new TH2D("h_2d_cc","", 16, 4, 8,50,0,50);
  h_2d_cc->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_cc->GetYaxis()->SetTitle("#of Tracklets");
  h_2d_cc->SetMarkerStyle(21);
  h_2d_cc->SetMarkerSize(.5);
  h_2d_cc->SetMarkerColor( kBlue );
  t2->Project("h_2d_cc","nTracklets:mass",cut_cc_pm);
  h_2d_cc->Sumw2();

  TF2 *f_2d_cc = new TF2("f_2d_cc",_func_2D_cc,4,8,0,50,6);
  f_2d_cc->SetParameter(0,-1.408);
  f_2d_cc->SetParameter(1,0.01969);
  f_2d_cc->SetParameter(2,9.708);
  f_2d_cc->SetParameter(3,2.588);
  f_2d_cc->SetParameter(4,10);
  f_2d_cc->SetParameter(5,4806);

 // h_2d_cc->Fit( f_2d_cc , "MRL");
  TFitResultPtr r_cc = h_2d_cc->Fit(f_2d_cc, "MRLS");
  f_2d_cc->Draw("same");
  gMinuit->mnmatu(1);
  f_p0_cc = f_2d_cc->GetParameter(0);
  f_p1_cc = f_2d_cc->GetParameter(1);
  f_p2_cc = f_2d_cc->GetParameter(2);
  f_p3_cc = f_2d_cc->GetParameter(3);
  f_p4_cc = f_2d_cc->GetParameter(4);
  f_p5_cc = f_2d_cc->GetParameter(5);

  h_2d_cc->SetDrawOption("lego2");
  h_2d_cc->GetXaxis()->SetTitleOffset(1.8);   h_2d_cc->GetYaxis()->SetTitleOffset(1.8);


  c_2d_tmp->cd(4);//->SetGrid();
  TH2D *h_2d_pr = new TH2D("h_2d_pr","", 17, 3.75, 8,50,0,50);
  h_2d_pr->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_pr->GetYaxis()->SetTitle("#of Tracklets");
  h_2d_pr->SetMarkerStyle(21);
  h_2d_pr->SetMarkerSize(.5);
  h_2d_pr->SetMarkerColor( kBlue );
  t4->Project("h_2d_pr","nTracklets:mass",cut_pr_pm);
  h_2d_pr->Sumw2();

  TF2 *f_2d_pr = new TF2("f_2d_pr",_func_2D_jp,3.75,8,0,50,6);
  f_2d_pr->SetParameter(0,-1.408);
  f_2d_pr->SetParameter(1,0.01969);
  f_2d_pr->SetParameter(2,9.708);
  f_2d_pr->SetParameter(3,2.588);
  f_2d_pr->SetParameter(4,10);
  f_2d_pr->SetParameter(5,4806);

  h_2d_pr->Fit( f_2d_pr , "MRL");
  TFitResultPtr r_pr = h_2d_pr->Fit(f_2d_pr, "MRLS");
  f_2d_pr->Draw("same");
  gMinuit->mnmatu(1);
  f_p0_pr = f_2d_pr->GetParameter(0);
  f_p1_pr = f_2d_pr->GetParameter(1);
  f_p2_pr = f_2d_pr->GetParameter(2);
  f_p3_pr = f_2d_pr->GetParameter(3);
  f_p4_pr = f_2d_pr->GetParameter(4);
  f_p5_pr = f_2d_pr->GetParameter(5);

  h_2d_pr->SetDrawOption("lego2");
  h_2d_pr->GetXaxis()->SetTitleOffset(1.8);   h_2d_pr->GetYaxis()->SetTitleOffset(1.8);


//-----------------------------------------------Error Matrix-----------------------------------------------

TMatrixDSym cov_dy = r_dy->GetCovarianceMatrix();  
TMatrixDSym cov_bb = r_bb->GetCovarianceMatrix();  
TMatrixDSym cov_cc = r_cc->GetCovarianceMatrix();

r_dy->Print("V");  
r_bb->Print("V");  
r_cc->Print("V");

cov_dy.Print();
cov_bb.Print();
cov_cc.Print();

Double_t a[2] = {4,0};
Double_t b[2] = {8,50};

cout<<"Integral of template for DY: "<<    f_2d_dy->Integral(4,8,0,50,1E-5)/0.25<<"\n" 
    <<"Error of the fit for DY    : "<<f_2d_dy->IntegralError(2,a,b,f_2d_dy->GetParameters(), cov_dy.GetMatrixArray())/0.25<<"\n";

cout<<"Integral of template for bb: "<<    f_2d_bb->Integral(4,8,0,50,1E-5)/0.25<<"\n" 
    <<"Error of the fit for bb    : "<<    f_2d_bb->IntegralError(2,a,b,f_2d_bb->GetParameters(), cov_bb.GetMatrixArray())/0.25<<"\n";

cout<<"Integral of template for cc: "<<    f_2d_cc->Integral(4,8,0,50,1E-5)/0.25<<"\n" 
    <<"Error of the fit for cc    : "<<    f_2d_cc->IntegralError(2,a,b,f_2d_cc->GetParameters(), cov_cc.GetMatrixArray())/0.25<<"\n";



//---------------------------------  Real Data Histograms -----------------------------------------

  TCanvas *c_rd_hits = new TCanvas(TString("c_rd_hits"), "Real Data Histograms",1600 * .7, 1600 * .7);
  c_rd_hits->Divide(3,1);
  c_rd_hits->cd(1);//->SetGrid();
  TH2D *h_2d_rd_pm = new TH2D("h_2d_rd_pm","Unlike Sign (+-) Distribution", 14, 4.5, 8,50,0,50);
  h_2d_rd_pm->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_rd_pm->GetYaxis()->SetTitle("Number of Tracklets");
  h_2d_rd_pm->SetMarkerStyle(21);
  h_2d_rd_pm->SetMarkerSize(.5); h_2d_rd_pm->GetXaxis()->SetTitleOffset(1.8);
  h_2d_rd_pm->SetMarkerColor( kBlue );h_2d_rd_pm->GetYaxis()->SetTitleOffset(1.8);
  t3->Project("h_2d_rd_pm","nTracklets:mass",cut_rd_pm);
  h_2d_rd_pm->Sumw2();
  h_2d_rd_pm->Draw("lego2");

  c_rd_hits->cd(2);//->SetGrid();
  TH2D *h_2d_rd_pp = new TH2D("h_2d_rd_pp","Like Sign (++) Distribution", 16, 4, 8,50,0,50);
  h_2d_rd_pp->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_rd_pp->GetYaxis()->SetTitle("Number of Tracklets");
  h_2d_rd_pp->SetMarkerStyle(21);
  h_2d_rd_pp->SetMarkerSize(.5);h_2d_rd_pp->GetXaxis()->SetTitleOffset(1.8);
  h_2d_rd_pp->SetMarkerColor( kBlue );h_2d_rd_pp->GetYaxis()->SetTitleOffset(1.8);
  t3->Project("h_2d_rd_pp","nTracklets:mass",cut_rd_pp);
  h_2d_rd_pp->Sumw2();
  h_2d_rd_pp->Draw("lego2");

  c_rd_hits->cd(3);//->SetGrid();
  TH2D *h_2d_rd_mm = new TH2D("h_2d_rd_mm","Like sign (--) Distribution", 16, 4, 8,50,0,50);
  h_2d_rd_mm->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_rd_mm->GetYaxis()->SetTitle("Number of Tracklets");
  h_2d_rd_mm->SetMarkerStyle(21);
  h_2d_rd_mm->SetMarkerSize(.5);h_2d_rd_mm->GetXaxis()->SetTitleOffset(1.8);
  h_2d_rd_mm->SetMarkerColor( kBlue );h_2d_rd_mm->GetYaxis()->SetTitleOffset(1.8);
  t3->Project("h_2d_rd_mm","nTracklets:mass",cut_rd_mm);
  h_2d_rd_mm->Sumw2();
  h_2d_rd_mm->Draw("lego2");


  TCanvas *c_final_fit = new TCanvas(TString("Mass_vs_trk7"), "Mass_vs_trk7",1600 * .7, 1600 * .7);
  c_final_fit->cd();

  TH2D *h_2d_final = new TH2D("h_2d_final","", 46, 4.5, 16,150,0,150);
  h_2d_final->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_final->GetYaxis()->SetTitle("Number of Tracklets");
  h_2d_final->SetMarkerStyle(21);
  h_2d_final->SetMarkerSize(.5);
  h_2d_final->SetMarkerColor( kBlue );
  h_2d_final->GetXaxis()->SetTitleOffset(1.8);
  h_2d_final->GetYaxis()->SetTitleOffset(1.8);

  int total_test = 0;


  for (int i=1; i< 15; i++){
        for(int j = 1; j < 51; j++){

                h_2d_final->GetBin(i,j);
                h_2d_final->SetBinContent(i,j, h_2d_rd_pm->GetBinContent(i,j));
                h_2d_final->SetBinError(i,j, h_2d_rd_pm->GetBinError(i,j));
}}

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){
                h_2d_final->SetBinContent(i+14,j+50, h_2d_rd_pp->GetBinContent(i,j));
                h_2d_final->SetBinError(i+14,j+50, h_2d_rd_pp->GetBinError(i,j));

                h_2d_final->SetBinContent(i+30,j+100, h_2d_rd_mm->GetBinContent(i,j));
                h_2d_final->SetBinError(i+30,j+100, h_2d_rd_mm->GetBinError(i,j));

    total_test++;   

 
        }

}

  h_2d_final->Draw("surf1");

  cout<<" total fills = "<< total_test <<endl;

  TF2 *f_2d_1 = new TF2("f_2d_1",_func_final,4.5,16,0,150,9);

  f_2d_1->SetParameter(0,0.0387407);//f_2d_1->SetParLimits(0,0.0,1);
  f_2d_1->SetParameter(1,0.196813); //f_2d_1->SetParLimits(1,0.1,1);
  f_2d_1->SetParameter(2,0.239846); //f_2d_1->SetParLimits(2,0.2,1);
  f_2d_1->SetParameter(3,264499);//f_2d_1->SetParLimits(3,1,1);
  f_2d_1->SetParameter(4,-2.59109); 
  f_2d_1->SetParameter(5,0.12214); 
  f_2d_1->SetParameter(6,7.15352);
  f_2d_1->SetParameter(7,231479); 
  //f_2d_1->FixParameter(8,51.08);
  f_2d_1->SetParameter(8,105.6);//f_2d_1->SetParLimits(8,100,200);
 // f_2d_1->SetParameter(9,2.63435);  

  //h_2d_final->Fit( f_2d_1 , "MRL");
  TFitResultPtr r_sm = h_2d_final->Fit(f_2d_1, "MRLS");

  p0_fnl  = f_2d_1->GetParameter(0);
  p1_fnl  = f_2d_1->GetParameter(1);
  p2_fnl  = f_2d_1->GetParameter(2);
  p3_fnl  = f_2d_1->GetParameter(3);
  p4_fnl  = f_2d_1->GetParameter(4);
  p5_fnl  = f_2d_1->GetParameter(5);
  p6_fnl  = f_2d_1->GetParameter(6);
  p7_fnl  = f_2d_1->GetParameter(7);
  p8_fnl  = f_2d_1->GetParameter(8);
  //p9_fnl  = f_2d_1->GetParameter(9);
  //p10_fnl  = f_2d_1->GetParameter(10);

  p0_fnl_err  = f_2d_1->GetParError(0);
  p1_fnl_err  = f_2d_1->GetParError(1);
  p2_fnl_err  = f_2d_1->GetParError(2);
  p3_fnl_err  = f_2d_1->GetParError(3);
  p4_fnl_err  = f_2d_1->GetParError(4);
  p5_fnl_err  = f_2d_1->GetParError(5);
  p6_fnl_err  = f_2d_1->GetParError(6);
  p7_fnl_err  = f_2d_1->GetParError(7);
  p8_fnl_err  = f_2d_1->GetParError(8);
 // p9_fnl_err  = f_2d_1->GetParError(9);

  h_2d_final->SetDrawOption("lego2");

  TMatrixDSym cov_sm = r_sm->GetCovarianceMatrix();    
  r_sm->Print("V");    

 //-------------------------- Draw Combinatorial background function ------------------------

  TCanvas *c_comb = new TCanvas("c_comb","combinatorial",0,0,1600,1200);
  c_comb->cd();//->SetGrid();

  TF2 *f_comb = new TF2("f_comb","  exp([0]*x[0] + [1]*x[0]*x[0])*TMath::Poisson(x[1],[2])",4,8,0,50); 
  f_comb->SetParameters(p4_fnl,p5_fnl,p6_fnl);
  //Double_t err_cmb = sqrt(  (1/p3_fnl)*(p3_fnl_err*p3_fnl_err)+ (p3_fnl/1)*(1*1)  );  
  f_comb->SetContour(48);
  f_comb->SetFillColor(45);
  f_comb->Draw("surf1");
  f_comb->GetHistogram()->SetTitle("Final function for Combinatorial Background");
  f_comb->GetHistogram()->GetXaxis()->SetTitleOffset(1.8);   f_comb->GetHistogram()->GetYaxis()->SetTitleOffset(1.8);
  f_comb->GetHistogram()->GetXaxis()->SetTitle("Mass (GeV)");
  f_comb->GetHistogram()->GetYaxis()->SetTitle("Number of Tracklets");

//------------------------- Evaluate functional count for each process----------------------------
   Double_t err_cmb = sqrt(  (p7_fnl/p3_fnl)*(p3_fnl_err*p3_fnl_err)+ (p3_fnl/p7_fnl)*(p7_fnl_err*p7_fnl_err)  ); 
   double in_d_y = f_2d_dy->Integral(4.50,8,0,50,1E-5);
   cout << "integral_for_dy   = "<< in_d_y*p0_fnl/0.25 << "  err = " << in_d_y*p0_fnl_err/.25 <<endl;
   cout << "integral_for_dy_of template   = "<< in_d_y/0.25<<endl;

   double in_d_b = f_2d_bb->Integral(4.50,8,0,50,1E-5);
   cout << "integral_for_bb   = "<< in_d_b*p1_fnl/0.25 << "  err = " << in_d_b*p1_fnl_err/.25 <<endl;

   double in_d_c = f_2d_cc->Integral(4.50,8,0,50,1E-5);
   cout << "integral_for_cc   = "<< in_d_c*p2_fnl/0.25 << "  err = " << in_d_c*p2_fnl_err/.25 <<endl;
   double in_d_m = f_comb->Integral(4.50,8,0,50,1E-5);
   cout << "integral_for_comb = "<< in_d_m*sqrt(p3_fnl*p7_fnl)*2/0.25 << " err = " << in_d_m*err_cmb/.25  <<endl;

   TF1 *fit_trk_prime2 = new TF1("fit_trk_prime2",_func_prime2,0,50,1); fit_trk_prime2->SetLineColor(2);
   fit_trk_prime2->SetParameter(0,p8_fnl);

   cout << "4 - 4.25 psi prime= "<<f_2d_pr->Integral(4.50,8,0,50,1E-5)*p8_fnl/0.25 << "  err = " <<f_2d_pr->Integral(4.50,8,0,50,1E-5)*p8_fnl_err/.25 <<endl;

   cout <<" total = "<<(f_2d_pr->Integral(4.50,8,0,50,1E-5)*p8_fnl+in_d_m*sqrt(p3_fnl*p7_fnl)*2+in_d_b*p1_fnl+in_d_c*p2_fnl+in_d_y*p0_fnl)/.25<<"integral is = "<< f_2d_1->Integral(4.50,8,0,50,1E-5)/.25<<endl;

//-------------------------Projections of functions and simulations--------------------------------

//------------------------- Fill Histograms with the templates Found ------------------------------

  TCanvas *c_fill_fun = new TCanvas("c_fill_fun","c_fill_fun",0,0,1600,1200);
  c_fill_fun->Divide(3,2);
  c_fill_fun->cd(1)->SetGrid();
  TH2D *_h2_dy_f_fill = new TH2D("_h2_dy_f_fill","", 16, 4, 8,50,0,50);
  _h2_dy_f_fill->GetXaxis()->SetTitle("Mass (GeV)");
  _h2_dy_f_fill->GetYaxis()->SetTitle("Tracklets");
  _h2_dy_f_fill->SetMarkerStyle(21);  _h2_dy_f_fill->SetLineColor( kRed );
  _h2_dy_f_fill->SetMarkerSize(.5);   _h2_dy_f_fill->Sumw2();
  _h2_dy_f_fill->SetMarkerColor( kRed );


  double e_bin_dy = 0;  double e_bin_bb = 0;  double e_bin_cc = 0;

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                _h2_dy_f_fill->GetBin(i,j);
                _h2_dy_f_fill->SetBinContent(i,j, f_2d_dy->Integral(4+(i-1)*0.25,4+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);
		//Double_t a_dy[2] = {4+(i-1)*0.25,j-1};
		//Double_t b_dy[2] = {4+(i)*0.25,j};
      		//e_bin_dy = f_2d_dy->IntegralError(2,a_dy,b_dy,f_2d_dy->GetParameters(), cov_dy.GetMatrixArray())/0.25;
		e_bin_dy = h_2d_dy->GetBinError(i,j);
                _h2_dy_f_fill->SetBinError(i,j,e_bin_dy);
}}
  _h2_dy_f_fill->Draw("lego2");

  cout<< "Template histogram integral for dy is = " << _h2_dy_f_fill->Integral(0,16,0,50)<< 
         "Simulation actual integral  for dy is = " << h_2d_dy->Integral(0,16,0,50)<<endl;

  c_fill_fun->cd(2)->SetGrid();
  TH2D *_h2_bb_f_fill = new TH2D("_h2_bb_f_fill","", 16, 4, 8,50,0,50);
  _h2_bb_f_fill->GetXaxis()->SetTitle("Mass (GeV)");
  _h2_bb_f_fill->GetYaxis()->SetTitle("Tracklets");
  _h2_bb_f_fill->SetMarkerStyle(21);  _h2_bb_f_fill->SetLineColor( kRed );
  _h2_bb_f_fill->SetMarkerSize(.5);   _h2_bb_f_fill->Sumw2();
  _h2_bb_f_fill->SetMarkerColor( kRed );

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                _h2_bb_f_fill->GetBin(i,j);
                _h2_bb_f_fill->SetBinContent(i,j, f_2d_bb->Integral(4+(i-1)*0.25,4+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);
		//Double_t a_bb[2] = {4+(i-1)*0.25,j-1};
		//Double_t b_bb[2] = {4+(i)*0.25,j};
      		//e_bin_bb = f_2d_bb->IntegralError(2,a_bb,b_bb,f_2d_bb->GetParameters(), cov_bb.GetMatrixArray())/0.25;
		e_bin_bb = h_2d_bb->GetBinError(i,j);
                _h2_bb_f_fill->SetBinError(i,j,e_bin_bb);
}}
  _h2_bb_f_fill->Draw("lego2");

  cout<< "Template histogram integral for bb is = " << _h2_bb_f_fill->Integral(0,16,0,50)<< 
         "Simulation actual integral  for bb is = "<< h_2d_bb->Integral(0,16,0,50)<<endl;

  c_fill_fun->cd(3)->SetGrid();
  TH2D *_h2_cc_f_fill = new TH2D("_h2_cc_f_fill","", 16, 4, 8,50,0,50);
  _h2_cc_f_fill->GetXaxis()->SetTitle("Mass (GeV)");
  _h2_cc_f_fill->GetYaxis()->SetTitle("Tracklets");
  _h2_cc_f_fill->SetMarkerStyle(21);  _h2_cc_f_fill->SetLineColor( kRed );
  _h2_cc_f_fill->SetMarkerSize(.5);   _h2_cc_f_fill->Sumw2();
  _h2_cc_f_fill->SetMarkerColor( kRed );

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                _h2_cc_f_fill->GetBin(i,j);
                _h2_cc_f_fill->SetBinContent(i,j, f_2d_cc->Integral(4+(i-1)*0.25,4+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);
		//Double_t a_cc[2] = {4+(i-1)*0.25,j-1};
		//Double_t b_cc[2] = {4+(i)*0.25,j};
      		//e_bin_cc = f_2d_cc->IntegralError(2,a_cc,b_cc,f_2d_cc->GetParameters(), cov_cc.GetMatrixArray())/0.25;
		e_bin_cc = h_2d_cc->GetBinError(i,j);
                _h2_cc_f_fill->SetBinError(i,j,e_bin_cc);
}}

  _h2_cc_f_fill->Draw("lego2");

  cout<< "Template histogram integral for cc is = " << _h2_cc_f_fill->Integral(0,16,0,50)<<
         "Simulation actual integral  for cc is = "<< h_2d_cc->Integral(0,16,0,50)<<endl;


  c_fill_fun->cd(4)->SetGrid();
  TH2D *_h2_cb_f_fill = new TH2D("_h2_cb_f_fill","", 16, 4, 8,50,0,50);
  _h2_cb_f_fill->GetXaxis()->SetTitle("Mass (GeV)");
  _h2_cb_f_fill->GetYaxis()->SetTitle("Tracklets");
  _h2_cb_f_fill->SetMarkerStyle(21);  _h2_cb_f_fill->SetLineColor( kRed );
  _h2_cb_f_fill->SetMarkerSize(.5);   _h2_cb_f_fill->Sumw2();
  _h2_cb_f_fill->SetMarkerColor( kRed );

  double e_bin_cb = 0; 

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                _h2_cb_f_fill->GetBin(i,j);
                _h2_cb_f_fill->SetBinContent(i,j, f_comb->Integral(4+(i-1)*0.25,4+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);

		e_bin_cb = sqrt(_h2_cb_f_fill->GetBinContent(i,j));
                _h2_cb_f_fill->SetBinError(i,j,e_bin_cb);
}}
  _h2_cb_f_fill->Draw("lego2");




  cout<< "Final Fitting determined CB integral  = " << _h2_cb_f_fill->Integral(0,16,0,50)*sqrt(p3_fnl*p7_fnl)*2<<endl;
 // cout<< "Final Fitting determined Prime Total  = " << fit_trk_prime->Integral(0,50) + fit_trk_prime2->Integral(0,50)<<endl;

  c_fill_fun->cd(5)->SetGrid();
  TH2D *_h2_pr_f_fill = new TH2D("_h2_pr_f_fill","", 16, 4, 8,50,0,50);
  _h2_pr_f_fill->GetXaxis()->SetTitle("Mass (GeV)");
  _h2_pr_f_fill->GetYaxis()->SetTitle("Tracklets");
  _h2_pr_f_fill->SetMarkerStyle(21);  _h2_pr_f_fill->SetLineColor( kRed );
  _h2_pr_f_fill->SetMarkerSize(.5);   _h2_pr_f_fill->Sumw2();
  _h2_pr_f_fill->SetMarkerColor( kRed );

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                _h2_pr_f_fill->GetBin(i,j);
                _h2_pr_f_fill->SetBinContent(i,j, f_2d_pr->Integral(4+(i-1)*0.25,4+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);
	
		double e_bin_pr = h_2d_pr->GetBinError(i+1,j);
                _h2_pr_f_fill->SetBinError(i,j,e_bin_pr);
}}

  _h2_pr_f_fill->Draw("lego2");


//----------------------------Projection to mass--------------------------------------
   TCanvas *c_project_fnc = new TCanvas("c_project_fnc", "c_project_fnc",0,0,1600,1200);
   //gStyle->SetOptStat(0);
   c_project_fnc->Divide(3,2,0,0);

   c_project_fnc->cd(1); //->SetGrid();
   TH1D * h_m_dy_pjx = new TH1D("h_m_dy_pjx", "",  16, 4, 8); h_m_dy_pjx->SetMarkerColor( kBlue );
   h_m_dy_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_dy_pjx->SetLineWidth(0.6);
   h_m_dy_pjx = h_2d_dy->ProjectionX("h_m_dy_pjx", 0, 50);
   h_m_dy_pjx->Draw("E");

   TH1D * _h_m_dy_pjx = new TH1D("_h_m_dy_pjx", "",  16, 4, 8); _h_m_dy_pjx->SetMarkerColor( kRed );
   _h_m_dy_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_dy_pjx->SetLineWidth(0.6);
   _h_m_dy_pjx = _h2_dy_f_fill->ProjectionX("_h_m_dy_pjx ", 0, 50);
   _h_m_dy_pjx->Draw("SAME L E1");

   TLegend *leg_dyx = new TLegend(0.65, 0.65, .85, .85);
   leg_dyx->SetHeader("Drell-Yan Simulations");
   leg_dyx->AddEntry(h_m_dy_pjx,"Simulated Data","lp");
   leg_dyx->AddEntry(_h_m_dy_pjx,"Function Determined","lp");
   leg_dyx->Draw();

   c_project_fnc->cd(4); //->SetGrid();
   TH1D * _h_m_poox_dy = new TH1D("_h_m_poox_dy", "",  16, 4, 8); _h_m_poox_dy->SetMarkerColor( kBlue );
   _h_m_poox_dy->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_poox_dy->SetLineWidth(0.6);
   _h_m_poox_dy->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); _h_m_poox_dy->SetFillColor(15);_h_m_poox_dy->GetYaxis()->CenterTitle();
for (int i=1; i< 17; i++){
                _h_m_poox_dy->GetBin(i);
		if((h_m_dy_pjx->GetBinContent(i)- _h_m_dy_pjx->GetBinContent(i))==0 || h_m_dy_pjx->GetBinError(i) ==0){
		_h_m_poox_dy->SetBinContent(i,0);
                _h_m_poox_dy->SetBinError(i,0);
}

if((h_m_dy_pjx->GetBinContent(i)- _h_m_dy_pjx->GetBinContent(i))!=0 && h_m_dy_pjx->GetBinError(i) !=0){
_h_m_poox_dy->SetBinContent(i,(h_m_dy_pjx->GetBinContent(i)- _h_m_dy_pjx->GetBinContent(i))/h_m_dy_pjx->GetBinError(i));
_h_m_poox_dy->SetBinError(i,0);}

}
   _h_m_poox_dy->Draw();


   c_project_fnc->cd(2); //->SetGrid();
   TH1D * h_m_bb_pjx = new TH1D("h_m_bb_pjx", "",  16, 4, 8); h_m_bb_pjx->SetMarkerColor( kBlue );
   h_m_bb_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_bb_pjx->SetLineWidth(0.6);
   h_m_bb_pjx = h_2d_bb->ProjectionX("h_m_bb_pjx ", 0, 50);
   h_m_bb_pjx->Draw(); 

   TH1D * _h_m_bb_pjx = new TH1D("_h_m_bb_pjx", "",  16, 4, 8); _h_m_bb_pjx->SetMarkerColor( kRed );
   _h_m_bb_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_bb_pjx->SetLineWidth(0.6);
   _h_m_bb_pjx = _h2_bb_f_fill->ProjectionX("_h_m_bb_pjx ", 0, 50);
   _h_m_bb_pjx->Draw("SAME L E1");  

   TLegend *leg_bbx = new TLegend(0.65, 0.65, .85, .85);
   leg_bbx->SetHeader("b#bar{b} Simulation");
   leg_bbx->AddEntry(h_m_bb_pjx,"Simulated Data","lp");
   leg_bbx->AddEntry(_h_m_bb_pjx,"Function Determined","lp");
   leg_bbx->Draw();

   c_project_fnc->cd(5); //->SetGrid();
   TH1D * _h_m_poox_bb = new TH1D("_h_m_poox_bb", "",  16, 4, 8); _h_m_poox_bb->SetMarkerColor( kBlue );
   _h_m_poox_bb->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_poox_bb->SetLineWidth(0.6);
   _h_m_poox_bb->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); _h_m_poox_bb->SetFillColor(15);_h_m_poox_bb->GetYaxis()->CenterTitle();
for (int i=1; i< 17; i++){
                _h_m_poox_bb->GetBin(i);
		if((h_m_bb_pjx->GetBinContent(i)- _h_m_bb_pjx->GetBinContent(i))==0 || h_m_bb_pjx->GetBinError(i) ==0){
		_h_m_poox_bb->SetBinContent(i,0);
                _h_m_poox_bb->SetBinError(i,0);
}

if((h_m_bb_pjx->GetBinContent(i)- _h_m_bb_pjx->GetBinContent(i))!=0 && h_m_bb_pjx->GetBinError(i) !=0){
_h_m_poox_bb->SetBinContent(i,(h_m_bb_pjx->GetBinContent(i)- _h_m_bb_pjx->GetBinContent(i))/h_m_bb_pjx->GetBinError(i));
_h_m_poox_bb->SetBinError(i,0);}

}
   _h_m_poox_bb->Draw();

   c_project_fnc->cd(3); //->SetGrid();
   TH1D * h_m_cc_pjx = new TH1D("h_m_cc_pjx", "",  16, 4, 8); h_m_cc_pjx->SetMarkerColor( kBlue );
   h_m_cc_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_cc_pjx->SetLineWidth(0.6);
   h_m_cc_pjx = h_2d_cc->ProjectionX("h_m_cc_pjx ", 0, 50);
   h_m_cc_pjx->Draw();   

   TH1D * _h_m_cc_pjx = new TH1D("_h_m_cc_pjx", "",  16, 4, 8); _h_m_cc_pjx->SetMarkerColor( kRed );
   _h_m_cc_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_cc_pjx->SetLineWidth(0.6);
   _h_m_cc_pjx = _h2_cc_f_fill->ProjectionX(" _h_m_cc_pjx", 0, 50);
   _h_m_cc_pjx->Draw("SAME L E1"); 

   TLegend *leg_ccx = new TLegend(0.65, 0.65, .85, .85);
   leg_ccx->SetHeader("c#bar{c} Simulation");
   leg_ccx->AddEntry(h_m_cc_pjx,"Simulated Data","lp");
   leg_ccx->AddEntry(_h_m_cc_pjx,"Function Determined","lp");
   leg_ccx->Draw();

   c_project_fnc->cd(6); //->SetGrid();
   TH1D * _h_m_poox_cc = new TH1D("_h_m_poox_cc", "",  16, 4, 8); _h_m_poox_cc->SetMarkerColor( kBlue );
   _h_m_poox_cc->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_poox_cc->SetLineWidth(0.6);
   _h_m_poox_cc->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); _h_m_poox_cc->SetFillColor(15);_h_m_poox_cc->GetYaxis()->CenterTitle();
for (int i=1; i< 17; i++){
                _h_m_poox_cc->GetBin(i);
		if((h_m_cc_pjx->GetBinContent(i)- _h_m_cc_pjx->GetBinContent(i))==0 || h_m_cc_pjx->GetBinError(i) ==0){
		_h_m_poox_cc->SetBinContent(i,0);
                _h_m_poox_cc->SetBinError(i,0);
}

if((h_m_cc_pjx->GetBinContent(i)- _h_m_cc_pjx->GetBinContent(i))!=0 && h_m_cc_pjx->GetBinError(i) !=0){
_h_m_poox_cc->SetBinContent(i,(h_m_cc_pjx->GetBinContent(i)- _h_m_cc_pjx->GetBinContent(i))/h_m_cc_pjx->GetBinError(i));
_h_m_poox_cc->SetBinError(i,0);}

}
   _h_m_poox_cc->Draw();

//-----------------------Projection to Tracklet--------------------------------------------------------


   TCanvas *c_project_fnc_2 = new TCanvas("c_project_fnc_2", "c_project_fnc_2",0,0,1600,1200);
   //gStyle->SetOptStat(0);
   c_project_fnc_2->Divide(3,2,0,0);

   c_project_fnc_2->cd(1); //->SetGrid();
   TH1D * h_m_dy_pjy = new TH1D("h_m_dy_pjy", "",  50, 0, 50); h_m_dy_pjy->SetMarkerColor( kBlue );
   h_m_dy_pjy->GetXaxis()->SetTitle("# of Tracklets"); h_m_dy_pjy->SetLineWidth(0.6);
   h_m_dy_pjy = h_2d_dy->ProjectionY(" h_m_dy_pjy", 0, 16);
   h_m_dy_pjy->Draw();  

   TH1D * _h_m_dy_pjy = new TH1D("_h_m_dy_pjy", "",  50, 0, 50); _h_m_dy_pjy->SetMarkerColor( kRed );
   _h_m_dy_pjy->GetXaxis()->SetTitle("# of Tracklets"); _h_m_dy_pjy->SetLineWidth(0.6);
   _h_m_dy_pjy = _h2_dy_f_fill->ProjectionY("_h_m_dy_pjy ", 0, 16);
   _h_m_dy_pjy->Draw("SAME L E1");
 
   TLegend *leg_dyy = new TLegend(0.65, 0.65, .85, .85);
   leg_dyy->SetHeader("Drell-Yan Simulations");
   leg_dyy->AddEntry(h_m_dy_pjy,"Simulated Data","lp");
   leg_dyy->AddEntry(_h_m_dy_pjy,"Function Determined","lp");
   leg_dyy->Draw();

   c_project_fnc_2->cd(2); //->SetGrid();
   TH1D * h_m_bb_pjy = new TH1D("h_m_bb_pjy", "",  50, 0, 50); h_m_bb_pjy->SetMarkerColor( kBlue );
   h_m_bb_pjy->GetXaxis()->SetTitle("# of Tracklets"); h_m_bb_pjy->SetLineWidth(0.6);
   h_m_bb_pjy = h_2d_bb->ProjectionY(" h_m_bb_pjy", 0, 16);
   h_m_bb_pjy->Draw();  

   TH1D * _h_m_bb_pjy = new TH1D("_h_m_bb_pjy", "",  50, 0, 50); _h_m_bb_pjy->SetMarkerColor( kRed );
   _h_m_bb_pjy->GetXaxis()->SetTitle("# of Tracklets"); _h_m_bb_pjy->SetLineWidth(0.6);
   _h_m_bb_pjy = _h2_bb_f_fill->ProjectionY("_h_m_bb_pjy ", 0, 16);
   _h_m_bb_pjy->Draw("SAME L E1");  

   TLegend *leg_bby = new TLegend(0.65, 0.65, .85, .85);
   leg_bby->SetHeader("b#bar{b} Simulation");
   leg_bby->AddEntry(h_m_bb_pjy,"Simulated Data","lp");
   leg_bby->AddEntry(_h_m_bb_pjy,"Function Determined","lp");
   leg_bby->Draw();

   c_project_fnc_2->cd(3); //->SetGrid();
   TH1D * h_m_cc_pjy = new TH1D("h_m_cc_pjy", "",  50, 0, 50); h_m_cc_pjy->SetMarkerColor( kBlue );
   h_m_cc_pjy->GetXaxis()->SetTitle("# of Tracklets"); h_m_cc_pjy->SetLineWidth(0.6);
   h_m_cc_pjy = h_2d_cc->ProjectionY(" h_m_cc_pjy", 0, 16);
   h_m_cc_pjy->Draw(); 

   TH1D * _h_m_cc_pjy = new TH1D("_h_m_cc_pjy", "",  50, 0, 50); _h_m_cc_pjy->SetMarkerColor( kRed );
   _h_m_cc_pjy->GetXaxis()->SetTitle("# of Tracklets"); _h_m_cc_pjy->SetLineWidth(0.6);
   _h_m_cc_pjy = _h2_cc_f_fill->ProjectionY(" ", 0, 16);
   _h_m_cc_pjy->Draw("SAME L E1");

   TLegend *leg_ccy = new TLegend(0.65, 0.65, .85, .85);
   leg_ccy->SetHeader("c#bar{c} Simulation");
   leg_ccy->AddEntry(h_m_cc_pjy,"Simulated Data","lp");
   leg_ccy->AddEntry(_h_m_cc_pjy,"Function Determined","lp");
   leg_ccy->Draw();

   c_project_fnc_2->cd(4); //->SetGrid();

   TH1D * _h_m_pooy_dy = new TH1D("_h_m_pooy_dy", "",  50, 0, 50); _h_m_pooy_dy->SetMarkerColor( kRed );
   _h_m_pooy_dy->GetXaxis()->SetTitle("# of Tracklets"); _h_m_pooy_dy->SetLineWidth(0.6);
   _h_m_pooy_dy->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); _h_m_pooy_dy->SetFillColor(15);_h_m_pooy_dy->GetYaxis()->CenterTitle();

for (int i=1; i< 51; i++){
                _h_m_pooy_dy->GetBin(i);
		if((h_m_dy_pjy->GetBinContent(i)- _h_m_dy_pjy->GetBinContent(i))==0 || h_m_dy_pjy->GetBinError(i) ==0){
		_h_m_pooy_dy->SetBinContent(i,0);
                _h_m_pooy_dy->SetBinError(i,0);}

if((h_m_dy_pjy->GetBinContent(i)- _h_m_dy_pjy->GetBinContent(i))!=0 && h_m_dy_pjy->GetBinError(i) !=0){
_h_m_pooy_dy->SetBinContent(i,(h_m_dy_pjy->GetBinContent(i)- _h_m_dy_pjy->GetBinContent(i))/h_m_dy_pjy->GetBinError(i));
_h_m_pooy_dy->SetBinError(i,0);}

}

_h_m_pooy_dy->Draw();

   c_project_fnc_2->cd(5); //->SetGrid();

   TH1D * _h_m_pooy_bb = new TH1D("_h_m_pooy_bb", "",  50, 0, 50); _h_m_pooy_bb->SetMarkerColor( kRed );
   _h_m_pooy_bb->GetXaxis()->SetTitle("# of Tracklets"); _h_m_pooy_bb->SetLineWidth(0.6);
   _h_m_pooy_bb->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); _h_m_pooy_bb->SetFillColor(15);_h_m_pooy_bb->GetYaxis()->CenterTitle();

for (int i=1; i< 51; i++){
                _h_m_pooy_bb->GetBin(i);
		if((h_m_bb_pjy->GetBinContent(i)- _h_m_bb_pjy->GetBinContent(i))==0 || h_m_bb_pjy->GetBinError(i) ==0){
		_h_m_pooy_bb->SetBinContent(i,0);
                _h_m_pooy_bb->SetBinError(i,0);}

if((h_m_bb_pjy->GetBinContent(i)- _h_m_bb_pjy->GetBinContent(i))!=0 && h_m_bb_pjy->GetBinError(i) !=0){
_h_m_pooy_bb->SetBinContent(i,(h_m_bb_pjy->GetBinContent(i)- _h_m_bb_pjy->GetBinContent(i))/h_m_bb_pjy->GetBinError(i));
_h_m_pooy_bb->SetBinError(i,0);}

}

_h_m_pooy_bb->Draw();

   c_project_fnc_2->cd(6); //->SetGrid();

   TH1D * _h_m_pooy_cc = new TH1D("_h_m_pooy_cc", "",  50, 0, 50); _h_m_pooy_cc->SetMarkerColor( kRed );
   _h_m_pooy_cc->GetXaxis()->SetTitle("# of Tracklets"); _h_m_pooy_cc->SetLineWidth(0.6);
   _h_m_pooy_cc->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); _h_m_pooy_cc->SetFillColor(15);_h_m_pooy_cc->GetYaxis()->CenterTitle();

for (int i=1; i< 51; i++){
                _h_m_pooy_cc->GetBin(i);
		if((h_m_cc_pjy->GetBinContent(i)- _h_m_cc_pjy->GetBinContent(i))==0 || h_m_cc_pjy->GetBinError(i) ==0){
		_h_m_pooy_cc->SetBinContent(i,0);
                _h_m_pooy_cc->SetBinError(i,0);}

if((h_m_cc_pjy->GetBinContent(i)- _h_m_cc_pjy->GetBinContent(i))!=0 && h_m_cc_pjy->GetBinError(i) !=0){
_h_m_pooy_cc->SetBinContent(i,(h_m_cc_pjy->GetBinContent(i)- _h_m_cc_pjy->GetBinContent(i))/h_m_cc_pjy->GetBinError(i));
_h_m_pooy_cc->SetBinError(i,0);}

}

_h_m_pooy_cc->Draw();


  h_m_dy_pjx->GetXaxis()->SetTitleOffset(1);  _h_m_dy_pjx->GetXaxis()->SetTitleOffset(1);
  h_m_dy_pjy->GetXaxis()->SetTitleOffset(1);  _h_m_dy_pjx->GetXaxis()->SetTitleOffset(1);
  h_m_bb_pjx->GetXaxis()->SetTitleOffset(1);  _h_m_bb_pjx->GetXaxis()->SetTitleOffset(1);
  h_m_bb_pjy->GetXaxis()->SetTitleOffset(1);  _h_m_bb_pjx->GetXaxis()->SetTitleOffset(1);
  h_m_cc_pjx->GetXaxis()->SetTitleOffset(1);  _h_m_cc_pjx->GetXaxis()->SetTitleOffset(1);
  h_m_cc_pjy->GetXaxis()->SetTitleOffset(1);  _h_m_cc_pjx->GetXaxis()->SetTitleOffset(1);

//------------------------------------------------act data fit projection---------------------

  TCanvas *_c_f_fit = new TCanvas("_c_f_fit","_c_f_fit",0,0,1600,1200);
  _c_f_fit->cd()->SetGrid();

  TH2D *_h2_sim_f_fil = new TH2D("_h2_sim_f_fil","", 46, 4.50, 16,150,0,150);
  _h2_sim_f_fil->GetXaxis()->SetTitle("Mass (GeV)");
  _h2_sim_f_fil->GetYaxis()->SetTitle("Tracklets");
  _h2_sim_f_fil->SetMarkerStyle(21);  _h2_sim_f_fil->SetLineColor( kRed );
  _h2_sim_f_fil->SetMarkerSize(.5);   _h2_sim_f_fil->Sumw2();
  _h2_sim_f_fil->SetMarkerColor( kRed );
  r_sm->Print("V"); 
  double e_bin_sim = 0;

  for (int i=1; i< 47; i++){
        for(int j = 1; j < 151; j++){

                _h2_sim_f_fil->GetBin(i,j);
                _h2_sim_f_fil->SetBinContent(i,j, f_2d_1->Integral(4.50+(i-1)*0.25,4.50+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);
                e_bin_sim = sqrt(f_2d_1->Integral(4.50+(i-1)*0.25,4.50+(i)*0.25,(j-1)+0.0,j,1E-5)/0.25);
                _h2_sim_f_fil->SetBinError(i,j,e_bin_sim);
}}
  _h2_sim_f_fil->Draw("lego2");

  cout<< "new histogram integral for data fit = " << _h2_sim_f_fil->Integral(0,14,0,50)<< " Actual Data integral= "
  << h_2d_final->Integral(0,14,0,50)<<"integral of function ="<< f_2d_1->Integral(4.50,8,0,50,1E-5)/.25<<endl;


  TCanvas *c_f_proj = new TCanvas(TString("Final fit Projections"), "Final fit Projections",1600 * .7, 1600 * .7);
  c_f_proj->Divide(2,1);

  c_f_proj->cd(1)->SetGrid();
  gPad->SetLogy();
   TH1D * _h_m_sim_pjx = new TH1D("_h_m_sim_pjx", "",  46, 4.50, 16); _h_m_sim_pjx->SetMarkerColor( kRed );
   _h_m_sim_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); _h_m_sim_pjx->SetLineWidth(0.6);
   _h_m_sim_pjx = _h2_sim_f_fil->ProjectionX("_h_m_sim_pjx ", 0, 150);
   _h_m_sim_pjx->Draw("E"); 

   TH1D * h_m_sim_pjx = new TH1D("h_m_sim_pjx", "",  46, 4.50, 16); h_m_sim_pjx->SetMarkerColor( kBlue );
   h_m_sim_pjx->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_sim_pjx->SetLineWidth(0.6);
   h_m_sim_pjx = h_2d_final->ProjectionX("h_m_sim_pjx", 0, 150);
   h_m_sim_pjx->Draw("Esames");  



   c_f_proj->cd(2)->SetGrid();
   TH1D * h_m_sim_pjy = new TH1D("h_m_sim_pjy", "", 150,0,150); h_m_sim_pjy->SetMarkerColor( kBlue );
   h_m_sim_pjy->GetXaxis()->SetTitle("# of Tracklets"); h_m_sim_pjy->SetLineWidth(0.6);
   h_m_sim_pjy = h_2d_final->ProjectionY(" h_m_sim_pjy", 0, 46);
   h_m_sim_pjy->Draw(); 

   TH1D * _h_m_sim_pjy = new TH1D("_h_m_sim_pjy", "", 150,0,150); _h_m_sim_pjy->SetMarkerColor( kRed );
   _h_m_sim_pjy->GetXaxis()->SetTitle("# of Tracklets"); _h_m_sim_pjy->SetLineWidth(0.6);
   _h_m_sim_pjy = _h2_sim_f_fil->ProjectionY("_h_m_sim_pjy ", 0, 46);
   _h_m_sim_pjy->Draw("sames"); 


//************************************************ Final Fitting Comparison **********************************


  TCanvas *c_long = new TCanvas(TString("c_long"), "c_long",1600 * .7, 1600 * .7);
  c_long->cd();
  c_long->Divide(2,1);

  c_long->cd(1)->SetGrid();
  TH2D *h_2d_final_dy = new TH2D("h_2d_final_dy","", 46, 4.5, 16,150,0,150);
  h_2d_final_dy->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_final_dy->GetYaxis()->SetTitle("Tracklets");
  h_2d_final_dy->SetMarkerStyle(21);
  h_2d_final_dy->SetMarkerSize(.5);
  h_2d_final_dy->SetMarkerColor( kBlue );
  h_2d_final_dy->GetXaxis()->SetTitleOffset(1.4);
  h_2d_final_dy->GetYaxis()->SetTitleOffset(1.4);

  for (int i=1; i< 15; i++){
        for(int j = 1; j < 51; j++){

h_2d_final_dy->GetBin(i,j);
h_2d_final_dy->SetBinContent(i,j, (f_2d_dy->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p0_fnl);
h_2d_final_dy->SetBinError(i,j, sqrt(f_2d_dy->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p0_fnl);}}

 
 for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                h_2d_final_dy->SetBinContent(i+14,j+50,0);
                h_2d_final_dy->SetBinError(i+14,j+50,0);

                h_2d_final_dy->SetBinContent(i+30,j+100,0);
                h_2d_final_dy->SetBinError(i+30,j+100,0);
        }
}


  TH2D *h_2d_final_bb = new TH2D("h_2d_final_bb","", 46, 4.5, 16,150,0,150);
  h_2d_final_bb->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_final_bb->GetYaxis()->SetTitle("Tracklets");
  h_2d_final_bb->SetMarkerStyle(21);
  h_2d_final_bb->SetMarkerSize(.5);
  h_2d_final_bb->SetMarkerColor( kBlue );
  h_2d_final_bb->GetXaxis()->SetTitleOffset(1.4);
  h_2d_final_bb->GetYaxis()->SetTitleOffset(1.4);

  for (int i=1; i< 15; i++){
        for(int j = 1; j < 51; j++){

h_2d_final_bb->GetBin(i,j);
h_2d_final_bb->SetBinContent(i,j, (f_2d_bb->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p1_fnl);
h_2d_final_bb->SetBinError(i,j, sqrt(f_2d_bb->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p1_fnl);}}

 
 for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

    h_2d_final_bb->SetBinContent(i+14,j+50, (f_2d_bb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*0.2566*p1_fnl);
    h_2d_final_bb->SetBinError(i+14,j+50,sqrt(f_2d_bb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*0.2566*p1_fnl);

    h_2d_final_bb->SetBinContent(i+30,j+100, (f_2d_bb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*0.2695*p1_fnl);
    h_2d_final_bb->SetBinError(i+30,j+100,sqrt(f_2d_bb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*0.2695*p1_fnl);
        }
}


  TH2D *h_2d_final_cc = new TH2D("h_2d_final_cc","", 46, 4.5, 16,150,0,150);
  h_2d_final_cc->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_final_cc->GetYaxis()->SetTitle("Tracklets");
  h_2d_final_cc->SetMarkerStyle(21);
  h_2d_final_cc->SetMarkerSize(.5);
  h_2d_final_cc->SetMarkerColor( kBlue );
  h_2d_final_cc->GetXaxis()->SetTitleOffset(1.4);
  h_2d_final_cc->GetYaxis()->SetTitleOffset(1.4);

  for (int i=1; i< 15; i++){
        for(int j = 1; j < 51; j++){

h_2d_final_cc->GetBin(i,j);
h_2d_final_cc->SetBinContent(i,j, (f_2d_cc->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p2_fnl);
h_2d_final_cc->SetBinError(i,j, sqrt(f_2d_cc->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p2_fnl);}}

 
 for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                h_2d_final_cc->SetBinContent(i+14,j+50,0);
                h_2d_final_cc->SetBinError(i+14,j+50,0);

                h_2d_final_cc->SetBinContent(i+30,j+100,0);
                h_2d_final_cc->SetBinError(i+30,j+100,0);
        }
}


  TH2D *h_2d_final_cb = new TH2D("h_2d_final_cb","", 46, 4.5, 16,150,0,150);
  h_2d_final_cb->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_final_cb->GetYaxis()->SetTitle("Tracklets");
  h_2d_final_cb->SetMarkerStyle(21);
  h_2d_final_cb->SetMarkerSize(.5);
  h_2d_final_cb->SetMarkerColor( kBlue );
  h_2d_final_cb->GetXaxis()->SetTitleOffset(1.4);
  h_2d_final_cb->GetYaxis()->SetTitleOffset(1.4);

  for (int i=1; i< 15; i++){
        for(int j = 1; j < 51; j++){

h_2d_final_cb->GetBin(i,j);
h_2d_final_cb->SetBinContent(i,j, (f_comb->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*sqrt(p3_fnl*p7_fnl)*2);
h_2d_final_cb->SetBinError(i,j, sqrt(f_comb->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*sqrt(p3_fnl*p7_fnl)*2);}}

 
 for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){

                h_2d_final_cb->SetBinContent(i+14,j+50, (f_comb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p3_fnl);
                h_2d_final_cb->SetBinError(i+14,j+50,sqrt(f_comb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p3_fnl);

                h_2d_final_cb->SetBinContent(i+30,j+100, (f_comb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p7_fnl);
                h_2d_final_cb->SetBinError(i+30,j+100,sqrt(f_comb->Integral(4+(i-1)*0.25,4 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p7_fnl);
        }
}



  TH2D *h_2d_final_pr = new TH2D("h_2d_final_pr","", 46, 4.5, 16,150,0,150);
  h_2d_final_pr->GetXaxis()->SetTitle("Mass (GeV)");
  h_2d_final_pr->GetYaxis()->SetTitle("Tracklets");
  h_2d_final_pr->SetMarkerStyle(21);
  h_2d_final_pr->SetMarkerSize(.5);
  h_2d_final_pr->SetMarkerColor( kBlue );
  h_2d_final_pr->GetXaxis()->SetTitleOffset(1.4);
  h_2d_final_pr->GetYaxis()->SetTitleOffset(1.4);

  for (int i=1; i< 15; i++){
        for(int j = 1; j < 51; j++){

                h_2d_final_pr->GetBin(i,j);
                h_2d_final_pr->SetBinContent(i,j, (f_2d_pr->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p8_fnl);
                h_2d_final_pr->SetBinError(i,j, sqrt(f_2d_pr->Integral(4.5+(i-1)*0.25,4.5 +(i)*0.25,(j-1)+0.0,j,1E-5)/0.25)*p8_fnl);}}

  for (int i=1; i< 17; i++){
        for(int j = 1; j < 51; j++){
                h_2d_final_pr->SetBinContent(i+14,j+50,0);
                h_2d_final_pr->SetBinError(i+14,j+50,0);

                h_2d_final_pr->SetBinContent(i+30,j+100,0);
                h_2d_final_pr->SetBinError(i+30,j+100,0);
        }
}


   h_m_sim_pjx->Draw();  
   _h_m_sim_pjx->Draw("sames");


   TH1D * h_m_sim_pjx_bb = new TH1D("h_m_sim_pjx_bb", "h_m_sim_pjx_bb",  46, 4.5, 16); h_m_sim_pjx_bb->SetMarkerColor( kRed );
   h_m_sim_pjx_bb->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_sim_pjx_bb->SetLineWidth(0.6);
   h_m_sim_pjx_bb = h_2d_final_bb->ProjectionX("h_m_sim_pjx_bb", 0, 150);h_m_sim_pjx_bb->SetLineColor( 2 );
   h_m_sim_pjx_bb->Draw("sames");

   TH1D * h_m_sim_pjx_dy = new TH1D("h_m_sim_pjx_dy", "h_m_sim_pjx_dy",  46, 4.5, 16); h_m_sim_pjx_dy->SetMarkerColor( kBlue );
   h_m_sim_pjx_dy->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_sim_pjx_dy->SetLineWidth(0.6);
   h_m_sim_pjx_dy = h_2d_final_dy->ProjectionX("h_m_sim_pjx_dy", 0, 150); h_m_sim_pjx_dy->SetLineColor( 4 );
   h_m_sim_pjx_dy->Draw("sames");


   TH1D * h_m_sim_pjx_cc = new TH1D("h_m_sim_pjx_cc", "h_m_sim_pjx_cc",  46, 4.5, 16); h_m_sim_pjx_cc->SetMarkerColor( kGreen );
   h_m_sim_pjx_cc->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_sim_pjx_cc->SetLineWidth(0.6);
   h_m_sim_pjx_cc = h_2d_final_cc->ProjectionX("h_m_sim_pjx_cc", 0, 150); h_m_sim_pjx_cc->SetLineColor( 6 );
   h_m_sim_pjx_cc->Draw("sames");


   TH1D * h_m_sim_pjx_cb = new TH1D("h_m_sim_pjx_cb", "h_m_sim_pjx_cb",  46, 4.5, 16); h_m_sim_pjx_cb->SetMarkerColor( kGreen );
   h_m_sim_pjx_cb->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_sim_pjx_cb->SetLineWidth(0.6);
   h_m_sim_pjx_cb = h_2d_final_cb->ProjectionX("h_m_sim_pjx_cb", 0, 150); h_m_sim_pjx_cb->SetLineColor( 3 );
   h_m_sim_pjx_cb->Draw("sames");

   TH1D * h_m_sim_pjx_pr = new TH1D("h_m_sim_pjx_pr", "h_m_sim_pjx_pr",  46, 4.5, 16); h_m_sim_pjx_pr->SetMarkerColor( kGreen );
   h_m_sim_pjx_pr->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_m_sim_pjx_pr->SetLineWidth(0.6);
   h_m_sim_pjx_pr = h_2d_final_pr->ProjectionX("h_m_sim_pjx_pr", 0, 150); h_m_sim_pjx_pr->SetLineColor( 8 );
   h_m_sim_pjx_pr->Draw("sames");

   c_long->cd(2)->SetGrid();

   h_m_sim_pjy->Draw();  
   _h_m_sim_pjy->Draw("sames");


   TH1D * h_m_sim_pjy_bb = new TH1D("h_m_sim_pjy_bb", "h_m_sim_pjy_bb",  150,0,150); h_m_sim_pjy_bb->SetMarkerColor( kRed );
   h_m_sim_pjy_bb->GetXaxis()->SetTitle("# of Tracklets"); h_m_sim_pjy_bb->SetLineWidth(0.6);
   h_m_sim_pjy_bb = h_2d_final_bb->ProjectionY("h_m_sim_pjy_bb", 0, 46);h_m_sim_pjy_bb->SetLineColor( 2 );
   h_m_sim_pjy_bb->Draw("sames");

   TH1D * h_m_sim_pjy_dy = new TH1D("h_m_sim_pjy_dy", "h_m_sim_pjy_dy",  150,0,150); h_m_sim_pjy_dy->SetMarkerColor( kBlue );
   h_m_sim_pjy_dy->GetXaxis()->SetTitle("# of Tracklets"); h_m_sim_pjy_dy->SetLineWidth(0.6);
   h_m_sim_pjy_dy = h_2d_final_dy->ProjectionY("h_m_sim_pjy_dy", 0, 46); h_m_sim_pjy_dy->SetLineColor( 4 );
   h_m_sim_pjy_dy->Draw("sames");


   TH1D * h_m_sim_pjy_cc = new TH1D("h_m_sim_pjy_cc", "h_m_sim_pjy_cc",  150,0,150); h_m_sim_pjy_cc->SetMarkerColor( kGreen );
   h_m_sim_pjy_cc->GetXaxis()->SetTitle("# of Tracklets"); h_m_sim_pjy_cc->SetLineWidth(0.6);
   h_m_sim_pjy_cc = h_2d_final_cc->ProjectionY("h_m_sim_pjy_cc", 0, 46); h_m_sim_pjy_cc->SetLineColor( 6 );
   h_m_sim_pjy_cc->Draw("sames");


   TH1D * h_m_sim_pjy_cb = new TH1D("h_m_sim_pjy_cb", "h_m_sim_pjy_cb",  150,0,150); h_m_sim_pjy_cb->SetMarkerColor( kGreen );
   h_m_sim_pjy_cb->GetXaxis()->SetTitle("# of Tracklets"); h_m_sim_pjy_cb->SetLineWidth(0.6);
   h_m_sim_pjy_cb = h_2d_final_cb->ProjectionY("h_m_sim_pjy_cb", 0, 46); h_m_sim_pjy_cb->SetLineColor( 3 );
   h_m_sim_pjy_cb->Draw("sames");

   TH1D * h_m_sim_pjy_pr = new TH1D("h_m_sim_pjy_pr", "h_m_sim_pjy_pr",  150,0,150); h_m_sim_pjy_pr->SetMarkerColor( kGreen );
   h_m_sim_pjy_pr->GetXaxis()->SetTitle("# of Tracklets"); h_m_sim_pjy_pr->SetLineWidth(0.6);
   h_m_sim_pjy_pr = h_2d_final_pr->ProjectionY("h_m_sim_pjy_pr", 0, 46); h_m_sim_pjy_pr->SetLineColor( 8 );
   h_m_sim_pjy_pr->Draw("sames");

//---------------------------------PHOO VALUE------------------------------------------------------
   TCanvas *c_phoo_l = new TCanvas(TString("Final phoo Projection"), "Final phoo Projection",1600 * .7, 1600 * .7);
   //gStyle->SetOptStat(0);
   c_phoo_l->Divide(2,2,0,0);

   c_phoo_l->cd(1);//->SetGrid();

   h_m_sim_pjx->Draw("E");  //h_m_sim_pjx->Scale(1/h_m_sim_pjx->Integral());
   _h_m_sim_pjx->Draw("Esames");  //_h_m_sim_pjx->Scale(1/_h_m_sim_pjx->Integral());

   h_m_sim_pjx_bb->Draw("sames");
   h_m_sim_pjx_dy->Draw("sames");
   h_m_sim_pjx_cc->Draw("sames");
   h_m_sim_pjx_cb->Draw("sames");
   h_m_sim_pjx_pr->Draw("sames");

   TLegend *leg_fnlx = new TLegend(0.65, 0.65, .85, .85);
   leg_fnlx->SetHeader("Run 13 Data Fit");
   leg_fnlx->AddEntry(h_m_sim_pjx,"Run13 Data","lp");
   leg_fnlx->AddEntry(_h_m_sim_pjx,"Function Determined","lp");
   leg_fnlx->AddEntry(h_m_sim_pjx_dy,"Drell-Yan Yield","l");
   leg_fnlx->AddEntry(h_m_sim_pjx_bb,"b#bar{b} Yield","l");
   leg_fnlx->AddEntry(h_m_sim_pjx_cc,"c#bar{c} Yield","l");
   leg_fnlx->AddEntry(h_m_sim_pjx_cb,"Combinatorial Yield","l");
   leg_fnlx->AddEntry(h_m_sim_pjx_pr,"#Psi^{/} Yield","l");
   leg_fnlx->Draw();


   c_phoo_l->cd(3);//->SetGrid();
   TH1D * h_phoo_sim_m = new TH1D("h_phoo_sim_m", "",  46, 4.5, 16); h_phoo_sim_m->SetMarkerColor( kBlue );
   h_phoo_sim_m->GetXaxis()->SetTitle("Inv. Mass (GeV)"); h_phoo_sim_m->SetLineWidth(0.6);

   h_phoo_sim_m->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); h_phoo_sim_m->SetFillColor(15);h_phoo_sim_m->GetYaxis()->CenterTitle();

	for (int i=1; i< 47; i++){
                h_phoo_sim_m->GetBin(i);
		if((h_m_sim_pjx->GetBinContent(i)- _h_m_sim_pjx->GetBinContent(i))==0 || h_m_sim_pjx->GetBinError(i) ==0){
		h_phoo_sim_m->SetBinContent(i,0);
                h_phoo_sim_m->SetBinError(i,0);
	}

	if((h_m_sim_pjx->GetBinContent(i)- _h_m_sim_pjx->GetBinContent(i))!=0 && h_m_sim_pjx->GetBinError(i) !=0){
	h_phoo_sim_m->SetBinContent(i,(h_m_sim_pjx->GetBinContent(i)- _h_m_sim_pjx->GetBinContent(i))/h_m_sim_pjx->GetBinError(i));
	h_phoo_sim_m->SetBinError(i,0);}
	}
   h_phoo_sim_m->Draw();

   c_phoo_l->cd(2);//->SetGrid();

   h_m_sim_pjy->Draw();
   _h_m_sim_pjy->Draw("sames"); 
   h_m_sim_pjy_bb->Draw("sames");
   h_m_sim_pjy_dy->Draw("sames");
   h_m_sim_pjy_cc->Draw("sames");
   h_m_sim_pjy_cb->Draw("sames");
   h_m_sim_pjy_pr->Draw("sames");

   TLegend *leg_fnly = new TLegend(0.65, 0.65, .85, .85);
   leg_fnly->SetHeader("Run 13 Data Fit");
   leg_fnly->AddEntry(h_m_sim_pjy,"Run13 Data","lp");
   leg_fnly->AddEntry(_h_m_sim_pjy,"Function Determined","lp");
   leg_fnly->AddEntry(h_m_sim_pjy_dy,"Drell-Yan Yield","l");
   leg_fnly->AddEntry(h_m_sim_pjy_bb,"b#bar{b} Yield","l");
   leg_fnly->AddEntry(h_m_sim_pjy_cc,"c#bar{c} Yield","l");
   leg_fnly->AddEntry(h_m_sim_pjy_cb,"Combinatorial Yield","l");
   leg_fnly->AddEntry(h_m_sim_pjy_pr,"#Psi^{/} Yield","l");
   leg_fnly->Draw();

   c_phoo_l->cd(4);//->SetGrid();
   TH1D * h_phoo_sim_t = new TH1D("h_phoo_sim_t", "", 150,0,150); h_phoo_sim_t->SetMarkerColor( kRed );
   h_phoo_sim_t->GetXaxis()->SetTitle("# of Tracklets"); h_phoo_sim_t->SetLineWidth(0.6);
   h_phoo_sim_t->GetYaxis()->SetTitle("#frac{Data - Function}{#Delta Data}"); h_phoo_sim_t->SetFillColor(15);h_phoo_sim_t->GetYaxis()->CenterTitle();

	for (int i=1; i< 151; i++){
                h_phoo_sim_t->GetBin(i);
		if((h_m_sim_pjy->GetBinContent(i)- _h_m_sim_pjy->GetBinContent(i))==0 || h_m_sim_pjy->GetBinError(i) ==0){
		h_phoo_sim_t->SetBinContent(i,0);
                h_phoo_sim_t->SetBinError(i,0);}

	if((h_m_sim_pjy->GetBinContent(i)- _h_m_sim_pjy->GetBinContent(i))!=0 && h_m_sim_pjy->GetBinError(i) !=0){
	h_phoo_sim_t->SetBinContent(i,(h_m_sim_pjy->GetBinContent(i)- _h_m_sim_pjy->GetBinContent(i))/h_m_sim_pjy->GetBinError(i));
	h_phoo_sim_t->SetBinError(i,0);}
	}
	h_phoo_sim_t->Draw();


//--------------------------- Finding the Fractions ----------------------------------------------------

int lo_binx = h_2d_final_dy->GetXaxis()->FindBin(5.00001);
int hi_binx = h_2d_final_dy->GetXaxis()->FindBin(7.99999);
int lo_biny = h_2d_final_dy->GetYaxis()->FindBin(0.00001);//5.00001,10.00001,15.00001
int hi_biny = h_2d_final_dy->GetYaxis()->FindBin(4.9999);//9.9999,14.9999,49.9999

cout <<" lo x = "<<lo_binx<<" hi x = "<<hi_binx<<" lo y = "<<lo_biny<<" hi y = "<<hi_biny<<endl;


double all_prossesors_os = h_2d_final_dy-> Integral(lo_binx,hi_binx,lo_biny,hi_biny) + h_2d_final_bb-> Integral(lo_binx,hi_binx,lo_biny,hi_biny) + h_2d_final_cc-> Integral(lo_binx,hi_binx,lo_biny,hi_biny) + h_2d_final_cb-> Integral(lo_binx,hi_binx,lo_biny,hi_biny) + h_2d_final_pr-> Integral(lo_binx,hi_binx,lo_biny,hi_biny);
double hf_fraction = (h_2d_final_bb-> Integral(lo_binx,hi_binx,lo_biny,hi_biny) + h_2d_final_cc-> Integral(lo_binx,hi_binx,lo_biny,hi_biny))/all_prossesors_os;
double cb_fraction = (h_2d_final_cb-> Integral(lo_binx,hi_binx,lo_biny,hi_biny))/all_prossesors_os;
double dy_fraction = (h_2d_final_dy-> Integral(lo_binx,hi_binx,lo_biny,hi_biny))/all_prossesors_os;
double prime_fraction = (h_2d_final_pr-> Integral(lo_binx,hi_binx,lo_biny,hi_biny))/all_prossesors_os;

double bb_fraction_num = (h_2d_final_bb-> Integral(lo_binx,hi_binx,lo_biny,hi_biny));
double cc_fraction_num = (h_2d_final_cc-> Integral(lo_binx,hi_binx,lo_biny,hi_biny));
double cb_fraction_num = (h_2d_final_cb-> Integral(lo_binx,hi_binx,lo_biny,hi_biny));
double dy_fraction_num = (h_2d_final_dy-> Integral(lo_binx,hi_binx,lo_biny,hi_biny));
double prime_fraction_num = (h_2d_final_pr-> Integral(lo_binx,hi_binx,lo_biny,hi_biny));
//----------------------- same sign fractions
cout<<""<<endl;
cout<<""<<endl;

int lo_binx_pp = h_2d_final_dy->GetXaxis()->FindBin(9.00001);
int hi_binx_pp = h_2d_final_dy->GetXaxis()->FindBin(11.99999);
int lo_biny_pp = h_2d_final_dy->GetYaxis()->FindBin(50.00001);//50.001,55.00001,60.0001,65.0001
int hi_biny_pp = h_2d_final_dy->GetYaxis()->FindBin(54.9999);//54.999,59.9999,64.999,99.999

double all_prossesors_ss_pp =  h_2d_final_bb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp) +  h_2d_final_cb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp);
double hf_fraction_pp = (h_2d_final_bb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp))/all_prossesors_ss_pp;
double cb_fraction_pp = (h_2d_final_cb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp))/all_prossesors_ss_pp;

double hf_fraction_pp_num = (h_2d_final_bb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp));
double cb_fraction_pp_num = (h_2d_final_cb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp));

//----------------------------------------------------------------
cout<<""<<endl;
cout<<""<<endl;

int lo_binx_mm = h_2d_final_dy->GetXaxis()->FindBin(13.00001);
int hi_binx_mm = h_2d_final_dy->GetXaxis()->FindBin(15.99999);
int lo_biny_mm = h_2d_final_dy->GetYaxis()->FindBin(100.00001);//100.001,105.00001,110.0001,115.001
int hi_biny_mm = h_2d_final_dy->GetYaxis()->FindBin(104.9999);//104.999,109.9999,114.999,149.000

double all_prossesors_ss_mm =  h_2d_final_bb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm) +  h_2d_final_cb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm);
double hf_fraction_mm = (h_2d_final_bb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm))/all_prossesors_ss_mm;
double cb_fraction_mm = (h_2d_final_cb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm))/all_prossesors_ss_mm;

double hf_fraction_mm_num = (h_2d_final_bb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm));
double cb_fraction_mm_num = (h_2d_final_cb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm));
cout << " BB Events in pp is   = " << hf_fraction_pp_num <<endl;
cout << " CB Events in pp is   = " << cb_fraction_pp_num <<endl;
cout << " total Events in pp is = " << all_prossesors_ss_pp <<endl;
cout << " BB Events in mm is   = " << hf_fraction_mm_num <<endl;
cout << " CB Events in mm is   = " << cb_fraction_mm_num <<endl;
cout << " total Events mm is = " << all_prossesors_ss_mm <<endl;
cout << " ------------------------------------- " <<endl;
cout << " ------------------------------------- " <<endl;
cout << " HF Fraction in pp is   = " << hf_fraction_pp <<endl;
cout << " CB Fraction in pp is   = " << cb_fraction_pp <<endl;
cout << " total of frac in pp is = " << hf_fraction_pp + cb_fraction_pp <<endl;
cout << " HF Fraction in mm is  = " << hf_fraction_mm <<endl;
cout << " CB Fraction in mm is   = " << cb_fraction_mm <<endl;
cout << " total of fra in mm is = " << hf_fraction_mm + cb_fraction_mm <<endl;
cout << " ------------------------------------- " <<endl;
cout << " ------------------------------------- " <<endl;
cout << " DY Fraction is = " << dy_fraction <<endl;
cout << " HF Fraction is = " << hf_fraction <<endl;
cout << " CB Fraction is = " << cb_fraction <<endl;
cout << " PR Fraction is = " << prime_fraction <<endl;
cout << " total fracs is = " << hf_fraction + cb_fraction + dy_fraction + prime_fraction<<endl;


cout << "all same siggn hf frac = "<< (h_2d_final_bb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm)+ h_2d_final_bb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp))/(all_prossesors_ss_pp+all_prossesors_ss_mm)<<endl;

cout << "all same siggn cb frac = "<< (h_2d_final_cb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm)+ h_2d_final_cb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp))/(all_prossesors_ss_pp+all_prossesors_ss_mm)<<endl;
cout << " ------------------------------------- " <<endl;
cout << " ------------------------------------- " <<endl;
cout << " DY Events = " << dy_fraction_num <<endl;
cout << " BB Events = " << bb_fraction_num <<endl;
cout << " CC Events = " << cc_fraction_num <<endl;
cout << " CB Events = " << cb_fraction_num <<endl;
cout << " PR Events = " << prime_fraction_num <<endl;
cout << " total Events = " << all_prossesors_os<<endl;

cout << "all same siggn BB = "<< (h_2d_final_bb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm)+ h_2d_final_bb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp))<<endl;

cout << "all same siggn CB = "<< (h_2d_final_cb-> Integral(lo_binx_mm,hi_binx_mm,lo_biny_mm,hi_biny_mm)+ h_2d_final_cb-> Integral(lo_binx_pp,hi_binx_pp,lo_biny_pp,hi_biny_pp))<<endl;

}
