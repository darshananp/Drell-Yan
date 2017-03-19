#include "twobins_N.h"

using namespace std;

void twobins_N(){

gStyle->SetOptFit(1111);

  TCanvas *c_bin1 = new TCanvas("c_bin1","Graph2DErrors example",0,0,600,600);
  c_bin1->cd();
  int np = 8;

  // Double_t rx_bin1[8] = {0.1836,0.2627,0.3813,0.3923,0.6490,0.7812,0.8272,0.9875};
  // Double_t ry_bin1[8] = {0.1856,0.7373,0.2126,0.6077,0.0658,0.2188,0.0038,0.0125};
   Double_t rx_bin1[8] = {0.2173,
0.2627,
0.4430,
0.3923,
0.7326,
0.7812,
0.8656,
0.9875
};
   Double_t ry_bin1[8] = {0.2402,
0.7373,
0.2702,
0.6077,
0.0808,
0.2188,
0.0043,
0.0125
};
   Double_t rz_bin1[8] = {-0.0535,0.9684,0.1618,-1.2413,-0.0521,0.9837,0.1323,1.3066};
   Double_t ez_bin1[8] = {0.7023,0.7592,0.5513,0.5671,0.6862,0.6498,0.9667,0.8413};
   Double_t ex[8] = {0,0,0,0,0,0,0,0};
   Double_t ey[8] = {0,0,0,0,0,0,0,0};

/*   Double_t ex[8] = {0.0284,
0.0344,
0.0580,
0.0513,
0.0958,
0.1022,
0.1132,
0.1292
};
   Double_t ey[8] = {0.0399,
0.1224,
0.0449,
0.1009,
0.0134,
0.0363,
0.0007,
0.0021
};*/

   TGraph2DErrors *dte_bin1 = new TGraph2DErrors(np, rx_bin1, ry_bin1, rz_bin1, ex, ey, ez_bin1);
   dte_bin1->SetTitle("TGraph2D with error bars: option \"ERR\"");
   dte_bin1->SetFillColor(29);
   dte_bin1->SetMarkerSize(0.8);
   dte_bin1->SetMarkerStyle(20);
   dte_bin1->GetXaxis()->SetTitle("HF fraction");
   dte_bin1->GetXaxis()->SetTitleOffset(1.8);
   dte_bin1->GetYaxis()->SetTitle("CB fraction");
   dte_bin1->GetYaxis()->SetTitleOffset(1.8);
   dte_bin1->SetMinimum(-2);
   dte_bin1->SetMaximum(2);
   dte_bin1->SetMarkerColor(kRed);
   dte_bin1->SetLineColor(kBlue-3);
   dte_bin1->SetLineWidth(2);

   dte_bin1->Draw("err P0");

  TF2 *f_2d_bin1 = new TF2("f_2d_bin1",_func_2D,0,1,0,1,3);
  f_2d_bin1->SetParameters(0.01,0.01,0.01);
  //dte_bin1->Fit( f_2d_bin1 , "MRS++");
  TFitResultPtr r_dte_bin1 = dte_bin1->Fit(f_2d_bin1 , "MRS");
  f_2d_bin1->Draw("surf1 same");
  TMatrixDSym cov_dte_bin1 = r_dte_bin1->GetCovarianceMatrix();  
  r_dte_bin1->Print("V");  
  cov_dte_bin1.Print();

TCanvas* c_bin1_p= new TCanvas("ProjCan1","The Projections2",1000,400);
 c_bin1_p->Divide(2,1);
 c_bin1_p->cd(1);
 dte_bin1->Project("x")->Draw();
 c_bin1_p->cd(2);
 dte_bin1->Project("y")->Draw();

//*****************************************************************
  TCanvas *c_bin2 = new TCanvas("c_bin2","Graph2DErrors example",0,0,600,600);
  c_bin2->cd();

   //Double_t rx_bin2[8] = {0.2394,0.4675,0.4937,0.6139,0.7486,0.8979,0.8623,0.9949};
   //Double_t ry_bin2[8] = {0.1357,0.5325,0.1544,0.3861,0.0425,0.1021,0.0022,0.0051};
   Double_t rx_bin2[8] = {0.2266,
0.4674,
0.4671,
0.6139,
0.7105,
0.8979,
0.8182,
0.9949
};
   Double_t ry_bin2[8] = {0.1357,
0.5326,
0.1544,
0.3861,
0.0425,
0.1021,
0.0022,
0.0051
};
   Double_t rz_bin2[8] = {0.2439,0.1221,-0.4497,0.7750,0.1082,-0.9367,-2.1512,-0.6006};
   Double_t ez_bin2[8] = {0.5116,0.9526,0.4756,0.8166,0.4995,0.9076,0.6631,1.0944};

/*   Double_t ex2[8] = {0.0296,
0.0612,
0.0611,
0.0803,
0.0930,
0.1175,
0.1070,
0.1302
};
   Double_t ey2[8] = {0.0225,
0.0884,
0.0256,
0.0641,
0.0071,
0.0170,
0.0004,
0.0008
};*/

   TGraph2DErrors *dte_bin2 = new TGraph2DErrors(np, rx_bin2, ry_bin2, rz_bin2, ex, ey, ez_bin2);
   dte_bin2->SetTitle("TGraph2D with error bars: option \"ERR\"");
   dte_bin2->SetFillColor(29);
   dte_bin2->SetMarkerSize(0.8);
   dte_bin2->SetMarkerStyle(20);
   dte_bin2->GetXaxis()->SetTitle("HF fraction");
   dte_bin2->GetXaxis()->SetTitleOffset(1.8);
   dte_bin2->GetYaxis()->SetTitle("CB fraction");
   dte_bin2->GetYaxis()->SetTitleOffset(1.8);
   dte_bin2->SetMinimum(-2);
   dte_bin2->SetMaximum(2);
   dte_bin2->SetMarkerColor(kRed);
   dte_bin2->SetLineColor(kBlue-3);
   dte_bin2->SetLineWidth(2);

   dte_bin2->Draw("err P0");

  TF2 *f_2d_bin2 = new TF2("f_2d_bin2",_func_2D,0,1,0,1,3);
  f_2d_bin2->SetParameters(0.01,0.01,0.01);
  //dte_bin2->Fit( f_2d_bin2 , "MRS++");
  TFitResultPtr r_dte_bin2 = dte_bin2->Fit(f_2d_bin2 , "MRS");
  f_2d_bin2->Draw("surf1 same");
  TMatrixDSym cov_dte_bin2 = r_dte_bin2->GetCovarianceMatrix();  
  r_dte_bin2->Print("V");  
  cov_dte_bin2.Print();

 TCanvas* c_bin2_p= new TCanvas("ProjCan2","The Projections2",1000,400);
 c_bin2_p->Divide(2,1);
 c_bin2_p->cd(1);
 dte_bin2->Project("x")->Draw();
 c_bin2_p->cd(2);
 dte_bin2->Project("y")->Draw();

cout <<"DY Asymmetry bin 1= "<< f_2d_bin1->GetParameter(0) <<"  +-  "<<f_2d_bin1->GetParError(0) <<endl;
cout <<"DY Asymmetry bin 2= "<< f_2d_bin2->GetParameter(0) <<"  +-  "<<f_2d_bin2->GetParError(0) <<endl;
cout <<"Final DY Asymmetry= "<< (f_2d_bin1->GetParameter(0)*f_2d_bin1->GetParError(0)+f_2d_bin2->GetParameter(0)*f_2d_bin2->GetParError(0))/ (f_2d_bin1->GetParError(0)+ f_2d_bin2->GetParError(0))<<"  +-  "<<sqrt(f_2d_bin1->GetParError(0)*f_2d_bin1->GetParError(0)+ f_2d_bin2->GetParError(0)*f_2d_bin2->GetParError(0))/2 <<endl;


//*****************************************************************************

   Double_t r_ms[2],r_dy[2],e_ms[2],e_dy[2];
   r_ms[0]= 4.725;
   r_ms[1]= 5.935;

   e_ms[0]= 0;
   e_ms[1]= 0;

   r_dy[0]= f_2d_bin1->GetParameter(1);
   r_dy[1]= f_2d_bin2->GetParameter(1);

   e_dy[0]= f_2d_bin1->GetParError(1);
   e_dy[1]= f_2d_bin2->GetParError(1);

   TCanvas *c = new TCanvas("c","Graph2DErrors example",0,0,600,600);
  // c->Divide(3,2);
   //c->cd()->SetGrid();
      TH1 *frame_1 = new TH1F("frame_1","",40,4,8);
      frame_1->SetMinimum(-2.0);
      frame_1->SetMaximum(2.0);
      frame_1->SetDirectory(0);
      frame_1->SetStats(0);
      frame_1->SetTitle("Heavy-flavor A_{LL} North Arm");
      frame_1->GetXaxis()->SetTitle("Mass (GeV)");
      frame_1->GetXaxis()->SetTickLength(0.1);
      frame_1->GetXaxis()->SetLabelSize(0.03);
      frame_1->GetYaxis()->SetTitle("A_{LL}");
      frame_1->GetYaxis()->SetLabelSize(0.03);
      frame_1->Draw();

   TGraphErrors *dte_1 = new TGraphErrors(2, r_ms, r_dy, e_ms, e_dy);
   dte_1->SetMarkerSize(0.8);
   dte_1->SetMarkerStyle(20);
   dte_1->SetMarkerColor(kRed);
   dte_1->SetLineColor(4);
   dte_1->SetMarkerSize(1);      dte_1->SetLineWidth(2);
   dte_1->Draw("P");
   dte_1->Fit("pol0");
}
