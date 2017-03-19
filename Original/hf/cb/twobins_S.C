#include "twobins_S.h"

using namespace std;

void twobins_S(){

gStyle->SetOptFit(1111);

  TCanvas *c_bin1 = new TCanvas("c_bin1","Graph2DErrors example",0,0,600,600);
  c_bin1->cd();
  int np = 8;

//   Double_t rx_bin1[8] = {0.2043,0.3073,0.4458,0.4926,0.7326,0.8666,0.8814,0.9943};
 //  Double_t ry_bin1[8] = {0.1866,0.6927,0.1877,0.5074,0.0464,0.1334,0.0021,0.0057};

   Double_t rx_bin1[8] = {0.2041,
0.3073,
0.4475,
0.4926,
0.7194,
0.8666,
0.8321,
0.9943
};
   Double_t ry_bin1[8] = {0.2168,
0.6927,
0.2170,
0.5074,
0.0522,
0.1334,
0.0022,
0.0057
};
   Double_t rz_bin1[8] = {-0.4910,0.8368,0.4613,-0.0590,-0.9298,-0.2163,0.3581,-0.1594};
   Double_t ez_bin1[8] = {0.6006,0.5816,0.4211,0.4363,0.4866,0.4458,0.6309,0.5802};
   Double_t ex[8] = {0,0,0,0,0,0,0,0};
   Double_t ey[8] = {0,0,0,0,0,0,0,0};
   /*Double_t ex[8] = {0.0164,
0.0247,
0.0360,
0.0396,
0.0579,
0.0697,
0.0670,
0.0800
};
   Double_t ey[8] = {0.0755,
0.2412,
0.0756,
0.1767,
0.0182,
0.0465,
0.0008,
0.0020
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

 //  Double_t rx_bin2[8] = {0.2232,0.4367,0.4929,0.6292,0.7627,0.9403,0.8778,0.9968};
   //Double_t ry_bin2[8] = {0.1369,0.5633,0.1389,0.3708,0.0323,0.0597,0.0014,0.0032};
   Double_t rx_bin2[8] = {0.2039,
0.4367,
0.4528,
0.6292,
0.7034,
0.9190,
0.8134,
0.9968
};
   Double_t ry_bin2[8] = {0.1369,
0.5633,
0.1389,
0.3708,
0.0323,
0.0810,
0.0014,
0.0032
};


   Double_t rz_bin2[8] = {0.1802,-0.3225,-0.0601,-0.6057,0.2512,0.6165,0.3878,-0.6838};
   Double_t ez_bin2[8] = {0.3903,0.8071,0.3039,0.5333,0.3555,0.6033,0.4766,0.6762};
/*
   Double_t ex2[8] = {0.0164,
0.0351,
0.0364,
0.0506,
0.0566,
0.0740,
0.0655,
0.0802
};
   Double_t ey2[8] = {0.0477,
0.1961,
0.0484,
0.1291,
0.0112,
0.0282,
0.0005,
0.0011
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
   r_ms[0]= 4.728;
   r_ms[1]= 5.944;

   e_ms[0]= 0;
   e_ms[1]= 0;

   r_dy[0]= f_2d_bin1->GetParameter(2);
   r_dy[1]= f_2d_bin2->GetParameter(2);

   e_dy[0]= f_2d_bin1->GetParError(2);
   e_dy[1]= f_2d_bin2->GetParError(2);

   TCanvas *c = new TCanvas("c","Graph2DErrors example",0,0,600,600);
  // c->Divide(3,2);
   //c->cd()->SetGrid();
      TH1 *frame_1 = new TH1F("frame_1","",40,4,8);
      frame_1->SetMinimum(-2.0);
      frame_1->SetMaximum(2.0);
      frame_1->SetDirectory(0);
      frame_1->SetStats(0);
      frame_1->SetTitle("Combinatorial A_{LL} South Arm");
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
