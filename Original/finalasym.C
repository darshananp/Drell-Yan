#include <TGraphErrors.h>
#include <TH2.h>
#include <TF2.h>
#include <cassert>
#include <fstream.h>
#include <string>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <stdio.h>   
#include <math.h>
#include "TLegend.h"
#include "TLatex.h"
using namespace std;

void finalasym(){

gStyle->SetOptFit(1111);

   Double_t r_ms[2] = {1,2};
   Double_t e_ms[2] = {0,0}; 
//   Double_t r_dy[2] = {0.169925,-0.0224567}; //older version
//   Double_t e_dy[2] = {0.461418,0.600426};//older version 
   Double_t r_dy[2] = {0.185607,-0.000216383}; //with out err
   Double_t e_dy[2] = {0.449279,0.644912}; //with out err

//   Double_t r_dy[2] = {0.190919,0.195237}; //with  err
//   Double_t e_dy[2] = {0.497779,0.740524}; //with  err

   TCanvas *c = new TCanvas("c","Graph2DErrors example",0,0,1200,1200);
   //c->cd()->SetGrid();
      TH1 *frame_1 = new TH1F("frame_1","",2,0.5,2.5);
      frame_1->SetMinimum(-1);
      frame_1->SetMaximum(1);
      frame_1->SetDirectory(0);

      frame_1->SetStats(0);
      frame_1->SetTitle("Final Drell-Yan A_{LL} measurement for mass case");
      frame_1->GetXaxis()->SetTitle("Arm (1=South,2=North)");
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
