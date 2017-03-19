  #include <cassert>
  #include <string>
  #include "Math/Minimizer.h"
  #include "Math/Factory.h"
  #include "Math/Functor.h"
  #include <stdio.h>   
  #include <math.h>
  #include <cmath>
  #include <stdlib.h>
  #include <errno.h>
  #include <algorithm>
  #include <sys/stat.h>
  #include <TMinuit.h>
  #include <TMath.h>
  #include <TProfile.h>
  #include <TFile.h>
  #include <TTree.h>
  #include <TGraph.h>
  #include <TGraphPainter.h>
  #include <TGraphErrors.h>
  #include <TH2.h>
  #include <TF2.h>
  #include "TROOT.h"
  #include "TCut.h"
  #include "TTree.h"
  #include "TStyle.h"
  #include "TSystem.h"
  #include "TCanvas.h"
  #include "Riostream.h"
  #include "TApplication.h"
  #include "TException.h"
  #include "TClass.h"
  #include "TClassTable.h"
  #include "TEnv.h"
  #include "TBrowser.h"
  #include "TString.h"
  #include "TObjString.h"
  #include "TError.h"
  #include "TF12.h"
  #include "TVirtualFitter.h"
  #include "TFitResultPtr.h"
  #include "TMatrixDSym.h"
  #include "TMatrixD.h"
  #include "TMatrix.h"
  #include "TFitResult.h"
  #include "TCanvas.h"
  #include "TLegend.h"
  #include "TLatex.h"
  using namespace std;

  TTree *t1 = NULL;  TTree *t2 = NULL;  TTree *t3 = NULL; TTree *t4 = NULL; 

  TCut cut_dy_pm;  TCut cut_bb_pm;  TCut cut_cc_pm; TCut cut_rd_pm;  TCut cut_rd_pp;  TCut cut_rd_mm;  
  TCut cut_bb_pp;  TCut cut_bb_mm;  TCut cut_pr_pm;

  Double_t p1_dy ,e1_dy,p2_dy ,e2_dy ,p3_dy ,e3_dy ,p4_dy ,e4_dy ,p5_dy ,e5_dy;
  Double_t p6_dy ,e6_dy,p7_dy ,e7_dy ,p8_dy ,e8_dy ,p9_dy ,e9_dy ,p10_dy ,e10_dy ,p11_dy ,e11_dy;
  Double_t p1_cc ,e1_cc,p2_cc ,e2_cc ,p3_cc ,e3_cc ,p4_cc ,e4_cc ,p5_cc ,e5_cc;
  Double_t p6_cc ,e6_cc,p7_cc ,e7_cc ,p8_cc ,e8_cc ,p9_cc ,e9_cc ,p10_cc ,e10_cc ,p11_cc ,e11_cc;
  Double_t p1_bb ,e1_bb,p2_bb ,e2_bb ,p3_bb ,e3_bb, p4_bb ,e4_bb ,p5_bb ,e5_bb;
  Double_t p6_bb ,e6_bb,p7_bb ,e7_bb ,p8_bb ,e8_bb ,p9_bb ,e9_bb ,p10_bb ,e10_bb ,p11_bb ,e11_bb;

  Double_t f_p0_dy,f_p1_dy,f_p2_dy,f_p3_dy,f_p4_dy,f_p5_dy,f_p6_dy,f_p7_dy,f_p8_dy,f_p9_dy,f_p10_dy;
  Double_t f_p0_bb,f_p1_bb,f_p2_bb,f_p3_bb,f_p4_bb,f_p5_bb,f_p6_bb,f_p7_bb,f_p8_bb,f_p9_bb,f_p10_bb;
  Double_t f_p0_cc,f_p1_cc,f_p2_cc,f_p3_cc,f_p4_cc,f_p5_cc,f_p6_cc,f_p7_cc,f_p8_cc,f_p9_cc,f_p10_cc;
  Double_t f_p0_pr,f_p1_pr,f_p2_pr,f_p3_pr,f_p4_pr,f_p5_pr,f_p6_pr,f_p7_pr,f_p8_pr,f_p9_pr,f_p10_pr;

  Double_t p0_fnl,p1_fnl,p2_fnl,p3_fnl,p4_fnl,p5_fnl,p6_fnl,p7_fnl,p8_fnl,p9_fnl,p10_fnl,p11_fnl;
  Double_t p0_fnl_err,p1_fnl_err,p2_fnl_err,p3_fnl_err,p4_fnl_err,p5_fnl_err,p6_fnl_err,
           p7_fnl_err,p8_fnl_err,p9_fnl_err,p10_fnl_err,p11_fnl_err;
  void plots_N();


//-----------------------2 - D function for mass ---------------------------------------
  Double_t _func_2D_dy(Double_t *x, Double_t *par) ;
  Double_t _func_2D_bb(Double_t *x, Double_t *par) ;
  Double_t _func_2D_cc(Double_t *x, Double_t *par) ;
  Double_t _func_mixed(Double_t *x, Double_t *par) ;
  Double_t _func_prime1(Double_t *x, Double_t *par) ;  Double_t _func_prime2(Double_t *x, Double_t *par) ;
//-----------------------end fit function  ----------------------------------------------
  Double_t _func_final(Double_t *x, Double_t *par) ;

//-----------------------1 - D function for mass ---------------------------------------
  Double_t _func_mass_dy(Double_t *x, Double_t *par) { return (par[0] + par[1]*x[0] + par[2]*x[0]*x[0])*exp(par[3]+par[4]*x[0]) ; }
  Double_t _func_mass_bb(Double_t *x, Double_t *par) { return exp(par[0] + par[1]*x[0]+ par[2]*x[0]*x[0]); }

//-----------------------1 - D function for tracklet -----------------------------------
  Double_t _func_trk_dy(Double_t *x, Double_t *par) {
  return (par[0]*TMath::Power((par[1]/par[2]),(x[0]/par[2]))*(TMath::Exp(-(par[1]/par[2])))/TMath::Gamma((x[0]/par[2])+1))*(par[3] + par[4]*x[0] + par[5]*x[0]*x[0]); }

//para0=height,para1=mean,para3=
  Double_t _func_trk(Double_t *x, Double_t *par) {
  return (par[0]*TMath::Power((par[1]/par[2]),(x[0]/par[2]))*(TMath::Exp(-(par[1]/par[2])))/TMath::Gamma((x[0]/par[2])+1))*(par[3] + par[4]*x[0]); }

  Double_t _func_2D_dy(Double_t *x, Double_t *par) { return exp(par[0]*x[0] + par[1]*x[0]*x[0])*(TMath::Power(((par[2]+ par[7]*x[0])/(par[3]+ par[8]*x[0])),(x[1]/(par[3]+ par[8]*x[0])))*(TMath::Exp(-((par[2]+ par[7]*x[0])/(par[3]+ par[8]*x[0]))))/TMath::Gamma((x[1]/(par[3]+ par[8]*x[0]))+1))*(par[4] + par[5]*x[1] + par[6]*x[1]*x[1]) ; }

  Double_t _func_2D_bb(Double_t *x, Double_t *par) { return exp(par[0]*x[0] + par[1]*x[0]*x[0])*(TMath::Power((par[2]/par[3]),(x[1]/par[3]))*(TMath::Exp(-(par[2]/par[3])))/TMath::Gamma((x[1]/par[3])+1))*(par[4] + par[5]*x[1] ) ; }

  Double_t _func_2D_cc(Double_t *x, Double_t *par) { return exp(par[0]*x[0] + par[1]*x[0]*x[0])*(TMath::Power((par[2]/par[3]),(x[1]/par[3]))*(TMath::Exp(-(par[2]/par[3])))/TMath::Gamma((x[1]/par[3])+1))*(par[4] + par[5]*x[1] ) ;}

  Double_t _func_2D_jp(Double_t *x, Double_t *par) { return exp(par[0]*x[0] + par[1]*x[0]*x[0])*(TMath::Power((par[2]/par[3]),(x[1]/par[3]))*(TMath::Exp(-(par[2]/par[3])))/TMath::Gamma((x[1]/par[3])+1))*(par[4] + par[5]*x[1] ) ;}


  Double_t _func_final(Double_t *x, Double_t *par) { 
 if (x[0]>4 && x[0]<8 && x[1]<50){
 return par[0]*exp(f_p0_dy*x[0] + f_p1_dy*x[0]*x[0])*(TMath::Power(((f_p2_dy+ f_p7_dy*x[0])/(f_p3_dy+ f_p8_dy*x[0])),(x[1]/(f_p3_dy+ f_p8_dy*x[0])))*(TMath::Exp(-((f_p2_dy+ f_p7_dy*x[0])/(f_p3_dy+ f_p8_dy*x[0]))))/TMath::Gamma((x[1]/(f_p3_dy+ f_p8_dy*x[0]))+1))*(f_p4_dy + f_p5_dy*x[1]+ f_p6_dy*x[1]*x[1]) + 

 par[1]*exp(f_p0_bb*x[0] + f_p1_bb*x[0]*x[0])*(TMath::Power((f_p2_bb/f_p3_bb),(x[1]/f_p3_bb))*(TMath::Exp(-(f_p2_bb/f_p3_bb)))/TMath::Gamma((x[1]/f_p3_bb)+1))*(f_p4_bb + f_p5_bb*x[1])+

 par[2]*exp(f_p0_cc*x[0] + f_p1_cc*x[0]*x[0])*(TMath::Power((f_p2_cc/f_p3_cc),(x[1]/f_p3_cc))*(TMath::Exp(-(f_p2_cc/f_p3_cc)))/TMath::Gamma((x[1]/f_p3_cc)+1))*(f_p4_cc + f_p5_cc*x[1]) + 
 
par[8]*exp(f_p0_pr*x[0] + f_p1_pr*x[0]*x[0])*(TMath::Power((f_p2_pr/f_p3_pr),(x[1]/f_p3_pr))*(TMath::Exp(-(f_p2_pr/f_p3_pr)))/TMath::Gamma((x[1]/f_p3_pr)+1))*(f_p4_pr + f_p5_pr*x[1]) +

 sqrt(par[3]*par[7])*2*exp(par[4]*x[0] + par[5]*x[0]*x[0])*TMath::Poisson(x[1],par[6]);}

  if (x[0]> 8 && x[0]<12 && x[1]>50 && x[1]<100){

  //return par[1]*0.2566*exp(f_p0_bb*(x[0]-4) + f_p1_bb*(x[0]-4)*(x[0]-4))*(TMath::Power((f_p2_bb/f_p3_bb),((x[1]-50)/f_p3_bb))*(TMath::Exp(-(f_p2_bb/f_p3_bb)))/TMath::Gamma(((x[1]-50)/f_p3_bb)+1))*(f_p4_bb + f_p5_bb*(x[1]-50))+ 

//north
  return par[1]*0.2571*exp(f_p0_bb*(x[0]-4) + f_p1_bb*(x[0]-4)*(x[0]-4))*(TMath::Power((f_p2_bb/f_p3_bb),((x[1]-50)/f_p3_bb))*(TMath::Exp(-(f_p2_bb/f_p3_bb)))/TMath::Gamma(((x[1]-50)/f_p3_bb)+1))*(f_p4_bb + f_p5_bb*(x[1]-50))+ 

 par[3]*exp(par[4]*(x[0]-4) + par[5]*(x[0]-4)*(x[0]-4))*TMath::Poisson((x[1]-50),par[6]); }

  if (x[0]> 12 && x[0]<16 && x[1]>100 && x[1]<150){

  //return par[1]*0.2695*exp(f_p0_bb*(x[0]-8) + f_p1_bb*(x[0]-8)*(x[0]-8))*(TMath::Power((f_p2_bb/f_p3_bb),((x[1]-100)/f_p3_bb))*(TMath::Exp(-(f_p2_bb/f_p3_bb)))/TMath::Gamma(((x[1]-100)/f_p3_bb)+1))*(f_p4_bb + f_p5_bb*(x[1]-100)) +
//north
  return par[1]*0.2528*exp(f_p0_bb*(x[0]-8) + f_p1_bb*(x[0]-8)*(x[0]-8))*(TMath::Power((f_p2_bb/f_p3_bb),((x[1]-100)/f_p3_bb))*(TMath::Exp(-(f_p2_bb/f_p3_bb)))/TMath::Gamma(((x[1]-100)/f_p3_bb)+1))*(f_p4_bb + f_p5_bb*(x[1]-100)) +

 par[7]*exp(par[4]*(x[0]-8) + par[5]*(x[0]-8)*(x[0]-8))*TMath::Poisson((x[1]-100),par[6]);}

  else {return 0;}
}


  Double_t _func_mixed(Double_t *x, Double_t *par) { 
 return par[0]*exp(f_p0_dy*x[0] + f_p1_dy*x[0]*x[0])*(TMath::Power(((f_p2_dy+ f_p7_dy*x[0])/(f_p3_dy+ f_p8_dy*x[0])),(x[1]/(f_p3_dy+ f_p8_dy*x[0])))*(TMath::Exp(-((f_p2_dy+ f_p7_dy*x[0])/(f_p3_dy+ f_p8_dy*x[0]))))/TMath::Gamma((x[1]/(f_p3_dy+ f_p8_dy*x[0]))+1))*(f_p4_dy + f_p5_dy*x[1]+ f_p6_dy*x[1]*x[1]) + 

 par[1]*exp(f_p0_bb*x[0] + f_p1_bb*x[0]*x[0])*(TMath::Power((f_p2_bb/f_p3_bb),(x[1]/f_p3_bb))*(TMath::Exp(-(f_p2_bb/f_p3_bb)))/TMath::Gamma((x[1]/f_p3_bb)+1))*(f_p4_bb + f_p5_bb*x[1])+

 par[2]*exp(f_p0_cc*x[0] + f_p1_cc*x[0]*x[0])*(TMath::Power((f_p2_cc/f_p3_cc),(x[1]/f_p3_cc))*(TMath::Exp(-(f_p2_cc/f_p3_cc)))/TMath::Gamma((x[1]/f_p3_cc)+1))*(f_p4_cc + f_p5_cc*x[1]);

}


  Double_t _func_prime1(Double_t *x, Double_t *par) { 
return (par[0]*0.1889*TMath::Power((7.368/1.634),(x[0]/1.634))*(TMath::Exp(-(7.368/1.634)))/TMath::Gamma((x[0]/1.634)+1))*(45.7 + 3.239*x[0] + 0.07618*x[0]*x[0]);
}

  Double_t _func_prime2(Double_t *x, Double_t *par) { 
return (par[0]*5.238*TMath::Power((9.438/1.272),(x[0]/1.272))*(TMath::Exp(-(9.438/1.272)))/TMath::Gamma((x[0]/1.272)+1))*(106.7 + -16.52*x[0] + 0.7532*x[0]*x[0]);
}
