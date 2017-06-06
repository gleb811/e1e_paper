TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;


TH1D *m_pip_p_bin_corr,*m_pip_pim_bin_corr,*m_pim_p_bin_corr;
TH1D *theta_p_bin_corr,*theta_pim_bin_corr,*theta_pip_bin_corr;
TH1D *alpha_p_bin_corr,*alpha_pim_bin_corr,*alpha_pip_bin_corr;

int i;

TH1D *h_w[12];

TCanvas *c = new TCanvas("c","c",900,300);


Float_t Q2_bin,W_bin[30];

Double_t w_err;

TLegend *leg;

TSpline5 *spline;
ostringstream qqq;


 
Double_t f_fit(Double_t *x, Double_t *par) {

return spline->Eval(x[0]);

}


TH1D *h_bin_corr_theta(TH1D *h) {

TH1D *h_out;

Double_t *x = new double[(h->GetNbinsX())-1];
Double_t *y = new double[(h->GetNbinsX())-1];

for (Int_t aa = 1; aa <=((h->GetNbinsX())-1); aa++) {
x[aa-1] = (h->GetBinCenter(aa)+ h->GetBinCenter(aa+1))/2.;
y[aa-1] = (h->GetBinContent(aa) + h->GetBinContent(aa+1))/2.;
};


TGraph *gr = new TGraphErrors(((h->GetNbinsX())-1),x,y);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(1);
h_out->SetMarkerColor(1);


qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");


f = new TF1("func",f_fit,h->GetBinLowEdge(1),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1),0);

tmp = new TH1D("tmp","tmp",100,h->GetBinLowEdge(0),h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1));
for (Int_t subbin=1; subbin<101;subbin++) {

tmp->SetBinContent(subbin,f->Eval(tmp->GetBinCenter(subbin)));
};
tmp->SetLineWidth(2);
tmp->Draw("same c");

for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

/*
for (Int_t subbin=0; subbin<100;subbin++) {

position = left + width*subbin/99.;
avrg = avrg + spline->Eval(position);

};
avrg = avrg/100.;
*/

avrg = (f->Integral(left,left+width))/width;



if ((h->GetBinContent(bin)) == 0.) avrg = 0;

//cout << m_pip_p->GetBinContent(bin) << "  " << avrg << endl;

if ((h->GetBinContent(bin)) != 0.) {
h_out->SetBinContent(bin,(h->GetBinContent(bin))*(f->Eval(h->GetBinCenter(bin)))/avrg);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));
} else {
h_out->SetBinContent(bin,0.);
h_out->SetBinError(bin,0.000001);
};

};




return h_out;


}; 





TH1D *h_alpha_sym(TH1D *h) {

TH1D *h_out;




Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName();
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(1);
h_out->SetMarkerColor(1);


qqq.str("");
qqq << h->GetName();
h_out->SetName(qqq.str().c_str());
qqq.str("");




for (Int_t bin=1; bin<=Int_t((h_out->GetNbinsX())/2.);bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

h_out->SetBinContent(bin,((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinContent(((h->GetNbinsX())-bin+1),((h->GetBinContent(bin))+(h->GetBinContent((h->GetNbinsX())-bin+1)))/2.);


h_out->SetBinError(bin,sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);

h_out->SetBinError(((h->GetNbinsX())-bin+1),sqrt((h->GetBinError(bin))*(h->GetBinError(bin))+(h->GetBinError((h->GetNbinsX())-bin+1))*(h->GetBinError((h->GetNbinsX())-bin+1)))/2.);


};

return h_out;


}; 



TH1D *h_bin_corr(TH1D *h) {

TH1D *h_out;


TGraph *gr = new TGraph(h);
spline = new TSpline5("spline",gr);
spline->SetLineWidth(2);

Float_t left;
Float_t width;
Float_t position;
Float_t avrg;

qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out = (TH1D*)h->Clone(qqq.str().c_str());
qqq.str("");

h_out->SetLineColor(1);
h_out->SetMarkerColor(1);


qqq.str("");
qqq << h->GetName() << "bin_corr";
h_out->SetName(qqq.str().c_str());
qqq.str("");




for (Int_t bin=1; bin<=(h_out->GetNbinsX());bin++) {

left = h->GetBinLowEdge(bin);
width = h->GetBinWidth(bin);
avrg = 0;

for (Int_t subbin=0; subbin<100;subbin++) {

position = left + width*subbin/99.;
avrg = avrg + spline->Eval(position);

};
avrg = avrg/100.;

if ((h->GetBinContent(bin)) == 0.) avrg = 0;

//cout << m_pip_p->GetBinContent(bin) << "  " << avrg << endl;

if ((h->GetBinContent(bin)) != 0.) {
h_out->SetBinContent(bin,(h->GetBinContent(bin))*(spline->Eval(h->GetBinCenter(bin)))/avrg);
h_out->SetBinError(bin,(h->GetBinError(bin))/(h->GetBinContent(bin))*(h_out->GetBinContent(bin)));
} else {
h_out->SetBinContent(bin,0.);
h_out->SetBinError(bin,0.000001);
};

};

return h_out;


}; 


void draw_1d_hist (Int_t canvas, TH1D *h, string title, string name, string ytitle, string xtitle, Int_t color, string draw_options, string distr_flag,Int_t i) {

c->cd(canvas);
c->cd(canvas)->SetBottomMargin(0.2);
c->cd(canvas)->SetTopMargin(0.1);
c->cd(canvas)->SetLeftMargin(0.23);
c->cd(canvas)->SetRightMargin(0.01);
gPad->SetFillStyle(0);

h->SetMarkerStyle(20);
h->SetMarkerColor(color);
h->SetLineColor(color);
h->SetOption("pX0");
h->SetTitle(title.c_str());
h->SetTitleSize(0.1);


h->SetName(name.c_str());


 h->GetYaxis()->SetTitle(ytitle.c_str());
 h ->GetXaxis()->SetTitle(xtitle.c_str());
 h ->GetXaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleSize(0.08);
 h->GetYaxis()->SetTitleOffset(1.3);
 h->GetXaxis()->SetLabelSize(0.08);
 h->GetXaxis()->SetNdivisions(6);
 h->GetYaxis()->SetLabelSize(0.07);
 h->GetYaxis()->SetNdivisions(5);


h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum())/2., "Y");
if ((name == "model_thetaP_")&&(W_bin[i]>1.65))h->SetAxisRange(0, h-> GetMaximum()+(h-> GetMaximum()), "Y");


};



Int_t get_max_w (Float_t Q2_bin) {
Int_t get_max_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.5))get_max_w = 21;
if ((Q2_bin> 0.5)&&(Q2_bin< 0.6))get_max_w = 20;
if ((Q2_bin> 0.6)&&(Q2_bin< 0.65))get_max_w = 19;
if ((Q2_bin> 0.65)&&(Q2_bin< 0.7))get_max_w = 18;
if ((Q2_bin> 0.7)&&(Q2_bin< 0.75))get_max_w = 17;
if ((Q2_bin> 0.75)&&(Q2_bin< 0.8))get_max_w = 16;
if ((Q2_bin> 0.8)&&(Q2_bin< 0.85))get_max_w = 14;
if ((Q2_bin> 0.85)&&(Q2_bin< 0.9))get_max_w = 13;
if ((Q2_bin> 0.9)&&(Q2_bin< 0.95))get_max_w = 12;
if ((Q2_bin> 0.95)&&(Q2_bin< 1.0))get_max_w = 10;
if ((Q2_bin> 1.0)&&(Q2_bin< 1.05))get_max_w = 9;
if ((Q2_bin> 1.05)&&(Q2_bin< 1.1))get_max_w = 7;
if ((Q2_bin> 1.1)&&(Q2_bin< 1.15))get_max_w = 6;
if ((Q2_bin> 1.15)&&(Q2_bin< 1.2))get_max_w = 4;
if ((Q2_bin> 1.2)&&(Q2_bin< 1.25))get_max_w = 2;
if ((Q2_bin> 1.25)&&(Q2_bin< 1.3))get_max_w = 1;
return get_max_w;
};

Int_t get_min_w (Float_t Q2_bin) {
Int_t get_min_w = 0;
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 12;
return get_min_w;
};

void read_data_eff(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p_bin_corr);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim_bin_corr);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p_bin_corr);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p_bin_corr);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim_bin_corr);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip_bin_corr);
qqq.str("");






qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p_bin_corr);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim_bin_corr);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip_bin_corr);
qqq.str("");



};


void read_data_rec(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_pim_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_p);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);
qqq.str("");



qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);
qqq.str("");






qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);
qqq.str("");





qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);
qqq.str("");




qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_bin_corr";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");



};




void draw_1d_canvas( Int_t i, Int_t qq2 ) {

c = new TCanvas("c","c",700,700);

c->Divide(3,3);

qqq << "Q^{2} = " << Q2_bin << " GeV^{2}" <<"; W = " << W_bin[i] <<" GeV" ;

Float_t max;

max = m_pip_p->GetMaximum();

if ((m_pip_pim->GetMaximum()) > max) max = m_pip_pim->GetMaximum();

if ((m_pim_p->GetMaximum()) > max) max = m_pim_p->GetMaximum();


for (Int_t kk=1; kk<4; kk++) { 
c->cd(kk)->SetBottomMargin(0.2);
c->cd(kk)->SetTopMargin(0.15);
c->cd(kk)->SetLeftMargin(0.23);
c->cd(kk)->SetRightMargin(0.01);
c->cd(kk)->SetFillStyle(4000);
};

for (Int_t kk=4; kk<10; kk++) { 
c->cd(kk)->SetBottomMargin(0.2);
c->cd(kk)->SetTopMargin(0.005);
c->cd(kk)->SetLeftMargin(0.23);
c->cd(kk)->SetRightMargin(0.01);
c->cd(kk)->SetFillStyle(4000);
};
gPad->SetFillStyle(4000);

c->cd(1);
m_pip_p->SetLineColor(1);
m_pip_p->SetMarkerColor(1);
m_pip_p->SetMinimum(0);
m_pip_p->SetMaximum(1.2*max);
m_pip_p->Draw("e1");
m_pip_p->Draw("e0 same");


c->cd(2);
m_pip_pim->SetLineColor(1);
m_pip_pim->SetMarkerColor(1);
m_pip_pim->SetMinimum(0);
m_pip_pim->SetMaximum(1.2*max);
m_pip_pim->Draw("e1");
m_pip_pim->Draw("e0 same");


c->cd(3);
m_pim_p->SetLineColor(1);
m_pim_p->SetMarkerColor(1);
m_pim_p->SetMinimum(0);
m_pim_p->SetMaximum(1.2*max);
m_pim_p->Draw("e1");
m_pim_p->Draw("e0 same");


max = theta_p->GetMaximum();

if ((theta_pip->GetMaximum()) > max) max = theta_pip->GetMaximum();

if ((theta_pim->GetMaximum()) > max) max = theta_pim->GetMaximum();

c->cd(4);
theta_p->SetLineColor(1);
theta_p->SetMarkerColor(1);
theta_p->SetMinimum(0);
theta_p->SetMaximum(1.2*max);
theta_p->GetXaxis()->SetTitle("#theta_{p'} (deg)");
theta_p->Draw("e1");


c->cd(5);
theta_pim->SetLineColor(1);
theta_pim->SetMarkerColor(1);
theta_pim->SetMinimum(0);
theta_pim->SetMaximum(1.2*max);
theta_pim->GetXaxis()->SetTitle("#theta_{#pi^{-}} (deg)");
theta_pim->Draw("e1");

c->cd(6);
theta_pip->SetLineColor(1);
theta_pip->SetMarkerColor(1);
theta_pip->SetMinimum(0);
theta_pip->SetMaximum(1.2*max);
theta_pip->GetXaxis()->SetTitle("#theta_{#pi^{+}} (deg)");
theta_pip->Draw("e1");

max = alpha_p->GetMaximum();

if ((alpha_pip->GetMaximum()) > max) max = alpha_pip->GetMaximum();

if ((alpha_pim->GetMaximum()) > max) max = alpha_pim->GetMaximum();


c->cd(7);
alpha_p->SetLineColor(1);
alpha_p->SetMarkerColor(1);
alpha_p->SetMinimum(0);
alpha_p->SetMaximum(1.2*max);
alpha_p->Draw("e1");


c->cd(8);
alpha_pim->SetLineColor(1);
alpha_pim->SetMarkerColor(1);
alpha_pim->SetMinimum(0);
alpha_pim->SetMaximum(1.2*max);
alpha_pim->GetXaxis()->SetTitle("#alpha_{#pi^{-}} (deg)");
alpha_pim->Draw("e1");


c->cd(9);
alpha_pip->SetLineColor(1);
alpha_pip->SetMarkerColor(1);
alpha_pip->SetMinimum(0);
alpha_pip->SetMaximum(1.2*max);
alpha_pip->Draw("e1");


c->cd();
TLatex latex;
latex.SetTextAlign(12);  //centered
latex.DrawLatex(.2,.9745,qqq.str().c_str());
qqq.str("");
c->Update();

qqq << "Q_" << Q2_bin << "_W_" << W_bin[i] <<".pdf";
c->Print(qqq.str().c_str());
qqq.str("");
//c->Delete();

};


void e1e_e1c() {
#include <TH2.h>
#include <TH1.h>
#include <TH3.h>
#include <THnSparse.h>
gStyle->SetTitleSize(0.07,"t");
gStyle->SetTitleY(1.01);
gStyle->SetOptStat(0);
gStyle->SetErrorX(0);
gErrorIgnoreLevel = kError;
gStyle->SetStatY(0.88); 





std::ostringstream file3;   

file3 << "e1c04505.dat";
// cout << file3.str() << "\n";




 float sigma4[8],error4[8],theta4[8],errorx4[8],error4_sum[8];
 
 

std::ifstream sect4(file3.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <8; i++)
{
  sect4 >> theta4[i] >> sigma4[i] >> error4[i];
  errorx4[i] = 0.0125;
  
  if (i == 0) error4_sum[i] = sqrt(error4[i]*error4[i]+0.5*0.5);
  if (i == 1) error4_sum[i] = sqrt(error4[i]*error4[i]+1.5*1.5);
  if (i == 2) error4_sum[i] = sqrt(error4[i]*error4[i]+2.*2.);
  if (i == 3) error4_sum[i] = sqrt(error4[i]*error4[i]+2.*2.); 
  if (i == 4) error4_sum[i] = sqrt(error4[i]*error4[i]+2.2*2.2);
  if (i == 5) error4_sum[i] = sqrt(error4[i]*error4[i]+2.3*2.3);
  if (i == 6) error4_sum[i] = sqrt(error4[i]*error4[i]+3.7*3.7);
  if (i == 7) error4_sum[i] = sqrt(error4[i]*error4[i]+5.*5.);    
  
//  sigma[i-20]=sigma[i-20]*2.0;

// cout << theta4[i] << "  " << sigma4[i] << "  " << error4[i] << "\n";

}
 sect4.close();   
 
   TGraphErrors *Tg0475 = new TGraphErrors(8,theta4,sigma4,errorx4,error4);
   
   Tg0475->SetMarkerStyle(33);
   Tg0475->SetMarkerColor(1);
   Tg0475->SetMarkerSize(1.5);
   Tg0475->SetLineColor(1); 
   Tg0475->SetMaximum(33.);
   Tg0475->GetXaxis()->SetLimits(1.3,1.9);
   Tg0475->SetMinimum(0.);
   Tg0475->GetXaxis()->SetTitle("W (GeV)");
   Tg0475->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg0475->GetXaxis()->SetTitleSize(0.09);
   Tg0475->GetYaxis()->SetTitleSize(0.09); 
   Tg0475->GetYaxis()->SetTitleOffset(0.85);
   Tg0475->GetXaxis()->SetTitleOffset(0.85);    
   Tg0475->GetXaxis()->SetLabelSize(0.07); 
   Tg0475->GetYaxis()->SetLabelSize(0.07);     
   


//std::ostringstream file3;   
file3.str("");

file3 << "e1c05055.dat";
// cout << file3.str() << "\n";




 float sigma5[8],error5[8],theta5[8],errorx5[8],error5_sum[8];
 
 

std::ifstream sect5(file3.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <8; i++)
{
  sect5 >> theta5[i] >> sigma5[i] >> error5[i];
  errorx5[i]= 0.0125;  
  if (i == 0) error5_sum[i] = sqrt(error5[i]*error5[i]+1.*1.);
  if (i == 1) error5_sum[i] = sqrt(error5[i]*error5[i]+1.5*1.5);
  if (i == 2) error5_sum[i] = sqrt(error5[i]*error5[i]+3.8*3.8);
  if (i == 3) error5_sum[i] = sqrt(error5[i]*error5[i]+3.7*3.7); 
  if (i == 4) error5_sum[i] = sqrt(error5[i]*error5[i]+1.8*1.8);
  if (i == 5) error5_sum[i] = sqrt(error5[i]*error5[i]+2.*2.);
  if (i == 6) error5_sum[i] = sqrt(error5[i]*error5[i]+2.2*2.2);
  if (i == 7) error5_sum[i] = sqrt(error5[i]*error5[i]+2.3*2.3);    
//  sigma[i-20]=sigma[i-20]*2.0;

// cout << theta5[i] << "  " << sigma5[i] << "  " << error5[i] << "\n";

}
 sect5.close();   
 
   TGraphErrors *Tg0525 = new TGraphErrors(8,theta5,sigma5,errorx5,error5);
   
   Tg0525->SetMarkerStyle(33);
   Tg0525->SetMarkerColor(1);
   Tg0525->SetMarkerSize(1.5);   
   Tg0525->SetLineColor(1); 
   Tg0525->SetMaximum(33.);
   Tg0525->GetXaxis()->SetLimits(1.3,1.9);
   Tg0525->SetMinimum(0.);
      Tg0525->GetXaxis()->SetTitle("W (GeV)");
   Tg0525->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg0525->GetXaxis()->SetTitleSize(0.09);
   Tg0525->GetYaxis()->SetTitleSize(0.09); 
   Tg0525->GetYaxis()->SetTitleOffset(0.85);
   Tg0525->GetXaxis()->SetTitleOffset(0.85); 
   Tg0525->GetXaxis()->SetLabelSize(0.07); 
   Tg0525->GetYaxis()->SetLabelSize(0.07);    
   
   
   
   
   std::ostringstream file6;   

file6 << "e1c05506.dat";
// cout << file6.str() << "\n";




 float sigma6[6],error6[6],theta6[6],errorx6[6],error6_sum[6];
 
 

std::ifstream sect6(file6.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <6; i++)
{
  sect6 >> theta6[i] >> sigma6[i] >> error6[i];
  errorx6[i] = 0.0125;  
  
  if (i == 0) error6_sum[i] = sqrt(error6[i]*error6[i]+0.5*0.5);
  if (i == 1) error6_sum[i] = sqrt(error6[i]*error6[i]+1.2*1.2);
  if (i == 2) error6_sum[i] = sqrt(error6[i]*error6[i]+1.7*1.7);
  if (i == 3) error6_sum[i] = sqrt(error6[i]*error6[i]+2.2*2.2); 
  if (i == 4) error6_sum[i] = sqrt(error6[i]*error6[i]+1.9*1.9);
  if (i == 5) error6_sum[i] = sqrt(error6[i]*error6[i]+1.8*1.8);

//  sigma[i-20]=sigma[i-20]*2.0;

// cout << theta6[i] << "  " << sigma6[i] << "  " << error6[i] << "\n";

}
 sect6.close();   
 
   TGraphErrors *Tg0575 = new TGraphErrors(6,theta6,sigma6,errorx6,error6);
   
   Tg0575->SetMarkerStyle(33);
   Tg0575->SetMarkerColor(1);
   Tg0575->SetMarkerSize(1.5);   
   Tg0575->SetLineColor(1); 
   Tg0575->SetMaximum(33.);
   Tg0575->GetXaxis()->SetLimits(1.3,1.9);
   Tg0575->SetMinimum(0.);
      Tg0575->GetXaxis()->SetTitle("W (GeV)");
   Tg0575->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg0575->GetXaxis()->SetTitleSize(0.09);
   Tg0575->GetXaxis()->SetTitleOffset(0.85);
   Tg0575->GetYaxis()->SetTitleSize(0.09); 
   Tg0575->GetYaxis()->SetTitleOffset(0.85); 
   Tg0575->GetXaxis()->SetLabelSize(0.07); 
   Tg0575->GetYaxis()->SetLabelSize(0.07);     

   
    
  



c->Divide(3,1);


//Define input files

TFile *file_cr_sec_pim = new TFile("out_cr_sec_all_top_final_bin_corr.root","READ");

TFile *file_max_eff = new TFile("out_cr_sec_all_top_final_bin_corr_eff.root","READ");

 for (Int_t qq2=1; qq2<4;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;
   double ax[21];
   double ay[21];
   double aexl[21];
   double aexh[21];
   double aeyl[21];
   double aeyh[21];
for (Int_t i=0; i<21;i++) { 
ax[i] = 1.3125+0.025*i;
ay[i]=0;
aexl[i]=0;
aexh[i]=0;
aeyl[i]=0;
aeyh[i]=0;

};  
   TGraphAsymmErrors* gae;
qqq.str("");
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
//h_w[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),(get_max_w(Q2_bin)-get_min_w(Q2_bin)),1.3+0.025*get_min_w(Q2_bin),1.3+0.025*get_max_w(Q2_bin));

h_w[qq2] = new TH1D(qqq.str().c_str(),qqq.str().c_str(),21,1.3,1.825);

qqq.str("");
//for (Int_t i=13; i<14;i++) {
for (Int_t i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {  


 W_bin[i] = 1.3125+0.025*i; 

/* 
Double_t nu_g = (W_bin[i]*W_bin[i] + Q2_bin - 0.938272*0.938272)/2./0.938272;

Double_t theta_el = acos(1.- Q2_bin/2.039/(2.039 - nu_g)/2.);

Double_t eps_t = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2_bin)*tan(theta_el/2.)*tan(theta_el/2.));

Double_t eps_l = Q2_bin*eps_t/nu_g/nu_g; 




Double_t theta_el_1 = acos(1.- Q2_bin/1.515/(1.515 - nu_g)/2.);

Double_t eps_t_1 = 1./(1.+ 2.*(1. + nu_g*nu_g/Q2_bin)*tan(theta_el_1/2.)*tan(theta_el_1/2.));

Double_t eps_l_1 = Q2_bin*eps_t_1/nu_g/nu_g;
*/
//cout << eps_l/eps_l_1 << endl;


read_data_rec(file_cr_sec_pim,i);
read_data_eff(file_max_eff,i);
m_pip_pim->Scale(m_pip_pim->GetBinWidth(1));
m_pip_p->Scale(m_pip_p->GetBinWidth(1));
m_pim_p->Scale(m_pim_p->GetBinWidth(1));
alpha_pim->Scale((alpha_pim->GetBinWidth(1))*3.1415/180.);
alpha_p->Scale((alpha_p->GetBinWidth(1))*3.1415/180.);
alpha_pip->Scale((alpha_pip->GetBinWidth(1))*3.1415/180.);

m_pip_pim_bin_corr->Scale(m_pip_pim_bin_corr->GetBinWidth(1));
m_pip_p_bin_corr->Scale(m_pip_p_bin_corr->GetBinWidth(1));
m_pim_p_bin_corr->Scale(m_pim_p_bin_corr->GetBinWidth(1));
alpha_pim_bin_corr->Scale((alpha_pim_bin_corr->GetBinWidth(1))*3.1415/180.);
alpha_p_bin_corr->Scale((alpha_p_bin_corr->GetBinWidth(1))*3.1415/180.);
alpha_pip_bin_corr->Scale((alpha_pip_bin_corr->GetBinWidth(1))*3.1415/180.);

TH1D *h_cos_th;
//FILLING OUT -d(cos) HISTOGRAM
Int_t n_theta_bins;
if ((i==0)||(i==1)) n_theta_bins = 6;
if ((i==2)||(i==3)) n_theta_bins = 8;
if (i>=4) n_theta_bins = 10; 
h_cos_th = new TH1D("h_cos_th","h_cos_th",n_theta_bins,0.,180.);  
Double_t temp;   
for (Int_t j=1; j<=n_theta_bins; j++) {
temp = cos((h_cos_th->GetBinLowEdge(j))*3.1415/180.)-cos(3.1415/180.*(h_cos_th->GetBinLowEdge(j)+h_cos_th->GetBinWidth(j)));
temp = 1./temp;
h_cos_th->SetBinContent(j,temp);
h_cos_th->SetBinError(j,0.);
};
theta_p->Divide(h_cos_th);
theta_pip->Divide(h_cos_th);
theta_pim->Divide(h_cos_th);

theta_p_bin_corr->Divide(h_cos_th);
theta_pip_bin_corr->Divide(h_cos_th);
theta_pim_bin_corr->Divide(h_cos_th);

Double_t w_int;
//w_int = (m_pip_pim->Integral())+(m_pip_p->Integral())+(m_pim_p->Integral());
w_int = (alpha_pim->Integral())+(alpha_p->Integral())+(alpha_pip->Integral());
//w_int = w_int + (theta_pim->Integral())+(theta_p->Integral())+(theta_pip->Integral());

w_int = w_int/3.;
//cout << w_int << "  " << theta_p->Integral() << endl;

Double_t w_int_eff;
//w_int = (m_pip_pim->Integral())+(m_pip_p->Integral())+(m_pim_p->Integral());
w_int_eff = (alpha_pim_bin_corr->Integral())+(alpha_p_bin_corr->Integral())+(alpha_pip_bin_corr->Integral());
//w_int = w_int + (theta_pim->Integral())+(theta_p->Integral())+(theta_pip->Integral());

w_int_eff = w_int_eff/3.;

h_w[qq2]->SetBinContent(i+1,w_int);
//cout << m_pip_pim->Integral() << endl;
Double_t err_total;
//m_pip_pim->IntegralAndError(1,m_pip_pim->GetNbinsX(),w_err);
//err_total = w_err*w_err;
//err_total = w_err;
//m_pip_p->IntegralAndError(1,m_pip_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;
//m_pim_p->IntegralAndError(1,m_pim_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;

alpha_pim->IntegralAndError(1,alpha_pim->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
err_total = w_err;
alpha_p->IntegralAndError(1,alpha_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
err_total = err_total + w_err;
alpha_pip->IntegralAndError(1,alpha_pip->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
err_total = err_total + w_err;

//theta_pim->IntegralAndError(1,theta_pim->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;
//theta_p->IntegralAndError(1,theta_p->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;
//theta_pip->IntegralAndError(1,theta_pip->GetNbinsX(),w_err);
//err_total = err_total + w_err*w_err;
//err_total = err_total + w_err;

//err_total = sqrt(err_total)/9.;
err_total = err_total/3.;


h_w[qq2]->SetBinError(i+1,err_total);
//draw_1d_canvas(i,qq2);
ax[i] = W_bin[i];
ay[i] = w_int;
aexl[i] = 0.0125;
aexh[i] = 0.0125;
aeyl[i] = w_int*0.05;
aeyh[i] = w_int*0.05;

Double_t sys_err;
/*
sys_err = ((w_int - m_pip_pim->Integral())*(w_int - m_pip_pim->Integral()));
sys_err = sys_err + ((w_int - m_pip_p->Integral())*(w_int - m_pip_p->Integral()));
sys_err = sys_err + ((w_int - m_pim_p->Integral())*(w_int - m_pim_p->Integral()));
*/
sys_err = ((w_int - alpha_p->Integral())*(w_int - alpha_p->Integral()));
sys_err = sys_err + ((w_int - alpha_pip->Integral())*(w_int - alpha_pip->Integral()));
sys_err = sys_err + ((w_int - alpha_pim->Integral())*(w_int - alpha_pim->Integral()));

/*
sys_err = sys_err + ((w_int - theta_p->Integral())*(w_int - theta_p->Integral()));
sys_err = sys_err + ((w_int - theta_pip->Integral())*(w_int - theta_pip->Integral()));
sys_err = sys_err + ((w_int - theta_pim->Integral())*(w_int - theta_pim->Integral()));
*/
sys_err = sys_err + (w_int*0.03)*(w_int*0.03);

aeyl[i] = aeyl[i]*aeyl[i]+sys_err;
aeyh[i] = aeyh[i]*aeyh[i]+sys_err;

aeyl[i] = aeyl[i]+err_total*err_total;
aeyh[i] = aeyh[i]+err_total*err_total;

aeyl[i] = sqrt(aeyl[i]);
aeyh[i] = sqrt(aeyh[i]);

};
c->cd();
//c->SetFillColor(kGray);
c->cd(qq2); 
//c->cd(qq2)->SetFillColor(kGray);
c->cd(qq2)->SetBottomMargin(0.16);
c->cd(qq2)->SetTopMargin(0.08);
c->cd(qq2)->SetLeftMargin(0.23);
c->cd(qq2)->SetRightMargin(0.01);
c->cd(qq2)->SetFillStyle(4000);
h_w[qq2] ->GetXaxis()->SetTitleSize(0.08);
h_w[qq2] ->GetXaxis()->SetRangeUser(1.2875,1.85);
 h_w[qq2]->GetYaxis()->SetTitleSize(0.08);
 h_w[qq2]->GetYaxis()->SetTitleOffset(1.3);
 h_w[qq2]->GetXaxis()->SetLabelSize(0.08);
// h_w[qq2]->GetXaxis()->SetNdivisions(6);
 h_w[qq2]->GetYaxis()->SetLabelSize(0.07);
// h_w[qq2]->GetYaxis()->SetNdivisions(5);
h_w[qq2]->SetMinimum(0.);
h_w[qq2]->SetMaximum(40.);
h_w[qq2]->GetXaxis()->SetTitle("W (GeV)");
h_w[qq2]->GetYaxis()->SetTitle("#sigma (#mubn)");
h_w[qq2]->SetMarkerStyle(20);
h_w[qq2]->SetMarkerSize(0.5);
h_w[qq2]->SetMarkerColor(1);
h_w[qq2]->Draw("e1P");
gae = new TGraphAsymmErrors(21, ax, ay, aexl, aexh, aeyl, aeyh);

TGraphAsymmErrors* gae1 = new TGraphAsymmErrors(8, theta4, sigma4, errorx4, errorx4, error4_sum, error4_sum);

TGraphAsymmErrors* gae2 = new TGraphAsymmErrors(8, theta5, sigma5, errorx5, errorx5, error5_sum, error5_sum);

TGraphAsymmErrors* gae3 = new TGraphAsymmErrors(6, theta6, sigma6, errorx6, errorx6, error6_sum, error6_sum);

qqq.str("");
qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
gae->SetTitle(qqq.str().c_str());
gae->GetXaxis()->SetTitleSize(0.08);
 gae->GetYaxis()->SetTitleSize(0.08);
 gae->GetYaxis()->SetTitleOffset(0.8);
 gae->GetXaxis()->SetLabelSize(0.08);
gae->GetXaxis()->SetTitle("W (GeV)");
gae->GetYaxis()->SetTitle("#sigma (#mubn)");   
   gae->SetFillColor(1);
   gae->SetFillStyle(3001);
   gae->SetMinimum(0.);
   gae->SetMaximum(40.);
   gae->GetXaxis()->SetRangeUser(1.2875,1.85);
   gae->Draw("a2");
   gae1->SetFillColor(kBlack);
   gae1->SetFillStyle(3002);
//   gae1->SetFillColorAlpha(kRed, 0.35);
   gae2->SetFillColor(kBlack);
   gae2->SetFillStyle(3002); 
   gae3->SetFillColor(kBlack);
   gae3->SetFillStyle(3002);      
if (qq2 == 1)   gae1->Draw("2");
if (qq2 == 2)   gae2->Draw("2");
if (qq2 == 3)   gae3->Draw("2");
   
h_w[qq2]->Draw("e1P same");

if (qq2 == 1) Tg0475->Draw("P same");
if (qq2 == 2) Tg0525->Draw("P same");
if (qq2 == 3) Tg0575->Draw("P same");
TLegend *leg_w_int = new TLegend(0.65,0.2,0.95,0.35); 
leg_w_int->SetFillStyle(0);
leg_w_int->AddEntry(Tg0475,"e1c data","P");
leg_w_int->AddEntry(h_w[qq2],"e1e data","P"); 
leg_w_int->Draw();
};

/*
leg_w_int = new TLegend(0.17,0.75,0.95,0.95); 
leg_w_int->SetNColumns(3);
leg_w_int->SetFillStyle(0);

TCanvas *c2 = new TCanvas("c2","c2",800.,600.);
c2->cd();
c2->SetFillColor(kGray);
c2->SetTopMargin(0.05);
c2->SetRightMargin(0.05);
c2->SetBottomMargin(0.17);
c2->SetLeftMargin(0.17); 
 qqq.str("");
 qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
 leg_w_int->AddEntry(h_w[0],qqq.str().c_str(),"p");
 qqq.str("");
 
h_w[0]->SetMaximum(40.);
h_w[0]->SetMarkerStyle(20);
h_w[0]->SetMarkerSize(1);
h_w[0]->SetTitle("");
h_w[0]->GetYaxis()->SetTitleOffset(0.9);
h_w[0]->Draw("e1P");
 for (Int_t qq2=1; qq2<12;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;
h_w[qq2]->SetTitle(""); 
if (qq2 < 6) {
h_w[qq2]->SetMarkerStyle(20);
h_w[qq2]->SetMarkerSize(1); 
h_w[qq2]->SetMarkerColor(qq2+1);
h_w[qq2]->SetLineColor(qq2+1);
} else {
h_w[qq2]->SetMarkerStyle(3);
h_w[qq2]->SetMarkerSize(1); 
h_w[qq2]->SetMarkerColor(qq2+1-6);
h_w[qq2]->SetLineColor(qq2+1-6);
}; 
h_w[qq2]->Draw("e1P same"); 
 qqq.str("");
 qqq << "Q^{2} = " << Q2_bin << " GeV^{2}";
 leg_w_int->AddEntry(h_w[qq2],qqq.str().c_str(),"p");
 qqq.str("");
 };
leg_w_int->Draw();
//file_out->Close();

*/
}; //end of main program

