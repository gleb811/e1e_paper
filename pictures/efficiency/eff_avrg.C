TH1D *m_pip_p,*m_pip_pim,*m_pim_p;
TH1D *theta_p,*theta_pim,*theta_pip;
TH1D *alpha_p,*alpha_pim,*alpha_pip;

Float_t Q2_bin,W_bin[30];

ostringstream qqq;
 
void eff_avrg() {
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
TCanvas *c = new TCanvas("c","c",700,700);

leg_w_int = new TLegend(0.11,0.7,0.89,0.89); 
leg_w_int->SetNColumns(3);
leg_w_int->SetFillStyle(0);




//Define input files

TFile *file_eff = new TFile("out_eff.root","READ");


 for (Int_t qq2=2; qq2<3;qq2++) {
 Q2_bin = 0.425 + 0.05*qq2;


for (Int_t i=13; i<14;i++) {
// for (i=get_min_w(Q2_bin); i<get_max_w(Q2_bin);i++) {
 W_bin[i] = 1.3125+0.025*i; 

read_eff(file_eff,i);

draw_1d_canvas(i,qq2);
//alpha_p->Draw();

};
};

c->Print("eff_plot.pdf");

}; //end of main program


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
//if (name == "h1prj_inv_m_pip_p_1_") leg->AddEntry(h,"empty cells from model","p");
//if (name == "h1prj_inv_m_pip_p_1_noempty") leg->AddEntry(h,"zero in empty cells","p");
//if (name == "model_pipP_1_") leg->AddEntry(h,"model","p");
h->Draw(draw_options.c_str());
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
if ((Q2_bin> 0.4)&&(Q2_bin< 0.45))get_min_w = 8;
return get_min_w;
};



void read_eff(TFile *file_eff,Int_t i) {

file_eff->cd();

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_inv_m_pip_pim_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pip_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_inv_m_pim_p_1_";
gDirectory->GetObject(qqq.str().c_str(),m_pim_p);
qqq.str("");


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_P_";
gDirectory->GetObject(qqq.str().c_str(),theta_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIm_";
gDirectory->GetObject(qqq.str().c_str(),theta_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_th_PIp_";
gDirectory->GetObject(qqq.str().c_str(),theta_pip);
qqq.str("");


qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h1prj_alpha_PIpPIm_pipf_";
gDirectory->GetObject(qqq.str().c_str(),alpha_p);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h2prj_alpha_PPIp_piPIm_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pim);
qqq.str("");

qqq.str("");
qqq << "q2_" << Q2_bin << "/w_" << W_bin[i] << "/h3prj_alpha_PPIm_piPIp_";
gDirectory->GetObject(qqq.str().c_str(),alpha_pip);
qqq.str("");

};




void draw_1d_canvas( Int_t i, Int_t qq2 ) {

 c->Divide(3,3); 
 
leg = new TLegend(0.3,0.965,0.99,1.0); 
leg->SetNColumns(3);
leg->SetFillStyle(0);

// pi+ proton mass model

qqq << "Q^{2} = " << Q2_bin << " GeV^{2}" <<", W = " << W_bin[i] <<" GeV" ;

//draw_1d_hist(1,model_1[i]->Projection(0,""),qqq.str(),"model_pipP_1_","d#sigma/dM, #mubn/GeV","m_{#pi^{+}p}, GeV",3,"e1","mass",i);
//qqq.str("");

 Float_t max;
 
 max = m_pip_p->GetMaximum();
 m_pip_p->SetMaximum(max);

draw_1d_hist(1,m_pip_p,qqq.str(),"h1prj_inv_m_pip_p_1_","avrg eff","m_{#pi^{+}p} (GeV)",1,"e1","mass",i);
qqq.str("");


m_pip_pim->SetMaximum(max);

draw_1d_hist(2,m_pip_pim,"","h1prj_inv_m_pip_pim_1_","avrg eff","m_{#pi^{+}#pi^{-}} (GeV)",1,"e1","mass",i);

m_pim_p->SetMaximum(max);

draw_1d_hist(3,m_pim_p,"","h3prj_inv_m_pim_p_1_","avrg eff","m_{#pi^{-}p} (GeV)",1,"e1","mass",i);

max = (theta_p->GetMaximum())/1.5;
theta_p->SetMaximum(max);

draw_1d_hist(4,theta_p,"","h1prj_th_P_","avrg eff","#theta_{p'} in c.m. (deg)",1,"e1","theta",i);

theta_pim->SetMaximum(max);

draw_1d_hist(5,theta_pim,"","h1prj_th_PIm_","avrg eff","#theta_{#pi^{-}} in c.m. (deg)",1,"e1","theta",i);


theta_pip->SetMaximum(max);

draw_1d_hist(6,theta_pip,"","h1prj_th_PIp_","avrg eff","#theta_{#pi^{+}} in c.m. (deg)",1,"e1","theta",i);

max = (alpha_p->GetMaximum())/1.5;
alpha_p->SetMaximum(max);

draw_1d_hist(7,alpha_p,"","h1prj_alpha_PIpPIm_pipf_","avrg eff","#alpha_{p'} (deg)",1,"e1","alpha",i);

alpha_pim->SetMaximum(max);

draw_1d_hist(8,alpha_pim,"","h2prj_alpha_PPIp_piPIm_","avrg eff","#alpha_{#pi^{-}} (deg)",1,"e1","alpha",i);

alpha_pip->SetMaximum(max);

draw_1d_hist(9,alpha_pip,"","h3prj_alpha_PPIm_piPIp_","avrg eff","#alpha_{#pi^{+}} (deg)",1,"e1","alpha",i);



 c->cd();
// leg->Draw();
 c->Update();




};
