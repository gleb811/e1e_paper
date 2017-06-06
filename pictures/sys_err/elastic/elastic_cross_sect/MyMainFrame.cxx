#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include <math.h>
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TGObject.h"
#include "TObject.h"
#include "TSystem.h"
#include "TMinuit.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGButtonGroup.h>
#include <RQ_OBJECT.h>
#include <TGNumberEntry.h>
#include <TGProgressBar.h>
#include <TGLabel.h>
#include <stdio.h>
#include <dlfcn.h>
#include "MyMainFrame.h"
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <TGFileDialog.h>
#include <GuiTypes.h>
#include <TGDoubleSlider.h>
#include <TGComboBox.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TError.h> 
#include <auto_ptr.h>
#ifndef __CINT__
#include <cstdlib>
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <RooLandau.h> 
#include <RooNumConvPdf.h>
#include <RooDataHist.h>
#include "RooBinning.h"
#include <sys/types.h>
#include <wait.h>
#include <unistd.h>
#include <cstring>
#include <getopt.h>
#include <cstdlib>



 using namespace std; 
 
#define _USE_MATH_DEFINES
 
extern "C"  { double elasrad_(float&, float&, float&, float&);    }




    void MyMainFrame::MainFrame(UChar_t flag, Float_t E_beam, Short_t nfiles, string inp_files[], Short_t nfiles_sim, string inp_files_sim[], string outfile_in) { 


	inpfile_inp = inp_files[0];
	outfile_inp = outfile_in;
	
	E0 = E_beam;

        n_files = nfiles;
        n_files_sim = nfiles_sim;
 
 
        file = new string[nfiles];
        file_sim = new string[nfiles_sim];
  
  
        file = inp_files;
        file_sim = inp_files_sim;	
	
	data_sim = flag;
	

DoDraw();


   
}

void MyMainFrame::DoDraw() {

 
 
       t20tot21(); 
   
}








    void MyMainFrame::t20tot21() { 
    

    ostringstream adc_num;
    ostringstream tdc_num;
    ostringstream ref_tdc;
    
    
	 
    Long64_t j;
    
 Double_t * adc_offset;
 adc_offset = new Double_t [12];
 Double_t * adc_cut;
 adc_cut = new Double_t [12]; 
 Short_t m;
 Short_t m_old=0;
 Long64_t i,nstart,nstop,n_incl,n_elast;
 
TH1I *hist_adc_off[12]; 

 TFile *finp;
 UInt_t block_curr = 0;
 UInt_t block_prev = 0;
 UInt_t block_total = 0; 
 Float_t Qfull = 0.;
 Float_t LiveTime,inclusive,elastic;
 Float_t P_EL_old,th_EL_old,W_old,Q2_old,P_EL_new,th_EL_new,W_new,Q2_new;
 
 
 
 
 
 
 
 
/*  	t21 = new TTree("t21","Tree 21"); 
  t21->SetDirectory(0);
  
*/  
	UChar_t         npart,segment,sector,indtype,pdhit;
	Short_t pmt_hit;
	Int_t PdHit_EL;
	Float_t * p; 
        p = new Float_t [20];
        Float_t  P_EL,ph_EL,th_EL,W,Q2,m_proton,NpheCC_EL,ECtot_EL;
	Float_t  x_EL,y_EL,z_EL;
	Float_t  ECin_EL,ECout_EL;
	UInt_t block;
        Long64_t gpart,k,last_i,last_k;
	Float_t  q_l,Qdiff,Qcurr,Qprev,Qtotal,deltaQ;
	Int_t sc_part_local;
	Float_t sc_pd_local,sc_sect_local;
	Float_t sc_x,sc_y,sc_z,a,b;
	Float_t nx,ny,nz,t,theta_cc,ECT,EIN,EOUT,par1,par2,th_min;
	Float_t sx,sy,sz,px,py,pz,ph_cc;
	bool cut_fiduch;
	m_proton = 0.93821;

	
/*	sx=-0.000784;
	sy=0.;
	sz=-0.00168;	
	
//	E0 = 2.039;
	
	m_proton = 0.93821;
 
t21->Branch("indtype",&indtype);  
t21->Branch("npart",&npart); 
t21->Branch("P_EL",&P_EL);
t21->Branch("block",&block_total);
t21->Branch("deltaQ",&deltaQ);
t21->Branch("n_incl",&n_incl);
t21->Branch("n_elast",&n_elast);
t21->Branch("ph_EL",&ph_EL);
t21->Branch("th_EL",&th_EL);
t21->Branch("W",&W);
t21->Branch("Q2",&Q2);
t21->Branch("NpheCC_EL",&NpheCC_EL);
t21->Branch("ECtot_EL",&ECtot_EL);
t21->Branch("ECin_EL",&ECin_EL);
t21->Branch("ECout_EL",&ECout_EL);
t21->Branch("x_EL",&x_EL);
t21->Branch("y_EL",&y_EL);
t21->Branch("z_EL",&z_EL);
t21->Branch("PdHit_EL",&pdhit);
t21->Branch("sc_x",&sc_x);
t21->Branch("sc_y",&sc_y);
t21->Branch("sc_z",&sc_z);
t21->Branch("pmt_hit",&pmt_hit);
t21->Branch("segment",&segment);
t21->Branch("theta_cc",&theta_cc);
t21->Branch("ph_cc",&ph_cc);
t21->Branch("sector",&sector); 
 
 
 
 
 
 
 
 
 
 
*/ 
 
 
 
 
 if (data_sim == 1) {
 indtype = 1;
 };

TH2F *hist_sector1 = new TH2F("cc_sector1","cc_sector1",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector2 = new TH2F("cc_sector2","cc_sector2",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector3 = new TH2F("cc_sector3","cc_sector3",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector4 = new TH2F("cc_sector4","cc_sector4",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector5 = new TH2F("cc_sector5","cc_sector5",200, -5., 55.,200., -25., 25.);
TH2F *hist_sector6 = new TH2F("cc_sector6","cc_sector6",200, -5., 55.,200., -25., 25.);


TH2F *hist_nphe_sector1 = new TH2F("cc_nphe_sector1","cc_nphe_sector1",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector2 = new TH2F("cc_nphe_sector2","cc_nphe_sector2",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector3 = new TH2F("cc_nphe_sector3","cc_nphe_sector3",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector4 = new TH2F("cc_nphe_sector4","cc_nphe_sector4",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector5 = new TH2F("cc_nphe_sector5","cc_nphe_sector5",200, -5., 55.,200., -25., 25.);
TH2F *hist_nphe_sector6 = new TH2F("cc_nphe_sector6","cc_nphe_sector6",200, -5., 55.,200., -25., 25.);


TH1F *nphe_sector1 = new TH1F("nphe_sector1","nphe_sector1",501, -1., 500.);
TH1F *nphe_sector1_after = new TH1F("nphe_sector1_after","nphe_sector1_after",501, -1., 500.);
TH1F *nphe_sector2 = new TH1F("nphe_sector2","nphe_sector2",501, -1., 500.);
TH1F *nphe_sector2_after = new TH1F("nphe_sector2_after","nphe_sector2_after",501, -1., 500.);
TH1F *nphe_sector3 = new TH1F("nphe_sector3","nphe_sector3",501, -1., 500.);
TH1F *nphe_sector3_after = new TH1F("nphe_sector3_after","nphe_sector3_after",501, -1., 500.);
TH1F *nphe_sector4 = new TH1F("nphe_sector4","nphe_sector4",501, -1., 500.);
TH1F *nphe_sector4_after = new TH1F("nphe_sector4_after","nphe_sector4_after",501, -1., 500.);
TH1F *nphe_sector5 = new TH1F("nphe_sector5","nphe_sector5",501, -1., 500.);
TH1F *nphe_sector5_after = new TH1F("nphe_sector5_after","nphe_sector5_after",501, -1., 500.);
TH1F *nphe_sector6 = new TH1F("nphe_sector6","nphe_sector6",501, -1., 500.);
TH1F *nphe_sector6_after = new TH1F("nphe_sector6_after","nphe_sector6_after",501, -1., 500.);


TH2C *hist_ltime = new TH2C("ltime","ltime",60000, -.5, 59999.5,500, 0., 1.1);
TH1C *hist_ltime_1d = new TH1C("ltime_1d","ltime_1d",1100, 0., 1.1);
TH2C *hist_n_incl = new TH2C("n_incl","n_incl",60000, -.5, 59999.5,400, 40000., 80000.);
TH1C *hist_n_incl_1d = new TH1C("n_incl_1d","n_incl_1d",90000, 0., 90000.);
TH2C *hist_n_elast = new TH2C("n_elast","n_elast",60000, -.5, 59999.5,400, 0., 40000.);
TH1C *hist_n_elast_1d = new TH1C("n_elast_1d","n_elast_1d",40000, 0., 40000.);


TH2F *hist_ein_etot_sector1 = new TH2F("ein_1","ein_1",100, 0., 2.1,100, 0., 1.1);
TH2F *hist_ein_etot_sector2 = new TH2F("ein_2","ein_2",100, 0., 2.1,100, 0., 1.1);
TH2F *hist_ein_etot_sector3 = new TH2F("ein_3","ein_3",100, 0., 2.1,100, 0., 1.1);
TH2F *hist_ein_etot_sector4 = new TH2F("ein_4","ein_4",100, 0., 2.1,100, 0., 1.1);
TH2F *hist_ein_etot_sector5 = new TH2F("ein_5","ein_5",100, 0., 2.1,100, 0., 1.1);
TH2F *hist_ein_etot_sector6 = new TH2F("ein_6","ein_6",100, 0., 2.1,100, 0., 1.1);

TH2F *hist_ein_eout_sector1 = new TH2F("ein_eout_1","ein_eout_1",100, 0., 0.5,100, 0., 0.5);
TH2F *hist_ein_eout_sector2 = new TH2F("ein_eout_2","ein_eout_2",100, 0., 0.5,100, 0., 0.5);
TH2F *hist_ein_eout_sector3 = new TH2F("ein_eout_3","ein_eout_3",100, 0., 0.5,100, 0., 0.5);
TH2F *hist_ein_eout_sector4 = new TH2F("ein_eout_4","ein_eout_4",100, 0., 0.5,100, 0., 0.5);
TH2F *hist_ein_eout_sector5 = new TH2F("ein_eout_5","ein_eout_5",100, 0., 0.5,100, 0., 0.5);
TH2F *hist_ein_eout_sector6 = new TH2F("ein_eout_6","ein_eout_6",100, 0., 0.5,100, 0., 0.5);


TH1I *hist_W_old = new TH1I("W_old","W_old",100, 0.8, 1.1);
TH1I *hist_W_new = new TH1I("W_new","W_new",100, 0.8, 1.1);



TH1F *elast_data_sector[6][31];
TH1F *elast_gen_sector[6][31];
TH1F *elast_rec_sector[6][31];
TH1F *elast_eff_sector[6][31];
TH1F *elast_xsect_sector[6][31];
TH1F *elast_ratio_sector[6][31];

TH1F *elast_model;
elast_model = new TH1F("elast_model","elast_model",35, 15.5, 50.5);
TH1F *elast_factor;
elast_factor = new TH1F("elast_factor","elast_factor",35, 15.5, 50.5);

float es,theta_d,tt,wcut;
es = 2.039;
theta_d = 25;
tt = 0.002246181491;
wcut = 1.02;

for (i=0; i<35; i++) {
theta_d = 16 + i;
cout  << "thetd_d = " << theta_d << "elasrar = "  << elasrad_(es,theta_d,tt,wcut) << "\n"; 
elast_model->Fill(theta_d, elasrad_(es,theta_d,tt,wcut));
};

for (i=1; i<=6; i++) {

for (m=0; m<=30; m++) {
adc_num.str("");
adc_num  << "sector" << i << "_phi" << m; 
elast_data_sector[i-1][m] = new TH1F(adc_num.str().c_str(),adc_num.str().c_str(),35, 15.5, 50.5);
};

for (m=0; m<=30; m++) {
adc_num.str("");
adc_num  << "gen_sector" << i << "_phi" << m; 
elast_gen_sector[i-1][m] = new TH1F(adc_num.str().c_str(),adc_num.str().c_str(),35, 15.5, 50.5);
};

for (m=0; m<=30; m++) {
adc_num.str("");
adc_num  << "rec_sector" << i << "_phi" << m; 
elast_rec_sector[i-1][m] = new TH1F(adc_num.str().c_str(),adc_num.str().c_str(),35, 15.5, 50.5);
};

};



  for (m=1; m<=n_files; m++) {
  
  finp = new TFile(file[m-1].c_str()); 
  
  cout << "Processing file" << m << "\n";
  
  
  
  

    
  TTree *t21 = (TTree*)finp->Get("t21");
  TBranch *br_deltaQ = t21->GetBranch("deltaQ"); 
  TBranch *br_Q2 = t21->GetBranch("Q2");   
  TBranch *br_n_incl = t21->GetBranch("n_incl");
  TBranch *br_n_elast = t21->GetBranch("n_elast");   
  TBranch *br_sector = t21->GetBranch("sector"); 
  TBranch *br_ph_cc = t21->GetBranch("ph_cc");
  TBranch *br_theta_cc = t21->GetBranch("theta_cc");      
  TBranch *br_p_el = t21->GetBranch("P_EL");  
  TBranch *br_ectot_el = t21->GetBranch("ECtot_EL"); 
  TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
  TBranch *br_ecin_el = t21->GetBranch("ECin_EL"); 
  TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL"); 
  TBranch *br_th_EL = t21->GetBranch("th_EL"); 
  TBranch *br_ph_EL = t21->GetBranch("ph_EL"); 
  TBranch *br_LiveTime = t21->GetBranch("LiveTime");
  TBranch *br_block =  t21->GetBranch("block");     
   
//   TBranch *br_q_l = t20->GetBranch("q_l"); 
//  TLeaf *adc1 = eventn->GetLeaf("adc1");

                         
//  cout << "E0 = " << E0 << "\n";
  
//  cout << "Processing file " << m << "  N enteries = " << br_npart->GetEntries() << "\n";


     
//   for (j=0; j<12; j++) {
//  adc_num << "adc" << j+1;
//  
//  adc[j] = eventn->GetLeaf(adc_num.str().c_str())->GetValue();

//  hist_adc_off[j]->Fill(eventn->GetLeaf(adc_num.str().c_str())->GetValue(),1);
//  adc_num.str("");
//  };h10tot21
  
//  hist_adc[3]->Fill(adc[3],1.);  
  

  
 
  
//  };



  
  
  
  // Creating and filling new tree
  



// Loop over old tree 

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
  
//  cout << M_PI << "\n";
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;

  br_deltaQ->GetEntry(i);
   br_Q2->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_block->GetEntry(i);  
  br_sector->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);   
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i); 
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();

  theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
  ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
          

 
   EIN = br_ecin_el->GetLeaf("ECin_EL")->GetValue();
       EOUT = br_ecout_el->GetLeaf("ECout_EL")->GetValue();  
      
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (EOUT+EIN)) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (EOUT+EIN);    
          };
   block_curr = br_block->GetLeaf("block")->GetValue(); 
   

   
//   if (block_curr == 861)cout << "block_curr = " << block_curr << "\n"; 
   
    if ((m_old != 0) && (m != m_old))block_total=block_total+block_prev;
    
    m_old=m;

   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();  
   if (block_curr != block_prev) {
//   cout << "block_total+block_curr = " << block_total+block_curr << " block_prev = " << block_prev << " Qfull = " << Qfull << " n_elast = " << br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue() << " m = " << m <<  "\n";
   hist_ltime->Fill(block_total+block_curr,LiveTime);
   hist_ltime_1d->Fill(LiveTime);
   hist_n_incl->Fill(block_total+block_curr,inclusive);
hist_n_incl_1d->Fill(inclusive);  
 hist_n_elast->Fill(block_total+block_curr,elastic);
hist_n_elast_1d->Fill(elastic); 

if ((LiveTime>0.89) && (LiveTime<0.94) && (inclusive>62000) && (inclusive<78000) && (elastic>20000) &&  (elastic<24000)) { 
   Qfull = Qfull + br_deltaQ->GetLeaf("deltaQ")->GetValue();
   };
   
//   if ((LiveTime < 0.) || (LiveTime > 1.1))cout << "LiveTime = " << LiveTime << "\n";	  
   
   block_prev = block_curr;
   };
	  
if ((LiveTime>0.89) && (LiveTime<0.94) && (inclusive>62000) && (inclusive<78000) && (elastic>20000) &&  (elastic<24000)) {

//cout << "P_EL = " << P_EL << " P_EL_corr = " << corrfunc.correct_pel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL) << "\n";
P_EL_new = corrfunc.correct_pel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
th_EL_new = corrfunc.correct_thel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
P_EL_old = P_EL;
P_EL = P_EL_new;
th_EL_old = th_EL;
th_EL = th_EL_new;

W_new = pow(float((E0+m_proton-P_EL)),2);
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*cos(ph_EL*M_PI/180.),2);	  
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*sin(ph_EL*M_PI/180.),2);	
W_new = W_new -pow(float(E0 - P_EL*cos(th_EL*M_PI/180.)),2); 
W_new = sqrt(W_new);

 
 
W = W_new;

  th_min=(9.5+17./(P_EL+0.2));
  par1=0.705+1.1*P_EL;
  par2=-63.5-30.*P_EL; 

   a=37.14*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1500./th_EL/th_EL));
     
   b=-37.14*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1500./th_EL/th_EL));  
 
  if ((P_EL > 0.461)  && (W<1.02)   && (W>0.85)) {

hist_W_new->Fill(W);  
hist_W_old->Fill(W);  
  
//  cout << "sector = " << br_sector->GetLeaf("sector")->GetValue() << "\n";
  switch (int(br_sector->GetLeaf("sector")->GetValue())) {
case 1 : 
if (((ECT/P_EL) > (-0.03606*P_EL*P_EL+0.1357*P_EL+0.08835)) && ((ECT/P_EL) < (0.008336*P_EL*P_EL-0.0301*P_EL)+0.38)) {

if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+360.) ||  (ph_EL < a)){

nphe_sector1->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

 hist_sector1->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 hist_nphe_sector1->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  
 
 
/*if (ph_EL <30.) { 

//cout  << "ph_EL = " << int((ph_EL-1)/2)   <<  "\n";
//elast_data_sector[0][ int((29-ph_EL)/2) ]->Fill(th_EL);
elast_data_sector[0][ int((ph_EL+31)/2) ]->Fill(th_EL);
} else {

//cout  << "ph_EL = " << int((ph_EL-300-1)/2)   <<  "\n";
elast_data_sector[0][ int((ph_EL-329)/2) ]->Fill(th_EL);
};
*/
//     if (th_EL >=20) {
//     cout  << "th_EL = " << int(th_EL -19) << "\n";
//1`      };
 
     };
     };
     };
     };
     };
     };
     };
     };
break; 
case 2 : 
if (((ECT/P_EL) > (-0.03616*P_EL*P_EL+0.1228*P_EL+0.1099)) && ((ECT/P_EL) < (0.001805*P_EL*P_EL-0.01003*P_EL)+0.419)) {

if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+60.) &&  (ph_EL < a+60.)){

nphe_sector2->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
 hist_sector2->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 

hist_nphe_sector2->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  

 


//cout  << "ph_EL = " << int((ph_EL-1)/2)   <<  "\n";
//elast_data_sector[1][ int((ph_EL-29)/2) ]->Fill(th_EL);


 };
 };
 };
 };
 };
 };
 };
 }; 
 
 
break; 
case 3 : 

if (((ECT/P_EL) > (-0.024*P_EL*P_EL+0.1014*P_EL+0.09593)) && ((ECT/P_EL) < (0.003918*P_EL*P_EL-0.01219*P_EL)+0.3616)) {


if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+120.) &&  (ph_EL < a+120.)){



nphe_sector3->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {



 hist_sector3->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 

hist_nphe_sector3->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue()); 


//elast_data_sector[2][ int((ph_EL-89)/2) ]->Fill(th_EL);
 
  };
  };
  };
  };
  };
  };
  };
  };  
  
break; 
case 4 : 

if (((ECT/P_EL) > (-0.03557*P_EL*P_EL+0.1322*P_EL+0.0951)) && ((ECT/P_EL) < (0.001583*P_EL*P_EL-0.008766*P_EL)+0.3724)) {




if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+180.) &&  (ph_EL < a+180.)){



nphe_sector4->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {


 hist_sector4->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 
hist_nphe_sector4->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  
 
//elast_data_sector[3][ int((ph_EL-149)/2) ]->Fill(th_EL); 


  };
  };
  };
  };
  };
  };
  };
  }; 
break;   
case 5 : 

if (((ECT/P_EL) > (-0.03453*P_EL*P_EL+0.132*P_EL+0.08782)) && ((ECT/P_EL) < (0.004611*P_EL*P_EL-0.01966*P_EL)+0.3846)) {



if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+240.) &&  (ph_EL < a+240.)){


nphe_sector5->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {


 hist_sector5->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 
hist_nphe_sector5->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());   
 
//elast_data_sector[4][ int((ph_EL-209)/2) ]->Fill(th_EL);  
 
  }; 
  };
  };
  };
  }; 
  };
  };
  };  
break;
case 6 : 

if (((ECT/P_EL) > (-0.033*P_EL*P_EL+0.1336*P_EL+0.08992)) && ((ECT/P_EL) < (0.002539*P_EL*P_EL-0.005779*P_EL)+0.33669)) {



if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+300.) &&  (ph_EL < a+300.)){

nphe_sector6->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

 hist_sector6->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 
 hist_nphe_sector6->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());   
 
//elast_data_sector[5][ int((ph_EL-269)/2) ]->Fill(th_EL);   
 
 
 };
 };
 };
 };
 };
 };
 };
 }; 
break;
  };
  };
    };
    };
    
  t21->Delete();

 finp->Close();    
    };
    
    
    
    
    
      TH2F *avrg_nphe_sector1 = (TH2F*)hist_nphe_sector1->Clone();
      avrg_nphe_sector1->SetTitle("avrg_nphe_sector1");
      avrg_nphe_sector1->SetName("avrg_nphe_sector1");
      avrg_nphe_sector1->Divide(hist_sector1);
     
      
      TH2F *avrg_nphe_sector2 = (TH2F*)hist_nphe_sector2->Clone();
      avrg_nphe_sector2->SetTitle("avrg_nphe_sector2");
      avrg_nphe_sector2->SetName("avrg_nphe_sector2");
      avrg_nphe_sector2->Divide(hist_sector2);
     
      TH2F *avrg_nphe_sector3 = (TH2F*)hist_nphe_sector3->Clone();
      avrg_nphe_sector3->SetTitle("avrg_nphe_sector3");
      avrg_nphe_sector3->SetName("avrg_nphe_sector3");
      avrg_nphe_sector3->Divide(hist_sector3);
      
      TH2F *avrg_nphe_sector4 = (TH2F*)hist_nphe_sector4->Clone();
      avrg_nphe_sector4->SetTitle("avrg_nphe_sector4");
      avrg_nphe_sector4->SetName("avrg_nphe_sector4");
      avrg_nphe_sector4->Divide(hist_sector4);
     
      TH2F *avrg_nphe_sector5 = (TH2F*)hist_nphe_sector5->Clone();
      avrg_nphe_sector5->SetTitle("avrg_nphe_sector5");
      avrg_nphe_sector5->SetName("avrg_nphe_sector5");
      avrg_nphe_sector5->Divide(hist_sector5);

      TH2F *avrg_nphe_sector6 = (TH2F*)hist_nphe_sector6->Clone();
      avrg_nphe_sector6->SetTitle("avrg_nphe_sector6");
      avrg_nphe_sector6->SetName("avrg_nphe_sector6");
      avrg_nphe_sector6->Divide(hist_sector6);

    
    
    
    
    
 //end of first loop
 
 
 
   for (m=1; m<=n_files; m++) {
  
  finp = new TFile(file[m-1].c_str()); 
  
  cout << "Reading file" << m << " second time" << "\n";
  
     
  TTree *t21 = (TTree*)finp->Get("t21");  
  TBranch *br_deltaQ = t21->GetBranch("deltaQ"); 
  TBranch *br_Q2 = t21->GetBranch("Q2");   
  TBranch *br_n_incl = t21->GetBranch("n_incl");
  TBranch *br_n_elast = t21->GetBranch("n_elast");   
  TBranch *br_LiveTime = t21->GetBranch("LiveTime");
  TBranch *br_block =  t21->GetBranch("block");    
  TBranch *br_sector = t21->GetBranch("sector"); 
  TBranch *br_ph_cc = t21->GetBranch("ph_cc");
  TBranch *br_theta_cc = t21->GetBranch("theta_cc");      
  TBranch *br_p_el = t21->GetBranch("P_EL");  
  TBranch *br_ectot_el = t21->GetBranch("ECtot_EL"); 
  TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
  TBranch *br_ecin_el = t21->GetBranch("ECin_EL"); 
  TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL"); 
  TBranch *br_th_EL = t21->GetBranch("th_EL"); 
  TBranch *br_ph_EL = t21->GetBranch("ph_EL");        

  Bool_t adc_cut_switch,tdc_cut_switch;
  
  Qdiff = 0.;
  Qcurr = 0.;
  Qprev = 0.;
  Qtotal = 0.;
  k = 0;
  block = 0;
  last_k = 0;
  nstart = 0;
  nstop = 0;
  n_incl = 0;
  n_elast = 0;
  
  
//  cout << M_PI << "\n";
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  
  Qprev = Qcurr;


  br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_sector->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);   
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i); 
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();
  
P_EL_new = corrfunc.correct_pel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
th_EL_new = corrfunc.correct_thel_e1_2039_2250_feb09(P_EL,th_EL,ph_EL);
P_EL_old = P_EL;
P_EL = P_EL_new;
th_EL_old = th_EL;
th_EL = th_EL_new;

W_new = pow(float((E0+m_proton-P_EL)),2);
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*cos(ph_EL*M_PI/180.),2);	  
W_new = W_new -pow(P_EL*sin(th_EL*M_PI/180.)*sin(ph_EL*M_PI/180.),2);	
W_new = W_new -pow(float(E0 - P_EL*cos(th_EL*M_PI/180.)),2); 
W_new = sqrt(W_new);

// cout  << " P_EL_old "  <<  P_EL_old  << " P_EL = " << P_EL << "\n";
 
W = W_new;  
  
  th_min=(9.5+17./(P_EL+0.2));
  par1=0.705+1.1*P_EL;
  par2=-63.5-30.*P_EL; 
  theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
  ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();
          
   a=37.14*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1500./th_EL/th_EL));
     
   b=-37.14*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1500./th_EL/th_EL));  
 
    EIN = br_ecin_el->GetLeaf("ECin_EL")->GetValue();
       EOUT = br_ecout_el->GetLeaf("ECout_EL")->GetValue();
      
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (EOUT+EIN)) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = EOUT+EIN;    
          };
	  

	  
   LiveTime=br_LiveTime->GetLeaf("LiveTime")->GetValue();
   inclusive = br_n_incl->GetLeaf("n_incl")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue();
   elastic =  br_n_elast->GetLeaf("n_elast")->GetValue()/br_deltaQ->GetLeaf("deltaQ")->GetValue(); 	  
	  
if ((LiveTime>0.89) && (LiveTime<0.94) && (inclusive>62000) && (inclusive<78000) && (elastic>20000) &&  (elastic<24000)) {	 
 


	  
	  
  if ((P_EL > 0.461)  && (W>0.85)  && (W<1.02)){
  
  
//  cout << "sector = " << br_sector->GetLeaf("sector")->GetValue() << "\n";
  switch (int(br_sector->GetLeaf("sector")->GetValue())) {
case 1 : 
if (((ECT/P_EL) > (-0.03606*P_EL*P_EL+0.1357*P_EL+0.08835)) && ((ECT/P_EL) < (0.008336*P_EL*P_EL-0.0301*P_EL)+0.38)) {

if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {


if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+360.) ||  (ph_EL < a)){

//cout  << "avrg = " << avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.))   << "\n";

if (avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
nphe_sector1_after->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.); 
hist_ein_etot_sector1->Fill(P_EL,EIN/ECT,1);
hist_ein_eout_sector1->Fill(EOUT/P_EL,EIN/P_EL,1);

if (ph_EL <30.) { 

//cout  << "ph_EL = " << int((ph_EL-1)/2)   <<  "\n";
//elast_data_sector[0][ int((29-ph_EL)/2) ]->Fill(th_EL);
elast_data_sector[0][ int((ph_EL+31)/2) ]->Fill(th_EL);
} else {

//cout  << "ph_EL = " << int((ph_EL-300-1)/2)   <<  "\n";
elast_data_sector[0][ int((ph_EL-329)/2) ]->Fill(th_EL);
};

/*if (P_EL>1.1) {
cout << "P_EL = " << P_EL << "ECin/ECtot = " << br_ecin_el->GetLeaf("ECin_EL")->GetValue()/ECT  <<"\n";
   };
   */
     };
     };
     };
     };
     };
     };
     };
     };
     };
break; 
case 2 : 
if (((ECT/P_EL) > (-0.03616*P_EL*P_EL+0.1228*P_EL+0.1099)) && ((ECT/P_EL) < (0.001805*P_EL*P_EL-0.01003*P_EL)+0.419)) {

if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+60.) &&  (ph_EL < a+60.)){

if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
  
if (avrg_nphe_sector2->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) >70.) {
nphe_sector2_after->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


hist_ein_etot_sector2->Fill(P_EL,EIN/ECT,1);
hist_ein_eout_sector2->Fill(EOUT/P_EL,EIN/P_EL,1);

elast_data_sector[1][ int((ph_EL-29)/2) ]->Fill(th_EL);

 };
 };
 };
 };
 };
 };
 };
 }; 
 };
 
break; 
case 3 : 

if (((ECT/P_EL) > (-0.024*P_EL*P_EL+0.1014*P_EL+0.09593)) && ((ECT/P_EL) < (0.003918*P_EL*P_EL-0.01219*P_EL)+0.3616)) {


if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+120.) &&  (ph_EL < a+120.)){
 
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
  
if (avrg_nphe_sector3->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
nphe_sector3_after->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);

hist_ein_etot_sector3->Fill(P_EL,EIN/ECT,1);
hist_ein_eout_sector3->Fill(EOUT/P_EL,EIN/P_EL,1);

elast_data_sector[2][ int((ph_EL-89)/2) ]->Fill(th_EL);

  };
  };
  };
  };
  };
  };
  }; 
  }; 
  };  
break; 
case 4 : 

if (((ECT/P_EL) > (-0.03557*P_EL*P_EL+0.1322*P_EL+0.0951)) && ((ECT/P_EL) < (0.001583*P_EL*P_EL-0.008766*P_EL)+0.3724)) {




if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+180.) &&  (ph_EL < a+180.)){

if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
  
if (avrg_nphe_sector4->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
nphe_sector4_after->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


hist_ein_etot_sector4->Fill(P_EL,EIN/ECT,1);
hist_ein_eout_sector4->Fill(EOUT/P_EL,EIN/P_EL,1);

elast_data_sector[3][ int((ph_EL-149)/2) ]->Fill(th_EL);
 
  };
  };
  };
  };
  };
  };
  };
  };
  };
break;  
case 5 : 

if (((ECT/P_EL) > (-0.03453*P_EL*P_EL+0.132*P_EL+0.08782)) && ((ECT/P_EL) < (0.004611*P_EL*P_EL-0.01966*P_EL)+0.3846)) {



if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+240.) &&  (ph_EL < a+240.)){

if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
  
if (avrg_nphe_sector5->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
nphe_sector5_after->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.); 
 
hist_ein_etot_sector5->Fill(P_EL,EIN/ECT,1);
hist_ein_eout_sector5->Fill(EOUT/P_EL,EIN/P_EL,1);
 
 elast_data_sector[4][int((ph_EL-209)/2) ]->Fill(th_EL);
 
 
  };
  }; 
  };
  };
  };
  }; 
  };
  };
  };
break;
case 6 : 

if (((ECT/P_EL) > (-0.033*P_EL*P_EL+0.1336*P_EL+0.08992)) && ((ECT/P_EL) < (0.002539*P_EL*P_EL-0.005779*P_EL)+0.33669)) {



if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+300.) &&  (ph_EL < a+300.)){

if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
  
if (avrg_nphe_sector6->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) >70.) {
nphe_sector6_after->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.); 

hist_ein_etot_sector6->Fill(P_EL,EIN/ECT,1);
hist_ein_eout_sector6->Fill(EOUT/P_EL,EIN/P_EL,1);

elast_data_sector[5][ int((ph_EL-269)/2) ]->Fill(th_EL);

 };
 };
 };
 };
 }; 
 };
 };
 };
 };
break;
  };
  };
    };
    };
    
  t21->Delete();

 finp->Close();    
    };
 
 
 
   for (m=1; m<=n_files_sim; m++) {
  
  finp = new TFile(file_sim[m-1].c_str());   
  
  TTree *t21 = (TTree*)finp->Get("t21");  
 
  TBranch *br_indtype= t21->GetBranch("indtype");  
  TBranch *br_Q2 = t21->GetBranch("Q2");   
  TBranch *br_deltaQ = t21->GetBranch("deltaQ"); 
  TBranch *br_n_incl = t21->GetBranch("n_incl");
  TBranch *br_n_elast = t21->GetBranch("n_elast");   
  TBranch *br_LiveTime = t21->GetBranch("LiveTime");
  TBranch *br_block =  t21->GetBranch("block");    
  TBranch *br_sector = t21->GetBranch("sector"); 
  TBranch *br_ph_cc = t21->GetBranch("ph_cc");
  TBranch *br_theta_cc = t21->GetBranch("theta_cc");      
  TBranch *br_p_el = t21->GetBranch("P_EL");  
  TBranch *br_ectot_el = t21->GetBranch("ECtot_EL"); 
  TBranch *br_ecout_el = t21->GetBranch("ECout_EL");
  TBranch *br_ecin_el = t21->GetBranch("ECin_EL"); 
  TBranch *br_nphe_el = t21->GetBranch("NpheCC_EL"); 
  TBranch *br_th_EL = t21->GetBranch("th_EL"); 
  TBranch *br_ph_EL = t21->GetBranch("ph_EL"); 
  TBranch *br_W = t21->GetBranch("W");          

  
  cout << "Reading file"  << file_sim[m-1].c_str() << " with Monte Carlo" << "\n";
  
  
  for (i=0; i<br_sector->GetEntries(); i++) { 
  

  br_indtype->GetEntry(i);
  br_Q2->GetEntry(i); 
  br_deltaQ->GetEntry(i);
  br_n_incl->GetEntry(i);
  br_n_elast->GetEntry(i);  
  br_LiveTime->GetEntry(i);
  br_sector->GetEntry(i);  
  br_ph_cc->GetEntry(i);
  br_theta_cc->GetEntry(i);  
  br_p_el->GetEntry(i);   
  br_ectot_el->GetEntry(i); 
  br_ecin_el->GetEntry(i); 
  br_ecout_el->GetEntry(i);
  br_nphe_el->GetEntry(i); 
  br_th_EL->GetEntry(i);
  br_ph_EL->GetEntry(i); 
  br_W->GetEntry(i);   
  P_EL = br_p_el->GetLeaf("P_EL")->GetValue();
  th_EL = br_th_EL->GetLeaf("th_EL")->GetValue();
  ph_EL = br_ph_EL->GetLeaf("ph_EL")->GetValue();  
  
  W =  br_W->GetLeaf("W")->GetValue();
  
  if ((br_indtype->GetLeaf("indtype")->GetValue() == 2)    && (W<1.02)   && (W>0.85)  && (th_EL>20)){
  
//    if ((W<1.02)   &&  (W>0.89)) {


//if (ph_EL<30.)cout  << "ph_EL = " << int((ph_EL-1)/2)   <<  "\n";
if ( br_sector->GetLeaf("sector")->GetValue() == 1) {
 if (ph_EL<30.) {
 elast_gen_sector[0][ int((ph_EL+31)/2) ]->Fill(th_EL,1.);
};

//cout  << "ph_EL = " << int((ph_EL-300-1)/2)   <<  "\n";
 if (ph_EL >330.) {
 elast_gen_sector[0][ int((ph_EL-329)/2) ]->Fill(th_EL,1.);
};
};


if ( br_sector->GetLeaf("sector")->GetValue() == 2)elast_gen_sector[1][ int((ph_EL-29)/2) ]->Fill(th_EL,1.);  
if ( br_sector->GetLeaf("sector")->GetValue() == 3)elast_gen_sector[2][ int((ph_EL-89)/2) ]->Fill(th_EL,1.);  
if ( br_sector->GetLeaf("sector")->GetValue() == 4)elast_gen_sector[3][ int((ph_EL-149)/2) ]->Fill(th_EL,1.);  
if ( br_sector->GetLeaf("sector")->GetValue() == 5)elast_gen_sector[4][ int((ph_EL-209)/2) ]->Fill(th_EL,1.);  
if ( br_sector->GetLeaf("sector")->GetValue() == 6)elast_gen_sector[5][ int((ph_EL-269)/2) ]->Fill(th_EL,1.);  



//  };
  
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  if ((br_indtype->GetLeaf("indtype")->GetValue() == 1) &&  (P_EL > 0.461)  && (W<1.02)   && (W>0.89)  &&  (br_Q2->GetLeaf("Q2")->GetValue()>0.00001) && (th_EL>20)) {
  
  
   EIN = br_ecin_el->GetLeaf("ECin_EL")->GetValue();
       EOUT = br_ecout_el->GetLeaf("ECout_EL")->GetValue();  
      
  if (br_ectot_el->GetLeaf("ECtot_EL")->GetValue() > (EOUT+EIN)) { 
          ECT  = br_ectot_el->GetLeaf("ECtot_EL")->GetValue();  
         } else {
          ECT  = (EOUT+EIN);    
          };  
	  
  th_min=(9.5+17./(P_EL+0.2));	  
  par1=0.705+1.1*P_EL;
  par2=-63.5-30.*P_EL; 
  theta_cc = br_theta_cc->GetLeaf("theta_cc")->GetValue();
  ph_cc = br_ph_cc->GetLeaf("ph_cc")->GetValue();  
          
   a=37.14*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1500./th_EL/th_EL));
     
   b=-37.14*pow((sin((th_EL-th_min)*0.01745)),(par1+par2/th_EL+1500./th_EL/th_EL));  	  
	  
	  
    switch (int(br_sector->GetLeaf("sector")->GetValue())) {
case 1 : 
if (((ECT/P_EL) > (-0.03283*P_EL*P_EL+0.117*P_EL+0.07099)) && ((ECT/P_EL) < (0.01303*P_EL*P_EL-0.05557*P_EL)+0.345)) {

if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+360.) ||  (ph_EL < a)){

//nphe_sector1->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1);
 
if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

// hist_sector1->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
// hist_nphe_sector1->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  
 
if (avrg_nphe_sector1->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
 
if (ph_EL <30.) { 

//cout  << "ph_EL = " << int((ph_EL-1)/2)   <<  "\n";
//elast_data_sector[0][ int((29-ph_EL)/2) ]->Fill(th_EL);
elast_rec_sector[0][ int((ph_EL+31)/2) ]->Fill(th_EL,1.);
} else {

//cout  << "ph_EL = " << int((ph_EL-300-1)/2)   <<  "\n";
elast_rec_sector[0][ int((ph_EL-329)/2) ]->Fill(th_EL,1.);
};
//     if (th_EL >=20) {
//     cout  << "th_EL = " << int(th_EL -19) << "\n";
//1`      };

    };
 
     };
     };
     };
     };
     };
     };
     };
     };
break; 
case 2 : 
if (((ECT/P_EL) > (-0.03283*P_EL*P_EL+0.117*P_EL+0.07099)) && ((ECT/P_EL) < (0.01303*P_EL*P_EL-0.05557*P_EL)+0.345)) {

if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+60.) &&  (ph_EL < a+60.)){

//nphe_sector2->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {
// hist_sector2->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 

//hist_nphe_sector2->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  

 if (avrg_nphe_sector2->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {


//cout  << "ph_EL = " << int((ph_EL-29)/2)   <<  "\n";
elast_rec_sector[1][ int((ph_EL-29)/2) ]->Fill(th_EL,1.);



};


 };
 };
 };
 };
 };
 };
 };
 }; 
 
 
break; 
case 3 : 

if (((ECT/P_EL) > (-0.03283*P_EL*P_EL+0.117*P_EL+0.07099)) && ((ECT/P_EL) < (0.01303*P_EL*P_EL-0.05557*P_EL)+0.345)) {


if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+120.) &&  (ph_EL < a+120.)){



//nphe_sector3->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {



// hist_sector3->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 

//hist_nphe_sector3->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue()); 


if (avrg_nphe_sector3->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {

elast_rec_sector[2][ int((ph_EL-89)/2) ]->Fill(th_EL,1.);


  };
 
  };
  };
  };
  };
  };
  };
  };
  };  
  
break; 
case 4 : 

if (((ECT/P_EL) > (-0.03283*P_EL*P_EL+0.117*P_EL+0.07099)) && ((ECT/P_EL) < (0.01303*P_EL*P_EL-0.05557*P_EL)+0.345)) {




if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+180.) &&  (ph_EL < a+180.)){



//nphe_sector4->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {


// hist_sector4->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 
//hist_nphe_sector4->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  

if (avrg_nphe_sector4->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
 
elast_rec_sector[3][ int((ph_EL-149)/2) ]->Fill(th_EL,1.); 


  };

  };
  };
  };
  };
  };
  };
  };
  }; 
break;   
case 5 : 

if (((ECT/P_EL) > (-0.03283*P_EL*P_EL+0.117*P_EL+0.07099)) && ((ECT/P_EL) < (0.01303*P_EL*P_EL-0.05557*P_EL)+0.345)) {



if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+240.) &&  (ph_EL < a+240.)){


//nphe_sector5->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {


// hist_sector5->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 
//hist_nphe_sector5->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());  

if (avrg_nphe_sector5->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) { 
 
elast_rec_sector[4][ int((ph_EL-209)/2) ]->Fill(th_EL,1.);  


  };
 
  }; 
  };
  };
  };
  }; 
  };
  };
  };  
break;
case 6 : 

if (((ECT/P_EL) > (-0.03283*P_EL*P_EL+0.117*P_EL+0.07099)) && ((ECT/P_EL) < (0.01303*P_EL*P_EL-0.05557*P_EL)+0.345)) {



if (th_EL > th_min) {
if (th_EL < 50.) {
if ((ph_EL > b+300.) &&  (ph_EL < a+300.)){

//nphe_sector6->Fill(br_nphe_el->GetLeaf("NpheCC_EL")->GetValue(),1.);


if  (theta_cc > 7.0+0.0032*ph_cc+0.0499*ph_cc*ph_cc) {
if ((pow((theta_cc-45.5)/34.5,2)+pow((ph_cc)/21.,2)) <= 1.) {
if ((pow((theta_cc-45.5)/1.75,2)+pow((ph_cc)/21.,2)) > 1.) {
if  (theta_cc < 45.) {

// hist_sector6->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),1.); 
 
// hist_nphe_sector6->Fill(br_theta_cc->GetLeaf("theta_cc")->GetValue(),br_ph_cc->GetLeaf("ph_cc")->GetValue(),br_nphe_el->GetLeaf("NpheCC_EL")->GetValue());   

if (avrg_nphe_sector6->GetBinContent(int((theta_cc+5.)*200./60.),int((ph_cc+30.)*200./60.)) > 70.) {
 
elast_rec_sector[5][ int((ph_EL-269)/2) ]->Fill(th_EL,1.);   


 };
 
 
 };
 };
 };
 };
 };
 };
 };
 }; 
break;
  };
  
  
  
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  };
  
  
  t21->Delete();
  
  };
   
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 //end of second loopS   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    TFile *outFile;
     
     outFile = new TFile(outfile_inp.c_str(),"recreate");

     for (i=1; i<=6; i++) { 
     
     adc_num.str("");
     adc_num << "sector" << i << "/data";
     outFile->mkdir(adc_num.str().c_str());
     outFile->cd(adc_num.str().c_str());  
     for (m=0; m<=30; m++) {
    elast_data_sector[i-1][m]->Write("", TObject::kOverwrite);
    };
    
    };
    
     TDirectory *eprst; 
     for (i=1; i<=6; i++) { 
     
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i ;
     eprst = outFile->GetDirectory(adc_num.str().c_str());
     eprst ->mkdir("gen");
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i << "/gen" ;     
     outFile->cd(adc_num.str().c_str());  
     for (m=0; m<=30; m++) {
    elast_gen_sector[i-1][m]->Write("", TObject::kOverwrite);
    };
    
    };
    
    
//     TDirectory *eprst; 
     for (i=1; i<=6; i++) { 
     
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i ;
     eprst = outFile->GetDirectory(adc_num.str().c_str());
     eprst ->mkdir("rec");
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i << "/rec" ;     
     outFile->cd(adc_num.str().c_str());  
     for (m=0; m<=30; m++) {
     elast_rec_sector[i-1][m]->Write("", TObject::kOverwrite);  
       adc_num.str("");   
       adc_num << "eff_ sector" << i << "phi_" << -30 + 2*m+1 ;  
     elast_eff_sector[i-1][m] = (TH1F*)elast_rec_sector[i-1][m]->Clone();
    elast_eff_sector[i-1][m] ->SetTitle(adc_num.str().c_str());
       adc_num.str("");   
       adc_num << "eff_ sector" << i << "phi_" << m ;             
    elast_eff_sector[i-1][m] ->SetName(adc_num.str().c_str());     
     elast_eff_sector[i-1][m]  ->Divide( elast_gen_sector[i-1][m]);
     
       adc_num.str("");   
       adc_num << "xsect_ sector" << i << "phi_" << -30 + 2*m+1 ;  
     elast_xsect_sector[i-1][m] = (TH1F*)elast_data_sector[i-1][m]->Clone();
    elast_xsect_sector[i-1][m] ->SetTitle(adc_num.str().c_str());
       adc_num.str("");   
       adc_num << "xsect_ sector" << i << "phi_" << m ;             
    elast_xsect_sector[i-1][m] ->SetName(adc_num.str().c_str());     
     elast_xsect_sector[i-1][m]  ->Divide( elast_eff_sector[i-1][m]);     
     
    };
    
    };    
    
    for (i=1; i<=6; i++) { 
     
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i ;
     eprst = outFile->GetDirectory(adc_num.str().c_str());
     eprst ->mkdir("eff");
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i << "/eff" ;     
     outFile->cd(adc_num.str().c_str());  
     for (m=0; m<=30; m++) {
      elast_eff_sector[i-1][m]->Write("", TObject::kOverwrite);
    };
    
    };
    
    

        
    
    
    
    outFile->cd();
    outFile->mkdir("model");
    outFile->cd("model");
    elast_model->Write("", TObject::kOverwrite);
    
     outFile->cd();
     
     
     for (i=0; i<35; i++) {
theta_d = 16 + i;
//cout  << "thetd_d = " << theta_d << "elasrar = "  << elasrad_(es,theta_d,tt,wcut) << "\n"; 
elast_factor->Fill(theta_d, ((cos((theta_d-0.5)*M_PI/180.)-cos((theta_d+0.5)*M_PI/180.))*Qfull*M_PI*2./180./1.893785246E-6));
//elast_factor->Fill(    theta_d, (    (    cos(0.5*M_PI/180.)    )      )    );
};

    outFile->cd();
    outFile->mkdir("factor");
    outFile->cd("factor");
    elast_factor->Write("", TObject::kOverwrite);
    
     outFile->cd();  
     
     
     
         for (i=1; i<=6; i++) { 
     
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i ;
     eprst = outFile->GetDirectory(adc_num.str().c_str());
     eprst ->mkdir("xsect");
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i << "/xsect" ;     
     outFile->cd(adc_num.str().c_str());  
     for (m=0; m<=30; m++) {
      elast_xsect_sector[i-1][m]->Divide(elast_factor);
      elast_xsect_sector[i-1][m]->Write("", TObject::kOverwrite);      
    };
    
    }; 
    
     
         for (i=1; i<=6; i++) { 
     
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i ;
     eprst = outFile->GetDirectory(adc_num.str().c_str());
     eprst ->mkdir("ratio");
     adc_num.str("");
     outFile->cd();
     adc_num << "sector" << i << "/ratio" ;     
     outFile->cd(adc_num.str().c_str());  
     for (m=0; m<=30; m++) {
            adc_num.str("");  
            adc_num << "ratio_sector" << i << "phi_" << -30 + 2*m+1 ;       
       elast_ratio_sector[i-1][m] = (TH1F*)elast_xsect_sector[i-1][m]->Clone();
       elast_ratio_sector[i-1][m] ->SetTitle(adc_num.str().c_str());
       adc_num.str("");   
         adc_num << "ratio_sector" << i << "phi_" << m ;      
    elast_ratio_sector[i-1][m] ->SetName(adc_num.str().c_str());     
    elast_ratio_sector[i-1][m]  ->Divide(elast_model);  
    
         for (k=1; k<=35; k++) {
      
      if (elast_data_sector[i-1][m] ->GetBinContent(k) > 0.) {
      
      elast_ratio_sector[i-1][m] ->SetBinError(k,elast_ratio_sector[i-1][m] ->GetBinContent(k)/sqrt(elast_data_sector[i-1][m] ->GetBinContent(k) + elast_rec_sector[i-1][m] ->GetBinContent(k)*abs(elast_gen_sector[i-1][m] ->GetBinContent(k)-elast_rec_sector[i-1][m] ->GetBinContent(k))/pow(elast_gen_sector[i-1][m] ->GetBinContent(k),3)));
     
     };
     
     };
     
     
     elast_ratio_sector[i-1][m]->Write("", TObject::kOverwrite);      
    };
    
    };     
    
       
    
       outFile->cd();    
     
  
    hist_ein_etot_sector1->Write("", TObject::kOverwrite);
    hist_ein_etot_sector2->Write("", TObject::kOverwrite);
    hist_ein_etot_sector3->Write("", TObject::kOverwrite);
    hist_ein_etot_sector4->Write("", TObject::kOverwrite);    
    hist_ein_etot_sector5->Write("", TObject::kOverwrite);
    hist_ein_etot_sector6->Write("", TObject::kOverwrite);    
    
    hist_ein_eout_sector1->Write("", TObject::kOverwrite);
    hist_ein_eout_sector2->Write("", TObject::kOverwrite);
    hist_ein_eout_sector3->Write("", TObject::kOverwrite);
    hist_ein_eout_sector4->Write("", TObject::kOverwrite);    
    hist_ein_eout_sector5->Write("", TObject::kOverwrite);
    hist_ein_eout_sector6->Write("", TObject::kOverwrite);        
     
    
    
    
      hist_W_new->Write("", TObject::kOverwrite); 
      hist_W_old->Write("", TObject::kOverwrite);
     
     
     
      hist_ltime->Write("", TObject::kOverwrite); 
      hist_n_incl->Write("", TObject::kOverwrite);
      hist_n_elast->Write("", TObject::kOverwrite); 
      
      hist_ltime_1d->Write("", TObject::kOverwrite); 
      hist_n_incl_1d->Write("", TObject::kOverwrite);
      hist_n_elast_1d->Write("", TObject::kOverwrite);       
           
      hist_sector1->Write("", TObject::kOverwrite);
      hist_nphe_sector1->Write("", TObject::kOverwrite);
      
      
     
      avrg_nphe_sector1->Write("", TObject::kOverwrite);
      
      hist_sector2->Write("", TObject::kOverwrite);
      hist_nphe_sector2->Write("", TObject::kOverwrite);
      
     
      avrg_nphe_sector2->Write("", TObject::kOverwrite);
            
      hist_sector3->Write("", TObject::kOverwrite);
      hist_nphe_sector3->Write("", TObject::kOverwrite); 
      
 

      avrg_nphe_sector3->Write("", TObject::kOverwrite);     
      
      
      hist_sector4->Write("", TObject::kOverwrite);
      hist_nphe_sector4->Write("", TObject::kOverwrite); 
      
      

      avrg_nphe_sector4->Write("", TObject::kOverwrite);     
         
        
      hist_sector5->Write("", TObject::kOverwrite);
      hist_nphe_sector5->Write("", TObject::kOverwrite); 
      

      avrg_nphe_sector5->Write("", TObject::kOverwrite);     
      
      
      
      
      hist_sector6->Write("", TObject::kOverwrite);  
      hist_nphe_sector6->Write("", TObject::kOverwrite);    
      

      avrg_nphe_sector6->Write("", TObject::kOverwrite);     
      
        
      nphe_sector1->Write("", TObject::kOverwrite); 
      nphe_sector1_after->Write("", TObject::kOverwrite);  
      nphe_sector2->Write("", TObject::kOverwrite); 
      nphe_sector2_after->Write("", TObject::kOverwrite);        
      nphe_sector3->Write("", TObject::kOverwrite); 
      nphe_sector3_after->Write("", TObject::kOverwrite);  
      nphe_sector4->Write("", TObject::kOverwrite); 
      nphe_sector4_after->Write("", TObject::kOverwrite);   
      nphe_sector5->Write("", TObject::kOverwrite); 
      nphe_sector5_after->Write("", TObject::kOverwrite);  
      nphe_sector6->Write("", TObject::kOverwrite); 
      nphe_sector6_after->Write("", TObject::kOverwrite);             
 
    outFile->Write();
    
      nphe_sector1->Delete();
      nphe_sector1_after->Delete();
      nphe_sector2->Delete();
      nphe_sector2_after->Delete();     
      nphe_sector3->Delete();
      nphe_sector3_after->Delete();
      nphe_sector4->Delete();
      nphe_sector4_after->Delete();      
      nphe_sector5->Delete();
      nphe_sector5_after->Delete();
      nphe_sector6->Delete();
      nphe_sector6_after->Delete();     
     
 	outFile->Close();
//	avrg_nphe_sector1->Delete();
	hist_sector1->Delete();
	hist_nphe_sector1->Delete();
//	avrg_nphe_sector1->Delete();
	hist_sector2->Delete();
	hist_nphe_sector2->Delete();
	hist_sector3->Delete();
	hist_nphe_sector3->Delete();
	hist_sector4->Delete();
	hist_nphe_sector4->Delete();
	hist_sector5->Delete();
	hist_nphe_sector5->Delete();
	hist_sector6->Delete();
	hist_nphe_sector6->Delete();
	
	cout << "Q total = " << Qfull << "\n";
	
    };

    
     

















