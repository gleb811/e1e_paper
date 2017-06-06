void vertex_sim() {
gStyle->SetCanvasBorderSize(50);
gStyle->SetCanvasBorderMode(2);
//gStyle->SetOptStat("e");
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetTitleSize(0.12,"t");
gStyle->SetTitleAlign(33); 
gStyle->SetTitleY(1.02);
gStyle->SetTitleX(0.7);
gStyle->SetStatY(0.9);                
gStyle->SetStatX(0.95);                
gStyle->SetStatW(0.65);                
gStyle->SetStatH(0.1);
gStyle->SetStatFontSize(0.07);
//gDirectory->cd("sector1");
Short_t i;
Short_t j;
ostringstream qqq;
ostringstream qqq1;
ostringstream qqq2;
ostringstream qqq3;
ostringstream qqq4;
ostringstream qqq5;
ostringstream qqq6;
ostringstream qqq7;
TGaxis::SetMaxDigits(3);
//TFile *MyFile = new TFile("../fiduch/test_proton_fiduch.root","READ");
TFile *MyFile = new TFile("test_zel_new.root","READ");
//TFile *MyFile_sim = new TFile("../fiduch/th_vs_p_p_sim/out.root","READ");
TFile *MyFile_sim = new TFile("test_zel_new.root","READ");

TH1 *h_z_el[6];
TH1 *h_z_el_sim[6];
TH1 *h_z_el_sim1[6];
TH1 *h_z_el_empty[6];

TCanvas *c = new TCanvas("c","c",900,450);
c->Divide(3,2);

for (i=0; i<6; i++) {

//cout << i << endl;
c->cd(i+1);
c->cd(i+1)->SetBottomMargin(0.15);
c->cd(i+1)->SetRightMargin(0.02);

qqq4 << "hist_z_el_" << i+1 << "_empty";
MyFile->GetObject(qqq4.str().c_str(),h_z_el_empty[i]);
h_z_el_empty[i]->Rebin(4);
qqq4.str("");
h_z_el_empty[i]->SetAxisRange(-4.,4.,"X");
h_z_el_empty[i]->SetLineColor(2);
h_z_el_empty[i]->Scale(-18.0839/334.603);


qqq4 << "hist_z_el_" << i+1;
MyFile->GetObject(qqq4.str().c_str(),h_z_el[i]);
h_z_el[i]->Rebin(4);
qqq4.str("");
qqq4 << "sector " << i+1;
h_z_el[i]->SetTitle(qqq4.str().c_str());
qqq4.str("");
h_z_el[i]->SetAxisRange(-4.,4.,"X");
h_z_el[i]->SetAxisRange(-4.,4.,"X");
h_z_el[i]->GetXaxis()->SetTitle("Z_{e'} (cm)");
h_z_el[i]->GetXaxis()->SetLabelSize(0.07);
h_z_el[i]->GetXaxis()->SetTitleSize(0.085);
h_z_el[i]->GetXaxis()->SetNdivisions(5);
h_z_el[i]->GetXaxis()->SetTitleOffset(0.8);
h_z_el[i]->GetYaxis()->SetLabelSize(0.07);
h_z_el[i]->GetYaxis()->SetNdivisions(5);
h_z_el[i]->GetXaxis()->SetNoExponent(1);
h_z_el[i]->Add(h_z_el_empty[i]);
h_z_el[i]->SetLineWidth(2);
h_z_el[i]->GetYaxis()->SetTitle("counts");
h_z_el[i]->GetYaxis()->SetTitleSize(0.085);
h_z_el[i]->GetYaxis()->SetTitleOffset(0.5);

h_z_el[i]->Draw(); 


h_z_el_empty[i]->SetLineWidth(2);
//h_z_el_empty[i]->Draw("same"); 

/*
TLine *left_cut = new TLine(-2.3,0.,-2.3,1000.);
left_cut->SetLineColor(kGreen);
left_cut->SetLineWidth(2);
left_cut->Draw("same");

TLine *right_cut = new TLine(2.,0.,2.,1000.);
right_cut->SetLineColor(kGreen);
right_cut->SetLineWidth(2);
right_cut->Draw("same");

*/
qqq4 << "hist_z_el_" << i+1 << "_sim_2";
MyFile_sim->GetObject(qqq4.str().c_str(),h_z_el_sim[i]);
h_z_el_sim[i]->Rebin(4);
qqq4.str("");

h_z_el_sim[i]->SetAxisRange(-4.,4.,"X");
h_z_el_sim[i]->Scale((h_z_el[i]->GetMaximum())/(h_z_el_sim[i]->GetMaximum()));
//h_z_el_sim[i]->Scale(1./(h_z_el_sim[i]->GetMaximum()));
h_z_el_sim[i]->SetLineColor(2);
h_z_el_sim[i]->SetTitle("");
h_z_el_sim[i]->SetLineWidth(2);
h_z_el_sim[i]->Draw("same"); 



qqq4 << "hist_z_el_" << i+1 << "_sim_1";
MyFile_sim->GetObject(qqq4.str().c_str(),h_z_el_sim1[i]);
h_z_el_sim1[i]->Rebin(4);
qqq4.str("");

h_z_el_sim1[i]->SetAxisRange(-4.,4.,"X");
h_z_el_sim1[i]->Scale((h_z_el[i]->GetMaximum())/(h_z_el_sim1[i]->GetMaximum()));
//h_z_el_sim[i]->Scale(1./(h_z_el_sim[i]->GetMaximum()));
h_z_el_sim1[i]->SetLineColor(3);
h_z_el_sim1[i]->SetTitle("");
h_z_el_sim1[i]->SetLineWidth(2);
h_z_el_sim1[i]->Draw("same"); 


};


/*c->cd(1)->SetRightMargin(0.01);
c->cd(1)->SetLeftMargin(0.15);
c->cd(1)->SetBottomMargin(0.2);
c->cd(1)->SetTopMargin(0.);
c->cd(2)->SetRightMargin(0.01);
c->cd(2)->SetLeftMargin(0.15);
c->cd(2)->SetBottomMargin(0.2);
c->cd(2)->SetTopMargin(0.1);
c->cd(3)->SetRightMargin(0.01);
c->cd(3)->SetLeftMargin(0.15);
c->cd(3)->SetBottomMargin(0.2);
c->cd(3)->SetTopMargin(0.1);
*/


//qqq4 << "sector6/ph_vs_th_6pe[" << j << "]";
//cout << qqq4.str().c_str() <<"\n";
//MyFile->GetObject("ltime",ltime);
//MyFile->GetObject("n_incl",n_incl);
//MyFile->GetObject("n_elast",n_elast);

/*
c->cd(1);
//ltime->Rebin(5.);
//ltime->Scale(1./5.);
ltime->SetTitle("");
ltime->SetMaximum(0.98);
ltime->SetMinimum(0.8);
ltime->GetXaxis()->SetLabelSize(0.09);
ltime->GetXaxis()->SetNoExponent(1);
ltime->GetXaxis()->SetTitle("block number");
ltime->GetXaxis()->SetTitleSize(0.12);
ltime->GetXaxis()->SetTitleOffset(0.66);
ltime->GetYaxis()->SetLabelSize(0.09);
ltime->GetYaxis()->SetNdivisions(5);
ltime->GetYaxis()->SetTitle("live time");
ltime->GetYaxis()->SetTitleSize(0.12);
ltime->GetYaxis()->SetTitleOffset(0.5);
ltime->Draw();

TLine *ltime_1 = new TLine(0.,0.925,35000.,0.925);
ltime_1->SetLineColor(kRed);
ltime_1->SetLineWidth(2);
ltime_1->Draw("same");

TLine *ltime_2 = new TLine(0.,0.89,35000.,0.89);
ltime_2->SetLineColor(kRed);
ltime_2->SetLineWidth(2);
ltime_2->Draw("same");

c->cd(2);
n_incl->SetTitle("");
n_incl->SetMaximum(70000.);
n_incl->SetMinimum(40000.);
n_incl->GetXaxis()->SetLabelSize(0.09);
n_incl->GetXaxis()->SetNoExponent(1);
n_incl->GetXaxis()->SetTitle("block number");
n_incl->GetXaxis()->SetTitleSize(0.12);
n_incl->GetXaxis()->SetTitleOffset(0.66);
n_incl->GetYaxis()->SetLabelSize(0.09);
n_incl->GetYaxis()->SetNdivisions(5);
n_incl->GetYaxis()->SetTitle("#frac{N_{inclusive events}}{Q_{block}} (#frac{1}{#muC})");
n_incl->GetYaxis()->SetTitleSize(0.12);
n_incl->GetYaxis()->SetTitleOffset(0.5);
n_incl->Draw();

TLine *incl_1 = new TLine(0.,60000.,35000.,60000.);
incl_1->SetLineColor(kRed);
incl_1->SetLineWidth(2);
incl_1->Draw("same");

TLine *incl_2 = new TLine(0.,48000.,35000.,48000.);
incl_2->SetLineColor(kRed);
incl_2->SetLineWidth(2);
incl_2->Draw("same");



c->cd(3);
n_elast->SetTitle("");
n_elast->SetMaximum(22000.);
n_elast->SetMinimum(15000.);
n_elast->GetXaxis()->SetLabelSize(0.09);
n_elast->GetXaxis()->SetNoExponent(1);
n_elast->GetXaxis()->SetTitle("block number");
n_elast->GetXaxis()->SetTitleSize(0.12);
n_elast->GetXaxis()->SetTitleOffset(0.66);
n_elast->GetYaxis()->SetLabelSize(0.09);
n_elast->GetYaxis()->SetNdivisions(5);
n_elast->GetYaxis()->SetTitle("#frac{N_{elastic events}}{Q_{block}} (#frac{1}{#muC})");
n_elast->GetYaxis()->SetTitleSize(0.12);
n_elast->GetYaxis()->SetTitleOffset(0.5);
TGaxis::SetMaxDigits(3);
n_elast->Draw();

TLine *elast_1 = new TLine(0.,20500.,35000.,20500.);
elast_1->SetLineColor(kRed);
elast_1->SetLineWidth(2);
elast_1->Draw("same");

TLine *elast_2 = new TLine(0.,16400.,35000.,16400.);
elast_2->SetLineColor(kRed);
elast_2->SetLineWidth(2);
elast_2->Draw("same");
*/
c->SaveAs("vertex_sim.pdf");


 };   
