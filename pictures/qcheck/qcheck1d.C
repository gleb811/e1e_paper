void qcheck1d() {
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

TFile *MyFile = new TFile("../fiducial_cuts/test_proton_fiduch.root","READ");

TH1 *ltime;
TH1 *n_incl;
TH1 *n_elast;

TCanvas *c = new TCanvas("c","c",300,600);
c->Divide(1,3);

c->cd(1)->SetRightMargin(0.13);
c->cd(1)->SetBottomMargin(0.4);
c->cd(1)->SetTopMargin(0.11);
c->cd(1)->SetLeftMargin(0.19);
c->cd(2)->SetRightMargin(0.13);
c->cd(2)->SetBottomMargin(0.21); 
c->cd(2)->SetTopMargin(0.1);
c->cd(2)->SetLeftMargin(0.12);
c->cd(3)->SetRightMargin(0.13);
c->cd(3)->SetBottomMargin(0.4);
c->cd(3)->SetTopMargin(0.1);
c->cd(3)->SetLeftMargin(0.12);
//qqq4 << "sector6/ph_vs_th_6pe[" << j << "]";
//cout << qqq4.str().c_str() <<"\n";
MyFile->GetObject("ltime_1d",ltime);
MyFile->GetObject("n_incl_1d",n_incl);
MyFile->GetObject("n_elast_1d",n_elast);

c->cd(1);
//ltime->Rebin(5.);
//ltime->Scale(1./5.);
ltime->SetTitle("");
ltime->SetMaximum(4800.);
//ltime->SetMinimum(0.8);
ltime->GetXaxis()->SetLabelSize(0.1);
ltime->GetXaxis()->SetNdivisions(5);
//ltime->GetXaxis()->SetNoExponent(1);
ltime->GetXaxis()->SetTitle("live time");
ltime->GetXaxis()->SetTitleSize(0.12);
ltime->GetXaxis()->SetTitleOffset(0.8);
ltime->SetAxisRange(0.85,0.98,"X");
ltime->GetYaxis()->SetLabelSize(0.1);
ltime->GetYaxis()->SetNdivisions(5);
ltime->GetYaxis()->SetTitle("number of blocks");
ltime->GetYaxis()->SetTitleSize(0.12);
ltime->GetYaxis()->SetTitleOffset(0.55);
ltime->SetLineColor(kBlack);
ltime->SetLineWidth(3);
ltime->Draw();

TLine *ltime_1 = new TLine(0.925,0.,0.925,4800.);
ltime_1->SetLineColor(kBlack);
ltime_1->SetLineWidth(3);
ltime_1->Draw("same");

TLine *ltime_2 = new TLine(0.89,0.,0.89,4800.);
ltime_2->SetLineColor(kBlack);
ltime_2->SetLineWidth(3);
ltime_2->Draw("same");

c->cd(2);
c->cd(2)->SetRightMargin(0.13);
c->cd(2)->SetRightMargin(0.13);
c->cd(2)->SetBottomMargin(0.4);
c->cd(2)->SetTopMargin(0.11);
c->cd(2)->SetLeftMargin(0.19);

n_incl->SetTitle("");
//n_incl->SetMaximum(70000.);
//n_incl->SetMinimum(40000.);
n_incl->GetXaxis()->SetLabelSize(0.1);
n_incl->GetXaxis()->SetNdivisions(5);
n_incl->GetXaxis()->SetNoExponent(0);
n_incl->SetAxisRange(40000.,70000.,"X");
n_incl->GetXaxis()->SetTitle("#frac{N_{inclusive events}}{Q_{block}} (#frac{1}{#muC})");
n_incl->GetXaxis()->SetTitleSize(0.12);
n_incl->GetXaxis()->SetTitleOffset(1.3);
n_incl->GetYaxis()->SetLabelSize(0.1);
n_incl->GetYaxis()->SetNdivisions(5);
n_incl->GetYaxis()->SetTitle("number of blocks");
n_incl->GetYaxis()->SetTitleSize(0.12);
n_incl->GetYaxis()->SetTitleOffset(0.55);
n_incl->SetLineColor(kBlack);
n_incl->SetLineWidth(3);
//n_incl->GetXaxis()->SetWmax(1.2);
n_incl->Draw();


TLine *incl_1 = new TLine(60000.,0.,60000.,4000.);
incl_1->SetLineColor(kBlack);
incl_1->SetLineWidth(3);
incl_1->Draw("same");

TLine *incl_2 = new TLine(48000.,0.,48000.,4000.);
incl_2->SetLineColor(kBlack);
incl_2->SetLineWidth(3);
incl_2->Draw("same");



c->cd(3);
c->cd(3)->SetRightMargin(0.13);
c->cd(3)->SetRightMargin(0.13);
c->cd(3)->SetBottomMargin(0.4);
c->cd(3)->SetTopMargin(0.11);
c->cd(3)->SetLeftMargin(0.19);

n_elast->SetTitle("");
n_elast->SetMaximum(3000.);
//n_elast->SetMinimum(15000.);
n_elast->GetXaxis()->SetLabelSize(0.1);
n_elast->GetXaxis()->SetNdivisions(5);
n_elast->GetXaxis()->SetNoExponent(0);
n_elast->SetAxisRange(12000.,25000.,"X");
n_elast->GetXaxis()->SetTitle("#frac{N_{elastic events}}{Q_{block}} (#frac{1}{#muC})");
n_elast->GetXaxis()->SetTitleSize(0.12);
n_elast->GetXaxis()->SetTitleOffset(1.3);
n_elast->GetYaxis()->SetLabelSize(0.1);
n_elast->GetYaxis()->SetNdivisions(5);
n_elast->GetYaxis()->SetTitle("number of blocks");
n_elast->GetYaxis()->SetTitleSize(0.12);
n_elast->GetYaxis()->SetTitleOffset(0.55);
n_elast->SetLineColor(kBlack);
n_elast->SetLineWidth(3);
TGaxis::SetMaxDigits(3);
n_elast->Draw();


TLine *elast_1 = new TLine(20500.,0.,20500.,3000.);
elast_1->SetLineColor(kBlack);
elast_1->SetLineWidth(3);
elast_1->Draw("same");

TLine *elast_2 = new TLine(16400.,0.,16400.,3000.);
elast_2->SetLineColor(kBlack);
elast_2->SetLineWidth(3);
elast_2->Draw("same");

c->SaveAs("qcheck1d.pdf");
 };   
