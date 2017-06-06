void vertex_empty() {
//gStyle->SetCanvasBorderSize(50);
//gStyle->SetCanvasBorderMode(2);
//gStyle->SetOptStat("e");
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetTitleSize(0.12,"t");
gStyle->SetTitleAlign(33); 
gStyle->SetTitleY(1.02);
/*
gStyle->SetTitleX(0.7);
gStyle->SetStatY(0.9);                
gStyle->SetStatX(0.95);                
gStyle->SetStatW(0.65);                
gStyle->SetStatH(0.1);
gStyle->SetStatFontSize(0.07);
*/
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

TCanvas *c = new TCanvas("c","c",750,450);
//c->Divide(3,2);

c->SetBottomMargin(0.16);

qqq4 << "hist_z_el_" << 1 << "_empty";
MyFile->GetObject(qqq4.str().c_str(),h_z_el_empty[0]);
h_z_el_empty[0]->Rebin(2);
qqq4.str("");

qqq4 << "hist_z_el_" << 1;
MyFile->GetObject(qqq4.str().c_str(),h_z_el[0]);
h_z_el[0]->Rebin(2);
qqq4.str("");

for (i=1; i<6; i++) {



qqq4 << "hist_z_el_" << i+1 << "_empty";
MyFile->GetObject(qqq4.str().c_str(),h_z_el_empty[i]);
h_z_el_empty[i]->Rebin(2);
qqq4.str("");



qqq4 << "hist_z_el_" << i+1;
MyFile->GetObject(qqq4.str().c_str(),h_z_el[i]);
h_z_el[i]->Rebin(2);
qqq4.str("");

qqq4 << "sector " << i+1;
h_z_el[i]->SetTitle(qqq4.str().c_str());
qqq4.str("");


h_z_el[0]->Add(h_z_el[i]); 



h_z_el_empty[0]->Add(h_z_el_empty[i]); 






};


h_z_el_empty[0]->SetAxisRange(-4.,4.,"X");
h_z_el_empty[0]->SetLineColor(kBlack);
h_z_el_empty[0]->SetLineStyle(2);
h_z_el_empty[0]->SetLineWidth(4);
h_z_el_empty[0]->Scale(20.0839/334.603);
qqq4 << "";

h_z_el[0]->SetTitle(qqq4.str().c_str());
qqq4.str("");
h_z_el[0]->SetAxisRange(-4.,4.,"X");
h_z_el[0]->SetAxisRange(-4.,4.,"X");
h_z_el[0]->GetXaxis()->SetTitle("Z_{e'} (cm)");
h_z_el[0]->GetXaxis()->SetLabelSize(0.07);
h_z_el[0]->GetXaxis()->SetTitleSize(0.085);
h_z_el[0]->GetXaxis()->SetNdivisions(5);
h_z_el[0]->GetXaxis()->SetTitleOffset(0.8);
h_z_el[0]->GetYaxis()->SetLabelSize(0.07);
h_z_el[0]->GetYaxis()->SetNdivisions(5);
h_z_el[0]->GetXaxis()->SetNoExponent(1);
h_z_el[0]->Draw();
h_z_el[0]->SetLineColor(kBlack);
h_z_el[0]->SetLineWidth(3);
h_z_el[0]->GetYaxis()->SetTitle("counts/Q_{FC} (#mu C^{-1})");
h_z_el[0]->GetYaxis()->SetTitleSize(0.085);
h_z_el[0]->GetYaxis()->SetTitleOffset(0.5);

gROOT->ForceStyle();
h_z_el_empty[0]->SetLineWidth(2);
h_z_el_empty[0]->Draw("hist same");

TLine *left_cut = new TLine(-2.3,0.,-2.3,3500.);
left_cut->SetLineColor(kBlack);
left_cut->SetLineWidth(3);
left_cut->Draw("same");

TLine *right_cut = new TLine(2.,0.,2.,3500.);
right_cut->SetLineColor(kBlack);
right_cut->SetLineWidth(3);
right_cut->Draw("same");



 TArrow *ar1 = new TArrow(-1.,1000.,-1.5,450.,0.02,"|>");
   ar1->SetLineWidth(3);
   ar1->Draw();
   
 TArrow *ar2 = new TArrow(0.,1000.,0.5,450.,0.02,"|>");
   ar2->SetLineWidth(3);
   ar2->Draw();   
   
  TArrow *ar3 = new TArrow(3.1,1000.,2.6,450.,0.02,"|>");
   ar3->SetLineWidth(3);
   ar3->Draw();   


   TText *windows = new TText(-.5,1450.,"Target");
   windows->SetTextAlign(21);
   windows->SetTextFont(43);
   windows->SetTextSize(30);
   windows->Draw();
   
   TText *windows1 = new TText(-.5,1100.,"windows");
   windows1->SetTextAlign(21);
   windows1->SetTextFont(43);
   windows1->SetTextSize(30);
   windows1->Draw();   
   
   
   TText *windows2 = new TText(3.,1450.,"Upstream");
   windows2->SetTextAlign(21);
   windows2->SetTextFont(43);
   windows2->SetTextSize(30);
   windows2->Draw();
   
   TText *windows3 = new TText(3.,1100.,"foil");
   windows3->SetTextAlign(21);
   windows3->SetTextFont(43);
   windows3->SetTextSize(30);
   windows3->Draw();   
      
   
   

//c->SaveAs("vertex_empty.pdf");


 };   
