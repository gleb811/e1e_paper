


double testfunc2(double *x, double *par) {
   double s;
   s = sqrt(1. - pow((x[0]-45.2)/34.5,2))*21.;
   return s;
}

double testfunc3(double *x, double *par) {
   double s;
   s = -1.*sqrt(1. - pow((x[0]-45.2)/34.5,2))*21.;
   return s;
}

double testfunc4(double *x, double *par) {
   double s;
   s = sqrt(1. - pow((x[0]-45.7)/1.75,2))*21.;
   return s;
}

double testfunc5(double *x, double *par) {
   double s;
   s = -1.*sqrt(1. - pow((x[0]-45.7)/1.75,2))*21.;
   return s;
}





void ratio_qqq() {

TFile *MyFile = new TFile("new_ratio.root","READ");
MyFile->ls();

TCanvas *c1 = new TCanvas("rose512", "DATA", 450, 300);

gStyle->SetOptLogz(0);
 Int_t palette[1];
   palette[0] = 1;
/*   palette[1] = 20;
   palette[2] = 23;
   palette[3] = 30;
   palette[4] = 32;
   */
//   gStyle->SetPalette(5,palette);
gStyle->SetPalette(1,palette);
//gStyle->SetPalette(0);
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
gStyle->SetTitleSize(0.09,"t");
//c1->Divide(2,3);

std::stringstream out;
std::ostringstream title_2d;


for (Int_t ii=1; ii<2; ii++)
{
out.str("");
out << "h30" << ii;

TH2F *h1 = (TH2F*)MyFile->Get(out.str().c_str());

h1->SetMinimum(0.8);
//h1->SetMaximum(0.99);

//c1->cd(ii);
c1->SetBottomMargin(0.14);
c1->SetTopMargin(0.1);
c1->SetLeftMargin(0.15);
c1->SetRightMargin(0.05);
title_2d.str("");
title_2d << "sectror " << ii;
h1->SetTitle(title_2d.str().c_str());
//h1->SetTitleSize(0.5);
h1->GetYaxis()->SetTitle("#varphi_{cc} (deg)");
h1->GetYaxis()->SetTitleSize(0.08); 
h1->GetYaxis()->SetTitleOffset(0.65);
h1->GetYaxis()->SetLabelSize(0.07);
h1->GetXaxis()->SetTitleSize(0.08); 
h1->GetXaxis()->SetLabelSize(0.07); 
h1->GetZaxis()->SetLabelSize(0.07);
h1->GetXaxis()->SetTitleOffset(0.8);
h1->GetXaxis()->SetTitle("#theta_{cc} (deg)"); 
h1->GetXaxis()->SetRange(50,180);
h1->Draw("col");
};

//   TF1 *f1 = new TF1("f1",testfunc,6.9,45,0);
//   f1->Draw("same");

//   TF1 *f2 = new TF1("f2",testfunc1,6.9,45,0);
//   f2->Draw("same");
   
   TF1 *f3 = new TF1("f3",testfunc2,10.5262,45.7,0);
   f3->SetLineColor(kBlack);
   f3->SetLineWidth(4);
   f3->Draw("same");  
   
   TF1 *f4 = new TF1("f4",testfunc3,10.5262,45.7,0);
   f4->SetLineColor(kBlack);
   f4->SetLineWidth(4);   
   f4->Draw("same");    
   
   TF1 *f5 = new TF1("f5",testfunc4,43.941708543,45.6,0);
   f5->SetLineColor(kBlack);
   f5->SetLineWidth(4);   
   f5->Draw("same");        

   TF1 *f6 = new TF1("f6",testfunc5,43.941708543,45.6,0);
   f6->SetLineColor(kBlack);
   f6->SetLineWidth(4);     
   f6->Draw("same");   



}
