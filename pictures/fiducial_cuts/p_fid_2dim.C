double ap(double *x, double *par) {
   double s;
   s =  24.*(1-exp(-1.*0.08*(x[0]-9.)));
       return s;
};

double bp(double *x, double *par) {
   double w; 
   w = -1.*25.*(1-exp(-1.*0.1*(x[0]-10.)));
          return w;
};





void p_fid_2dim () {

gStyle ->SetOptLogz(1);
gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.1,"t");
gStyle->SetStatY(1.);                
gStyle->SetStatX(0.9);                
gStyle->SetStatW(0.26);                
gStyle->SetStatH(0.11);
gStyle->SetStatFontSize(0.15);


UInt_t Number = 5;
  Double_t Red[Number];
  Double_t Green[Number];
  Double_t Blue[Number];
  Double_t Stops[Number];

  double dcol = 1/double(Number);
  double grey = 1;
  
    for (int j = 0; j < Number; j++){  

    // ...... Define color with RGB equal to : gray, gray, gray .......

    Stops[j]=double(j)/double(Number-1);
    Red[j]=grey;
    Blue[j]=grey;
    Green[j]=grey;

    grey = grey - dcol;
  }

  Int_t nb=100;

  TColor::CreateGradientColorTable(Number,Stops,Red,Green,Blue,nb);





TCanvas *c_arr[5];
TH2 *hist;
Short_t i, j;
ostringstream qqq;
ostringstream qqq1;
ostringstream qqq2;
ostringstream qqq3;
ostringstream qqq7;

TFile *MyFile = new TFile("test_proton_fiduch.root","READ");

for (i=2; i<3; i++) {
qqq1.str("");
qqq1 << "c_arr" << i+1;
c_arr[i] = new TCanvas (qqq1.str().c_str(),qqq1.str().c_str(),0,0,500,300);
c_arr[i]->Divide(2,3);//0.048,0.

for (j=0; j<1; j++) {
c_arr[i]->SetBottomMargin(0.25);
c_arr[i]->SetTopMargin(0.17);
c_arr[i]->SetLeftMargin(0.14);
c_arr[i]->SetRightMargin(0.13);
c_arr[i]->cd();

TH2D *fr = new TH2D("fr","", 1000,0.,110.,1000,-30.,30.); 
fr->GetYaxis()->SetLabelSize(0.1);
fr->GetZaxis()->SetLabelSize(0.1);
fr->GetYaxis()->SetNdivisions(4);
fr->GetXaxis()->SetNdivisions(6);
fr->GetXaxis()->SetLabelSize(0.1);
fr->GetXaxis()->SetTitle("#theta_{p'} (deg)");
fr->GetYaxis()->SetTitle("#varphi_{p'} (deg)");
fr->GetYaxis()->SetTitleOffset(0.4);
fr->GetZaxis()->SetLabelOffset(-0.005);
fr->GetXaxis()->SetTitleOffset(0.7);
fr->GetXaxis()->SetTitleSize(0.14);
fr->GetYaxis()->SetTitleSize(0.14);
fr->SetMinimum(1.);
fr->SetMaximum(1000.);
qqq3.str("");
qqq3 << 200+200*i << " < P_{p'}< " << 200+200*(i+1) << " MeV";
fr->SetTitle(qqq3.str().c_str());
/* 
TPad *null=new TPad("null", "null", 0, 0, 1, 1);
null->SetFillStyle(0);
null->SetFrameFillStyle(0);
null->Draw();
null->cd();
*/
fr->Draw("cont3");

qqq.str("");
//qqq << "sector1_p_fid/ph_vs_th_p_1_w";
qqq << "sector" << j+1 << "_p_fid/ph_th_p_" << j+1 << "["  << i << "]";
cout << qqq.str().c_str() << "\n";
MyFile->GetObject(qqq.str().c_str(),hist); 
gStyle->SetTitleFont(32);
hist->SetMinimum(1.);
hist->SetMaximum(1000.);

 TCutG *cutg = new TCutG("cutg",5);
   cutg->SetPoint(0,0.,-29.5);
   cutg->SetPoint(1,109.,-29.5);
   cutg->SetPoint(2,109.,29.5);
   cutg->SetPoint(3,0., 29.5);
   cutg->SetPoint(4,0.,-29.5);
hist->Draw("colz [cutg] same");
hist->GetYaxis()->SetLabelSize(0.1);
hist->GetZaxis()->SetLabelSize(0.1);
hist->GetYaxis()->SetNdivisions(4);
hist->GetXaxis()->SetLabelSize(0.1);
hist->GetXaxis()->SetTitle("#theta (deg)");
hist->GetYaxis()->SetTitle("#varphi (deg)");
hist->GetYaxis()->SetTitleOffset(0.4);
hist->GetZaxis()->SetLabelOffset(-0.005);
hist->GetXaxis()->SetTitleOffset(0.7);
hist->GetXaxis()->SetTitleSize(0.15);
hist->GetYaxis()->SetTitleSize(0.15);
qqq3.str("");
qqq3 << 200+200*i << " < P_{p'}< " << 200+200*(i+1) << " MeV";
//qqq3 << 0.3*(60+6.666666*i)<< "< theta < " << 0.3*(60+6.66666*(i+1));
hist->SetTitle(qqq3.str().c_str());
TF1 *f1 = new TF1("f1",ap,9.,100.,1);
 //f1->SetParameter(0,p);
//f1->SetParameter(1,1.);
  f1->SetLineColor(kBlack);
  f1->SetLineWidth(3);
   f1->Draw("same");
   
TF1 *f2 = new TF1("f2",bp,9.,100.,1);
 //f1->SetParameter(0,p);
//f1->SetParameter(1,1.);
  f2->SetLineColor(kBlack);
  f2->SetLineWidth(3);
   f2->Draw("same");   
   
   
   
 }; 
 
 /*
  c_arr[i]->cd();
newpad = new TPad("newpad","a transparent pad",0.,0.,1.,1.);
 newpad->SetFillStyle(4000);
 newpad->Draw();
  newpad->cd();
 TLatex tex;
 qqq2.str("");
qqq2 << 200+200*i << " < P_{p}< " << 200+200*(i+1) << " MeV";
 tex.SetTextSize(0.04);
// tex.DrawLatex(0.45,0.96,qqq2.str().c_str());
 
 qqq7 << "p_fid_2dim_" << i << ".pdf"; 
c_arr[i]->SaveAs(qqq7.str().c_str());
qqq7.str("");
*/
/*
qqq7 << "p_fid_2dim_" << i << ".eps"; 
c_arr[i]->SaveAs(qqq7.str().c_str());
qqq7.str("");
*/
};
};

