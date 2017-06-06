

void ec() {
gROOT->LoadMacro("DynamicExec.C");
gStyle->SetPalette(kGreyScale);


UInt_t Number = 50;
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









TFile *MyFile = new TFile("ec_data.root","READ");
MyFile->ls();
gStyle->SetOptLogz(1);
//gStyle->SetPalette(1);
gStyle->SetOptStat(0);
gStyle->SetOptFit(1);
gStyle->SetTitleSize(0.09,"t"); 
TCanvas *c1 = new TCanvas("rose512", "DATA", 550, 500);
//c2 = new TCanvas("rose514", "Lower age", 800, 800);
//c2->SetBorderMode(0);
//c3 = new TCanvas("rose515", "Upper age", 800, 800);
//c3->SetBorderMode(0);
c1->Divide(1,2);
//c2->Divide(2,3);
//c3->Divide(2,3);
int ii;
std::ostringstream title_2d;

for (ii=1; ii<2; ii++)
{
std::string histName;
std::stringstream out;
out << "h" << ii*100+ii*10+ii;
histName = out.str();
cout << histName << "\n";
TH2F *h1 = (TH2F*)MyFile->Get(histName.c_str());
c1->cd(ii);
c1->cd(ii)->SetBottomMargin(0.14);
c1->cd(ii)->SetTopMargin(0.1);
c1->cd(ii)->SetLeftMargin(0.13);
title_2d.str("");
title_2d << "sectror " << ii << " (data)";
//TH2F *h2;
//h2 = eprst(h1);
h1->SetTitle(title_2d.str().c_str());
//h1->SetTitleSize(0.5);

h1->GetYaxis()->SetTitle("E_{tot}/P_{e'} (dimensionless)");
h1->GetYaxis()->SetTitleSize(0.08); 
h1->GetYaxis()->SetTitleOffset(0.65);
h1->GetYaxis()->SetLabelSize(0.07);
h1->GetXaxis()->SetTitleSize(0.08); 
h1->GetXaxis()->SetLabelSize(0.07); 
h1->GetXaxis()->SetTitleOffset(0.8);
h1->GetXaxis()->SetTitle("P_{e'} (GeV)"); 
h1->GetXaxis()->SetTitleFont(22);
h1->GetYaxis()->SetTitleFont(22);
h1->SetTitleFont(22);
h1->GetXaxis()->SetNdivisions(414);
h1->GetYaxis()->SetNdivisions(414);
h1->GetXaxis()->SetLabelFont(22);
h1->GetYaxis()->SetLabelFont(22);

 TCutG *cutg = new TCutG("cutg",5);
   cutg->SetPoint(0,0.,0.);
   cutg->SetPoint(1,1.98,0.);
   cutg->SetPoint(2,1.98,0.5);
   cutg->SetPoint(3,0., 0.5);
   cutg->SetPoint(4,0.,0.);

h1->Draw("colz [cutg]");
   TLine *line_up1 = new TLine(0.,0.5,2.,0.5);
   line_up1->Draw("same");

//h1->SetTitleSize(0.5);
//c1->AddExec("dynamic","DynamicExec()");

int a = h1->GetNbinsX();
float xnew[65];
float ynew[65];
float ynew1[65];
float err[65];
cout << a << "\n";

TH1F *hhh[65];
short i;
for (i=25; i<a-10; i++)
{
hhh[i-25] = (TH1F*) h1->ProjectionY(" ",i,i);
xnew[i-25] = h1->GetXaxis()->GetBinCenter(i);
hhh[i-25]->Fit("gaus","0QW");
TF1 *fit =hhh[i-25]->GetFunction("gaus");
Double_t p1 = fit->GetParameter(1);
Double_t p2 = fit->GetParameter(2);
//cout << xnew[i-21] << "\n"; 
c1->Update();
ynew[i-25]=p1-3.*p2;
ynew1[i-25]=p1+3.*p2;
//cin.get();
};

std::ostringstream title;
title << "Sector " << ii;

//c2->cd(ii);
TGraphErrors *Tg = new TGraphErrors(65,xnew,ynew,err,err);
 Tg->SetMarkerStyle(kFullCircle);
// Tg->Draw("AP");
//c2->Update();
Tg->Fit("pol2","W");
Tg->SetTitle(title.str().c_str());
//c2->Update();
TF1 *fitpol =Tg->GetFunction("pol2");
Double_t p11 = fitpol->GetParameter(0);
Double_t p12 = fitpol->GetParameter(1);
Double_t p13 = fitpol->GetParameter(2);
TF1 *lower = new TF1("lower","pol2",0.461,2.0);
lower->SetParameters(p11,p12,p13);
c1->cd(ii);
lower->SetLineColor(1);
lower->SetLineWidth(3);
lower->Draw("same");




//c3->cd(ii);



TGraphErrors *Tg1 = new TGraphErrors(65,xnew,ynew1,err,err);
 Tg1->SetMarkerStyle(kFullCircle);
// Tg1->Draw("AP");
//c3->Update();
Tg1->Fit("pol2","W"); 
Tg1->SetTitle(title.str().c_str());
//c3->Update();
TF1 *fitpo2 =Tg1->GetFunction("pol2");
Double_t p21 = fitpo2->GetParameter(0);
Double_t p22 = fitpo2->GetParameter(1);
Double_t p23 = fitpo2->GetParameter(2);
TF1 *upper = new TF1("upper","pol2",0.461,2.0);
upper->SetParameters(p21,p22,p23);
c1->cd(ii);
upper->SetLineColor(1);
upper->SetLineWidth(3);
upper->Draw("same");

TLine *line = new TLine(0.461,lower->Eval(0.461),0.461,upper->Eval(0.461));
line->SetLineColor(1);
line->SetLineWidth(3);
line->Draw("same");

c1->Update();
//c2->Update();
//c3->Update();

TFile *MyFile1 = new TFile("ec_sim_new.root","READ");


out.str("");
out << "h" << ii*100+ii*10+ii;
TH2F *h11 = (TH2F*)MyFile1->Get(out.str().c_str());
c1->cd(ii+1);
c1->cd(ii+1)->SetBottomMargin(0.14);
c1->cd(ii+1)->SetTopMargin(0.1);
c1->cd(ii+1)->SetLeftMargin(0.13);

title_2d.str("");
title_2d << "sectror " << ii << "  (Monte Carlo)";
//TH2F *h2;
//h2 = eprst(h1);


//h1->SetTitleSize(0.5);

h11->SetTitle(title_2d.str().c_str());
h11->GetYaxis()->SetTitle("E_{tot}/P_{e'} (dimensionless)");
h11->GetYaxis()->SetTitleSize(0.08); 
h11->GetYaxis()->SetTitleOffset(0.65);
h11->GetYaxis()->SetLabelSize(0.07);
h11->GetXaxis()->SetTitleSize(0.08); 
h11->GetXaxis()->SetLabelSize(0.07); 
h11->GetXaxis()->SetTitleOffset(0.8);
h11->GetXaxis()->SetTitle("P_{e'} (GeV)"); 
h11->SetTitleFont(h1->GetTitleFont());
h11->GetXaxis()->SetTitleFont(22);
h11->GetYaxis()->SetTitleFont(22);
h11->GetXaxis()->SetNdivisions(414);
h11->GetYaxis()->SetNdivisions(414);
h11->GetXaxis()->SetLabelFont(22);
h11->GetYaxis()->SetLabelFont(22);

h11->GetXaxis()->SetRange(h11->GetXaxis()->FindBin(h1->GetXaxis()->GetXmin()),h11->GetXaxis()->FindBin(h1->GetXaxis()->GetXmax())-1);
h11->GetYaxis()->SetRange(h11->GetYaxis()->FindBin(h1->GetYaxis()->GetXmin()),h11->GetYaxis()->FindBin(h1->GetYaxis()->GetXmax()));
h11->Draw("colz [cutg]");
//h1->SetTitleSize(0.5);
//c1->AddExec("dynamic","DynamicExec()");
TF1 *upper1 = new TF1("upper1","pol2",0.461,1.5);
upper1->SetParameters(0.345,-0.05557,0.01303);
upper1->SetLineColor(1);
upper1->SetLineWidth(3);
upper1->Draw("same");
TF1 *lower1 = new TF1("lower1","pol2",0.461,1.5);
lower1->SetParameters(0.07099,0.117,-0.03283);
lower1->SetLineColor(1);
lower1->SetLineWidth(3);
lower1->Draw("same");
TLine *line1 = new TLine(0.461,lower1->Eval(0.461),0.461,upper1->Eval(0.461));
line1->SetLineColor(1);
line1->SetLineWidth(3);
line1->Draw("same");

//h11->Draw("sameaxis"); 
//gPad->RedrawAxis("g");

   TGaxis *axis = new TGaxis();
   Double_t xmin,ymin,xmax,ymax;
   gPad->GetRangeAxis(xmin,ymin,xmax,ymax);
//  redraw x axis
   TLine *line_up2 = new TLine(0.,0.5,2.,0.5);
   line_up2->Draw("same");



c1->Update();

};
};

TH2F *eprst(TH2F *h) {
  TH2F *h2 = (TH2F*)h->Clone();
  Int_t nby = h2->GetYaxis()->GetNbins();
  for (Int_t bx = 1; bx < (h2->GetXaxis()->GetNbins() + 2); bx++) {
    h2->GetXaxis()->SetRange(bx,bx);
    TH1 *h1 = h2->ProjectionY();
    Double_t max = h1->GetMaximum();
    for (Int_t by = 1; by < nby; by++) {
      Double_t bc = h2->GetBinContent(bx,by);
      if (bc > 0) {
        Int_t global_bin = h2->GetBin(bx,by);
        h2->SetBinContent(global_bin,bc/max);
      };
    };
  };
 h2->GetXaxis()->SetRange(0,h2->GetXaxis()->GetNbins()+1);
  return h2;
};



















