Float_t fida(double *x, double *par) {
   Float_t s;
   double th_min1;
     th_min1=(9.5+17./(par[0]+0.2));
     double par1,par2;
      par1=0.85+1.1*par[0];
     par2=-62.8-30.*par[0]; 
     s = par[1]*37.3*pow((sin((x[0]-th_min1)*0.01745)),(par1+par2/x[0]+1525./x[0]/x[0]))-par[1]*0;
   return s;
}


void fid_new2() {

//gStyle->SetCanvasBorderSize(50);
//gStyle->SetCanvasBorderMode(1);
gStyle->SetTitleSize(0.12,"t"); 
gStyle->SetStatY(0.9);                
gStyle->SetStatX(0.7);                
gStyle->SetStatW(0.2);                
gStyle->SetStatH(0.2);

gStyle ->SetPalette(1);
gStyle->SetOptStat(0);
//gStyle->SetOptStat("e");
gStyle ->SetStatFontSize(0.09);

gStyle ->SetOptLogz(1);
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



//TFile *MyFile = new TFile("../../event_selection/b_vs_p/test3.root","READ");
TFile *MyFile = new TFile("test_proton_fiduch.root","READ");
//TCanvas *c1 = new TCanvas("Y vs X", "Y vs X", 1000, 800);
//c1->Divide(2,3);

int ii, j;
TCanvas *c_arr[17];
ostringstream qqq2;
ostringstream qqq3;
ostringstream qqq4;
Float_t p,th_min;

for (j=8; j<9; j++){
qqq2.str("");
qqq2 << "c_arr" << j+1;
c_arr[j] = new TCanvas (qqq2.str().c_str(),qqq2.str().c_str(),0,0,500,300);
//c_arr[j]->Divide(2,3);//0.048,0.
c_arr[j]->SetBottomMargin(0.25);
c_arr[j]->SetTopMargin(0.17);
c_arr[j]->SetLeftMargin(0.13);
c_arr[j]->SetRightMargin(0.14);
c_arr[j]->cd();  

TH2D *fr = new TH2D("fr","", 1000,20.,70.,1000,-30.,30.); 
fr->GetYaxis()->SetLabelSize(0.1);
fr->GetZaxis()->SetLabelSize(0.1);
fr->GetYaxis()->SetNdivisions(4);
fr->GetXaxis()->SetNdivisions(6);
fr->GetXaxis()->SetLabelSize(0.1);
fr->GetXaxis()->SetTitle("#theta_{e'} (deg)");
fr->GetYaxis()->SetTitle("#varphi_{e'} (deg)");
fr->GetYaxis()->SetTitleOffset(0.4);
fr->GetZaxis()->SetLabelOffset(-0.005);
fr->GetXaxis()->SetTitleOffset(0.7);
fr->GetXaxis()->SetTitleSize(0.14);
fr->GetYaxis()->SetTitleSize(0.14);
fr->SetMinimum(1.);
fr->SetMaximum(1000.);
qqq3.str("");
qqq3 <<  400+80*j << " < P_{e'}< " << 400+80*(j+1) << " MeV";;
fr->SetTitle(qqq3.str().c_str());
/* 
TPad *null=new TPad("null", "null", 0, 0, 1, 1);
null->SetFillStyle(0);
null->SetFrameFillStyle(0);
null->Draw();
null->cd();
*/
fr->Draw("cont3");

 
 for (ii=1; ii<2; ii++) {
 


std::string histName;
std::stringstream out;
out << "sector" << ii << "/ph_vs_th_" << ii << "pe[" << j << "]";

histName = out.str();
cout << histName << "\n";
TH2F *h1 = (TH2F*)MyFile->Get(histName.c_str());
gStyle->SetTitleFont(32);
h1->SetMinimum(1.);
h1->SetMaximum(1000.);
//h1->SetAxisRange(20., 70.,"X");
 TCutG *cutg = new TCutG("cutg",5);
   cutg->SetPoint(0,20,-29);
   cutg->SetPoint(1,69.7,-29);
   cutg->SetPoint(2,69.7, 29);
   cutg->SetPoint(3,20, 29);
   cutg->SetPoint(4,20,-29);
h1->Draw("colz [cutg] same"); 
h1->GetYaxis()->SetLabelSize(0.1);
h1->GetZaxis()->SetLabelSize(0.1);
h1->GetYaxis()->SetNdivisions(4);
h1->GetXaxis()->SetLabelSize(0.1);
h1->GetXaxis()->SetTitle("#theta (deg)");
h1->GetYaxis()->SetTitle("#varphi (deg)");
h1->GetYaxis()->SetTitleOffset(0.4);
h1->GetZaxis()->SetLabelOffset(-0.005);
h1->GetXaxis()->SetTitleOffset(0.7);
h1->GetXaxis()->SetTitleSize(0.15);
h1->GetYaxis()->SetTitleSize(0.15);
//h1->SetAxisRange(20., 70.,"X");
qqq3.str("");
qqq3 <<  400+80*j << " < P_{e'}< " << 400+80*(j+1) << " MeV";;
h1->SetTitle(qqq3.str().c_str());
p = 0.4+0.08*(j+1) - 0.04;
th_min=(9.5+17./(p+0.2))+1.;
cout << th_min << endl;
 TF1 *f1 = new TF1("f1",fida,th_min,50.,2);
 f1->SetParameter(0,p);
 f1->SetParameter(1,1.);
  f1->SetLineColor(kBlack);
  f1->SetLineWidth(3);
   f1->Draw("same");
   
 TF1 *f2 = new TF1("f2",fida,th_min,50.,2);
 f2->SetParameter(0,p);
 f2->SetParameter(1,-1.); 
 f2->SetLineColor(kBlack);
 f2->SetLineWidth(3);

   f2->Draw("same");
 
    
TLine *line_l = new TLine(th_min+0.1,f2->Eval(th_min+0.1),th_min+0.1,f1->Eval(th_min+0.1));
line_l->SetLineColor(kBlack);
line_l->SetLineWidth(3);
line_l->Draw("same");

TLine *line_r = new TLine(49.9,f2->Eval(50.+0.1),49.9,f1->Eval(50.+0.1));
line_r->SetLineColor(kBlack);
line_r->SetLineWidth(3);
line_r->Draw("same");
 
};
gPad->RedrawAxis();
/*
//c_arr[j]->cd();
//newpad = new TPad("newpad","a transparent pad",0.,0.,1.,1.);
// newpad->SetFillStyle(4000);
// newpad->Draw();
//  newpad->cd();
 TLatex tex;
 qqq4.str("");
qqq4 << 400+80*j << " < P_{el}< " << 400+80*(j+1) << " MeV";
 tex.SetTextSize(0.03);
// tex.DrawLatex(0.4,0.95,qqq4.str().c_str());

 qqq4.str("");
qqq4 << "el_fid_2dim_" << j << ".pdf"; 
c_arr[j]->SaveAs(qqq4.str().c_str());
qqq4.str("");
*/
}; 
};


