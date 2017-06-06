#include "poissonf.C"

void ph_el(){

//gROOT->LoadMacro("poissonf.C");

Bool_t left_right = 0;

gStyle->SetOptStat(0);
gStyle ->SetOptFit(0);
gStyle->SetStatY(0.9);                
gStyle->SetStatX(0.995);                

gStyle->SetStatW(0.25);                
gStyle->SetStatH(0.1);

gStyle->SetStatFontSize(0.15);
gStyle->SetTitleSize(0.1,"t");

gStyle->SetTitleAlign(33); 
gStyle->SetTitleX(.65);
gStyle->SetTitleY(1.02);

TCanvas *c_arr[6];
TCanvas *c1_arr[6];
Short_t i, j;
TH1 *hist,*hist1;
Double_t p0,p1,p2,p3,p4,p5;

Double_t fract[3][6][18];
/*Double_t ledg[6][15] =[ [22.,24.,25.,26.,28.,30.,30.,32.,26.,30.,25.,27.,36.,29.,28.]
[25.,16.,17.,22.,23.,21.,22.,24.,21.,30.,30.,30.,31.,29.,29.]
[20.,21.,24.,23.,19.,23.,31.,28.,27.,28.,27.,30.,32.,37.,37.]
[24.,20.,24.,19.,28.,26.,24.,27.,25.,29.,22.,23.,34.,37.,27.]
[22.,27.,25.,28.,29.,23.,22.,28.,35.,29.,51.,29.,39.,44.,36.]
[28.,28.,22.,22.,28.,70.,70.,41.,27.,34.,34.,29.,37.,38.,36.]];*/

Double_t ledg1[15]= {23.,24.,25.,23.,25.,27.,25.,34.,26.,30.,28.,27.,34.,31.,28.};
Double_t ledg2[15]= {30.,25.,25.,22.,30.,21.,22.,24.,23.,25.,31.,27.,27.,29.,29.};
Double_t ledg3[15]= {23.,25.,24.,23.,23.,23.,31.,28.,27.,28.,31.,32.,32.,32.,30.};
Double_t ledg4[15]= {30.,20.,24.,1.,25.,26.,24.,27.,25.,27.,22.,20.,32.,39.,30.};
Double_t ledg5[15]= {31.,30.,25.,28.,27.,23.,22.,26.,31.,26.,47.,25.,35.,42.,36.};
Double_t ledg6[15]= {28.,28.,22.,28.,27.,26.,26.,28.,26.,29.,28.,27.,31.,36.,34.};


Double_t redg1[15]= {22.,20.,24.,25.,22.,24.,26.,26.,32.,26.,23.,30.,28.,25.,25.};
Double_t redg2[15]= {28.,20.,22.,30.,25.,23.,30.,33.,23.,24.,27.,35.,30.,33.,36.};
Double_t redg3[15]= {23.,26.,28.,23.,23.,25.,24.,27.,25.,27.,26.,30.,33.,29.,36.};
Double_t redg4[15]= {1.,1.,1.,1.,20.,20.,22.,26.,1.,27.,20.,23.,27.,25.,25.};
Double_t redg5[15]= {27.,26.,25.,51.,38.,32.,29.,30.,31.,36.,35.,27.,38.,41.,31.};
Double_t redg6[15]= {29.,21.,26.,33.,29.,26.,29.,32.,36.,33.,26.,32.,50.,35.,36.};

TFile *MyFile = new TFile("data_no_nphe.root","READ");

ostringstream qqq;
for (i=3; i<4; i++) {

qqq.str("");
qqq << "c_arr" << i+1;
c_arr[i] = new TCanvas (qqq.str().c_str(),qqq.str().c_str(),0,0,400,300);
//c_arr[i]->Divide(3,5);//0.048,0.
for (j=9; j<10; j++) {
//c_arr[i]->cd(j-1);
c_arr[i]->SetBottomMargin(0.2);
c_arr[i]->SetTopMargin(0.11);
c_arr[i]->SetLeftMargin(0.12);
c_arr[i]->SetRightMargin(0.1);
qqq.str("");
//qqq << "sector1_p_fid/ph_vs_th_p_1_w";
if (left_right) {
qqq << "ph_el_sector_" << i+1 << "/left/sector_" << i+1 << "_left_"  << j+1;
} else {
qqq << "ph_el_sector_" << i+1 << "/right/sector_" << i+1 << "_right_"  << j+1;
};
//cout << qqq.str().c_str() << "\n";
MyFile->GetObject(qqq.str().c_str(),hist); 
qqq.str("");
qqq << "PMT " << j+1;
hist->SetTitle(qqq.str().c_str());
hist->GetXaxis()->SetRange(1,300);
hist->GetYaxis()->SetNoExponent(0);
hist->SetLineColor(kBlack);
hist->SetLineWidth(2);
hist->Draw();
hist->GetYaxis()->SetLabelSize(0.1);
hist->GetYaxis()->SetTitleOffset(1.1);
hist->GetYaxis()->SetNdivisions(4);
hist->GetXaxis()->SetNdivisions(4);
hist->GetXaxis()->SetLabelSize(0.1);
//hist->GetXaxis()->SetTitle("#theta, deg");
hist->GetXaxis()->SetTitle("N_{ph. el.} x 10");
hist->GetXaxis()->SetTitleSize(0.1);
hist->GetYaxis()->SetTitle("Counts");
hist->GetYaxis()->SetTitleSize(0.1);
hist->GetYaxis()->SetTitleOffset(0.5);
hist->GetXaxis()->SetTitleOffset(0.9);
TGaxis::SetMaxDigits(3);

TF1 *func1;

//cout << ledg[3][j-2] << "\n";

if (left_right) {
if (i == 0) func1 = new TF1("fit1",poissonf,ledg1[j-2],400,3);
if (i == 1) func1 = new TF1("fit1",poissonf,ledg2[j-2],400,3);
if (i == 2) func1 = new TF1("fit1",poissonf,ledg3[j-2],400,3);
if (i == 3) func1 = new TF1("fit1",poissonf,ledg4[j-2],400,3);
if (i == 4) func1 = new TF1("fit1",poissonf,ledg5[j-2],400,3);
if (i == 5) func1 = new TF1("fit1",poissonf,ledg6[j-2],400,3);
} else {
if (i == 0) func1 = new TF1("fit1",poissonf,redg1[j-2],400,3);
if (i == 1) func1 = new TF1("fit1",poissonf,redg2[j-2],400,3);
if (i == 2) func1 = new TF1("fit1",poissonf,redg3[j-2],400,3);
if (i == 3) func1 = new TF1("fit1",poissonf,redg4[j-2],400,3);
if (i == 4) func1 = new TF1("fit1",poissonf,redg5[j-2],400,3);
if (i == 5) func1 = new TF1("fit1",poissonf,redg6[j-2],400,3);
};


//func1->SetParameters(100000.,1.,20);
func1->SetParameter(0,1000.);
func1->SetParameter(1,1.);
func1->SetParameter(2,20.);

///func1->SetParLimits(0,99999.9,100000.1);
//func1->SetParLimits(1,1.,1.);
//func1->SetParLimits(2,20.,20.);
hist->Fit("fit1","QMNR");
func1->SetFillColor(kBlack);
func1->SetFillStyle(3436);
func1->SetLineColor(kBlack);
func1->Draw("lf same");

p0 = func1->GetParameter(0);
p1 = func1->GetParameter(1);
p2 = func1->GetParameter(2);

TF1 *f1;

if (left_right) {
f1 = new TF1("f1",poissonf,0,ledg4[j-2],3);
} else {
f1 = new TF1("f1",poissonf,0,redg4[j-2],3);
};

f1->SetParameter(0,p0);
f1->SetParameter(1,p1);
f1->SetParameter(2,p2);
  f1->SetLineColor(kBlack);
  f1->SetLineWidth(2);
f1->SetFillColor(kBlack);
f1->SetFillStyle(1001);  
   f1->Draw("lf same");

if (left_right) {   
if (i == 0)    fract[0][i][j] = f1->Integral(0,500.)/f1->Integral(ledg1[j-2],500.);
if (i == 1)    fract[0][i][j] = f1->Integral(0,500.)/f1->Integral(ledg2[j-2],500.);
if (i == 2)    fract[0][i][j] = f1->Integral(0,500.)/f1->Integral(ledg3[j-2],500.);
if (i == 3)    fract[0][i][j] =f1->Integral(0,500.)/f1-> Integral(ledg4[j-2],500.);
if (i == 4)    fract[0][i][j] = f1->Integral(0,500.)/f1->Integral(ledg5[j-2],500.);
if (i == 5)    fract[0][i][j] = f1->Integral(0,500.)/f1->Integral(ledg6[j-2],500.);
} else { 
if (i == 0)    fract[2][i][j] = f1->Integral(0,500.)/f1->Integral(redg1[j-2],500.);
if (i == 1)    fract[2][i][j] = f1->Integral(0,500.)/f1->Integral(redg2[j-2],500.);
if (i == 2)    fract[2][i][j] = f1->Integral(0,500.)/f1->Integral(redg3[j-2],500.);
if (i == 3)    fract[2][i][j] =f1->Integral(0,500.)/f1-> Integral(redg4[j-2],500.);
if (i == 4)    fract[2][i][j] = f1->Integral(0,500.)/f1->Integral(redg5[j-2],500.);
if (i == 5)    fract[2][i][j] = f1->Integral(0,500.)/f1->Integral(redg6[j-2],500.);
};

//cout << fract[i][j-2] << "\n";

//cout << p1 << "   " <<p2 <<"   " <<p3<< "\n";

TLine *line1;
if (left_right) { 
if (i == 0) line1 = new TLine(ledg1[j-2],0.,ledg1[j-2],hist->GetMaximum()); 
 if (i == 1) line1 = new TLine(ledg2[j-2],0.,ledg2[j-2],hist->GetMaximum());        
if (i == 2) line1 = new TLine(ledg3[j-2],0.,ledg3[j-2],hist->GetMaximum()); 	
if (i == 3) line1 = new TLine(ledg4[j-2],0.,ledg4[j-2],hist->GetMaximum()); 	
if (i == 4) line1 = new TLine(ledg5[j-2],0.,ledg5[j-2],hist->GetMaximum()); 	
if (i == 5) line1 = new TLine(ledg6[j-2],0.,ledg6[j-2],hist->GetMaximum()); 	} else { 
if (i == 0) line1 = new TLine(redg1[j-2],0.,redg1[j-2],hist->GetMaximum()); 
 if (i == 1) line1 = new TLine(redg2[j-2],0.,redg2[j-2],hist->GetMaximum());        
if (i == 2) line1 = new TLine(redg3[j-2],0.,redg3[j-2],hist->GetMaximum()); 	
if (i == 3) line1 = new TLine(redg4[j-2],0.,redg4[j-2],hist->GetMaximum()); 	
if (i == 4) line1 = new TLine(redg5[j-2],0.,redg5[j-2],hist->GetMaximum()); 	
if (i == 5) line1 = new TLine(redg6[j-2],0.,redg6[j-2],hist->GetMaximum()); 
};
	line1->SetLineColor(1);
        line1->SetLineWidth(3);
    line1->Draw("same");
};
};

for (i=0; i<0; i++) {

qqq.str("");
qqq << "c1_arr" << i+1;
c1_arr[i] = new TCanvas (qqq.str().c_str(),qqq.str().c_str(),0,0,1000,1000);
c1_arr[i]->Divide(3,5);//0.048,0.

for (j=2; j<17; j++) {

c1_arr[i]->cd(j-1);
c1_arr[i]->cd(j-1)->SetBottomMargin(0.15);
c1_arr[i]->cd(j-1)->SetTopMargin(0.15);
c1_arr[i]->cd(j-1)->SetLeftMargin(0.2);
qqq.str("");
qqq << "ph_el_sector_" << i+1 << "/left/sector_" << i+1 << "_left_"  << j+1;
MyFile->GetObject(qqq.str().c_str(),hist1); 

hist1->GetYaxis()->SetNoExponent(1);
hist1->SetLineColor(kBlack);
hist1->Draw();
if (j ==4) hist1->SetAxisRange(1,500,"X");
//if (j ==4) hist1->SetAxisRange(0,hist1->GetMaximum(),"Y");

hist1->GetYaxis()->SetLabelSize(0.1);
hist1->GetYaxis()->SetNdivisions(4);
hist1->GetXaxis()->SetNdivisions(6);
hist1->GetXaxis()->SetLabelSize(0.1);
hist1->GetXaxis()->SetTitle("Number of ph_el");
hist1->GetXaxis()->SetTitleSize(0.1);

TF1 *func2;

if (left_right) {
if (i == 0) func2 = new TF1("fit2",poissonf,ledg1[j-2],500,3);
if (i == 1) func2 = new TF1("fit2",poissonf,ledg2[j-2],500,3);
if (i == 2) func2 = new TF1("fit2",poissonf,ledg3[j-2],500,3);
if (i == 3) func2 = new TF1("fit2",poissonf,ledg4[j-2],500,3);
if (i == 4) func2 = new TF1("fit2",poissonf,ledg5[j-2],500,3);
if (i == 5) func2 = new TF1("fit2",poissonf,ledg6[j-2],500,3);
} else {
if (i == 0) func2 = new TF1("fit2",poissonf,redg1[j-2],500,3);
if (i == 1) func2 = new TF1("fit2",poissonf,redg2[j-2],500,3);
if (i == 2) func2 = new TF1("fit2",poissonf,redg3[j-2],500,3);
if (i == 3) func2 = new TF1("fit2",poissonf,redg4[j-2],500,3);
if (i == 4) func2 = new TF1("fit2",poissonf,redg5[j-2],500,3);
if (i == 5) func2 = new TF1("fit2",poissonf,redg6[j-2],500,3);
};

func2->SetParameter(0,100000.);
func2->SetParameter(1,1.);
func2->SetParameter(2,20.);

hist1->Fit("fit2","QMNR");

func2->Draw("same");

p3 = func2->GetParameter(0);
p4 = func2->GetParameter(1);
p5 = func2->GetParameter(2);

TF1 *f2 = new TF1("f2",poissonf,0,500,3);
f2->SetParameter(0,p3);
f2->SetParameter(1,p4);
f2->SetParameter(2,p5);
  f2->SetLineColor(kBlack);
  f2->SetLineWidth(2);
   f2->Draw("same");

if (left_right) {   
if (i == 0)    fract[0][i][j] = f2->Integral(0,500.)/f2->Integral(ledg1[j-2],500.);
if (i == 1)    fract[0][i][j] = f2->Integral(0,500.)/f2->Integral(ledg2[j-2],500.);
if (i == 2)    fract[0][i][j] = f2->Integral(0,500.)/f2->Integral(ledg3[j-2],500.);
if (i == 3)    fract[0][i][j] = f2->Integral(0,500.)/f2-> Integral(ledg4[j-2],500.);
if (i == 4)    fract[0][i][j] = f2->Integral(0,500.)/f2->Integral(ledg5[j-2],500.);
if (i == 5)    fract[0][i][j] = f2->Integral(0,500.)/f2->Integral(ledg6[j-2],500);
} else {
if (i == 0)    fract[2][i][j] = f2->Integral(0,500.)/f2->Integral(redg1[j-2],500.);
if (i == 1)    fract[2][i][j] = f2->Integral(0,500.)/f2->Integral(redg2[j-2],500.);
if (i == 2)    fract[2][i][j] = f2->Integral(0,500.)/f2->Integral(redg3[j-2],500.);
if (i == 3)    fract[2][i][j] = f2->Integral(0,500.)/f2-> Integral(redg4[j-2],500.);
if (i == 4)    fract[2][i][j] = f2->Integral(0,500.)/f2->Integral(redg5[j-2],500.);
if (i == 5)    fract[2][i][j] = f2->Integral(0,500.)/f2->Integral(redg6[j-2],500);
};


if (i == 0)    fract[1][i][j] = 1.;
if (i == 1)    fract[1][i][j] = 1.;
if (i == 2)    fract[1][i][j] = 1.;
if (i == 3)    fract[1][i][j] = 1.;
if (i == 4)    fract[1][i][j] = 1.;
if (i == 5)    fract[1][i][j] = 1.;


//cout << fract[1][i][j-2] << "\n";

//cout << p1 << "   " <<p2 <<"   " <<p3<< "\n";

TLine *line2;

if (i == 0) line2 = new TLine(redg1[j-2],0.,redg1[j-2],hist1->GetMaximum()); 
if (i == 1) line2 = new TLine(redg2[j-2],0.,redg2[j-2],hist1->GetMaximum());        
if (i == 2) line2 = new TLine(redg3[j-2],0.,redg3[j-2],hist1->GetMaximum()); 	
if (i == 3) line2 = new TLine(redg4[j-2],0.,redg4[j-2],hist1->GetMaximum()); 	
if (i == 4) line2 = new TLine(redg5[j-2],0.,redg5[j-2],hist1->GetMaximum()); 	
if (i == 5) line2 = new TLine(redg6[j-2],0.,redg6[j-2],hist1->GetMaximum()); 	
	line2->SetLineColor(2);
       line2->SetLineWidth(3);
  line2->Draw("same");



};
};
std::ofstream ofs ("phel_integr_fract.txt", std::ofstream::out);
for(k=0; k<3; k++){
for(i=0; i<6; i++){
fract[k][i][0] = 0.;
fract[k][i][1] = 0.;
fract[k][i][17] = 0.;
for(j=0; j<18; j++){
ofs << k << "," << i << "," <<j << "," << fract[k][i][j]<< "\n";
cout << k << "," << i << "," <<j << "," << fract[k][i][j]<< "\n";
};
};
};


  ofs << "lorem ipsum";

  ofs.close();
};

/*

double poissonf( double*x, double*par){

return s = par[0]*TMath::Poisson(x[0]/par[2],par[1]);
} ;
*/


