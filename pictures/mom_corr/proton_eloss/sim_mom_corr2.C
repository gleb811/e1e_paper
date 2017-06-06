void sim_mom_corr2() {

gStyle->SetOptStat(0);

TFile *MyFile = new TFile("sim_proton_mom_corr.root","READ");
TF1 *myfit = new TF1("myfit","gaus", -0.035, 0.035);
TF1 *myfit_pol2 = new TF1("myfit_pol2","pol5", 0.195, 2.0);
TF1 *myfit_pol2_th = new TF1("myfit_pol2_th","pol2", 12., 120.);

TH1F *h_p0 = new TH1F("P0","P0",18,12.,120.);
TH1F *h_p1 = new TH1F("P1","P1",18,12.,120.);
TH1F *h_p2 = new TH1F("P2","P2",18,12.,120.);
TH1F *h_p3 = new TH1F("P3","P3",18,12.,120.);
TH1F *h_p4 = new TH1F("P4","P4",18,12.,120.);
TH1F *h_p5 = new TH1F("P5","P5",18,12.,120.);

const Int_t npar = 18;
   Double_t eloss2d_params[npar];
   TF2* f_eloss2d = new TF2("f_eloss2d",&eloss2d,0.,120.,0.25,1.9, npar);
   TH2F *h_eloss2d = new TH2F("h_elos2d","h_elos2d",40,0.,120.,40,0.25,1.9);
   
TCanvas *c_corr2d = new TCanvas("c_corr2d","c_corr2d",0,0,600,600);


TCanvas *c_param = new TCanvas("parameters","parameters",0,0,1200,200);
c_param->Divide(6,1);

TH1 *h[20][20];
ostringstream qqq;
TCanvas *c[20];
TH1 *h_p[20];
c_p = new TCanvas("c_p","c_p",0,0,500,400);
c_p->Divide(5,4);

for (Int_t i=2; i<=19; i++){

qqq.str("");
qqq << "c_" << i;
c[i] = new TCanvas(qqq.str().c_str(),qqq.str().c_str(),0,0,500,400);
c[i]->Divide(5,4);

qqq.str("");
qqq << "th_" << i;

/*if (i==0) {
h_p[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),7,0.95,2.0);
} elseif (i==1) {
h_p[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),8,0.8,2.0);
} else {
*/
h_p[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),14,0.195,2.0);
//};

for (Int_t j=1; j<=19; j++){


c[i]->cd(j+1);



qqq.str("");
qqq << "Pdiff_" << j << "_" << i;

MyFile->GetObject(qqq.str().c_str(),h[i][j]);

qqq.str("");

h[i][j]->Rebin(4);
h[i][j]->Draw();
h[i][j]->Fit("myfit","Q");
//TF1 *fit = hist->GetFunction();

h_p[i]->Fill((0.1 + 0.095*j),myfit->GetParameter(1));
//cout << (21.5 + 3.*j) << " " << myfit->GetParameter(1) << endl;

};

c_p->cd(i+1);
h_p[i]->Draw();
h_p[i]->Fit("myfit_pol2","Q");

h_p0->Fill((6.*i+3.),myfit_pol2->GetParameter(0));
h_p1->Fill((6.*i+3.),myfit_pol2->GetParameter(1));
h_p2->Fill((6.*i+3.),myfit_pol2->GetParameter(2));
h_p3->Fill((6.*i+3.),myfit_pol2->GetParameter(3));
h_p4->Fill((6.*i+3.),myfit_pol2->GetParameter(4));
h_p5->Fill((6.*i+3.),myfit_pol2->GetParameter(5));

};
c_param->cd(1);
h_p0->Draw();
h_p0->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

eloss2d_params[0] = myfit_pol2_th->GetParameter(0);
eloss2d_params[1] = myfit_pol2_th->GetParameter(1);
eloss2d_params[2] = myfit_pol2_th->GetParameter(2);

c_param->cd(2);
h_p1->Draw();
h_p1->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

eloss2d_params[3] = myfit_pol2_th->GetParameter(0);
eloss2d_params[4] = myfit_pol2_th->GetParameter(1);
eloss2d_params[5] = myfit_pol2_th->GetParameter(2);


c_param->cd(3);
h_p2->Draw();
h_p2->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

eloss2d_params[6] = myfit_pol2_th->GetParameter(0);
eloss2d_params[7] = myfit_pol2_th->GetParameter(1);
eloss2d_params[8] = myfit_pol2_th->GetParameter(2);


c_param->cd(4);
h_p3->Draw();
h_p3->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

eloss2d_params[9] = myfit_pol2_th->GetParameter(0);
eloss2d_params[10] = myfit_pol2_th->GetParameter(1);
eloss2d_params[11] = myfit_pol2_th->GetParameter(2);


c_param->cd(5);
h_p4->Draw();
h_p4->Fit("myfit_pol2_th","Q");


cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

eloss2d_params[12] = myfit_pol2_th->GetParameter(0);
eloss2d_params[13] = myfit_pol2_th->GetParameter(1);
eloss2d_params[14] = myfit_pol2_th->GetParameter(2);

c_param->cd(6);
h_p5->Draw();
h_p5->Fit("myfit_pol2_th","Q");


cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

eloss2d_params[15] = myfit_pol2_th->GetParameter(0);
eloss2d_params[16] = myfit_pol2_th->GetParameter(1);
eloss2d_params[17] = myfit_pol2_th->GetParameter(2);


f_eloss2d->SetParameters(eloss2d_params);

c_corr2d->cd();
gPad->SetTheta(10.); // default is 30
gPad->SetPhi(65.); // default is 30
gPad->Update();



h_eloss2d->SetTitle("");
h_eloss2d->GetZaxis()->SetTitle("Corr. factor (%)");
h_eloss2d->GetYaxis()->SetTitle("P_{p} (GeV)");
h_eloss2d->GetXaxis()->SetTitle("#theta_{p} (deg)");
h_eloss2d->GetZaxis()->SetTitleOffset(1.2);
h_eloss2d->GetYaxis()->SetTitleOffset(2.1);
h_eloss2d->GetXaxis()->SetTitleOffset(2.1);
h_eloss2d->SetLineColor(1);

h_eloss2d->FillRandom("f_eloss2d",100000000);
h_eloss2d->Scale(1./100000000.*40.*40.);
h_eloss2d->Draw("LEGO2Z");
//f_eloss2d->Draw("LEGO2Z");

c_corr2d->SetRightMargin(0.13);
c_corr2d->Update();
//TPaletteAxis *palette = (TPaletteAxis*)f_eloss2d->GetObject("palette");
TPaletteAxis *palette = (TPaletteAxis*)h_eloss2d->GetListOfFunctions()->FindObject("palette");

   palette->SetX1NDC(0.9);
   palette->SetX2NDC(0.95);
   palette->SetY1NDC(0.2);
   palette->SetY2NDC(0.8);
   gPad->Modified();
   gPad->Update();
   
c_corr2d-> SaveAs("eloss2d.pdf");  

//cout << h_eloss2d->GetListOfFunctions() << endl;
};

Double_t eloss2d(Double_t *x, Double_t *par) {
   for (Int_t i=0; i<18; i++) {
   cout << "par[" << i << "] = " << par[i] << ";\n";
   };
   Double_t p0 = par[0] + x[0]*par[1] + x[0]*x[0]*par[2];
   Double_t p1 = par[3] + x[0]*par[4] + x[0]*x[0]*par[5];
   Double_t p2 = par[6] + x[0]*par[7] + x[0]*x[0]*par[8];
   Double_t p3 = par[9] + x[0]*par[10] + x[0]*x[0]*par[11];   
   Double_t p4 = par[12] + x[0]*par[13] + x[0]*x[0]*par[14];  
   Double_t p5 = par[15] + x[0]*par[16] + x[0]*x[0]*par[17]; 
   
   return -1.*(p0+p1*x[1]+p2*x[1]*x[1]+p3*x[1]*x[1]*x[1]+p4*x[1]*x[1]*x[1]*x[1]+p5*x[1]*x[1]*x[1]*x[1]*x[1])/x[1]*100.;
   }
