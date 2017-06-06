void sim_mom_corr2() {


TFile *MyFile = new TFile("sim_mom_corr.root","READ");
TF1 *myfit = new TF1("myfit","gaus", -0.035, 0.035);
TF1 *myfit_pol2 = new TF1("myfit_pol2","pol2", 0.5, 2.0);
TF1 *myfit_pol2_th = new TF1("myfit_pol2_th","pol2", 20., 50.);

TH1F *h_p0 = new TH1F("P0","P0",7,29.,50.);
TH1F *h_p1 = new TH1F("P1","P1",7,29.,50.);
TH1F *h_p2 = new TH1F("P2","P2",7,29.,50.);

TCanvas *c_param = new TCanvas("parameters","parameters",0,0,600,200);
c_param->Divide(3,1);

TH1 *h[10][10];
ostringstream qqq;
TCanvas *c[10];
TH1 *h_p[10];
c_p = new TCanvas("c_p","c_p",0,0,500,500);
c_p->Divide(4,3);

for (Int_t i=0; i<10; i++){

qqq.str("");
qqq << "c_" << i;
//c[i] = new TCanvas(qqq.str().c_str(),qqq.str().c_str(),0,0,500,500);
//c[i]->Divide(4,3);

qqq.str("");
qqq << "th_" << i;

if (i==0) {
h_p[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),7,0.95,2.0);
} elseif (i==1) {
h_p[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),8,0.8,2.0);
} else {
h_p[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),10,0.5,2.0);
};

for (Int_t j=0; j<10; j++){


//c[i]->cd(j+1);



qqq.str("");
qqq << "Pdiff_" << j << "_" << i;

MyFile->GetObject(qqq.str().c_str(),h[i][j]);

qqq.str("");

h[i][j]->Rebin(4);
//h[i][j]->Draw();
h[i][j]->Fit("myfit","QN");
//TF1 *fit = hist->GetFunction();

h_p[i]->Fill((0.575 + 0.15*j),myfit->GetParameter(1));
//cout << (21.5 + 3.*j) << " " << myfit->GetParameter(1) << endl;

};

c_p->cd(i+1);
h_p[i]->Draw();
h_p[i]->Fit("myfit_pol2","Q");

h_p0->Fill((21.5+3.*i),myfit_pol2->GetParameter(0));
h_p1->Fill((21.5+3.*i),myfit_pol2->GetParameter(1));
h_p2->Fill((21.5+3.*i),myfit_pol2->GetParameter(2));

};
c_param->cd(1);
h_p0->Draw();
h_p0->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

c_param->cd(2);
h_p1->Draw();
h_p1->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;

c_param->cd(3);
h_p2->Draw();
h_p2->Fit("myfit_pol2_th","Q");

cout << myfit_pol2_th->GetParameter(0) << "  " << myfit_pol2_th->GetParameter(1) << " " << myfit_pol2_th->GetParameter(2) << endl;
};
