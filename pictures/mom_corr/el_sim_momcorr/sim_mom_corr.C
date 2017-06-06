void sim_mom_corr() {


TFile *MyFile = new TFile("sim_mom_corr.root","READ");
TF1 *myfit = new TF1("myfit","gaus", -0.035, 0.035);

TH1 *h[10][10];
ostringstream qqq;
TCanvas *c[10];
TH1 *h_th[10];
c_th = new TCanvas("c_th","c_th",0,0,500,500);
c_th->Divide(4,3);

for (Int_t i=0; i<10; i++){

qqq.str("");
qqq << "c_" << i;
c[i] = new TCanvas(qqq.str().c_str(),qqq.str().c_str(),0,0,500,500);
c[i]->Divide(4,3);

qqq.str("");
qqq << "th_" << i;
h_th[i] = new TH1F(qqq.str().c_str(),qqq.str().c_str(),10,20.,50.);

for (Int_t j=0; j<10; j++){


c[i]->cd(j+1);



qqq.str("");
qqq << "Pdiff_" << i << "_" << j;

MyFile->GetObject(qqq.str().c_str(),h[i][j]);

qqq.str("");

h[i][j]->Rebin(4);
h[i][j]->Draw();
h[i][j]->Fit("myfit","Q");
//TF1 *fit = hist->GetFunction();

h_th[i]->Fill((21.5 + 3.*j),myfit->GetParameter(1));
//cout << (21.5 + 3.*j) << " " << myfit->GetParameter(1) << endl;

};

c_th->cd(i+1);
h_th[i]->Draw();
h_th[i]->Fit("pol1","Q");

};

};
