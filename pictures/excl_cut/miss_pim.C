void  miss_pim() {

gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.12,"t");
gStyle->SetTitleY(1.02);
gStyle->SetTitleX(0.65);
//gStyle->SetTitleAlign(33); 
TGaxis::SetMaxDigits(2);

ostringstream qqq;

TFile *MyFile_data = new TFile("out_ind4_miss_data.root","READ");

//pim_miss_radcorr_i

TFile *MyFile_sim = new TFile("out_ind4_sim_low_w.root","READ");


TFile *MyFile_ind4 = new TFile("out_ind4_miss_high_w.root","READ");

TCanvas *c = new TCanvas ("c","c",0,0,200,600);
//pim_miss_radcorr_sim_i




TH1 *h_data[20];
TH1 *h_data_pip[20];
TH1 *h_data_pr[20];
TH1 *h_data_0[20];
TH1 *h_sim[20];
TH1 *h_sim_pip[20];
TH1 *h_sim_pr[20];
TH1 *h_sim_0[20];
TH1 *h_bckgr[20];
TH1 *h_bckgr_pip[20];
TH1 *h_bckgr_pr[20];
TH1 *h_bckgr_0[20];

c->Divide(1,4);

Int_t i;

for (Int_t j=3; j<4; j++){

i = j*5;

cout << i << endl;

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pim_mis_radcorr_" << i;

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data[j]);

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pip_mis_radcorr_" << i;

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_pip[j]);

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pr_mis_radcorr_" << i;

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_pr[j]);


qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "0_mis_radcorr_" << i;

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_0[j]);



qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pim_mis_radcorr_sim_" << i;

//cout << qqq.str() << endl;
if(j<4)MyFile_sim->GetObject(qqq.str().c_str(),h_sim[j]);
if(j>=3) {
MyFile_ind4->GetObject(qqq.str().c_str(),h_sim[j]);

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pim_mis_radcorr_sim_bckg_" << i;
MyFile_ind4->GetObject(qqq.str().c_str(),h_bckgr[j]);
};


qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pip_mis_radcorr_sim_" << i;

//cout << qqq.str() << endl;
if(j<4)MyFile_sim->GetObject(qqq.str().c_str(),h_sim_pip[j]);
if(j>=3) {
MyFile_ind4->GetObject(qqq.str().c_str(),h_sim_pip[j]);

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pip_mis_radcorr_sim_bckg_" << i;
MyFile_ind4->GetObject(qqq.str().c_str(),h_bckgr_pip[j]);
};


qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pr_mis_radcorr_sim_" << i;

//cout << qqq.str() << endl;
if(j<4)MyFile_sim->GetObject(qqq.str().c_str(),h_sim_pr[j]);
if(j>=3) {
MyFile_ind4->GetObject(qqq.str().c_str(),h_sim_pr[j]);

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "pr_mis_radcorr_sim_bckg_" << i;
MyFile_ind4->GetObject(qqq.str().c_str(),h_bckgr_pr[j]);
};


qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "0_mis_radcorr_sim_" << i;

//cout << qqq.str() << endl;
if(j<4)MyFile_sim->GetObject(qqq.str().c_str(),h_sim_0[j]);
if(j>=3) {
MyFile_ind4->GetObject(qqq.str().c_str(),h_sim_0[j]);

qqq.str("");

qqq << "q2_0.475/" << "w_1." << 3125+250*i << "/" << "0_mis_radcorr_sim_bckg_" << i;
MyFile_ind4->GetObject(qqq.str().c_str(),h_bckgr_0[j]);
};

//MyFile_sim->GetObject(qqq.str().c_str(),h_sim[j]);

c->cd(1);
c->cd(1)->SetTopMargin(0.19);
c->cd(1)->SetBottomMargin(0.32);
c->cd(1)->SetRightMargin(0.01);
c->cd(1)->SetLeftMargin(0.25);


//h_data->Draw();
h_sim[j]->SetLineColor(1);
h_sim[j]->SetLineStyle(2);
h_data[j]->Rebin(2);
h_sim[j]->Rebin(2);
if(j>=3) {
h_bckgr[j]->Rebin(2);
if(j==4)h_bckgr[j]->Scale(7.2);
if(j==3)h_bckgr[j]->Scale(10.5);
h_sim[j]->Add(h_bckgr[j]);
};

h_sim[j]->Scale(h_data[j]->GetMaximum()/(h_sim[j]->GetMaximum()));
//h_sim[j]->Scale(h_data[j]->Integral(1,200)/(h_sim[j]->Integral(1,200)));


h_sim[j]->SetLineWidth(1);
h_data[j]->SetLineWidth(1);
h_data[j]->SetLineColor(1);
h_sim[j]->SetLineColor(1);

qqq.str("");

qqq << "1." << 300+25*i << " < W < 1." << 325+25*i << " GeV";

h_sim[j]->SetTitle(qqq.str().c_str());
h_sim[j]->SetTitle("");

h_sim[j]->GetXaxis()->SetTitle("M_{#pi^{-}}^{2} (GeV^{2})");
h_sim[j]->GetXaxis()->SetLabelSize(0.14);
h_sim[j]->GetXaxis()->SetNdivisions(5);
h_sim[j]->GetXaxis()->SetTitleSize(0.14);
h_sim[j]->GetXaxis()->SetNoExponent(1);

h_sim[j]->GetYaxis()->SetTitle("counts");
h_sim[j]->GetYaxis()->SetLabelSize(0.14);
h_sim[j]->GetYaxis()->SetNdivisions(5);
h_sim[j]->GetYaxis()->SetTitleSize(0.14);
h_sim[j]->GetYaxis()->SetTitleOffset(0.95);




//h_sim[j]->SetMarkerStyle(20);

h_sim[j]->SetLineStyle(1);
//h_sim[j]->SetMarkerStyle(20);
//h_sim[j]->SetFillColor(1);
//h_sim[j]->SetFillStyle(3844);
h_data[j]->SetMarkerStyle(29);
h_data[j]->SetMarkerSize(0.6);
h_sim[j]->Draw("HIST");
h_data[j]->Draw("P same");
TArrow *ar1 = new TArrow(0.06,(h_sim[j]->GetMaximum())/2.,0.06,0.,0.02,"|>");
ar1->SetLineWidth(2);
ar1->SetLineColor(kBlack);
ar1->SetFillColor(kBlack);
ar1->Draw();

TArrow *ar2 = new TArrow(-0.04,(h_sim[j]->GetMaximum())/2.,-0.04,0.,0.02,"|>");
ar2->SetLineWidth(2);
ar2->SetLineColor(kBlack);
ar2->SetFillColor(kBlack);
ar2->Draw();
//cout << h_data->GetNbinsX() << endl;


// pip


c->cd(2);
c->cd(2)->SetTopMargin(0.17);
c->cd(2)->SetBottomMargin(0.32);
c->cd(2)->SetRightMargin(0.01);
c->cd(2)->SetLeftMargin(0.25);

h_data_pip[j]->Rebin(2);

h_sim_pip[j]->SetLineColor(1);
h_sim_pip[j]->Rebin(2);
if(j>=3) {
h_bckgr_pip[j]->Rebin(2);
if(j==4)h_bckgr_pip[j]->Scale(7.2);
if(j==3)h_bckgr_pip[j]->Scale(10.5);
h_sim_pip[j]->Add(h_bckgr_pip[j]);
};
h_sim_pip[j]->Scale(h_data_pip[j]->GetMaximum()/(h_sim_pip[j]->GetMaximum()));
//h_sim[j]->Scale(h_data[j]->Integral(1,200)/(h_sim[j]->Integral(1,200)));
h_sim_pip[j]->SetLineWidth(1);
h_data_pip[j]->SetLineWidth(1);

qqq.str("");

qqq << "1." << 300+25*i << " < W < 1." << 325+25*i << " GeV";

h_sim_pip[j]->SetTitle(qqq.str().c_str());
h_sim_pip[j]->SetTitle("");

h_sim_pip[j]->GetXaxis()->SetTitle("M_{#pi^{+}}^{2} (GeV^{2})");
h_sim_pip[j]->GetXaxis()->SetLabelSize(0.14);
h_sim_pip[j]->GetXaxis()->SetNdivisions(5);
h_sim_pip[j]->GetXaxis()->SetTitleSize(0.14);
h_sim_pip[j]->GetXaxis()->SetNoExponent(1);

h_sim_pip[j]->GetYaxis()->SetTitle("counts");
h_sim_pip[j]->GetYaxis()->SetLabelSize(0.14);
h_sim_pip[j]->GetYaxis()->SetNdivisions(5);
h_sim_pip[j]->GetYaxis()->SetTitleSize(0.14);
h_sim_pip[j]->GetYaxis()->SetTitleOffset(0.95);



h_sim_pip[j]->SetLineStyle(1);

h_data_pip[j]->SetMarkerStyle(29);
h_data_pip[j]->SetMarkerSize(0.6);
h_sim_pip[j]->Draw("HIST");
h_data_pip[j]->Draw("P same");

qqq.str("");


TArrow *ar3 = new TArrow(0.06,(h_data_pip[j]->GetMaximum())/2.,0.06,0.,0.02,"|>");
ar3->SetLineWidth(2);
ar3->SetLineColor(kBlack);
ar3->SetFillColor(kBlack);
ar3->Draw();

TArrow *ar4 = new TArrow(-0.04,(h_data_pip[j]->GetMaximum())/2.,-0.04,0.,0.02,"|>");
ar4->SetLineWidth(2);
ar4->SetLineColor(kBlack);
ar4->SetFillColor(kBlack);
ar4->Draw();



// proton

c->cd(3);
c->cd(3)->SetTopMargin(0.17);
c->cd(3)->SetBottomMargin(0.32);
c->cd(3)->SetRightMargin(0.01);
c->cd(3)->SetLeftMargin(0.25);

h_data_pr[j]->Rebin(2);

h_sim_pr[j]->SetLineColor(1);
h_sim_pr[j]->Rebin(2);
if(j>=3) {
h_bckgr_pr[j]->Rebin(2);
if(j==4)h_bckgr_pr[j]->Scale(7.2);
if(j==3)h_bckgr_pr[j]->Scale(10.5);
h_sim_pr[j]->Add(h_bckgr_pr[j]);
};
h_sim_pr[j]->Scale(h_data_pr[j]->GetMaximum()/(h_sim_pr[j]->GetMaximum()));
//h_sim[j]->Scale(h_data[j]->Integral(1,200)/(h_sim[j]->Integral(1,200)));
h_sim_pr[j]->SetLineWidth(1);
h_data_pr[j]->SetLineWidth(1);

qqq.str("");

qqq << "1." << 300+25*i << " < W < 1." << 325+25*i << " GeV";

h_sim_pr[j]->SetTitle(qqq.str().c_str());
h_sim_pr[j]->SetTitle("");

h_sim_pr[j]->GetXaxis()->SetTitle("M_{p'}^{2} (GeV^{2})");
h_sim_pr[j]->GetXaxis()->SetLabelSize(0.14);
h_sim_pr[j]->GetXaxis()->SetNdivisions(2);
h_sim_pr[j]->GetXaxis()->SetTitleSize(0.14);
h_sim_pr[j]->GetXaxis()->SetNoExponent(1);

h_sim_pr[j]->GetYaxis()->SetTitle("counts");
h_sim_pr[j]->GetYaxis()->SetLabelSize(0.14);
h_sim_pr[j]->GetYaxis()->SetNdivisions(5);
h_sim_pr[j]->GetYaxis()->SetTitleSize(0.14);
h_sim_pr[j]->GetYaxis()->SetTitleOffset(0.955);



h_sim_pr[j]->SetLineStyle(1);

h_data_pr[j]->SetMarkerStyle(29);
h_data_pr[j]->SetMarkerSize(0.6);
h_sim_pr[j]->Draw("HIST");
h_data_pr[j]->Draw("P same");



TArrow *ar5 = new TArrow(0.882-0.04,(h_data_pr[j]->GetMaximum())/2.,0.882-0.04,0.,0.02,"|>");
ar5->SetLineWidth(2);
ar5->SetLineColor(kBlack);
ar5->SetFillColor(kBlack);
ar5->Draw();

TArrow *ar6 = new TArrow(0.882+0.04,(h_data_pr[j]->GetMaximum())/2.,0.882+0.04,0.,0.02,"|>");
ar6->SetLineWidth(2);
ar6->SetLineColor(kBlack);
ar6->SetFillColor(kBlack);
ar6->Draw();




// 0


c->cd(4);
c->cd(4)->SetTopMargin(0.17);
c->cd(4)->SetBottomMargin(0.32);
c->cd(4)->SetRightMargin(0.01);
c->cd(4)->SetLeftMargin(0.25);

//h_data_0[j]->Rebin(2);

h_sim_0[j]->SetLineColor(1);
//h_sim_0[j]->Rebin(2);
if(j>=3) {
//h_bckgr_0[j]->Rebin(2);
if(j==4)h_bckgr_0[j]->Scale(7.2);
if(j==3)h_bckgr_0[j]->Scale(10.5);
h_sim_0[j]->Add(h_bckgr_0[j]);
};
h_sim_0[j]->Scale(h_data_0[j]->GetMaximum()/(h_sim_0[j]->GetMaximum()));
//h_sim[j]->Scale(h_data[j]->Integral(1,200)/(h_sim[j]->Integral(1,200)));
h_sim_0[j]->SetLineWidth(1);
h_data_0[j]->SetLineWidth(1);

qqq.str("");

qqq << "1." << 300+25*i << " < W < 1." << 325+25*i << " GeV";

h_sim_0[j]->SetTitle(qqq.str().c_str());
h_sim_0[j]->SetTitle("");

h_sim_0[j]->GetXaxis()->SetTitle("M_{0}^{2} (GeV^{2})");
h_sim_0[j]->GetXaxis()->SetLabelSize(0.14);
h_sim_0[j]->GetXaxis()->SetNdivisions(5);
h_sim_0[j]->GetXaxis()->SetTitleSize(0.14);
h_sim_0[j]->GetXaxis()->SetNoExponent(1);
h_sim_0[j]->GetXaxis()->SetRangeUser(-0.025,0.025);
h_data_0[j]->GetXaxis()->SetRangeUser(-0.025,0.025);


h_sim_0[j]->GetYaxis()->SetTitle("counts");
h_sim_0[j]->GetYaxis()->SetLabelSize(0.14);
h_sim_0[j]->GetYaxis()->SetNdivisions(5);
h_sim_0[j]->GetYaxis()->SetTitleSize(0.14);
h_sim_0[j]->GetYaxis()->SetTitleOffset(0.95);


h_sim_0[j]->SetLineStyle(1);

h_data_0[j]->SetMarkerStyle(29);
h_data_0[j]->SetMarkerSize(0.6);
h_sim_0[j]->Draw("HIST");
h_data_0[j]->Draw("P same");





TArrow *ar7 = new TArrow(-0.01,(h_data_0[j]->GetMaximum())/2.,-0.01,0.,0.02,"|>");
ar7->SetLineWidth(2);
ar7->SetLineColor(kBlack);
ar7->SetFillColor(kBlack);
ar7->Draw();

TArrow *ar8 = new TArrow(0.01,(h_data_0[j]->GetMaximum())/2.,0.01,0.,0.02,"|>");
ar8->SetLineWidth(2);
ar8->SetLineColor(kBlack);
ar8->SetFillColor(kBlack);
ar8->Draw();



};


c->Print("miss_mass.pdf");

}
