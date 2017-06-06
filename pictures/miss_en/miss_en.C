void  miss_en() {

gStyle->SetOptStat(0);
gStyle->SetTitleSize(0.12,"t");
gStyle->SetTitleY(1.02);
gStyle->SetTitleX(0.65);
//gStyle->SetTitleAlign(33); 
TGaxis::SetMaxDigits(2);

ostringstream qqq;

TFile *MyFile_data = new TFile("data_new_time_corr.root","READ");

//pim_miss_radcorr_i

TFile *MyFile_sim = new TFile("out_miss_en.root","READ");


//TFile *MyFile_ind4 = new TFile("out_ind4_miss_high_w.root","READ");

TCanvas *c = new TCanvas ("c","c",0,0,800,300);
//pim_miss_radcorr_sim_i




TH1 *h_data_pim;
TH1 *h_data_pip;
TH1 *h_data_p;
TH1 *h_data_0;

TH1 *h_sim_pim;
TH1 *h_sim_pip;
TH1 *h_sim_p;
TH1 *h_sim_0;

/*
TH1 *h_sim[20];
TH1 *h_sim_pip[20];
TH1 *h_sim_pr[20];
TH1 *h_sim_0[20];
TH1 *h_bckgr[20];
TH1 *h_bckgr_pip[20];
TH1 *h_bckgr_pr[20];
TH1 *h_bckgr_0[20];
*/

c->Divide(4,1);


qqq.str("");

qqq << "hist_miss_en_0";

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_0);

qqq.str("");

qqq.str("");

qqq << "h_miss_en_0_sim";

//cout << qqq.str() << endl;
MyFile_sim->GetObject(qqq.str().c_str(),h_sim_0);

qqq.str("");

qqq << "hist_PIm_miss_en";

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_pim);

qqq.str("");



qqq << "h_PIm_miss_en_sim";

//cout << qqq.str() << endl;
MyFile_sim->GetObject(qqq.str().c_str(),h_sim_pim);

qqq.str("");

//h_data_pim->Rebin(4);
//h_sim_pim->Rebin(4);


qqq << "hist_PIp_miss_en";

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_pip);

qqq.str("");



qqq << "h_PIp_miss_en_sim";

//cout << qqq.str() << endl;
MyFile_sim->GetObject(qqq.str().c_str(),h_sim_pip);

qqq.str("");

//h_data_pip->Rebin(4);
//h_sim_pip->Rebin(4);


qqq << "hist_P_miss_en";

//cout << qqq.str() << endl;
MyFile_data->GetObject(qqq.str().c_str(),h_data_p);

qqq.str("");



qqq << "hist_P_miss_en_sim";

//cout << qqq.str() << endl;
MyFile_sim->GetObject(qqq.str().c_str(),h_sim_p);

qqq.str("");

//h_data_p->Rebin(4);
//h_sim_p->Rebin(4);


c->cd(1)->SetBottomMargin(0.13);
c->cd(1)->SetLeftMargin(0.15);
c->cd(1);

TGaxis::SetMaxDigits(4);
h_data_0->GetXaxis()->SetRangeUser(-0.3,0.5);
h_data_0->GetXaxis()->SetNdivisions(5);
h_data_0->SetLineColor(1);
h_data_0->SetLineWidth(2);
h_data_0->GetXaxis()->SetLabelSize(0.07);
h_data_0->GetYaxis()->SetLabelSize(0.07);
h_data_0->GetXaxis()->SetTitle("E_{0} (GeV)");
h_data_0->GetYaxis()->SetTitle("counts");
h_data_0->GetXaxis()->SetTitleSize(0.08);
h_data_0->GetXaxis()->SetTitleOffset(0.8);
h_data_0->GetYaxis()->SetTitleSize(0.08);
h_data_0->GetYaxis()->SetTitleOffset(0.95);
h_data_0->SetTitle("");
h_data_0->SetMaximum((h_data_0->GetMaximum())*1.5);
h_data_0->Draw();



h_sim_0->Scale((h_data_0->Integral())/(h_sim_0->Integral()));
h_sim_0->SetLineColor(1);
h_sim_0->SetLineWidth(2);
//h_sim_0->Draw("same");

TLine *ar_0 = new TLine(-0.05,(h_data_0->GetMaximum())/2.,-0.05,0.);
ar_0->SetLineColor(1);
ar_0->SetLineWidth(2);
//ar_0->SetFillColor(kGreen+2);
ar_0->Draw();




c->cd(2)->SetBottomMargin(0.13);
c->cd(2)->SetLeftMargin(0.15);
c->cd(2);

TGaxis::SetMaxDigits(4);
h_data_pim->GetXaxis()->SetRangeUser(-0.15,1.5);
h_data_pim->GetXaxis()->SetNdivisions(5);
h_data_pim->SetLineColor(1);
h_data_pim->SetLineWidth(2);
h_data_pim->GetXaxis()->SetLabelSize(0.07);
h_data_pim->GetYaxis()->SetLabelSize(0.07);
h_data_pim->GetXaxis()->SetTitle("E_{#pi^{-}} (GeV)");
h_data_pim->GetYaxis()->SetTitle("counts");
h_data_pim->GetXaxis()->SetTitleSize(0.08);
h_data_pim->GetXaxis()->SetTitleOffset(0.8);
h_data_pim->GetYaxis()->SetTitleSize(0.08);
h_data_pim->GetYaxis()->SetTitleOffset(0.95);
h_data_pim->SetTitle("");
h_data_pim->SetMaximum((h_data_pim->GetMaximum())*1.5);
h_data_pim->Draw();



h_sim_pim->Scale((h_data_pim->Integral())/(h_sim_pim->Integral()));
h_sim_pim->SetLineColor(1);
h_sim_pim->SetLineWidth(2);
//h_sim_pim->Draw("same");

TLine *ar_pim = new TLine(-0.05+0.13957,(h_data_pim->GetMaximum())/2.,-0.05+0.13957,0.);
ar_pim->SetLineColor(1);
ar_pim->SetLineWidth(2);
//ar_pim->SetFillColor(kGreen+2);
ar_pim->Draw();





c->cd(3)->SetBottomMargin(0.13);
c->cd(3)->SetLeftMargin(0.15);
c->cd(3);

TGaxis::SetMaxDigits(3);
h_data_pip->GetXaxis()->SetRangeUser(-0.15,1.5);
h_data_pip->GetXaxis()->SetNdivisions(5);
h_data_pip->SetLineColor(1);
h_data_pip->SetLineWidth(2);
h_data_pip->GetXaxis()->SetLabelSize(0.07);
h_data_pip->GetYaxis()->SetLabelSize(0.07);
h_data_pip->GetXaxis()->SetTitle("E_{#pi^{+}} (GeV)");
h_data_pip->GetYaxis()->SetTitle("counts");
h_data_pip->GetXaxis()->SetTitleSize(0.08);
h_data_pip->GetXaxis()->SetTitleOffset(0.8);
h_data_pip->GetYaxis()->SetTitleSize(0.08);
h_data_pip->GetYaxis()->SetTitleOffset(0.95);
h_data_pip->SetTitle("");
h_data_pip->SetMaximum((h_data_pip->GetMaximum())*1.5);
h_data_pip->Draw();



h_sim_pip->Scale((h_data_pip->Integral())/(h_sim_pip->Integral()));
h_sim_pip->SetLineColor(1);
h_sim_pip->SetLineWidth(2);
//h_sim_pip->Draw("same");

TLine *ar_pip = new TLine(-0.05+0.13957,(h_data_pip->GetMaximum())/2.,-0.05+0.13957,0.);
ar_pip->SetLineColor(1);
ar_pip->SetLineWidth(2);
//ar_pip->SetFillColor(kGreen+2);
ar_pip->Draw();










c->cd(4)->SetBottomMargin(0.13);
c->cd(4)->SetLeftMargin(0.15);
c->cd(4);

TGaxis::SetMaxDigits(3);
h_data_p->GetXaxis()->SetRangeUser(0.8,2.);
h_data_p->GetXaxis()->SetNdivisions(5);
h_data_p->SetLineColor(1);
h_data_p->SetLineWidth(2);
h_data_p->GetXaxis()->SetLabelSize(0.07);
h_data_p->GetYaxis()->SetLabelSize(0.07);
h_data_p->GetXaxis()->SetTitle("E_{p'} (GeV)");
h_data_p->GetYaxis()->SetTitle("counts");
h_data_p->GetXaxis()->SetTitleSize(0.08);
h_data_p->GetXaxis()->SetTitleOffset(0.8);
h_data_p->GetYaxis()->SetTitleSize(0.08);
h_data_p->GetYaxis()->SetTitleOffset(0.95);
h_data_p->SetTitle("");
h_data_p->SetMaximum((h_data_p->GetMaximum())*1.5);
h_data_p->Draw();



h_sim_p->Scale((h_data_p->Integral())/(h_sim_p->Integral()));
h_sim_p->SetLineColor(2);
h_sim_p->SetLineWidth(2);
//h_sim_p->Draw("same");

TLine *ar_p = new TLine(-0.05+0.938272,(h_data_p->GetMaximum())/2.,-0.05+0.938272,0.);
ar_p->SetLineColor(1);
ar_p->SetLineWidth(2);
//ar_p->SetFillColor(kGreen+2);
ar_p->Draw();











}
