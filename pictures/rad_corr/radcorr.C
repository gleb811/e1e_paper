#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>



void radcorr() {
gStyle->SetOptStat(0);
Int_t i;
Float_t R;

TH1F *avrg = new TH1F("avrg","avrg",22,1.3,1.85);
TH1F *avrg_unit = new TH1F("avrg_unit","avrg_unit",22,1.3,1.85);


TFile *file = new TFile("no_rad.root");


TH2F* q2vsw_norad = (TH2F*)file->Get("Q2vsW"); 
TH2F* q2vsw_norad_1 = (TH2F*)file->Get("Q2vsW_1"); 
Float_t err;
for (Int_t i=1; i<=q2vsw_norad_1->GetNbinsY(); i++){
if (i==1){
q2vsw_norad_1->SetBinContent(1,i,0);
q2vsw_norad_1->SetBinContent(2,i,0);
q2vsw_norad_1->SetBinContent(3,i,0);
q2vsw_norad_1->SetBinContent(4,i,0);
q2vsw_norad_1->SetBinContent(5,i,0);
q2vsw_norad_1->SetBinContent(6,i,0);
q2vsw_norad_1->SetBinContent(7,i,0);
q2vsw_norad_1->SetBinContent(8,i,0);
q2vsw_norad_1->SetBinContent(22,i,0);
};
if (i==2)q2vsw_norad_1->SetBinContent(22,i,0);
if (i==3){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
};
if (i==4){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
};
if (i==5){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
};
if (i==6){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
};
if (i==7){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
q2vsw_norad_1->SetBinContent(18,i,0);
};
if (i==8){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
q2vsw_norad_1->SetBinContent(18,i,0);
q2vsw_norad_1->SetBinContent(17,i,0);
};
if (i==9){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
q2vsw_norad_1->SetBinContent(18,i,0);
q2vsw_norad_1->SetBinContent(17,i,0);
q2vsw_norad_1->SetBinContent(16,i,0);
q2vsw_norad_1->SetBinContent(15,i,0);
};
if (i==10){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
q2vsw_norad_1->SetBinContent(18,i,0);
q2vsw_norad_1->SetBinContent(17,i,0);
q2vsw_norad_1->SetBinContent(16,i,0);
q2vsw_norad_1->SetBinContent(15,i,0);
q2vsw_norad_1->SetBinContent(14,i,0);
};
if (i==11){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
q2vsw_norad_1->SetBinContent(18,i,0);
q2vsw_norad_1->SetBinContent(17,i,0);
q2vsw_norad_1->SetBinContent(16,i,0);
q2vsw_norad_1->SetBinContent(15,i,0);
q2vsw_norad_1->SetBinContent(14,i,0);
q2vsw_norad_1->SetBinContent(13,i,0);
};
if (i==12){
q2vsw_norad_1->SetBinContent(22,i,0);
q2vsw_norad_1->SetBinContent(21,i,0);
q2vsw_norad_1->SetBinContent(20,i,0);
q2vsw_norad_1->SetBinContent(19,i,0);
q2vsw_norad_1->SetBinContent(18,i,0);
q2vsw_norad_1->SetBinContent(17,i,0);
q2vsw_norad_1->SetBinContent(16,i,0);
q2vsw_norad_1->SetBinContent(15,i,0);
q2vsw_norad_1->SetBinContent(14,i,0);
q2vsw_norad_1->SetBinContent(13,i,0);
q2vsw_norad_1->SetBinContent(12,i,0);
q2vsw_norad_1->SetBinContent(11,i,0);
};

};

for (Int_t i=1; i<=q2vsw_norad->GetNbinsX(); i++){
for (Int_t j=1; j<=q2vsw_norad->GetNbinsY(); j++){
if ((q2vsw_norad_1->GetBinContent(i,j)) > 0) {
err = 1./(sqrt(q2vsw_norad_1->GetBinContent(i,j)));
err = err*(q2vsw_norad->GetBinContent(i,j));
q2vsw_norad->SetBinError(i,j,err);
} else {
q2vsw_norad->SetBinError(i,j,0.);
q2vsw_norad->SetBinContent(i,j,0.);
};


};
};


TFile *file1 = new TFile("rad.root");


//TH2F* q2vsw_rad = (TH2F*)file1->Get("h105");
 TH2F* q2vsw_rad = (TH2F*)file1->Get("Q2vsW");
 TH2F* q2vsw_rad_1 = (TH2F*)file1->Get("Q2vsW_1");
for (Int_t i=1; i<=q2vsw_rad->GetNbinsX(); i++){
for (Int_t j=1; j<=q2vsw_rad->GetNbinsY(); j++){

if ((q2vsw_rad_1->GetBinContent(i,j)) > 20) {
err = 1./(sqrt(q2vsw_rad_1->GetBinContent(i,j)));
err = err*(q2vsw_rad->GetBinContent(i,j));
q2vsw_rad->SetBinError(i,j,err);
} else {
q2vsw_rad->SetBinError(i,j,0.);
q2vsw_rad->SetBinContent(i,j,0.);
};


};
};



//q2vsw_rad->Scale(1./(5250000));
//q2vsw_rad->Scale(1./(6237.*250000.));

q2vsw_norad->Divide(q2vsw_rad);

//q2vsw_norad->Draw("LEGO2");

TCanvas *c = new TCanvas("c","c",500,500);
TCanvas *c_avrg = new TCanvas("c_avrg","c_avrg",500,500);
c->cd();
//c->cd()->SetBottomMargin(0.2);
c->cd()->SetLeftMargin(0.14);
c->cd()->SetRightMargin(0.);
c->cd()->SetTopMargin(0.05);
TH1 *odn = q2vsw_norad->ProjectionX("odn",2,2);

odn->SetMarkerStyle(20);
//odn->Scale(0.1);
odn->SetMinimum(0.85);
odn->SetMaximum(1.35);
odn->SetLineColor(kBlack);
odn->SetMarkerColor(kBlack);
odn->SetTitle(" ");

odn->GetYaxis()->SetLabelSize(0.04);
odn->GetXaxis()->SetLabelSize(0.04);
odn->GetXaxis()->SetTitle("W (GeV)");
odn->GetYaxis()->SetTitle("1/R");
odn->GetYaxis()->SetTitleOffset(1.);
odn->GetXaxis()->SetTitleOffset(0.7);
odn->GetXaxis()->SetNdivisions(10);
odn->GetYaxis()->SetNdivisions(15);
odn->GetXaxis()->SetTitleSize(0.06);
odn->GetYaxis()->SetTitleSize(0.06);


odn->Draw("E1PX0");



//TH1 *odnq = q2vsw_rad->ProjectionX("odnq",2,2);
//odnq->SetMarkerStyle(20);
//odnq->SetMarkerColor(kRed);
//odnq->Draw("E1PX0 same");

TH1 *odn11 = q2vsw_norad->ProjectionX("odn11",1,1);
odn11->SetMarkerStyle(20);
odn11->SetLineColor(11);
odn11->SetMarkerColor(11);
odn11->Draw("E1PX0 same");

TH1 *odn12 = q2vsw_norad->ProjectionX("odn12",3,3);
odn12->SetMarkerStyle(20);
odn12->SetLineColor(13);
odn12->SetMarkerColor(13);
odn12->Draw("E1PX0 same");

TH1 *odn13 = q2vsw_norad->ProjectionX("odn13",5,5);
odn13->SetMarkerStyle(20);
odn13->SetLineColor(25);
odn13->SetMarkerColor(25);
odn13->Draw("E1PX0 same");

TH1 *odn14 = q2vsw_norad->ProjectionX("odn14",7,7);
odn14->SetMarkerStyle(20);
odn14->SetLineColor(35);
odn14->SetMarkerColor(35);
odn14->Draw("E1PX0 same");

TH1 *odn15 = q2vsw_norad->ProjectionX("odn15",9,9);
odn15->SetMarkerStyle(20);
odn15->SetLineColor(52);
odn15->SetMarkerColor(52);
odn15->Draw("E1PX0 same");


TH1 *odn2 = q2vsw_norad->ProjectionX("odn2",4,4);
odn2->SetMarkerStyle(20);
odn2->SetLineColor(kBlue);
odn2->SetMarkerColor(kBlue);
odn2->Draw("E1PX0 same");


TH1 *odn3 = q2vsw_norad->ProjectionX("odn3",6,6);
odn3->SetMarkerStyle(20);
odn3->SetLineColor(kGreen);
odn3->SetMarkerColor(kGreen);
odn3->Draw("E1PX0 same");

TH1 *odn4 = q2vsw_norad->ProjectionX("odn4",8,8);
odn4->SetMarkerStyle(20);
odn4->SetLineColor(kRed);
odn4->SetMarkerColor(kRed);
odn4->Draw("E1PX0 same");


TH1 *odn5 = q2vsw_norad->ProjectionX("odn5",10,10);
odn5->SetMarkerStyle(20);
odn5->SetLineColor(6);
odn5->SetMarkerColor(6);
odn5->Draw("E1PX0 same");


TH1 *odn6 = q2vsw_norad->ProjectionX("odn6",10,10);
odn6->SetMarkerStyle(20);
odn6->SetLineColor(9);
odn6->SetMarkerColor(9);
//odn6->Draw("E1PX0 same");


leg = new TLegend(0.3,0.8,0.99,0.95);
leg->SetNColumns(2);
leg->AddEntry(odn11,"0.4 GeV^{2} < Q^{2} < 0.45 GeV^{2} ","p");
leg->AddEntry(odn,"0.45 GeV^{2} < Q^{2} < 0.5 GeV^{2} ","p");
leg->AddEntry(odn12,"0.5 GeV^{2} < Q^{2} < 0.55 GeV^{2} ","p");
leg->AddEntry(odn2,"0.55 GeV^{2} < Q^{2} < 0.6 GeV^{2} ","p");
leg->AddEntry(odn13,"0.6 GeV^{2} < Q^{2} < 0.65 GeV^{2} ","p");
leg->AddEntry(odn3,"0.65 GeV^{2} < Q^{2} < 0.7 GeV^{2} ","p");
leg->AddEntry(odn14,"0.7 GeV^{2} < Q^{2} < 0.75 GeV^{2} ","p");
leg->AddEntry(odn4,"0.75 GeV^{2} < Q^{2} <  0.8 GeV^{2} ","p");
leg->AddEntry(odn15,"0.8 GeV^{2} < Q^{2} < 0.85 GeV^{2} ","p");
leg->AddEntry(odn5,"0.85 GeV^{2} < Q^{2} < 0.9 GeV^{2} ","p");
leg->Draw();
for (Int_t j=1; j<=22; j++){
avrg->SetBinError(j,0.);
avrg->SetBinContent(j,0.);
avrg_unit->SetBinError(j,0.);
avrg_unit->SetBinContent(j,0.);
};
for (Int_t i=1; i<=q2vsw_norad->GetNbinsX(); i++){
for (Int_t j=1; j<=q2vsw_norad->GetNbinsY(); j++){

if ((q2vsw_norad_1->GetBinContent(i,j)) > 0) {
err = 1./(sqrt(q2vsw_norad_1->GetBinContent(i,j)));
err = err*(q2vsw_norad->GetBinContent(i,j));
q2vsw_norad->SetBinError(i,j,err);
avrg->SetBinContent(i,avrg->GetBinContent(i)+q2vsw_norad->GetBinContent(i,j));
avrg->SetBinError(i,avrg->GetBinError(i)+(q2vsw_norad->GetBinError(i,j))*(q2vsw_norad->GetBinError(i,j)));
avrg_unit->SetBinContent(i,avrg_unit->GetBinContent(i)+1.);
};


};
};

for (Int_t j=1; j<=22; j++){
avrg->SetBinError(j,sqrt(avrg->GetBinError(j)));
avrg_unit->SetBinError(j,0.);
};

avrg->Divide(avrg_unit);
c_avrg->cd();
c_avrg->cd()->SetLeftMargin(0.14);
c_avrg->cd()->SetRightMargin(0.);
c_avrg->cd()->SetTopMargin(0.05);

for (Int_t j=1; j<=22; j++){
cout << avrg->GetBinCenter(j) << "   " << avrg->GetBinContent(j) << "   " << avrg->GetBinError(j)<< endl;
};

avrg->SetMarkerStyle(20);
//odn->Scale(0.1);
avrg->SetMinimum(0.85);
avrg->SetMaximum(1.35);
avrg->SetLineColor(kBlack);
avrg->SetMarkerColor(kBlack);
avrg->SetTitle(" ");

avrg->GetYaxis()->SetLabelSize(0.04);
avrg->GetXaxis()->SetLabelSize(0.04);
avrg->GetXaxis()->SetTitle("W (GeV)");
avrg->GetYaxis()->SetTitle("1/R");
avrg->GetYaxis()->SetTitleOffset(1.);
avrg->GetXaxis()->SetTitleOffset(0.7);
avrg->GetXaxis()->SetNdivisions(10);
avrg->GetYaxis()->SetNdivisions(15);
avrg->GetXaxis()->SetTitleSize(0.06);
avrg->GetYaxis()->SetTitleSize(0.06);




avrg->SetMarkerStyle(20);
avrg->SetMinimum(0.85);
avrg->SetMaximum(1.35);
avrg->Draw("E1PX0");

//for (Int_t j=1; j<=25; j++){


//};
//odn->Draw("P");

c->SaveAs("rad_corr_all_q2.pdf");
c_avrg->SaveAs("rad_corr_avrg.pdf");


};
