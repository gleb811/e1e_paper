void plot_ratio() {

short ii,kk;
Int_t nbins[35];
Float_t  xsect[35];
Float_t err[35];
TCanvas *c2 = new TCanvas("all sectors", "all sectors", 1000, 600);
c2->Divide(3,2);

for (Short_t sector=1; sector<=6; sector++) {

for (kk=0; kk<=34; kk++) {
nbins[kk] = 0;
xsect[kk] = 0.;
err[kk] = 0.;
};

gStyle->SetPadBottomMargin(0.2);
gStyle->SetPadLeftMargin(0.2);
gStyle->SetOptStat(0);

TFile *MyFile = new TFile("out_ltime.root","READ");

TCanvas *c1 = new TCanvas("elast ration", "elast ratio", 1000, 800);
c1->Divide(5,6);
std::string histName;
std::stringstream out;

out << "sector" << sector;

c1->SetTitle(out.str().c_str());
out.str("");

out << "sector" << sector << "/ratio";

MyFile->cd(out.str().c_str());
TLine *line[30];
TLine *lineup[30];
TLine *linedown[30];
 for (ii=1; ii<31; ii++) {
 c1->cd(ii);
 out.str("");
out << "sector" << sector <<"/ratio/ratio_sector" << sector << "phi_" << ii;
histName = out.str();
cout << histName << "\n";
TH1F *h1 = (TH1F*)MyFile->Get(histName.c_str());
line[ii-1] = new TLine(16.,1.,50.,1.);
lineup[ii-1] = new TLine(16.,1.10,50.,1.10);
linedown[ii-1] = new TLine(16.,0.9,50.,0.9);
h1->SetMaximum(1.4);
h1->SetMinimum(0.6);
h1->SetLabelSize(0.1,"x");
h1->SetLabelSize(0.1,"y");
h1->GetXaxis()->SetTitle("#theta_{e}, deg");
h1->GetXaxis()->SetTitleSize(0.1);
h1->GetXaxis()->SetTitleOffset(0.85);
h1->GetYaxis()->SetTitle("d#sigma/d#theta_{e}, #mubn/rad");
h1->GetYaxis()->SetTitleSize(0.1);
h1->GetYaxis()->SetTitleOffset(0.75);

//cout  << h1->GetNbinsX() << "\n";

for (kk=1; kk<=h1->GetNbinsX(); kk++) {

//cout << h1->GetBinContent(kk) << "\n";

if (( h1->GetBinContent(kk) < 1.05) && ( h1->GetBinContent(kk) > 0.95)) {
//if (( h1->GetBinContent(kk) < 2.) && ( h1->GetBinContent(kk) > 0.)) {
nbins[kk-1] = nbins[kk-1]+1;
xsect[kk-1] = xsect[kk-1] + h1->GetBinContent(kk);
err[kk-1] =err[kk-1] + pow(((h1->GetBinError(kk))/(h1->GetBinContent(kk))),2);
};

};

 out.str("");
out << ii*2-31 <<  " < #phi_{e} < " << ii*2-29 << ", deq";
h1->SetTitle(out.str().c_str());
gStyle->SetTitleSize(0.1,"t");
out.str("");
h1->Draw("e");
line[ii-1]->SetLineColor(2);
line[ii-1]->SetLineWidth(3);
line[ii-1]->Draw("same");
lineup[ii-1]->SetLineColor(3);
lineup[ii-1]->SetLineWidth(3);
lineup[ii-1]->Draw("same");
linedown[ii-1]->SetLineColor(3);
linedown[ii-1]->SetLineWidth(3);
linedown[ii-1]->Draw("same");


}; 

TH1F *h2;
h2 = new TH1F("avrg xsect","avrg xsect",35, 15.5, 50.5);
for (kk=0; kk<=34; kk++) {

if (nbins[kk] != 0) {
cout << xsect[kk]/nbins[kk] << "     " << (xsect[kk]/nbins[kk])*(sqrt(err[kk])/nbins[kk])  << "\n"; 
h2->SetBinContent(kk+1, xsect[kk]/nbins[kk]);
h2->SetBinError(kk+1, (xsect[kk]/nbins[kk])*(sqrt(err[kk])/nbins[kk]));
} else {
cout <<"0." << "\n";
};
 
};

out.str("");
out << "sector " << sector;
//TCanvas *c2 = new TCanvas(out.str().c_str(), out.str().c_str(), 500, 400);
c2->cd(sector);
c2->cd(sector)->SetBottomMargin(0.2);
c2->cd(sector)->SetLeftMargin(0.25);
c2->cd(sector)->SetRightMargin(0.01);
c2->cd(sector)->SetTopMargin(0.02);
h2->SetTitle(out.str().c_str());
h2->SetMaximum(1.1);
h2->SetMinimum(0.9);
h2->SetLabelSize(0.06,"x");
h2->SetLabelSize(0.06,"y");
h2->GetXaxis()->SetTitle("#theta_{e'} (deg)");
h2->GetXaxis()->SetTitleSize(0.08);
h2->GetXaxis()->SetTitleOffset(0.85);

h2->GetYaxis()->SetTitle("#frac{(d#sigma/d#theta_{e'})_{measured}}{(d#sigma/d#theta_{e'})_{parameterization}}");
h2->GetYaxis()->SetTitleSize(0.06);
h2->GetYaxis()->SetTitleOffset(1.75);
h2->Draw("e");

TLine *line_avrg;
TLine *lineup_avrg];
TLine *linedown_avrg;
line_avrg = new TLine(16.,1.,50.,1.);
lineup_avrg = new TLine(16.,1.050,50.,1.05);
linedown_avrg = new TLine(16.,0.95,50.,0.95);
line_avrg->SetLineColor(2);
line_avrg->SetLineWidth(3);
line_avrg->Draw("same");
lineup_avrg->SetLineColor(3);
lineup_avrg->SetLineWidth(3);
lineup_avrg->Draw("same");
linedown_avrg->SetLineColor(3);
linedown_avrg->SetLineWidth(3);
linedown_avrg->Draw("same");

};

c2->Print("elast_ratio_6_sectors.pdf");

};
