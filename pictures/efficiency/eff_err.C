void eff_err() {

TH2F *h_eff_err;

TFile *file_eff_err = new TFile("eff_err.root","READ");

file_eff_err->cd();

gDirectory->GetObject("h_eff_err",h_eff_err);

TCanvas *c = new TCanvas("c","c",500,500);

c->cd();

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




gStyle->SetOptStat(0);
h_eff_err->SetTitle("");

c->cd()->SetRightMargin(0.12);
c->cd()->SetBottomMargin(0.12);
c->cd()->SetTopMargin(0.01);
c->cd()->SetLogz();

h_eff_err->GetXaxis()->SetTitle("efficiency");
h_eff_err->GetXaxis()->SetTitleSize(0.05);
h_eff_err->GetXaxis()->SetTitleOffset(0.8);
h_eff_err->GetYaxis()->SetTitle("relative efficiency error");
h_eff_err->GetYaxis()->SetTitleOffset(0.8);
h_eff_err->GetYaxis()->SetTitleSize(0.05);
h_eff_err->SetMaximum(50);

h_eff_err->Draw("colz");

TPaletteAxis *palette = (TPaletteAxis*)h_eff_err->GetListOfFunctions()->FindObject("palette");
   palette->SetX1NDC(0.9);
   palette->SetX2NDC(0.94);
   palette->SetY1NDC(0.12);
   palette->SetY2NDC(0.99);

TLine *line1 = new TLine(0.001,0.3,0.2,0.3);
line1->SetLineColor(1);
line1->SetLineWidth(3);
line1->Draw();
   c->Modified();
   c->Update();


c->Print("eff_err.pdf");




};
