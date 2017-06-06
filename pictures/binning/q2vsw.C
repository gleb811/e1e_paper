void q2vsw() {

gStyle->SetPalette(kGreyScale);


UInt_t Number = 50;
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









TFile *MyFile = new TFile("q2vsw.root","READ");
MyFile->ls();
gStyle->SetOptLogz(1);
//gStyle->SetPalette(8);
gStyle->SetOptStat(0);
gStyle->SetOptFit(0);
gStyle->SetOptTitle(0);
gStyle->SetTitleBorderSize(0);
 
TCanvas *c1 = new TCanvas("Q^{2} vs W", "Q^{2} vs W", 800, 500);
c1->SetFillColor(0);
c1->SetBorderMode(0);	
c1->cd(1)->SetBottomMargin(0.22);
c1->cd(1)->SetTopMargin(0.05);
c1->cd(1)->SetLeftMargin(0.13);
c1->cd(1)->SetRightMargin(0.12);
c1->cd(1)->SetFrameBorderMode(0);

//c1->SetTopMargin(-20.);
c1->SetGrid(1,1);
gStyle->SetGridStyle(1);
gStyle->SetGridWidth(2);

TH2F *h1 = (TH2F*)MyFile->Get("h111");
h1->Draw("colz");

h1->GetXaxis()->SetNdivisions(-528);
h1->GetYaxis()->SetNdivisions(-520);
h1->GetXaxis()->SetLabelSize(0.);
h1->GetYaxis()->SetLabelSize(0.045);
h1->GetXaxis()->SetTitle("W (GeV)");
h1->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
h1->GetXaxis()->SetTitleSize(0.08);
h1->GetXaxis()->SetTitleOffset(1.2);
h1->GetYaxis()->SetTitleSize(0.08);
h1->GetYaxis()->SetTitleOffset(0.8);
h1->GetZaxis()->SetLabelSize(0.06);

   const Int_t nx = 29;
   char *months[nx] = {"1.2","1.225","1.25","1.275","1.3","1.325","1.35","1.375","1.4","1.425","1.45","1.475","1.5","1.525","1.55","1.575","1.6","1.625","1.65","1.675","1.7","1.725","1.75","1.775","1.8","1.825","1.85","1.875","1.9"};
		      
   Float_t x, y;
   y = 0.28;
   TText t;
   t.SetTextAngle(60);
   t.SetTextSize(0.045);
   t.SetTextAlign(33);
   for (i=0;i<nx;i++) {
      cout << "y=" << y << "\n";
      x = (1.2+0.025*i);
      t.DrawText(x,y,months[i]);
   }

/*		      
TBox *ripani = new TBox(1.4,0.5,1.8,0.55);
ripani->SetFillStyle(3325);
//ripani->SetFillColor(0);
ripani->SetLineWidth(3);
ripani->SetLineColor(0);
ripani->Draw("L");
*/

TLine *line1 = new TLine(1.3,1.0,1.55,1.0);
line1->SetLineWidth(3);
line1->SetLineColor(0);
line1->Draw();

TLine *line2 = new TLine(1.7,0.8,1.7,0.75);
line2->SetLineWidth(3);
line2->SetLineColor(0);
line2->Draw();

TLine *line3 = new TLine(1.7,0.75,1.725,0.75);
line3->SetLineWidth(3);
line3->SetLineColor(0);
line3->Draw();

TLine *line4 = new TLine(1.725,0.75,1.725,0.7);
line4->SetLineWidth(3);
line4->SetLineColor(0);
line4->Draw();

TLine *line5 = new TLine(1.725,0.7,1.75,0.7);
line5->SetLineWidth(3);
line5->SetLineColor(0);
line5->Draw();

TLine *line6 = new TLine(1.75,0.7,1.75,0.65);
line6->SetLineWidth(3);
line6->SetLineColor(0);
line6->Draw();

TLine *line7 = new TLine(1.75,0.65,1.775,0.65);
line7->SetLineWidth(3);
line7->SetLineColor(0);
line7->Draw();

TLine *line8 = new TLine(1.775,0.65,1.775,0.6);
line8->SetLineWidth(3);
line8->SetLineColor(0);
line8->Draw();


TLine *line9 = new TLine(1.775,0.6,1.8,0.6);
line9->SetLineWidth(3);
line9->SetLineColor(0);
line9->Draw();


TLine *line10 = new TLine(1.8,0.6,1.8,0.5);
line10->SetLineWidth(3);
line10->SetLineColor(0);
line10->Draw();

TLine *line10_1 = new TLine(1.8,0.5,1.825,0.5);
line10_1->SetLineWidth(3);
line10_1->SetLineColor(0);
line10_1->Draw();


TLine *line10_2 = new TLine(1.825,0.5,1.825,0.4);
line10_2->SetLineWidth(3);
line10_2->SetLineColor(0);
line10_2->Draw();


TLine *line11 = new TLine(1.825,0.4,1.6,0.4);
line11->SetLineWidth(3);
line11->SetLineColor(0);
line11->Draw();

TLine *line12 = new TLine(1.6,0.4,1.6,0.45);
line12->SetLineWidth(3);
line12->SetLineColor(0);
line12->Draw();

TLine *line13 = new TLine(1.6,0.45,1.3,0.45);
line13->SetLineWidth(3);
line13->SetLineColor(0);
line13->Draw();

TLine *line14 = new TLine(1.3,0.45,1.3,1.0);
line14->SetLineWidth(3);
line14->SetLineColor(0);
line14->Draw();

TLine *line15 = new TLine(1.55,1.0,1.55,0.95);
line15->SetLineWidth(3);
line15->SetLineColor(0);
line15->Draw();

TLine *line16 = new TLine(1.55,0.95,1.6,0.95);
line16->SetLineWidth(3);
line16->SetLineColor(0);
line16->Draw();

TLine *line17 = new TLine(1.6,0.95,1.6,0.9);
line17->SetLineWidth(3);
line17->SetLineColor(0);
line17->Draw();


TLine *line18 = new TLine(1.6,0.9,1.625,0.9);
line18->SetLineWidth(3);
line18->SetLineColor(0);
line18->Draw();

TLine *line19 = new TLine(1.625,0.9,1.625,0.85);
line19->SetLineWidth(3);
line19->SetLineColor(0);
line19->Draw();

TLine *line191 = new TLine(1.625,0.85,1.65,0.85);
line191->SetLineWidth(3);
line191->SetLineColor(0);
line191->Draw();

TLine *line20 = new TLine(1.65,0.85,1.65,0.8);
line20->SetLineWidth(3);
line20->SetLineColor(0);
line20->Draw();

TLine *line21 = new TLine(1.65,0.8,1.7,0.8);
line21->SetLineWidth(3);
line21->SetLineColor(0);
line21->Draw();


/*
TBox *ya1 = new TBox(1.3,0.45,1.45,0.6);
ya1->SetFillStyle(3354);
//ripani->SetFillColor(0);
ya1->SetLineWidth(3);
ya1->SetLineColor(0);
ya1->Draw("L");


TBox *ya2 = new TBox(1.45,0.45,1.5,0.55);
ya2->SetFillStyle(3354);
//ripani->SetFillColor(0);
ya2->SetLineWidth(3);
ya2->SetLineColor(0);
ya2->Draw("L");


TBox *ya3 = new TBox(1.5,0.4,1.525,0.45);
ya3->SetFillStyle(3354);
//ripani->SetFillColor(0);
ya3->SetLineWidth(3);
ya3->SetLineColor(0);
ya3->Draw("L");
*/

c1->Print("q2vsw.pdf");





}
