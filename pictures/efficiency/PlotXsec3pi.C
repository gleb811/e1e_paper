
#include <string>

using namespace std;

// exp data 
#include "saphir.h"
#include "abbhhm.h"



//_________________________ MAIN ____________________________________

void PlotXsec3pi(){

void TH1::SetDrawMark(int lineCol=kBlack, int markType=kDot, 
                      char *dOpt="", double msize=0.8)
{  
  const std::string op="P";
  this->SetLineColor(lineCol); 
  this->SetMarkerColor(lineCol); 
  this->SetMarkerStyle(markType);
  this->SetMarkerSize(msize);
  const std::string opt = op+dOpt;
  this->DrawCopy( opt.c_str());
  gPad->GetCanvas()->Update();
  return;
};


  gStyle->SetErrorX(0.0001);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(0.5);
  TCanvas *can = new TCanvas("can","PLOT",50,50,300,300); //

  can->cd()->SetRightMargin(0.001);
  can->cd()->SetTopMargin(0.02);
  
  TH1F *hSap1 = GetSaphir1Hist("hSap1");
  TH1F *hSap2 = GetSaphir2Hist("hSap2");
  TH1F *hAbb2 = GetAbbhhmHist2("hAbb2");
  TH1F *hAbb1 = GetAbbhhmHist1("hAbb1");
  TH1F *hSap4 = GetSaphir1Hist_3pi("hSap4");  //abbhhm ??
  
 
  TH2F *fr = new TH2F("fr","", 10,1.3,1.85,100,0.,100.); 
  fr->GetXaxis()->SetTitle("W (GeV)");
  fr->GetYaxis()->SetTitle("#sigma (#mubarn)");
  fr->GetYaxis()->SetTitleOffset(1.3); 
  fr->Draw();
  hSap1->SetMarkerSize(0.5); 
  hSap1->SetDrawMark(kRed,   kFullCircle, "e1 same");
  hSap2->SetDrawMark(kBlack,   kFullCircle, "e1 same");
  hAbb2->SetDrawMark(kGreen, kFullCircle, "e1 same");
  hAbb1->SetDrawMark(kRed+2, kFullCircle,  "e1 same");
  hSap4->SetDrawMark(kBlue,  kFullCircle, "e1 same");

  can->Print("2pi_3pi.pdf");
  // SAPHIR   Eur. Phys. J. A 23, 317 (2005)
  // ABBHHM   Phys. Rev. 175, 1669 (1968)
  // ABBHHM   Phys. Rev. 188, 2060 (1969)  
  
  
  
}
