void PlotXsec(){
 gStyle->SetTitleX(0.5);
 gStyle->SetTitleY(0.87);
 gStyle->SetTitleAlign(23);
 gStyle->SetOptTitle(0); 
 gStyle->SetTitleBorderSize(0);
  
 cout << "\nPlotXsec started.\n";
 
  TCanvas *can = new TCanvas("can","PLqOT",50,50,450,450); //
  can->SetFillColor(0);
  //can->SetWindowSize(1000,700);
  //can->SetCanvasSize(950,650);
  can->Divide(2,2);
  can->SetBorderMode(0);
  can->cd(1);
  can->cd(1)->SetBorderMode(0);
  can->cd(1)->SetFrameBorderMode(0);
can->cd(1)->SetBottomMargin(0.17);
can->cd(1)->SetTopMargin(0.03);
can->cd(1)->SetLeftMargin(0.15);
can->cd(1)->SetRightMargin(0.05);
can->cd(1)->SetFrameBorderMode(0);  
int j,i;

std::ostringstream file;
std::ostringstream title;






 

  TMarker *mark = new TMarker();
   
;


std::ostringstream file3;   

file3 << "e1c.dat";
 cout << file3.str() << "\n";



int i;
 float sigma4[30],error4[30],theta4[30],errorx4[30];
 
 

std::ifstream sect4(file3.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 9; i++)
{
  sect4 >> theta4[i] >> sigma4[i] >> error4[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta4[i] << "  " << sigma4[i] << "  " << error4[i] << "\n";

}
 sect4.close();   
 
   TGraphErrors *Tg3 = new TGraphErrors(9,theta4,sigma4,errorx4,error4);
   Tg3->SetMarkerStyle(33);
   Tg3->SetMarkerColor(1);
   Tg3->SetMarkerSize(1.5);
   Tg3->SetLineColor(1); 
   Tg3->SetMaximum(33.);
   Tg3->GetXaxis()->SetLimits(1.3,1.9);
   Tg3->SetMinimum(0.);
   
   Tg3->SetTitle(title.str().c_str());
   Tg3->Draw("AP");
   Tg3->GetXaxis()->SetTitle("W (GeV)");
   Tg3->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg3->GetXaxis()->SetTitleSize(0.09);
   Tg3->GetYaxis()->SetTitleSize(0.09); 
   Tg3->GetYaxis()->SetTitleOffset(0.85);
   Tg3->GetXaxis()->SetTitleOffset(0.85);   
    Tg3->GetXaxis()->SetLabelSize(0.07);  
    Tg3->GetYaxis()->SetLabelSize(0.07);  
            
  TMarker *mark1 = new TMarker();
   
   mark1->SetMarkerStyle(33);
   mark1->SetMarkerColor(1);
   mark1->SetMarkerSize(1.5);
   mark1->DrawMarker(1.48,9.);

   TLatex t1;
//   t.SetTextAlign(33);
   t1.SetTextSize(0.06);
   t1.DrawLatex(1.52,8.2,"e1c Q^{2} = 0.425 GeV^{2}");
   
   
  TMarker *mark2 = new TMarker();
   
   mark2->SetMarkerStyle(24);
   mark2->SetMarkerColor(1);
   mark2->SetMarkerSize(0.9);
   mark2->DrawMarker(1.48,6.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.06);
   t2.DrawLatex(1.52,5.2,"e1e Q^{2} = 0.425 GeV^{2}");   
   
   
   

















title << " ";
file << "/home/gleb/e1e/e1e/2picross-sect/q204045/jopa";
// cout << file.str() << "\n";



int i;
 float sigma[30],error[30],theta[30],errorx[30];
 
 

std::ifstream sect1(file.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect1 >> theta[i] >> sigma[i] >> error[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta[i] << "  " << sigma[i] << "  " << error[i] << "\n";

}
 sect1.close();


   

   
  TGraphErrors *Tg = new TGraphErrors(20,theta,sigma,errorx,error);
   
   Tg->SetMarkerStyle(24);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg->SetTitle(title.str().c_str());
   Tg->Draw("PS");
   Tg->SetMarkerColor(1);
   Tg->SetLineColor(1);
   Tg->SetMarkerSize(0.9);
   Tg->GetXaxis()->SetTitle("W (GeV)");
   Tg->GetXaxis()->SetTitleSize(0.06);
   Tg->GetXaxis()->SetTitleOffset(0.7);
 
    
std::ostringstream file1;   
file1 << "/home/gleb/e1e/e1e/2picross-sect/q206507/jopa";
 cout << file.str() << "\n";



int i;
 float sigma2[30],error2[30],theta2[30],errorx2[30];
 
 

std::ifstream sect2(file1.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect2 >> theta2[i] >> sigma2[i] >> error2[i];
  

 cout << theta2[i] << "  " << sigma2[i] << "  " << error2[i] << "\n";
// sigma2[i]=sigma2[i]*2.0;

}
 sect2.close();   
 

   


   
/*   TLine *line1 = new TLine(17.,1.05,52.,1.05);
   line1->SetLineColor(kRed);
   line1->SetLineWidth(3); 
   line1->Draw(); 
   TLine *line1 = new TLine(17.,0.95,52.,0.95);
   line1->SetLineWidth(3); 
   line1->SetLineColor(kRed);
   line1->Draw();
   
*/
   
//   delete line1;
       

file.str("");
title.str("");














































can->cd(2);
can->cd(2)->SetBorderMode(0);
can->cd(2)->SetFrameBorderMode(0);
can->cd(2)->SetBottomMargin(0.17);
can->cd(2)->SetTopMargin(0.03);
can->cd(2)->SetLeftMargin(0.15);
can->cd(2)->SetRightMargin(0.05);

  
std::ostringstream file;
std::ostringstream title;





/* std::ostringstream file2;   
file2 << "ripani.dat";
 cout << file2.str() << "\n";




 float sigma3[30],error3[30],theta3[30],errorx3[30];
 
 

std::ifstream sect3(file2.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect3 >> theta3[i] >> sigma3[i] >> error3[i] >> errorx3[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta3[i] << "  " << sigma3[i] << "  " << error3[i] << "\n";

}
 sect3.close();   
 
   TGraphErrors *Tg2 = new TGraphErrors(20,theta3,sigma3,errorx3,error3);
   
   Tg2->SetMarkerStyle(kFullCircle);
   Tg2->SetMarkerColor(2);
   Tg2->SetLineColor(2); 
   Tg2->SetMaximum(33.);
   Tg2->GetXaxis()->SetLimits(1.3,1.9);
   Tg2->SetMinimum(0.);
   
   Tg2->SetTitle(title.str().c_str());
   Tg2->Draw("AP");  */


std::ostringstream file3;   

file3 << "e1c04505.dat";
 cout << file3.str() << "\n";



int i;
 float sigma4[30],error4[30],theta4[30],errorx4[30];
 
 

std::ifstream sect4(file3.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 9; i++)
{
  sect4 >> theta4[i] >> sigma4[i] >> error4[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta4[i] << "  " << sigma4[i] << "  " << error4[i] << "\n";

}
 sect4.close();   
 
   TGraphErrors *Tg3 = new TGraphErrors(8,theta4,sigma4,errorx4,error4);
   
   Tg3->SetMarkerStyle(33);
   Tg3->SetMarkerColor(1);
   Tg3->SetMarkerSize(1.5);
   Tg3->SetLineColor(1); 
   Tg3->SetMaximum(33.);
   Tg3->GetXaxis()->SetLimits(1.3,1.9);
   Tg3->SetMinimum(0.);
   Tg3->GetXaxis()->SetTitle("W (GeV)");
   Tg3->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg3->GetXaxis()->SetTitleSize(0.09);
   Tg3->GetYaxis()->SetTitleSize(0.09); 
   Tg3->GetYaxis()->SetTitleOffset(0.85);
   Tg3->GetXaxis()->SetTitleOffset(0.85);    
   Tg3->GetXaxis()->SetLabelSize(0.07); 
   Tg3->GetYaxis()->SetLabelSize(0.07);     
   Tg3->SetTitle(title.str().c_str());
   Tg3->Draw("AP");
   
   
   
   
  
  TMarker *mark = new TMarker();
   
/*   mark->SetMarkerStyle(kFullCircle);
   mark->SetMarkerColor(2);
   mark->DrawMarker(1.8,30.);

   TLatex t;
//   t.SetTextAlign(33);
   t.SetTextSize(0.02);
   t.DrawLatex(1.808,29.7,"M. Ripani 0.65 GeV^{2}");

  TMarker *mark1 = new TMarker(); */

   
   mark1->SetMarkerStyle(33);
   mark1->SetMarkerColor(1);
   mark1->SetMarkerSize(1.5);
   mark1->DrawMarker(1.48,9.);

   TLatex t1;
//   t.SetTextAlign(33);
   t1.SetTextSize(0.06);
   t1.DrawLatex(1.52,8.2,"e1c Q^{2} = 0.475 GeV^{2}");
   
   
  TMarker *mark2 = new TMarker();
   
   mark2->SetMarkerStyle(24);
   mark2->SetMarkerColor(1);
   mark2->SetMarkerSize(0.9);
   mark2->DrawMarker(1.48,6.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.06);
   t2.DrawLatex(1.52,5.2,"e1e Q^{2} = 0.475 GeV^{2}");   
   
   
/*  
  TMarker *mark3 = new TMarker();
   
   mark3->SetMarkerStyle(kFullCircle);
   mark3->SetMarkerColor(1);
   mark3->DrawMarker(1.8,24.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.02);
   t2.DrawLatex(1.808,23.7,"e1e 0.675 GeV^{2}"); 



*/












title << " ";
file << "/home/gleb/e1e/e1e/2picross-sect/q204505/jopa";
// cout << file.str() << "\n";



int i;
 float sigma[30],error[30],theta[30],errorx[30];
 
 

std::ifstream sect1(file.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect1 >> theta[i] >> sigma[i] >> error[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta[i] << "  " << sigma[i] << "  " << error[i] << "\n";

}
 sect1.close();


   

   
  TGraphErrors *Tg = new TGraphErrors(20,theta,sigma,errorx,error);
   
   Tg->SetMarkerStyle(24);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg->SetTitle(title.str().c_str());
   Tg->Draw("PS");
   Tg->SetMarkerColor(1);
   Tg->SetLineColor(1);
   Tg->SetMarkerSize(0.9);
   Tg->GetXaxis()->SetTitle("W (GeV)");
   Tg->GetXaxis()->SetTitleSize(0.06);
   Tg->GetXaxis()->SetTitleOffset(0.7);

/*   
std::ostringstream file1;   
file1 << "/home/gleb/e1e/e1e/2picross-sect/q206507/jopa";
 cout << file.str() << "\n";



int i;
 float sigma2[30],error2[30],theta2[30],errorx2[30];
 
 

std::ifstream sect2(file1.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect2 >> theta2[i] >> sigma2[i] >> error2[i];
  

 cout << theta2[i] << "  " << sigma2[i] << "  " << error2[i] << "\n";
// sigma2[i]=sigma2[i]*2.0;

}
 sect2.close();   
 
   TGraphErrors *Tg1 = new TGraphErrors(20,theta2,sigma2,errorx2,error2);
   
   Tg1->SetMarkerStyle(kFullCircle);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg1->SetTitle(title.str().c_str());
   Tg1->Draw("PS");
   

*/
   
/*   TLine *line1 = new TLine(17.,1.05,52.,1.05);
   line1->SetLineColor(kRed);
   line1->SetLineWidth(3); 
   line1->Draw(); 
   TLine *line1 = new TLine(17.,0.95,52.,0.95);
   line1->SetLineWidth(3); 
   line1->SetLineColor(kRed);
   line1->Draw();
   
*/
   
//   delete line1;
       

file.str("");
title.str("");





















































































can->cd(3);
can->cd(3)->SetBorderMode(0);
can->cd(3)->SetFrameBorderMode(0);
can->cd(3)->SetBottomMargin(0.17);
can->cd(3)->SetTopMargin(0.03);
can->cd(3)->SetLeftMargin(0.15);
can->cd(3)->SetRightMargin(0.05);  
  
  
std::ostringstream file;
std::ostringstream title;





/* std::ostringstream file2;   
file2 << "ripani.dat";
 cout << file2.str() << "\n";




 float sigma3[30],error3[30],theta3[30],errorx3[30];
 
 

std::ifstream sect3(file2.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect3 >> theta3[i] >> sigma3[i] >> error3[i] >> errorx3[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta3[i] << "  " << sigma3[i] << "  " << error3[i] << "\n";

}
 sect3.close();   
 
   TGraphErrors *Tg2 = new TGraphErrors(20,theta3,sigma3,errorx3,error3);
   
   Tg2->SetMarkerStyle(kFullCircle);
   Tg2->SetMarkerColor(2);
   Tg2->SetLineColor(2); 
   Tg2->SetMaximum(33.);
   Tg2->GetXaxis()->SetLimits(1.3,1.9);
   Tg2->SetMinimum(0.);
   
   Tg2->SetTitle(title.str().c_str());
   Tg2->Draw("AP");  */


std::ostringstream file3;   

file3 << "e1c05055.dat";
 cout << file3.str() << "\n";



int i;
 float sigma4[30],error4[30],theta4[30],errorx4[30];
 
 

std::ifstream sect4(file3.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 9; i++)
{
  sect4 >> theta4[i] >> sigma4[i] >> error4[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta4[i] << "  " << sigma4[i] << "  " << error4[i] << "\n";

}
 sect4.close();   
 
   TGraphErrors *Tg3 = new TGraphErrors(8,theta4,sigma4,errorx4,error4);
   
   Tg3->SetMarkerStyle(33);
   Tg3->SetMarkerColor(1);
   Tg3->SetMarkerSize(1.5);   
   Tg3->SetLineColor(1); 
   Tg3->SetMaximum(33.);
   Tg3->GetXaxis()->SetLimits(1.3,1.9);
   Tg3->SetMinimum(0.);
      Tg3->GetXaxis()->SetTitle("W (GeV)");
   Tg3->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg3->GetXaxis()->SetTitleSize(0.09);
   Tg3->GetYaxis()->SetTitleSize(0.09); 
   Tg3->GetYaxis()->SetTitleOffset(0.85);
   Tg3->GetXaxis()->SetTitleOffset(0.85); 
   Tg3->GetXaxis()->SetLabelSize(0.07); 
   Tg3->GetYaxis()->SetLabelSize(0.07);     
   Tg3->SetTitle(title.str().c_str());
   Tg3->Draw("AP");
   
   
   
   
  
  TMarker *mark = new TMarker();
   
/*   mark->SetMarkerStyle(kFullCircle);
   mark->SetMarkerColor(2);
   mark->DrawMarker(1.8,30.);

   TLatex t;
//   t.SetTextAlign(33);
   t.SetTextSize(0.02);
   t.DrawLatex(1.808,29.7,"M. Ripani 0.65 GeV^{2}");

  TMarker *mark1 = new TMarker(); */
   
   mark1->SetMarkerStyle(33);
   mark1->SetMarkerColor(1);
   mark1-> SetMarkerSize(1.5);  
   mark1->DrawMarker(1.48,9.);

   TLatex t1;
//   t.SetTextAlign(33);
   t1.SetTextSize(0.06);
   t1.DrawLatex(1.52,8.2,"e1c Q^{2} = 0.525 GeV^{2}");
   
   
  TMarker *mark2 = new TMarker();
   
   mark2->SetMarkerStyle(24);
   mark2->SetMarkerColor(1);
   mark2->SetMarkerSize(0.9);
   mark2->DrawMarker(1.48,6.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.06);
   t2.DrawLatex(1.52,5.2,"e1e Q^{2} = 0.525 GeV^{2}");   
   
   
/*  
  TMarker *mark3 = new TMarker();
   
   mark3->SetMarkerStyle(kFullCircle);
   mark3->SetMarkerColor(1);
   mark3->DrawMarker(1.8,24.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.02);
   t2.DrawLatex(1.808,23.7,"e1e 0.675 GeV^{2}"); 



*/












title << " ";
file << "/home/gleb/e1e/e1e/2picross-sect/q205055/jopa";
// cout << file.str() << "\n";



int i;
 float sigma[30],error[30],theta[30],errorx[30];
 
 

std::ifstream sect1(file.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect1 >> theta[i] >> sigma[i] >> error[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta[i] << "  " << sigma[i] << "  " << error[i] << "\n";

}
 sect1.close();


   

   
  TGraphErrors *Tg = new TGraphErrors(20,theta,sigma,errorx,error);
   
   Tg->SetMarkerStyle(24);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg->SetTitle(title.str().c_str());
   Tg->Draw("PS");
   Tg->SetMarkerColor(1);
   Tg->SetMarkerSize(0.9);   
   Tg->SetLineColor(1);
   Tg->GetXaxis()->SetTitle("W (GeV)");
   Tg->GetXaxis()->SetTitleSize(0.06);
   Tg->GetXaxis()->SetTitleOffset(0.7);

/*   
std::ostringstream file1;   
file1 << "/home/gleb/e1e/e1e/2picross-sect/q206507/jopa";
 cout << file.str() << "\n";



int i;
 float sigma2[30],error2[30],theta2[30],errorx2[30];
 
 

std::ifstream sect2(file1.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect2 >> theta2[i] >> sigma2[i] >> error2[i];
  

 cout << theta2[i] << "  " << sigma2[i] << "  " << error2[i] << "\n";
// sigma2[i]=sigma2[i]*2.0;

}
 sect2.close();   
 
   TGraphErrors *Tg1 = new TGraphErrors(20,theta2,sigma2,errorx2,error2);
   
   Tg1->SetMarkerStyle(kFullCircle);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg1->SetTitle(title.str().c_str());
   Tg1->Draw("PS");
   

*/
   
/*   TLine *line1 = new TLine(17.,1.05,52.,1.05);
   line1->SetLineColor(kRed);
   line1->SetLineWidth(3); 
   line1->Draw(); 
   TLine *line1 = new TLine(17.,0.95,52.,0.95);
   line1->SetLineWidth(3); 
   line1->SetLineColor(kRed);
   line1->Draw();
   
*/
   
//   delete line1;
       

file.str("");
title.str("");

























































































































can->cd(4);
can->cd(4)->SetBorderMode(0);
can->cd(4)->SetFrameBorderMode(0);
can->cd(4)->SetBottomMargin(0.17);
can->cd(4)->SetTopMargin(0.03);
can->cd(4)->SetLeftMargin(0.15);
can->cd(4)->SetRightMargin(0.05);  
  
std::ostringstream file;
std::ostringstream title;





/* std::ostringstream file2;   
file2 << "ripani.dat";
 cout << file2.str() << "\n";




 float sigma3[30],error3[30],theta3[30],errorx3[30];
 
 

std::ifstream sect3(file2.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect3 >> theta3[i] >> sigma3[i] >> error3[i] >> errorx3[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta3[i] << "  " << sigma3[i] << "  " << error3[i] << "\n";

}
 sect3.close();   
 
   TGraphErrors *Tg2 = new TGraphErrors(20,theta3,sigma3,errorx3,error3);
   
   Tg2->SetMarkerStyle(kFullCircle);
   Tg2->SetMarkerColor(2);
   Tg2->SetLineColor(2); 
   Tg2->SetMaximum(33.);
   Tg2->GetXaxis()->SetLimits(1.3,1.9);
   Tg2->SetMinimum(0.);
   
   Tg2->SetTitle(title.str().c_str());
   Tg2->Draw("AP");  */


std::ostringstream file3;   

file3 << "e1c05506.dat";
 cout << file3.str() << "\n";



int i;
 float sigma4[30],error4[30],theta4[30],errorx4[30];
 
 

std::ifstream sect4(file3.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 9; i++)
{
  sect4 >> theta4[i] >> sigma4[i] >> error4[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta4[i] << "  " << sigma4[i] << "  " << error4[i] << "\n";

}
 sect4.close();   
 
   TGraphErrors *Tg3 = new TGraphErrors(6,theta4,sigma4,errorx4,error4);
   
   Tg3->SetMarkerStyle(33);
   Tg3->SetMarkerColor(1);
   Tg3->SetMarkerSize(1.5);   
   Tg3->SetLineColor(1); 
   Tg3->SetMaximum(33.);
   Tg3->GetXaxis()->SetLimits(1.3,1.9);
   Tg3->SetMinimum(0.);
      Tg3->GetXaxis()->SetTitle("W (GeV)");
   Tg3->GetYaxis()->SetTitle("#sigma (#mubn)");
   Tg3->GetXaxis()->SetTitleSize(0.09);
   Tg3->GetXaxis()->SetTitleOffset(0.85);
   Tg3->GetYaxis()->SetTitleSize(0.09); 
   Tg3->GetYaxis()->SetTitleOffset(0.85); 
   Tg3->GetXaxis()->SetLabelSize(0.07); 
   Tg3->GetYaxis()->SetLabelSize(0.07);     
   Tg3->SetTitle(title.str().c_str());
   Tg3->Draw("AP");
   
   
   
   
  
  TMarker *mark = new TMarker();
   
/*   mark->SetMarkerStyle(kFullCircle);
   mark->SetMarkerColor(2);
   mark->DrawMarker(1.8,30.);

   TLatex t;
//   t.SetTextAlign(33);
   t.SetTextSize(0.02);
   t.DrawLatex(1.808,29.7,"M. Ripani 0.65 GeV^{2}");

  TMarker *mark1 = new TMarker(); */
   
   mark1->SetMarkerStyle(33);
   mark1->SetMarkerSize(1.5);   
   mark1->SetMarkerColor(1);
   mark1->DrawMarker(1.48,9);

   TLatex t1;
//   t.SetTextAlign(33);
   t1.SetTextSize(0.06);
   t1.DrawLatex(1.52,8.2,"e1c Q^{2} = 0.575 GeV^{2}");
   
   
  TMarker *mark2 = new TMarker();
   
   mark2->SetMarkerStyle(24);
   mark2->SetMarkerColor(1);
   mark2->SetMarkerSize(0.9);
   mark2->DrawMarker(1.48,6.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.06);
   t2.DrawLatex(1.52,5.2,"e1e Q^{2} = 0.575 GeV^{2}");   
   
   
/*  
  TMarker *mark3 = new TMarker();
   
   mark3->SetMarkerStyle(kFullCircle);
   mark3->SetMarkerColor(1);
   mark3->DrawMarker(1.8,24.);

   TLatex t2;
//   t.SetTextAlign(33);
   t2.SetTextSize(0.02);
   t2.DrawLatex(1.808,23.7,"e1e 0.675 GeV^{2}"); 



*/












title << " ";
file << "/home/gleb/e1e/e1e/2picross-sect/q205506/jopa";
// cout << file.str() << "\n";



int i;
 float sigma[30],error[30],theta[30],errorx[30];
 
 

std::ifstream sect1(file.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect1 >> theta[i] >> sigma[i] >> error[i];
  
//  sigma[i-20]=sigma[i-20]*2.0;

 cout << theta[i] << "  " << sigma[i] << "  " << error[i] << "\n";

}
 sect1.close();


   

   
  TGraphErrors *Tg = new TGraphErrors(20,theta,sigma,errorx,error);
   
   Tg->SetMarkerStyle(24);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg->SetTitle(title.str().c_str());
   Tg->Draw("PS");
   Tg->SetMarkerColor(1);
   Tg->SetLineColor(1);
      Tg->SetMarkerSize(0.9);
   Tg->GetXaxis()->SetTitle("W (GeV)");
   Tg->GetXaxis()->SetTitleSize(0.06);
   Tg->GetXaxis()->SetTitleOffset(0.7);

/*   
std::ostringstream file1;   
file1 << "/home/gleb/e1e/e1e/2picross-sect/q206507/jopa";
 cout << file.str() << "\n";



int i;
 float sigma2[30],error2[30],theta2[30],errorx2[30];
 
 

std::ifstream sect2(file1.str().c_str());
//  std::ifstream sect1("sect1.dat");

for (i= 0; i <= 19; i++)
{
  sect2 >> theta2[i] >> sigma2[i] >> error2[i];
  

 cout << theta2[i] << "  " << sigma2[i] << "  " << error2[i] << "\n";
// sigma2[i]=sigma2[i]*2.0;

}
 sect2.close();   
 
   TGraphErrors *Tg1 = new TGraphErrors(20,theta2,sigma2,errorx2,error2);
   
   Tg1->SetMarkerStyle(kFullCircle);
//   Tg->SetMaximum(1.2);
//   Tg->SetMinimum(0.8);
   Tg1->SetTitle(title.str().c_str());
   Tg1->Draw("PS");
   

*/
   
/*   TLine *line1 = new TLine(17.,1.05,52.,1.05);
   line1->SetLineColor(kRed);
   line1->SetLineWidth(3); 
   line1->Draw(); 
   TLine *line1 = new TLine(17.,0.95,52.,0.95);
   line1->SetLineWidth(3); 
   line1->SetLineColor(kRed);
   line1->Draw();
   
*/
   
//   delete line1;
       

file.str("");
title.str("");




 return;
}
