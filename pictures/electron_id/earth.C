


double testfunc(double *x, double *par) {
   double s;
   s = (-1.*0.0032+sqrt(0.0032*0.0032-4.*0.0499*(7.-x[0])))/2./0.0499;
   return s;
}

double testfunc1(double *x, double *par) {
   double s;
   s = (-1.*0.0032-sqrt(0.0032*0.0032-4.*0.0499*(7.-x[0])))/2./0.0499;
   return s;
}


double testfunc2(double *x, double *par) {
   double s;
   s = sqrt(1. - pow((x[0]-45.5)/34.5,2))*21.;
   return s;
}

double testfunc3(double *x, double *par) {
   double s;
   s = -1.*sqrt(1. - pow((x[0]-45.5)/34.5,2))*21.;
   return s;
}

double testfunc4(double *x, double *par) {
   double s;
   s = sqrt(1. - pow((x[0]-45.5)/1.75,2))*21.;
   return s;
}

double testfunc5(double *x, double *par) {
   double s;
   s = -1.*sqrt(1. - pow((x[0]-45.5)/1.75,2))*21.;
   return s;
}



void earth() {

gStyle->SetOptStat(0);

TH2D *earth1 = new TH2D("earth1", "earth1", 100, -5., 50., 100, -30., 30.);

double phi=0., cos_theta=0., theta=0.;

for(int i=0; i<100; i++){
phi=gRandom->Uniform(-30., 30.);
theta=7.0+0.0032*phi+0.0499*phi*phi;
earth1->Fill(theta, phi);
}
earth1->SetMarkerStyle(20);
TCanvas *c1 = new TCanvas("c1", "c1",1);

  earth1->Draw("L");


   TF1 *f1 = new TF1("f1",testfunc,6.9,45,0);
   f1->Draw("same");

   TF1 *f2 = new TF1("f2",testfunc1,6.9,45,0);
   f2->Draw("same");
   
   TF1 *f3 = new TF1("f3",testfunc2,10.9,45.6,0);
   f3->Draw("same");  
   
   TF1 *f4 = new TF1("f4",testfunc3,10.9,45.6,0);
   f4->Draw("same");    
   
   TF1 *f5 = new TF1("f5",testfunc4,43.741708543,45.4,0);
   f5->Draw("same");        

   TF1 *f6 = new TF1("f6",testfunc5,43.741708543,45.4,0);
   f6->Draw("same");   
}
