

TH1F *GetSaphir1Hist(const char* hName)
{
const double mp = 0.93827;
const int np = 19;
double e1[np] = 
  { 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 
    1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40 
  };
double e2[np] = 
  { 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 
    1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45 
  };
double sig[np] = 
  { 33.85, 55.66, 70.36, 74.98, 75.42, 72.85, 72.52, 74.22, 75.23, 
    72.60, 72.79, 72.59, 64.33, 60.81, 57.63, 54.62, 55.88, 54.94, 
    53.49
  };
double err[np] = 
  { 1.78, 1.84, 1.77, 1.84, 1.90, 1.91, 1.96, 2.09, 2.30, 2.55, 
    2.53, 2.56, 2.56, 2.40, 2.39, 2.21, 2.68, 2.39, 2.53 
  };

double wlow[np+1];
for(int i=0; i<np; i++)
  wlow[i]  = sqrt( mp*(mp+2.*e1[i]) );
  wlow[np] = wlow[np-1] + (wlow[np-1]-wlow[np-2]);
   
TH1F *hist = new TH1F(hName, "saphir34569 tot", np, wlow);
for(int i=1; i<=np; i++) {
  hist->SetBinContent(i, sig[i-1]);
  hist->SetBinError(i, err[i-1]);
}

return hist;
}; // end GetSaphir1Hist()



TH1F *GetSaphir2Hist(const char* hName)
{
const double mp = 0.93827;
const int np = 16;
double e1[np] = 
  { 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 
    2.00, 2.10, 2.20, 2.30, 2.40, 2.50 
  };
double e2[np] = 
  { 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 
    2.10, 2.20, 2.30, 2.40, 2.50, 2.60
  };
double sig[np] = 
  { 74.23, 62.24, 57.09, 54.00, 52.56, 50.06, 47.27, 45.08, 
    42.33, 39.41, 36.63, 32.98, 31.62, 30.14, 29.08, 28.90
  };
double err[np] = 
  { 3.24, 2.17, 2.54, 1.50, 1.52, 2.03, 2.12, 1.92, 1.90, 
    2.01, 1.54, 2.98, 1.74, 1.83, 1.49, 0.29
  };

double wlow[np+1];
for(int i=0; i<np; i++)
  {wlow[i]  = sqrt( mp*(mp+2.*e1[i]) );  }
   wlow[np] = wlow[np-1] + (wlow[np-1]-wlow[np-2]);

TH1F *hist = new TH1F(hName, "saphir2554 tot", np, wlow);
for(int i=1; i<=np; i++) {
  hist->SetBinContent(i, sig[i-1]);
  hist->SetBinError(i, err[i-1]);
}

return hist;
}; // end GetSaphir2Hist()




TH1F *GetAbbhhmHist2(const char* hName)
{
const double mp = 0.93827;
const int np = 33;
double eg[np] = 
  { 0.92, 0.97, 1.02, 1.07, 1.11, 1.17, 1.22, 1.27, 1.32, 1.37, 1.42,  
    1.47, 1.52, 1.57, 1.65, 1.74, 1.84, 1.95, 2.04, 2.14, 2.25, 2.34,  
    2.45, 2.54, 2.69, 2.89, 3.11, 3.36, 3.62, 3.86, 4.24, 4.73, 5.39  
  };

double sig[np] = 
  { 73.77, 73.20, 66.29, 68.97, 63.92, 57.45, 55.34, 53.36, 55.66, 
    56.43, 52.08, 52.91, 51.12, 51.25, 46.77, 45.87, 42.09, 38.32, 
    37.04, 39.47, 31.92, 30.77, 33.77, 30.00, 30.64, 27.50, 25.58, 
    24.49, 23.15, 24.69, 21.23, 22.00, 20.85
  };
double err[np] = 
  { 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 
    2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 
    2.4, 2.38, 2.36, 2.34, 2.32, 2.31, 
    2.3, 2.28, 2.25, 2.25, 2.2
  };


double wlow[np+1];
double elow;
  
for(int i=1; i<np; i++) {
  elow = eg[i] - (eg[i]-eg[i-1])/2.;
  wlow[i]  = sqrt( mp*(mp+2.*elow) ); }
  elow = eg[0] - (eg[1]-eg[0])/2.;
  wlow[0]  = sqrt( mp*(mp+2.*elow) );
  elow = eg[np-1] + (eg[np-1]-eg[np-2])/2.;
  wlow[np]  = sqrt( mp*(mp+2.*elow) );


TH1F *hist = new TH1F(hName, "saphir5542 tot", np, wlow);
for(int i=1; i<=np; i++) {
  hist->SetBinContent(i, sig[i-1]);
  hist->SetBinError(i, err[i-1]);
}

return hist;
};







TH1F *GetSaphir1Hist_3pi(const char* hName)
{
const double mp = 0.93827;
const int np = 20;
double eg[np] = 
  {  0.75, 0.85, 0.93, 1.05, 1.15, 1.24, 1.35, 1.45, 1.54, 
     1.65, 1.74, 1.85, 1.94, 2.24, 2.75, 3.23, 3.74, 4.23, 
     4.74, 5.29 
  };
double sig[np] = 
  { 3.60, 4.62, 6.06, 9.85, 15.0, 21.1, 23.0, 26.8, 27.6, 
    32.5, 39.1, 37.7, 42.4, 42.4, 44.9, 41.2, 35.8, 32.3, 
    22.4, 17.8 
  };
double err[np] = 
  { 1.1, 1.3, 1.35, 1.35, 1.4, 1.42, 1.5, 
    1.6, 1.6, 1.8,   1.9, 2.0, 2.2, 2.4, 
    2.44, 2.5, 2.5, 2.5, 2.5, 2.5, 
  };


double wlow[np+1];
double elow;
for(int i=1; i<np; i++) {
  elow = eg[i] - (eg[i]-eg[i-1])/2.;
  wlow[i]  = sqrt( mp*(mp+2.*elow) ); }
  elow = eg[0] - (eg[1]-eg[0])/2.;
  wlow[0]  = sqrt( mp*(mp+2.*elow) );
  elow = eg[np-1] + (eg[np-1]-eg[np-2])/2.;
  wlow[np] = sqrt( mp*(mp+2.*elow) );


TH1F *hist = new TH1F(hName, "saphir44444990 tot", np, wlow);
for(int i=1; i<=np; i++) {
  hist->SetBinContent(i, sig[i]);
  hist->SetBinError(i, err[i]);
}

return hist;
}; // end GetSaphir2Hist()






