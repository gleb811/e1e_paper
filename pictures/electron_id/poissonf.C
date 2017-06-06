double poissonf( double*x, double*par){

return par[0]*TMath::Poisson(x[0]/par[2],par[1]);
} ;
