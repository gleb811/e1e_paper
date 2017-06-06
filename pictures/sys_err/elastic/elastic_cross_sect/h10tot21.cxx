#include "TROOT.h"
#include "TFile.h"
#include "TLine.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include "TMacro.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TText.h"
#include "TStyle.h"
#include "TObject.h"
#include <TRint.h>
#include <stdio.h>
#include <dlfcn.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <iostream>
#include <stdio.h>
#include "MyMainFrame.h"


 using namespace std;      
     
     
     







int main(int argc, char** argv) {

// create an empty vector of strings
    vector<string> args;
    // copy program arguments into vector
    UChar_t flag;
    Int_t i;
    Int_t setnum_inp = 0;
    Int_t npoints_inp = 0;
    bool langau_inp = false;
    bool outtree_inp = false;
    string inpfile_inp = "clas_036516.A15.nt10.root";
    string outfile_inp = "out.root";
    Short_t number_of_files;
    Short_t number_of_files_sim;
    string qqq;
    string* file = NULL;
    string* file_sim = NULL;
    Float_t E0;
    
     getline (cin,qqq);
    qqq = qqq.substr(0, qqq.find(" ",0));    
    flag = atoi(qqq.c_str());
    if (flag == 1) {
    cout << "flag = 1 - data" << "\n";
    };
    if (flag == 2) {
    cout << "flag = 2 - simulation" << "\n";
    };
    
    
    
    getline (cin,qqq);
    qqq = qqq.substr(0, qqq.find(" ",0));
    E0 = atof(qqq.c_str());
    
    cout << "beam energy is " << E0 << " GeV" << "\n";
    
    getline (cin,qqq);
    qqq = qqq.substr(0, qqq.find(" ",0));    
    number_of_files = atoi(qqq.c_str());
    cout << "number of files for analysis = " << number_of_files << "\n";
    
    file = new string[number_of_files];
    
    
    for (i=1;i<=number_of_files; i++) {
    getline (cin,file[i-1]);
    file[i-1] = file[i-1].substr(0, file[i-1].find(" ",0));
    cout << "file " << i << " is " << file[i-1].c_str() << "\n"; 
    };
     
     
     
    getline (cin,qqq);
    qqq = qqq.substr(0, qqq.find(" ",0));    
    number_of_files_sim = atoi(qqq.c_str());
    cout << "number of files for simulation = " << number_of_files_sim << "\n";
    
    file_sim = new string[number_of_files_sim];
    
    
    for (i=1;i<=number_of_files_sim; i++) {
    getline (cin,file_sim[i-1]);
    file_sim[i-1] = file_sim[i-1].substr(0, file_sim[i-1].find(" ",0));
    cout << "file " << i << " for simulation is " << file_sim[i-1].c_str() << "\n"; 
    };     
     
     
     
     
     
     
     
     
     
     
     getline (cin,qqq);
    qqq = qqq.substr(0, qqq.find(" ",0));
         
    outfile_inp = qqq;
    
    cout << "output file is " << outfile_inp << "\n";
  
    TApplication *theApp = new TApplication("App", &argc, argv);
    MyMainFrame jopa;
    jopa.MainFrame(flag,E0,number_of_files,file,number_of_files_sim,file_sim,outfile_inp);
    
 //   TCanvas *cresolv = new TCanvas("cresolv", "Resolution", 500, 800);
//    TGMainFrame *fMain = new TGMainFrame(gClient->GetRoot(),200,200);
//    fMain->SetWindowName("Simple Example");
 //   TRootEmbeddedCanvas *fEcanvas = new TRootEmbeddedCanvas("Ecanvas",fMain,200,200); 
//    fMain->AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX| kLHintsExpandY,10,10,10,1));
   
 
//   cout << "jopa" << "\n";
   
//   theApp.Run();
  


  delete [] file;
  file = NULL;
  
  return 0;
}
