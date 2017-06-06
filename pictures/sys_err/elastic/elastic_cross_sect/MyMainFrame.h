#include <TGFrame.h>
#include <RQ_OBJECT.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <TGLabel.h>
#include <TGProgressBar.h>
#include <TGComboBox.h>
#include <TH1.h>
#include <TTree.h>
#include <TMacro.h>
#include <TMinuit.h>
#include "macro.h"
#include <RooLandau.h> 
#include <RooNumConvPdf.h>
#include <RooDataHist.h>
#include <RooBinning.h>
#include <RooAbsPdf.h>
#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooCmdArg.h>
#include <RooGaussian.h>
#include <RooMsgService.h>
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooProdPdf.h"
#include "RooCmdArg.h"
#include "RooMsgService.h"
#include "mom_corr.h"


        using namespace RooFit;
        using namespace std;

class MyMainFrame {
        RQ_OBJECT("MyMainFrame")
    private:
        mom_corr corrfunc;
        macro macros;
        TGComboBox *fCombo;
        TGComboBox *fNumPoints;	
	TGCheckButton *estat;
	TGCheckButton *inc_langau;	
        TGHProgressBar *fHProg1;
	TGLabel *label_ref_pmt;
	TGNumberEntry *ref_pmt;
        TGTextButton *go;
	TGTextButton *input_file;
	TGTextButton *output_file;
        TGMainFrame *fMain;
	TGTextButton *exit;
	TGHorizontalFrame *hframe;
        TRootEmbeddedCanvas *fEcanvas;
//	TGNumberEntry *n_points;
        TGLabel *label_n_points;
	TString path;	
	Double_t bar_length[62];
	const char* npoints_char[27];
	short nfiles;
	std::string _current_file;
	std::string _out_file;
	std::string _out_dir;
	TGLabel *labelq;
	TGLabel *labelqout;
	TGLabel *labellength;
	Bool_t file_status;
	TTree *t21;
        public:
	Int_t setnum_inp;
	Int_t npoints_inp;
	bool langau_inp;
	bool outtree_inp;
	string inpfile_inp;
	string outfile_inp;
	Float_t E0;
	string* file;
	Short_t n_files;
	string* file_sim;
	Short_t n_files_sim;	
	UChar_t data_sim;
        
        void MainFrame(UChar_t flag,Float_t E_beam,Short_t nfiles,string inp_files[],Short_t nfiles_sim,string inp_files_sim[], string outfile_in);
        void DoDraw();
	void t20tot21();
};
