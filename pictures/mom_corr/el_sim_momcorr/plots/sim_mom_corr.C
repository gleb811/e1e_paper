void sim_mom_corr() {

#define _USE_MATH_DEFINES

gStyle->SetOptStat(0);
gStyle->SetOptFit(1);

TFile *MyFile = new TFile("gpp194testnew31.root","READ");
TCanvas *c = new TCanvas("c","c",800,400);
TTree *t10 = (TTree*)MyFile->Get("h10");
 TBranch *br_p = t10->GetBranch("p");
 TBranch *br_pxpart = t10->GetBranch("pxpart"); 
 TBranch *br_pypart = t10->GetBranch("pypart"); 
 TBranch *br_pzpart = t10->GetBranch("pzpart"); 
 TBranch *br_id = t10->GetBranch("id");  
 TBranch *br_pidpart = t10->GetBranch("pidpart"); 
 TBranch *br_gpart = t10->GetBranch("gpart");  
 TBranch *br_stat = t10->GetBranch("stat");
 TBranch *br_cx = t10->GetBranch("cx");
 TBranch *br_cy = t10->GetBranch("cy"); 
 TBranch *br_cz = t10->GetBranch("cz"); 
 
 Float_t m_proton = 0.93821;
 Float_t P,pxpart,pypart,pzpart,ph_EL,th_EL,W_rec,Q2_rec,W_gen,Q2_gen,P_gen,th_gen;
 
 Float_t p0,p1,p2;
 
 Int_t id,pidpart,gpart,stat;
 
 TH1F *h = new TH1F("h","h",100,-0.05,0.05);
 TH1F *h1 = new TH1F("W_rec","W_rec",100,-0.05,0.05); 
 TH1F *h2 = new TH1F("dpoverp","dpoverp",100,-0.05,0.05);
 
   for (Int_t i=0; i<br_pxpart->GetEntries(); i++) { 
   br_p->GetEntry(i);
   br_pxpart->GetEntry(i);
   br_pypart->GetEntry(i);
   br_pzpart->GetEntry(i);
   br_id->GetEntry(i);
   br_pidpart->GetEntry(i);
   br_gpart->GetEntry(i);
   br_stat->GetEntry(i);
   br_cx->GetEntry(i);
   br_cy->GetEntry(i);   
   br_cz->GetEntry(i);      
    
    P = br_p->GetLeaf("p")->GetValue(0);
    pxpart = br_pxpart->GetLeaf("pxpart")->GetValue(0);
    pypart = br_pypart->GetLeaf("pypart")->GetValue(0); 
    pzpart = br_pzpart->GetLeaf("pzpart")->GetValue(0); 
    id = br_id->GetLeaf("id")->GetValue(0); 
    pidpart = br_pidpart->GetLeaf("pidpart")->GetValue(); 
    gpart = br_gpart->GetLeaf("gpart")->GetValue(0);
    stat = br_stat->GetLeaf("stat")->GetValue(0); 
    
    if (br_cx->GetLeaf("cx")->GetValue(0) != 0.) {
     ph_EL = (180./M_PI)*atan(br_cy->GetLeaf("cy")->GetValue(0)/br_cx->GetLeaf("cx")->GetValue(0));
     } else {
     if(br_cy->GetLeaf("cy")->GetValue(0) > 0.) ph_EL = 90.;
     if(br_cy->GetLeaf("cy")->GetValue(0) < 0.) ph_EL = 270.;
      };
     if ((br_cx->GetLeaf("cx")->GetValue(0) < 0.) && (br_cy->GetLeaf("cy")->GetValue(0) > 0)) ph_EL = ph_EL+180.;
     if ((br_cx->GetLeaf("cx")->GetValue(0) < 0.) && (br_cy->GetLeaf("cy")->GetValue(0) < 0)) ph_EL = ph_EL+180.;
     if ((br_cx->GetLeaf("cx")->GetValue(0) > 0.) && (br_cy->GetLeaf("cy")->GetValue(0) < 0)) ph_EL = ph_EL+360.;				   

    th_EL = (180./M_PI)*acos(br_cz->GetLeaf("cz")->GetValue(0));			    
   

    Q2_rec = 4.*2.039*P*(sin(th_EL*M_PI/2./180.))*(sin(th_EL*M_PI/2./180.));
    W_rec = m_proton*m_proton+2.*m_proton*(2.039-P)-Q2_rec;
    
    P_gen = sqrt(pxpart*pxpart+pypart*pypart+pzpart*pzpart);
    th_gen = (180./M_PI)*acos(pzpart/P_gen);
    Q2_gen = 4.*2.039*P_gen*(sin(th_gen*M_PI/2./180.))*(sin(th_gen*M_PI/2./180.));
    W_gen = m_proton*m_proton+2.*m_proton*(2.039-P_gen)-Q2_gen;
   
   
    if (W_rec > 0.) { 
     W_rec = sqrt(W_rec); } else {
     W_rec = 1000.;
     };
     
   if (W_gen > 0.) { 
     W_gen = sqrt(W_gen); } else {
     W_gen = 0.;
     };  
        
    if ((id == 11) && (pidpart == 3) && (gpart > 0) && (stat > 0) ) {
    if (sqrt(pxpart*pxpart+pypart*pypart+pzpart*pzpart) < 1.7) {
    h->Fill((sqrt(pxpart*pxpart+pypart*pypart+pzpart*pzpart)-P),1.);
   
    p0 = -0.00937397+th_EL*0.000462166+th_EL*th_EL*(-5.20498e-06);
    p1 = 0.021494+th_EL*(-0.00121626)+th_EL*th_EL*(1.39892e-05);
    p2 = -0.0110274+th_EL*0.000586353+th_EL*th_EL*(-6.90819e-06);
    
  //  cout <<  << endl;
    P = P-p0-p1*P-p2*P*P; 
 h2->Fill((sqrt(pxpart*pxpart+pypart*pypart+pzpart*pzpart)-P),1.);   //h2->Fill((sqrt(pxpart*pxpart+pypart*pypart+pzpart*pzpart)-P)/(sqrt(pxpart*pxpart+pypart*pypart+pzpart*pzpart)),1.);
    h1->Fill(W_gen-W_rec);
    };
    };
   // cout << id << endl;

    
    
   // cout << W_rec << endl;
   
   };
     c->Divide(2,1);
     c->cd(1);
     c->cd(1)->SetBottomMargin(0.11); 
     h->GetYaxis()->SetTitle("counts");
     h->GetYaxis()->SetTitleSize(0.05);
     h->GetYaxis()->SetTitleOffset(0.95);
     h->GetXaxis()->SetTitleSize(0.042);
     h->GetXaxis()->SetTitleOffset(1.08); 
     h->GetXaxis()->SetTitle("(P_{e,generated}-P_{e,reconstructed}) (GeV)");
     h->SetTitle("");
    h->Fit("gaus","","",-0.003,0.007);
    h->Draw();  
    
    
    c->cd(2); 
    c->cd(2)->SetBottomMargin(0.11); 
    h2->SetTitle("");
    h2->Draw(); 
    h2->GetYaxis()->SetTitle("counts");
    h2->GetYaxis()->SetTitleSize(0.05);
     h2->GetYaxis()->SetTitleOffset(0.95);
     h2->GetXaxis()->SetTitleSize(0.042);
     h2->GetXaxis()->SetTitleOffset(1.08); 
     h2->GetXaxis()->SetTitle("(P_{e,generated}-P_{e,reconstructed}) (GeV)");
     h2->Fit("gaus","","",-0.004,0.004); 
     
     c->SaveAs("el_sim_mom_corr.pdf");
};
