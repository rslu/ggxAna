#include <vector>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
using namespace std;
#include <iostream>

#include "untuplizer.h"
#include "PhotonSelections.h"

Int_t MINITREE=1;
Double_t XS=0.;

Double_t deltaPhi(Double_t phi1, Double_t phi2) {

  Double_t dPhi = phi1 - phi2;
  if (dPhi > TMath::Pi()) dPhi -= 2.*TMath::Pi();
  if (dPhi < -TMath::Pi()) dPhi += 2.*TMath::Pi();

  return dPhi;
}

Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2) {

  Double_t dEta, dPhi ;
  dEta = eta1 - eta2;
  dPhi = deltaPhi(phi1, phi2);

  return sqrt(dEta*dEta+dPhi*dPhi);
}

void xAna(Int_t dataset = 0) {

  Char_t fname[200];
  XS=1.;

  if (dataset == 0) {
    sprintf(fname, "./job_phys14_gjet_pt20to40_20bx25.root");//2.49480900000000000e+06
    //XS=1.; //154500 pb, filter eff 0.001776
    XS=0.3043;
  }
  if (dataset == 1){
    sprintf(fname, "./job_phys14_gjet_pt40_20bx25.root"); //2.501057e+06 event processed
    XS = 1.;//154500.*0.001776/(17180.*0.038); //17180pb filter eff 0.038 //weight down to match pt20to40
  }

  vector <string> pathes;
  pathes.push_back(fname);
  TreeReader data(pathes);
  printf(" reading file %s \n", fname);


  TFile *fout_;
  TTree *outtree_;

  TH1F *hdpt = new TH1F("hdpt","dpt", 100, 0., 1.);
  TH1F *hdR = new TH1F("hdR","dR", 100, 0., 2.);
  TH1F *hmcCalIso = new TH1F("hmcCalIso","mcCalIso", 100, 0., 20.);
  
  TH1F *h_npho = new TH1F("h_npho","n pho", 20, 0., 20);
  TH1F *h_truepho = new TH1F("h_truepho","true pho", 10, 0., 10.);
  TH1F *h_phoEt = new TH1F("h_phoEt","pho Et", 200, 0., 1000);
  TH1F *h_jetPt = new TH1F("h_jetPt","jet Pt", 200, 0., 1000);
  TH1F *h_npj = new TH1F("h_npj","n pho jet comb", 20, 0., 20.);
  TH1F *h_pjmass = new TH1F("h_pjmass","inv mass of pho jet", 1000, 0., 10000);


  Char_t oname[200];
  sprintf(oname, "skim.root");
//   if (dataset == 0) sprintf(oname, "output_job_spring14_gjet_pt20to40_s14.root"); 
//   if (dataset == 1) sprintf(oname, "output_job_spring14_gjet_pt40_s14.root");	 
//   if (dataset == 2) sprintf(oname, "output_job_spring14_gjet_pt20to40_20bx25.root");
//   if (dataset == 3) sprintf(oname, "output_job_spring14_gjet_pt40_20bx25.root");  
//    if (dataset == 4) sprintf(oname, "output_job_spring14_RSGravToGG_m1000_s14.root");
//    if (dataset == 5) sprintf(oname, "output_job_spring14_RSGravToGG_m5000_s14.root");
//    if (dataset == 6) sprintf(oname, "output_job_spring14_QCD_pt1000to1400.root");
//    if (dataset == 7) sprintf(oname, "output_job_spring14_QCD_pt1400to1800.root");


  fout_ = new TFile(oname,"recreate");
  outtree_ = new TTree("t", "mini tree");

  Float_t effArea[3][7] = { //[Ch,Nu,Pho][iPhof_eta]
    { 0.012 , 0.010 , 0.014 , 0.012 , 0.016 , 0.020 , 0.012 } ,
    { 0.030 , 0.057 , 0.039 , 0.015 , 0.024 , 0.039 , 0.072 } ,
    { 0.148 , 0.130 , 0.112 , 0.216 , 0.262 , 0.260 , 0.266 } 
  } ;
  
  Float_t pthat_, mcPt_, mcEta_, mcPhi_, recoPt, recoEta, recoPhi, recoSCEta, r9, puwei;
  Int_t   isMatched, isMatchedEle, idLoose, idMedium, idTight, nVtx, eleVeto, nPU;
  Float_t HoverE, sieie, chIso, phoIso, nhIso, chIsoRaw, phoIsoRaw, nhIsoRaw, chWorstIso, rho;
  Float_t sieip, sipip, e1x3, e3x3, e2x2, e2x5, e5x5, rawE, scEtaWidth, scPhiWidth, esRR, esEn, mva;
  Float_t sieie_2012, sipip_2012, sieip_2012, e1x3_2012, e3x3_2012, e2x2_2012, e2x5_2012, e5x5_2012;

  Float_t mcCalIso04, mcTrkIso04;
  Float_t xsweight=XS;

  outtree_->Branch("pthat",         &pthat_,        "pthat/F");
  outtree_->Branch("mcPt",         &mcPt_,        "mcPt/F");
  outtree_->Branch("mcEta",        &mcEta_,       "mcEta/F");
  outtree_->Branch("mcPhi",        &mcPhi_,       "mcPhi/F");
  outtree_->Branch("recoPt",       &recoPt,       "recoPt/F");
  outtree_->Branch("recoEta",      &recoEta,      "recoEta/F");
  outtree_->Branch("recoPhi",      &recoPhi,      "recoPhi/F");
  outtree_->Branch("recoSCEta",    &recoSCEta,    "recoSCEta/F");
  outtree_->Branch("r9",           &r9,           "r9/F");
  outtree_->Branch("isMatched",    &isMatched,    "isMatched/I");
  outtree_->Branch("isMatchedEle", &isMatchedEle, "isMatchedEle/I");
  outtree_->Branch("idLoose",      &idLoose,      "idLoose/I");
  outtree_->Branch("idMedium",     &idMedium,     "idMedium/I");
  outtree_->Branch("idTight",      &idTight,      "idTight/I");
  outtree_->Branch("nVtx",         &nVtx,         "nVtx/I");
  outtree_->Branch("nPU",          &nPU,          "nPU/I");
  //outtree_->Branch("puwei",        &puwei,        "puwei/F");
  outtree_->Branch("eleVeto",      &eleVeto,      "eleVeto/I");
  outtree_->Branch("HoverE",       &HoverE,       "HoverE/F");
  outtree_->Branch("sieie",        &sieie,        "sieie/F");
  outtree_->Branch("sieip",        &sieip,        "sieip/F");
  outtree_->Branch("sipip",        &sipip,        "sipip/F");
  outtree_->Branch("chIso",        &chIso,        "chIso/F");
  outtree_->Branch("phoIso",       &phoIso,       "phoIso/F");
  outtree_->Branch("nhIso",        &nhIso,        "nhIso/F");
  outtree_->Branch("chIsoRaw",     &chIsoRaw,     "chIsoRaw/F");
  outtree_->Branch("chWorstRaw",   &chWorstIso,   "chWorstIso/F");
  outtree_->Branch("phoIsoRaw",    &phoIsoRaw,    "phoIsoRaw/F");
  outtree_->Branch("nhIsoRaw",     &nhIsoRaw,     "nhIsoRaw/F");
  outtree_->Branch("rho",          &rho,          "rho/F"); 
  outtree_->Branch("e1x3",         &e1x3,         "e1x3/F");
  //outtree_->Branch("e3x3",         &e3x3,         "e3x3/F");
  outtree_->Branch("e2x2",         &e2x2,         "e2x2/F");
  outtree_->Branch("e2x5",         &e2x5,         "e2x5/F");
  outtree_->Branch("e5x5",         &e5x5,         "e5x5/F");
  outtree_->Branch("rawE",         &rawE,         "rawE/F");
  outtree_->Branch("scEtaWidth",   &scEtaWidth,   "scEtaWidth/F");
  outtree_->Branch("scPhiWidth",   &scPhiWidth,   "scPhiWidth/F");
  outtree_->Branch("esRR",         &esRR,         "esRR/F");   
  outtree_->Branch("esEn",         &esEn,         "esEn/F");   
  outtree_->Branch("mva",          &mva,          "mva/F");  

  outtree_->Branch("sieie_2012",        &sieie_2012,        "sieie_2012/F");
  outtree_->Branch("sieip_2012",        &sieip_2012,        "sieip_2012/F");
  outtree_->Branch("sipip_2012",        &sipip_2012,        "sipip_2012/F");
  outtree_->Branch("e1x3_2012",        &e1x3_2012,        "e1x3_2012/F");
  //outtree_->Branch("e3x3_2012",        &e3x3_2012,        "e3x3_2012/F");
  outtree_->Branch("e2x2_2012",        &e2x2_2012,        "e2x2_2012/F");
  outtree_->Branch("e2x5_2012",        &e2x5_2012,        "e2x5_2012/F");
  outtree_->Branch("e5x5_2012",        &e5x5_2012,        "e5x5_2012/F");

  outtree_->Branch("mcCalIso04",  &mcCalIso04, "mcCalIso04/F");
  outtree_->Branch("mcTrkIso04",  &mcTrkIso04, "mcTrkIso04/F");
  outtree_->Branch("xsweight",  &xsweight, "xsweight/F");

  // event loop
  //  for (Long64_t ev = 0; ev <200; ev++) {
  for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++) {
    //for (Long64_t ev = 0; ev < data.GetEntriesFast()/10.; ev++) {

    if (ev % 100000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());
    
    data.GetEntry(ev);

    Float_t pthat     = data.GetFloat("pthat");
    Int_t    nMC       = data.GetInt("nMC");
    Int_t*   mcPID     = data.GetPtrInt("mcPID");
    Int_t*   mcMomPID  = data.GetPtrInt("mcMomPID");
    Int_t*   mcGMomPID = data.GetPtrInt("mcGMomPID");
    Float_t* mcPt      = data.GetPtrFloat("mcPt");
    Float_t* mcEta     = data.GetPtrFloat("mcEta");
    Float_t* mcPhi     = data.GetPtrFloat("mcPhi");
    Float_t* mcE       = data.GetPtrFloat("mcE");
    
    nVtx = data.GetInt("nVtx");
    
    Int_t    run     = data.GetInt("run");
    Long64_t event   = data.GetLong64("event");  
    Int_t    nPUInfo = data.GetInt("nPUInfo");
    Int_t*   puBX    = data.GetPtrInt("puBX");
    Float_t* puTrue  = data.GetPtrFloat("puTrue");
    //puwei = (Float_t) puCalc.GetWeight(run, puTrue[1]);
    
    for (Int_t i=0; i<nPUInfo; ++i) {
      if (puBX[i] == 0) nPU = puTrue[i];
    }

    Int_t    nPho     = data.GetInt("nPho");
    Float_t* phoE   = data.GetPtrFloat("phoE");
    Float_t* phoEta   = data.GetPtrFloat("phoEta");
    Float_t* phoPhi   = data.GetPtrFloat("phoPhi");
    Float_t* phoEt    = data.GetPtrFloat("phoEt");
    Float_t* phoR9    = data.GetPtrFloat("phoR9");
    Float_t* phoSCEta = data.GetPtrFloat("phoSCEta");
    Float_t* phoSCPhi = data.GetPtrFloat("phoSCPhi");
    
    Int_t*   phoEleVeto          = data.GetPtrInt("phoEleVeto");
    Float_t* phoSigmaIEtaIEta    = data.GetPtrFloat("phoSigmaIEtaIEta");
    Float_t* phoHoverE           = data.GetPtrFloat("phoHoverE");
    Float_t* phoPFChIso          = data.GetPtrFloat("phoPFChIso");
    Float_t* phoPFNeuIso         = data.GetPtrFloat("phoPFNeuIso");
    Float_t* phoPFPhoIso         = data.GetPtrFloat("phoPFPhoIso");

    rho                 = data.GetFloat("rho");

    Float_t* phoSigmaIEtaIPhi  = data.GetPtrFloat("phoSigmaIEtaIPhi");
    Float_t* phoSigmaIPhiIPhi  = data.GetPtrFloat("phoSigmaIPhiIPhi");
    Float_t* phoE1x3           = data.GetPtrFloat("phoE1x3");
    //Float_t* phoE3x3           = data.GetPtrFloat("phoE3x3");
    Float_t* phoE2x2           = data.GetPtrFloat("phoE2x2");
    Float_t* phoE5x5           = data.GetPtrFloat("phoE5x5");
    Float_t* phoE2x5Max        = data.GetPtrFloat("phoE2x5Max");
    Float_t* phoSCRawE         = data.GetPtrFloat("phoSCRawE");
    Float_t* phoSCEtaWidth     = data.GetPtrFloat("phoSCEtaWidth");
    Float_t* phoSCPhiWidth     = data.GetPtrFloat("phoSCPhiWidth");
    Float_t* phoESEn           = data.GetPtrFloat("phoESEn");
    Float_t* phoESEffSigmaRR   = data.GetPtrFloat("phoESEffSigmaRR");
    Float_t* phoPFChWorstIso   = data.GetPtrFloat("phoPFChWorstIso");

    Float_t* phoSigmaIEtaIEta_2012  = data.GetPtrFloat("phoSigmaIEtaIEta_2012");
    Float_t* phoSigmaIEtaIPhi_2012  = data.GetPtrFloat("phoSigmaIEtaIPhi_2012");
    Float_t* phoSigmaIPhiIPhi_2012  = data.GetPtrFloat("phoSigmaIPhiIPhi_2012");

    Float_t* phoE1x3_2012           = data.GetPtrFloat("phoE1x3_2012");
    //Float_t* phoE3x3_2012           = data.GetPtrFloat("phoE3x3_2012");
    Float_t* phoE2x2_2012            = data.GetPtrFloat("phoE2x2_2012");
    Float_t* phoE5x5_2012            = data.GetPtrFloat("phoE5x5_2012");
    Float_t* phoE2x5Max_2012        = data.GetPtrFloat("phoE2x5Max_2012");

    Float_t* mcCalIsoDR04 = data.GetPtrFloat("mcCalIsoDR04");
    Float_t* mcTrkIsoDR04 = data.GetPtrFloat("mcTrkIsoDR04");

    vector<int> acc_phol;
    select_photon(0, data, acc_phol);
    vector<int> acc_phom;
    select_photon(1, data, acc_phom);
    vector<int> acc_phot;
    select_photon(2, data, acc_phot);

    vector<int> match;
    vector<int> match_ele;
    vector<float> mcpt;
    vector<float> mceta;
    vector<float> mcphi;

    h_npho->Fill(nPho);
    int nmatch=0;
    int nfake=0;
//     printf("-----------------------------------------------------------------------\n");
//     printf("event %d, npho %d, nMC %d\n", event, nPho, nMC);
    for (Int_t i=0; i<nPho; ++i) {
      isMatched    = -99;
      isMatchedEle = -1;
      //  printf("pho Et %.2f, eta %.2f, phi %.2f  \n", phoEt[i], phoEta[i], phoPhi[i]);
      for (Int_t k=0; k<nMC; ++k) {	
	//if (mcPID[k] == 22 && (mcMomPID[k] <= 22 || mcMomPID[k] == 5100039)) {
	if (mcPID[k] == 22 && TMath::Abs(mcMomPID[k]) <= 22 ) {
	  //  printf("mc  pt %.2f, eta %.2f, phi %.2f \n", mcPt[k], mcEta[k], mcPhi[k]);
	  float dr = deltaR(phoEta[i], phoPhi[i], mcEta[k], mcPhi[k]);
	  float dpt = fabs((phoEt[i] - mcPt[k])/mcPt[k]);
	  hdR->Fill(dr);
	  hdpt->Fill(dpt);
	  hmcCalIso->Fill(mcCalIsoDR04[k]);
	  //if (dr < 0.2 && dpt < 0.2 && mcCalIsoDR04[k] < 5.){
	  //  printf("  MCparticle %d, dr %.2f, dpt %.2f \n", k, dr, dpt);
	  if (dr < 0.2 && dpt < 0.2){
	    isMatched = 1;
 	    mcPt_  = mcPt[k];
 	    mcEta_ = mcEta[k];
 	    mcPhi_ = mcPhi[k];
	  
	  }
	}	
	
	if (fabs(mcPID[k]) == 11 && mcMomPID[k] == 23) {
	  if (deltaR(phoEta[i], phoPhi[i], mcEta[k], mcPhi[k]) < 0.2) {
	    if (fabs((phoEt[i] -mcPt[k])/mcPt[k]) < 0.2) {
	      isMatchedEle = 1;
	    }
	  }
	}	
      }
      mcpt.push_back(mcPt_);
      mceta.push_back(mcEta_);
      mcphi.push_back(mcPhi_);
      match.push_back(isMatched);
      match_ele.push_back(isMatchedEle);
      if(isMatched==1) nmatch++;

//       if(isMatched==1) printf("pho Et %.2f is  true photon\n", phoEt[i]);      
//       else printf("pho Et %.2f is                               fake photon\n", phoEt[i]);      

      //if(phoEt[i] < 150.) continue;
      //       if(isMatched==1 && TMath::Abs(phoSCEta[i])<=1.4442) h_phoEt->Fill(phoEt[i]);
      //       else if(isMatched!=1) h_jetPt->Fill(phoEt[i]);
      
    }
    //h_truepho->Fill((float)nmatch/(float)nPho);
    h_truepho->Fill(nmatch+0.001);
    //ask for only one mc true photon
    //if(nmatch < 1) continue;
//     printf(" nmatch %d, nfake %d \n", nmatch, nfake);

    int npj=0;

    for (Int_t i=0; i<nPho; ++i) {

      //isMatched    = -99;
      //isMatchedEle = -1;
      idLoose      = -1;
      idMedium     = -1;
      idTight      = -1;


      isMatched = match[i];
      isMatchedEle = match_ele[i];
      
      //if (isMatched == 1)
      //cout<<isMatched<<" "<<isMatchedEle<<" "<<phoGenMomPID[i]<<endl;
      pthat_    = pthat;
      mcPt_     = mcpt[i];
      mcEta_    = mceta[i];
      mcPhi_    = mcphi[i];

      recoPt    = phoEt[i];
      recoEta   = phoEta[i];
      recoPhi   = phoPhi[i];
      recoSCEta = phoSCEta[i];
      r9        = phoR9[i];
      eleVeto   = phoEleVeto[i];
      HoverE    = phoHoverE[i];
      sieie     = phoSigmaIEtaIEta[i];
      
      Int_t i_effArea = 0 ; // effective area for pile up correction for DR04 combine rel. Iso
      if      ( fabs(phoSCEta[i]) < 1.0                                        ) i_effArea = 0 ;
      else if ( fabs(phoSCEta[i]) >= 1.0   && fabs(phoSCEta[i]) < 1.479  ) i_effArea = 1 ;
      else if ( fabs(phoSCEta[i]) >= 1.479 && fabs(phoSCEta[i]) < 2.0    ) i_effArea = 2 ;
      else if ( fabs(phoSCEta[i]) >= 2.0   && fabs(phoSCEta[i]) < 2.2    ) i_effArea = 3 ;
      else if ( fabs(phoSCEta[i]) >= 2.2   && fabs(phoSCEta[i]) < 2.3    ) i_effArea = 4 ;
      else if ( fabs(phoSCEta[i]) >= 2.3   && fabs(phoSCEta[i]) < 2.4    ) i_effArea = 5 ;
      else if ( fabs(phoSCEta[i]) >= 2.4                                       ) i_effArea = 6 ;
      
      chIso  = TMath::Max( float(0.) , phoPFChIso[i] - ( effArea[0][i_effArea] * rho ) ) ;
      phoIso = TMath::Max( float(0.) , phoPFPhoIso[i] - ( effArea[2][i_effArea] * rho ) ) ;
      nhIso  = TMath::Max( float(0.) , phoPFNeuIso[i] - ( effArea[1][i_effArea] * rho ) ) ;
      chIsoRaw   = phoPFChIso[i];
      phoIsoRaw  = phoPFPhoIso[i];
      nhIsoRaw   = phoPFNeuIso[i];

      sieip      = phoSigmaIEtaIPhi[i];
      sipip      = phoSigmaIPhiIPhi[i];
      e1x3       = phoE1x3[i];
      //e3x3       = phoE3x3[i];
      e2x2       = phoE2x2[i];
      e2x5       = phoE2x5Max[i];
      e5x5       = phoE5x5[i];
      rawE       = phoSCRawE[i];
      scEtaWidth = phoSCEtaWidth[i];
      scPhiWidth = phoSCPhiWidth[i];
      esRR       = phoESEffSigmaRR[i];
      esEn       = phoESEn[i];
      chWorstIso = phoPFChWorstIso[i];

      sieie_2012     = phoSigmaIEtaIEta_2012[i];
      sieip_2012     = phoSigmaIEtaIPhi_2012[i];
      sipip_2012     = phoSigmaIPhiIPhi_2012[i];

      e1x3_2012       = phoE1x3_2012[i];
      //e3x3_2012       = phoE3x3_2012[i];
      e2x2_2012       = phoE2x2_2012[i];
      e2x5_2012       = phoE2x5Max_2012[i];
      e5x5_2012       = phoE5x5_2012[i];

      mcCalIso04 = mcCalIsoDR04[i];
      mcTrkIso04 = mcTrkIsoDR04[i];

      mva = select_photon_mva(data, i);

      for (Int_t j=0; j<acc_phol.size(); ++j)
	if (i == acc_phol[j]) idLoose = 1;
      
      for (Int_t j=0; j<acc_phom.size(); ++j)
	if (i == acc_phom[j]) idMedium = 1;
      
      for (Int_t j=0; j<acc_phot.size(); ++j)
	if (i == acc_phot[j]) idTight = 1;
      
      if(MINITREE==1) {
	if ((dataset == 4 || dataset==5)) {
  	  if(isMatched==1)	  outtree_->Fill();
	}else if ((dataset == 6 || dataset==7)) {
  	  if(isMatched!=1)	  outtree_->Fill();
	}else outtree_->Fill();
	
      }

      if(isMatched!=1) continue;

      //loop jets for photon+jet invarient mass
      for(int j=0; j<nPho; j++){	
	if(i==j) continue;
	if(match[j]==1) continue;
	if(phoEt[j]<150.) continue;
	
	TLorentzVector phoP4;// = new TLorentzVector();
	TLorentzVector jetP4;// = new TLorentzVector();
	TLorentzVector pjP4;//= new TLorentzVector();
	
	phoP4.SetPtEtaPhiM(phoEt[i], phoEta[i],phoPhi[i], 0);
	jetP4.SetPtEtaPhiM(phoEt[j], phoEta[j],phoPhi[j], phoE[j]);
	pjP4 = phoP4;
	pjP4 += jetP4;
	h_phoEt->Fill(phoEt[i]);
	h_jetPt->Fill(jetP4.Pt());
	h_pjmass->Fill(pjP4.M());
	npj++;
      }
      

    }
    
    h_npj->Fill(npj);

  }
  
  fout_->cd();
  outtree_->Write();
  hdR->Write();
  hdpt->Write();
  hmcCalIso->Write();

  h_truepho->Write();
  h_npho->Write();
  h_pjmass->Write();

  h_npj->Write();
  h_phoEt->Write();
  h_jetPt->Write();

  fout_->Close();

  
  fprintf(stderr, "Processed all events\n");

}
