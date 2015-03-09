#define Ana_gj_cxx
#include "Ana_gj.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <fstream>
#include <iostream>

#define MCTYPE 1
#define EBEECATE 1

#define combIso_bit 0
#define sieie_bit   1
#define HoverE_bit  2
#define Pixel_bit   3
#define ecalIso_bit  4
#define hcalIso_bit  5
#define trkIso_bit  6

#define JES 0 //+1, -1 sigma or 0 for Jet Energy Scale
#define JER 0 //+1, -1 Jet energy resolution
#define PHOER 0 //+1, -1, 0 Photon energy resolution

void Ana_gj::Loop() {
  
  TProfile* prof_phoEt_run = new TProfile("prof_phoEt_run","profile of photon Et vs run number", 1000, 160000, 170000);

  TH1F *h_nVtx = new TH1F("h_nVtx","nVtx distribution", 30, 0., 30);
  TH1F *h_nVtx_w = new TH1F("h_nVtx_w","weighted nVtx distribution", 30, 0., 30);
  TH1F *h_nPU = new TH1F("h_nPU","nPU distribution", 30, 0., 30);

  TH1F *h_njet = new TH1F("h_njet","number of jet",20, 0., 20);
  TH1F *h_jetPt = new TH1F("h_jetPt","pt of jet", 100, 0., 1000);
  
  TH1F *h_phoEt_all = new TH1F("h_phoEt_all","Et of photon", 100, 0., 1000);
  TH1F *h_phoEta = new TH1F("h_phoEta","Eta of photon",100, -3., 3.);
  TH1F *h_phoPhi = new TH1F("h_phoPhi","Phi of photon",100, -3.5, 3.5);
  TH1F *h_gj_dR = new TH1F("h_gj_dR","dR of gj",100, 0., 10.);
  TH1F *h_gj_dEta = new TH1F("h_gj_dEta","dEta of gj",100, 0., 5.);
  TH1F *h_gj_dPhi = new TH1F("h_gj_dPhi","dPhi of gj",100, 0., TMath::Pi());


  TH1F *h_npho[EBEECATE][MCTYPE];
  TH1F *h_phoEt[EBEECATE][MCTYPE];
  TH1F *h_combIso[EBEECATE][MCTYPE];
  TH1F *h_sieie[EBEECATE][MCTYPE];
  TH2F *h_phoEt_combIso[EBEECATE][MCTYPE];
  TH2F *h_phoEt_ecalIso[EBEECATE][MCTYPE];
  TH2F *h_phoEt_hcalIso[EBEECATE][MCTYPE];
  TH2F *h_phoEt_trkIso[EBEECATE][MCTYPE];
  TH2F *h_phoEt_sieie[EBEECATE][MCTYPE];
//   TH1F *h_combIso_jet[EBEECATE][MCTYPE];
  TH2F *h_mgj_ptg[EBEECATE][MCTYPE];
//   TH2F *h_mgj_ptj[EBEECATE][MCTYPE];
  TH2F *h_combIso_sieie[EBEECATE][MCTYPE];
  TH3F *h_combIso_sieie_phoEt[EBEECATE][MCTYPE];
  TH3F *h_combIso_sieie_mgj[EBEECATE][MCTYPE];
//   TH2F *h_trkIso_HoverE[EBEECATE][MCTYPE];
//   TH3F *h_trkIso_HoverE_phoEt[EBEECATE][MCTYPE];
//   TH3F *h_trkIso_HoverE_mgj[EBEECATE][MCTYPE];

  int category;
  float gjmass;
  float mcwei;

  TTree *tt = new TTree("tt","minitree");
  tt->Branch("category", &category, "categroy/I");
  tt->Branch("gjmass", &gjmass, "gjmass/F");
  tt->Branch("mcwei", &mcwei, "mcwei/F");


  char EBEE_cat[2][32] = {
    "EB", "EE",
  };
  
  char MCTYPE_cat[2][32] = {
    "phojet", "dijet",
  };

  char buffer[128];
  
  for(int cate=0; cate<EBEECATE; cate++){
    for(int type=0; type<MCTYPE; type++){
      sprintf(buffer,"h_%s_npho_%s",EBEE_cat[cate],MCTYPE_cat[type]);
      h_npho[cate][type] = new TH1F(buffer,buffer,  10, 0., 10.);
      h_npho[cate][type]->Sumw2();
      sprintf(buffer,"h_%s_phoEt_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_phoEt[cate][type] = new TH1F(buffer,buffer,  100, 0., 1000);
      h_phoEt[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_combIso_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_combIso[cate][type] = new TH1F(buffer,buffer,  50, -5., 20.);
      h_combIso[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_sieie_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
//       h_sieie[cate][type] = new TH1F(buffer,buffer,  80, 0., 0.04);
      h_sieie[cate][type] = new TH1F(buffer,buffer,  300, 0.008, 0.01);
      h_sieie[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_phoEt_combIso_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_phoEt_combIso[cate][type] = new TH2F(buffer,buffer, 100, 0.,1000., 50, -5., 20.);
      h_phoEt_combIso[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_phoEt_ecalIso_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_phoEt_ecalIso[cate][type] = new TH2F(buffer,buffer, 100, 0.,1000., 50, -5., 20.);
      h_phoEt_ecalIso[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_phoEt_hcalIso_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_phoEt_hcalIso[cate][type] = new TH2F(buffer,buffer, 100, 0.,1000., 50, -5., 20.);
      h_phoEt_hcalIso[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_phoEt_trkIso_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_phoEt_trkIso[cate][type] = new TH2F(buffer,buffer, 100, 0.,1000., 50, -5., 20.);
      h_phoEt_trkIso[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_phoEt_sieie_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_phoEt_sieie[cate][type] = new TH2F(buffer,buffer, 100, 0.,1000., 96, 0., 0.024);
      h_phoEt_sieie[cate][type]->Sumw2(); 
//       sprintf(buffer,"h_%s_combIso_jet_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
//       h_combIso_jet[cate][type] = new TH1F(buffer,buffer,  50, -5., 20.);
//       h_combIso_jet[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_mgj_ptg_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_mgj_ptg[cate][type] = new TH2F(buffer,buffer, 350, 0., 3500., 100, 0., 1000);
      h_mgj_ptg[cate][type]->Sumw2(); 
//       sprintf(buffer,"h_%s_mgj_ptj_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
//       h_mgj_ptj[cate][type] = new TH2F(buffer,buffer, 350, 0., 3500., 100, 0., 1000);
//       h_mgj_ptj[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_combIso_sieie_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_combIso_sieie[cate][type] = new TH2F(buffer,buffer, 50, -5., 20., 80, 0., 0.04);
      h_combIso_sieie[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_combIso_sieie_phoEt_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_combIso_sieie_phoEt[cate][type] = new TH3F(buffer,buffer, 50, -5., 20., 80, 0., 0.04, 100, 0., 1000.);
      h_combIso_sieie_phoEt[cate][type]->Sumw2(); 
      sprintf(buffer,"h_%s_combIso_sieie_mgj_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
      h_combIso_sieie_mgj[cate][type] = new TH3F(buffer,buffer, 50, -5., 20., 80, 0., 0.04, 350, 0., 3500.);
      h_combIso_sieie_mgj[cate][type]->Sumw2(); 
//       sprintf(buffer,"h_%s_trkIso_HoverE_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
//       h_trkIso_HoverE[cate][type] = new TH2F(buffer,buffer, 50, -5., 20., 100, 0., 0.1);  
//       h_trkIso_HoverE[cate][type]->Sumw2(); 
//       sprintf(buffer,"h_%s_trkIso_HoverE_phoEt_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
//       h_trkIso_HoverE_phoEt[cate][type] = new TH3F(buffer,buffer, 50, -5., 20., 100, 0., 0.1, 100, 0., 1000.);
//       h_trkIso_HoverE_phoEt[cate][type]->Sumw2(); 
//       sprintf(buffer,"h_%s_trkIso_HoverE_mgj_%s",EBEE_cat[cate],MCTYPE_cat[type]); 
//       h_trkIso_HoverE_mgj[cate][type] = new TH3F(buffer,buffer, 50, -5., 20., 100, 0., 0.1, 350, 0., 3500.);
//       h_trkIso_HoverE_mgj[cate][type]->Sumw2(); 

    }
  }

  // sumw2*h_nVtx->Write();
  h_nVtx_w->Sumw2();
  h_nPU->Sumw2();
  h_njet->Sumw2();
  h_jetPt->Sumw2();
  h_phoEt_all->Sumw2();
  h_phoEta->Sumw2();
  h_phoPhi->Sumw2();
  h_gj_dR->Sumw2();
  h_gj_dEta->Sumw2();
  h_gj_dPhi->Sumw2();

  TRandom3 *rnd = new TRandom3();
  rnd->SetSeed(0);

  //load JES table
  printf(" ----  Load JEC configuration ---- \n");
  const int NETAJEC = 32;
  const int NPTJEC = 39; //40 for Spring file
  double etajec[NETAJEC];
  double ptjec[NETAJEC][NPTJEC][3];
//   std::ifstream jecfile("Spring10DataV2_Uncertainty_AK5PF.txt");
  std::ifstream jecfile("GR_R_42_V19_AK5PF_Uncertainty.txt");
  for(int neta = 0; neta < NETAJEC; neta++) {
    double temp1, temp2;
    jecfile>>etajec[neta]>>temp1>>temp2;
    for(int npt = 0; npt < NPTJEC; npt++) {
      jecfile>>ptjec[neta][npt][0]>>ptjec[neta][npt][1]>>ptjec[neta][npt][2];
    }    
  }    
//   for(int neta = 0; neta < NETAJEC; neta++) {
//     for(int npt = 0; npt < NPTJEC; npt++) {
//       printf("JEC Eta %i, Pt %i, JEC %f, %f, %f, %f \n",
// 	     neta, npt, etajec[neta], ptjec[neta][npt][0], ptjec[neta][npt][1],
// 	     ptjec[neta][npt][2]);
//     }    
//   }    

  // load for PU weight

  TFile *fdata = new TFile("Pileup_2011_EPS_8_jul.root");
  TFile *fmc = new TFile("ana_gj_out_FNAME_MASK.root");

  TH1F *h_nPU_data = (TH1F*)fdata->Get("pileup");
  TH1F *h_nPU_mc = (TH1F*)fmc->Get("h_nPU");

  h_nPU_data->Scale(1./h_nPU_data->Integral());
  h_nPU_mc->Scale(1./h_nPU_mc->Integral());

  float weight[30];
  for(int ii=0; ii<30; ii++) weight[ii]=0.;
  for(int ibin=1; ibin<=h_nPU_data->GetNbinsX(); ibin++){
    if(h_nPU_mc->GetBinContent(ibin) > 0) {
      weight[ibin-1] = h_nPU_data->GetBinContent(ibin)/h_nPU_mc->GetBinContent(ibin);
    }
  }

  for(int ii=0; ii<30; ii++){
    printf("bin %d, weight %f \n", ii, weight[ii]);
  }
  fdata->Close();
  fmc->Close();

  int count1=0;
  int count2=0;
  int count3=0;


//   Float_t rhofracbad = 0.52, rhofrac = 0.17, combIso = 0, trkIso = 0;
//   Float_t rhofracbad = 0.52, rhofrac = 0.157/0.6, combIso = 0, trkIso = 0;
  float  rhofrac = 0.157/0.6;
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   printf("%d entries to be loaded \n", nentries);
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries; ++jentry) {
//     for (Long64_t jentry=0; jentry<20000; ++jentry) {
//      for (Long64_t jentry=0; jentry<100; ++jentry) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
     
     int hltbit=0;
     if(isData==1){
       //photon75_CaloIDVL_IsoLv* 
       if(run>=160431 && run<160954  && HLTIndex[29]>0){
	 if( HLT[HLTIndex[29]]==1 ) hltbit=1;
       }
       //photon75_CaloIDVL_v* 160954 - 163870
       if(run>=160954 && run <= 163870 && HLTIndex[30]>0){
	 if( HLT[HLTIndex[30]]==1 ) hltbit=2;
       }
       //photon75_CaloIDVL_IsoLv* 
       if(run>=165088 && run <= 165633 && HLTIndex[29]>0){
	 if( HLT[HLTIndex[29]]==1 ) hltbit=1;
       }
       //photon90_CaloIDVL_v* 
       if(run>165633 && run <= 167913 && HLTIndex[15]>0){
	 if( HLT[HLTIndex[15]]==1 ) hltbit=3;
       }
       //photon135
       if(run>=167039 && HLTIndex[31]>0){
	 if( HLT[HLTIndex[31]]==1 ) hltbit=4;
       }
     }

      if (isData == 1  &&  hltbit == 0 ) continue;
//       if (isData == 1  &&  hltbit == 1 ) continue; //for 70_noiso hlt
//       if (isData == 1  &&  hltbit == 2 ) continue;
//       if (isData == 1  &&  hltbit == 3 ) continue;
//       if (isData == 1  &&  hltbit == 4 ) continue;

     float ww = 1.;
     //--mc_mask if(isData!=1) {
     //--mc_mask   for(int iv=0; iv<3; iv++){
//--mc_mask 	 if(puBX[iv]==0 && nPU[iv]<30){
//--mc_mask 	   ww = weight[nPU[iv]];     
//--mc_mask 	   h_nPU->Fill(nPU[iv]);
//--mc_mask 	   h_nVtx->Fill(nVtx);
//--mc_mask 	   h_nVtx_w->Fill(nVtx,ww);
//--mc_mask 	 }
      //--mc_mask  }
     //--mc_mask }else if(isData==1) {
       h_nVtx->Fill(nVtx);
       h_nVtx_w->Fill(nVtx,ww);
     //--mc_mask }

       float Lint=1000.;
       //MC_job_gjet_80to120.root ww *= Lint*447.2/2046637.;
       //MC_job_gjet_120to170.root ww *= Lint *84.17/2.088216e+06;
       //MC_job_gjet_170to300.root ww *= Lint *22.64/1.869161e+06;
       //MC_job_gjet_300to470.root ww *= Lint *1.493/1.07688e+06;
       //MC_job_gjet_470to800.root ww *= Lint *0.1323/1.987212e+06;
       //MC_job_gjet_800to1000.root ww *= Lint *0.003481/4.050576e+06;

       //MC_job_qcd_80to120.root ww *= Lint *7.843e+05/6.589956e+06;
       //MC_job_qcd_120to170.root ww *= Lint*1.151e+05/6.127528e+06;
       //MC_job_qcd_170to300.root ww *= Lint*2.426e+04/5.72016e+06;
       //MC_job_qcd_300to470.root ww *= Lint *1168./6.432669e+06;
       //MC_job_qcd_470to600.root ww *= Lint *70.22/3.888286e+06;
       //MC_job_qcd_600to800.root ww *= Lint *15.55/4.240592e+06;
  


     int npho_EB=0;     int npho_EE=0;
     int njet=0;
     int MC_type=0;
     int EBEE_type=0;
     for(int ipho=0; ipho<nPho; ipho++){
       int photonid_3iso = PhotonID_3Iso(ipho);
       if(photonid_3iso != 0) continue; //skip if not passing photon ID
       h_phoEt_all->Fill(phoEt[ipho],ww);
       h_phoEta->Fill(phoEta[ipho],ww);
       h_phoPhi->Fill(phoPhi[ipho],ww);
     }

     for(int ipho=0; ipho<nPho; ipho++){
       //--mc_mask if(TMath::Abs(phoGenMomPID[ipho])>22) MC_type=1; else MC_type=0;

       //for pho+jet       
       //phojet_mask if(TMath::Abs(phoGenMomPID[ipho])>22) continue;
       //for qcd
       //qcd_mask if(TMath::Abs(phoGenMomPID[ipho])<=22) continue;

       if(TMath::Abs(phoSCEta[ipho]) >2.5) continue;
       if(TMath::Abs(phoSCEta[ipho]) >1.4442 && TMath::Abs(phoSCEta[ipho]) <1.556) continue;
       //temp cut
       if(TMath::Abs(phoEta[ipho])<0.9 || TMath::Abs(phoEta[ipho])>1.4) continue;


       int pho_eta_bin=0;
       if( TMath::Abs(phoSCEta[ipho]) > 1.5) {
	 pho_eta_bin=1;
	 EBEE_type=1;
       }else {
	 pho_eta_bin=0;
	 EBEE_type=0;
       }	 


       if(EBEECATE<2 && pho_eta_bin==1) continue;

       int cate = phocat(ipho);

       //-------------------------------------------------------------
       //photon Et correction
       Float_t Econst1 = 0; 
       Float_t Escale1 = 0;
       Float_t Esmear1= 1;
       Float_t Eshift = 0;
       if(isData == 1) {  // data
	  if(run >= 160431 && run <= 167913){
	    Eshift = Eoffset_jul5[cate];  // 160404 - 163869
          }else if(run >= 170249 && run <= 172619){
	    Eshift =  Eoffset_aug5[cate];  // 160404 - 163869
          }else if(run >= 172620 && run <= 172802){
	    Eshift = Eoffset_pr[cate];  // 160404 - 163869
          }
	  Escale1 = 0 - Eoffset[cate];
	  // Escale1 =   0 - (Eoffset[CiC_1] + sysEoffset[CiC_1]);	  
// 	  Escale1=0;     //!!!!!!!  no scale applied
       } 
       else{  // MC
	 EresPlus[cate]=Eres[cate]+sysEres[cate];
	 EresMinus[cate]=Eres[cate]-sysEres[cate];
	 // cout <<" Eres[cate]  with sys "<<EresPlus[cate]  <<endl;	 
	 Escale1 =  sysEoffset[cate];
	 
	 if(PHOER==1)  Econst1 = EresPlus[cate];
	 if(PHOER==-1) Econst1 = EresMinus[cate];
	 if(PHOER==0)  Econst1 = Eres[cate];
	 Float_t EscaleSmear = 1;
	 Esmear1= rnd->Gaus(EscaleSmear,Econst1);
       }
       float tmp = phoEt[ipho]/phoE[ipho];
       phoE[ipho] = (Escale1 + Esmear1) * (phoE[ipho] - Eshift);
       phoEt[ipho] = tmp*phoE[ipho];
       //end of photon Et correction-------------------------------------       


       if(phoEt[ipho]<80) continue;
//        if(phoEt[ipho]<150) continue;

       TLorentzVector g_P4;
       g_P4.SetPtEtaPhiM(phoEt[ipho], phoEta[ipho], phoPhi[ipho], 0.);

       int photonid = PhotonID(ipho);
       int photonid_3iso = PhotonID_3Iso(ipho);
       float combIso =   phoEcalIsoDR04[ipho]+phoHcalIsoDR04[ipho]+phoTrkIsoHollowDR04[ipho];       
       float sieie = phoSigmaIEtaIEta[ipho];

       float ecalIso = phoEcalIsoDR04[ipho]-rho25*0.110;
       float hcalIso = phoHcalIsoDR04[ipho]-rho25*0.037;
       float trkIso  = phoTrkIsoHollowDR04[ipho]-rho25*0.010;

       rhofrac = 0.157/0.6;
       if(pho_eta_bin==1){
	 rhofrac=0.181/0.6;
	 ecalIso = phoEcalIsoDR04[ipho]-rho25*0.054;     
	 hcalIso = phoHcalIsoDR04[ipho]-rho25*0.108;     
	 trkIso  = phoTrkIsoHollowDR04[ipho]-rho25*0.019;
       }
       combIso -= rhofrac*rho25;
       
       if(isData!=1) {
	 if(pho_eta_bin==0) {
	   combIso -= 0.09;
	   sieie -= 0.00003;
	 }else {
	   combIso -= 0.08;
 	   sieie -= 0.0003;
	 }
       }

       if( ( photonid & ~( 1<<combIso_bit | 1<<sieie_bit) )==0) {	 
	 h_combIso_sieie[EBEE_type][MC_type]->Fill(combIso,sieie,ww);
	 h_combIso_sieie_phoEt[EBEE_type][MC_type]->Fill(combIso,sieie, phoEt[ipho],ww);
	 
	 for(int ijet=0; ijet<nJet; ijet++){	   
	   if(jetPt[ijet]<30.) continue;
	   if(jetID(ijet)==0) continue;	   
	   TLorentzVector j_P4;
	   j_P4.SetPtEtaPhiE(jetPt[ijet], jetEta[ijet], jetPhi[ijet], jetEn[ijet]);

	   if(g_P4.DeltaR(j_P4)<0.4) {
	     continue;
	   }	   
// 	   h_jetPt->Fill(jetPt[ijet],ww);
	   TLorentzVector gj_P4 = g_P4;
	   gj_P4 += j_P4;
	   h_combIso_sieie_mgj[EBEE_type][MC_type]->Fill(combIso,sieie, gj_P4.M(),ww);
	   break; //break after leading jet (with ID)
	 }	 
       }

       if( isData==0 || (isData==1 && (hltbit==2 || hltbit==4)) ) {
	 if( ( photonid_3iso & ~(1<<ecalIso_bit) )==0 ) {	 
	   h_phoEt_ecalIso[EBEE_type][MC_type]->Fill(phoEt[ipho],ecalIso,ww);
	 }       
	 if( ( photonid_3iso & ~(1<<hcalIso_bit) )==0 ) {	 
	   h_phoEt_hcalIso[EBEE_type][MC_type]->Fill(phoEt[ipho],hcalIso,ww);
	 }
	 if( ( photonid_3iso & ~(1<<trkIso_bit) )==0 ) {	 
	   h_phoEt_trkIso[EBEE_type][MC_type]->Fill(phoEt[ipho],trkIso,ww);
	 }
       }

       if( ( photonid_3iso & ~( 1<<trkIso_bit | 1<<HoverE_bit) )==0) {	 
//  	   h_trkIso_HoverE[EBEE_type][MC_type]->Fill(trkIso,phoHoverE[ipho],ww);
//  	   h_trkIso_HoverE_phoEt[EBEE_type][MC_type]->Fill(trkIso,phoHoverE[ipho], phoEt[ipho],ww);
	 
	 for(int ijet=0; ijet<nJet; ijet++){	   
	   if(jetPt[ijet]<30.) continue;
	   if(jetID(ijet)==0) continue;	   
	   TLorentzVector j_P4;
	   j_P4.SetPtEtaPhiE(jetPt[ijet], jetEta[ijet], jetPhi[ijet], jetEn[ijet]);
	   if(g_P4.DeltaR(j_P4)<0.4) {
	     continue;
	   }	   
// 	   h_jetPt->Fill(jetPt[ijet],ww);
	   TLorentzVector gj_P4 = g_P4;
	   gj_P4 += j_P4;	   
// 	   h_trkIso_HoverE_mgj[EBEE_type][MC_type]->Fill(trkIso,phoHoverE[ipho], gj_P4.M(),ww);
	   break; //break after leading jet (with ID)
	 }	 
       }



//        if( (photonid & ~( 1<<sieie_bit )) == 0) {
       if( (photonid_3iso & ~( 1<<sieie_bit )) == 0) {
 	 if( isData==0 || 
	     (isData==1 && hltbit==1 && phoEt[ipho]>80) ||
	     (isData==1 && hltbit==2 && phoEt[ipho]>80) ||
	     (isData==1 && hltbit==3 && phoEt[ipho]>95) ||
 	     (isData==1 && hltbit==4 && phoEt[ipho]>140) ){
	   h_sieie[EBEE_type][MC_type]->Fill(sieie,ww);
	   h_phoEt_sieie[EBEE_type][MC_type]->Fill(phoEt[ipho],sieie,ww);
 	 }
       }

       if( (photonid & ~( 1<<combIso_bit )) == 0) {
	 if( isData==0 || (isData==1 && (hltbit==1 || hltbit==3)) ) {
	   h_combIso[EBEE_type][MC_type]->Fill(combIso,ww);
	   h_phoEt_combIso[EBEE_type][MC_type]->Fill(phoEt[ipho],combIso,ww);
	   
	   for(int ijet=0; ijet<nJet; ijet++){
	     if(jetPt[ijet]<30.) continue;
	     if(jetID(ijet)==0) continue;
	     TLorentzVector j_P4;
	     j_P4.SetPtEtaPhiE(jetPt[ijet], jetEta[ijet], jetPhi[ijet], jetEn[ijet]);
	     
	     // 	   if(jetPt[ijet]>100. && g_P4.DeltaR(j_P4)<0.1) {	     
	     // 	     h_combIso_jet[EBEE_type][MC_type]->Fill(combIso,ww);
	     // 	   }
	   }
	 }
       }

       if(photonid_3iso != 0) continue; //skip if not passing photon ID

       prof_phoEt_run->Fill(run,phoEt[ipho],ww);
       h_phoEt[EBEE_type][MC_type]->Fill(phoEt[ipho],ww);

       //count photon candidate
       if(pho_eta_bin==0) npho_EB++; else npho_EE++;
       //count jet candidate
       if((npho_EB+npho_EE)==1) {
	 for(int ijet=0; ijet<nJet; ijet++){
	   if(jetPt[ijet]<30.) continue;
	   if(jetID(ijet)==0) continue;
	   TLorentzVector j_P4;
	   j_P4.SetPtEtaPhiE(jetPt[ijet], jetEta[ijet], jetPhi[ijet], jetEn[ijet]);
	   if(g_P4.DeltaR(j_P4)>0.4) {
	     njet++;	   
	   }
	 }
       }

       //loop over jet and for photon+jet candidate
       for(int ijet=0; ijet<nJet; ijet++){

	 if(jetID(ijet)!=0) continue;

	 if(JER==1) jetPt[ijet] *= rnd->Gaus(1.,0.1);
	 //for JES	 
	 int neta = 0;
	 for(neta = 0; neta < NETAJEC; neta++) {
	   if(etajec[neta] > jetEta[ijet]) break;
	 }
	 neta = neta-1;
	 int npt = 0;
	 for(npt = 0; npt < NPTJEC; npt++) {
	   if(ptjec[neta][npt][0] > jetPt[ijet]) break;
	 }
	 npt = npt-1;
	 double jecuncertainty = ptjec[neta][npt][1];              
	 if(JES==1) {
// 	   printf(" --- move JES to +1 sigma \n");
	   jetPt[ijet] *= (1.+jecuncertainty);
	   jetEn[ijet] *= (1.+jecuncertainty);
	 }else if(JES==-1) {
// 	   printf(" --- move JES to -1 sigma \n");
	   jetPt[ijet] *= (1.-jecuncertainty);
	   jetEn[ijet] *= (1.-jecuncertainty);
	 }
	 //end of JES


	 if(jetPt[ijet]<30.) continue;
	 TLorentzVector j_P4;
	 j_P4.SetPtEtaPhiE(jetPt[ijet], jetEta[ijet], jetPhi[ijet], jetEn[ijet]);

	 h_gj_dR->Fill(g_P4.DeltaR(j_P4),ww);
	 if(g_P4.DeltaR(j_P4)<0.4) {
	   continue;
	 }
	 h_gj_dEta->Fill(TMath::Abs(g_P4.Eta()-j_P4.Eta()),ww);
	 h_gj_dPhi->Fill(g_P4.DeltaPhi(j_P4),ww);

	 h_jetPt->Fill(jetPt[ijet],ww);
	 TLorentzVector gj_P4 = g_P4;
	 gj_P4 += j_P4;
	 h_mgj_ptg[EBEE_type][MC_type]->Fill(gj_P4.M(), phoEt[ipho],ww);
// 	 h_mgj_ptj[EBEE_type][MC_type]->Fill(gj_P4.M(), jetPt[ipho],ww);

	 category = cate;
	 gjmass = gj_P4.M();
	 mcwei = ww;
	 tt->Fill();

	 break; //break after leading jet (with ID)
       }

       break; //break after leading photon (with ID)
     }

     for(int icat=0; icat<EBEECATE; icat++) {
       int nn = npho_EB;
       if(icat ==1 )  nn=npho_EE;
       h_npho[icat][MC_type]->Fill(nn,ww);
     }
     h_njet->Fill(njet,ww);


     
   }

   TFile *f = new TFile("ana_gj_out_FNAME_MASK.root","recreate");
   prof_phoEt_run->Write();

   h_nPU->Write();
   h_nVtx->Write();
   h_nVtx_w->Write();
   h_njet->Write();
   h_jetPt->Write();
   h_phoEt_all->Write();
   h_phoEta->Write();
   h_phoPhi->Write();
   h_gj_dR->Write();
   h_gj_dEta->Write();
   h_gj_dPhi->Write();

   tt->Write();


  for(int cate=0; cate<EBEECATE; cate++){
    for(int type=0; type<MCTYPE; type++){
      h_npho[cate][type]->Write();
      h_phoEt[cate][type]->Write(); 
      h_combIso[cate][type]->Write(); 
      h_sieie[cate][type]->Write(); 
      h_phoEt_combIso[cate][type]->Write(); 
      h_phoEt_ecalIso[cate][type]->Write(); 
      h_phoEt_hcalIso[cate][type]->Write(); 
      h_phoEt_trkIso[cate][type]->Write(); 
      h_phoEt_sieie[cate][type]->Write(); 
//       h_combIso_jet[cate][type]->Write(); 
      h_mgj_ptg[cate][type]->Write(); 
//       h_mgj_ptj[cate][type]->Write(); 
      h_combIso_sieie[cate][type]->Write(); 
      h_combIso_sieie_phoEt[cate][type]->Write(); 
      h_combIso_sieie_mgj[cate][type]->Write(); 
//       h_trkIso_HoverE[cate][type]->Write(); 
//       h_trkIso_HoverE_phoEt[cate][type]->Write(); 
//       h_trkIso_HoverE_mgj[cate][type]->Write(); 
    }
  }


   f->Close();

  delete rnd;
}



Int_t Ana_gj::PhotonID(int ipho) {

  Int_t passID = 0;
  if(phoSigmaIEtaIEta[ipho]<0.002) return 127;
  if(phoSigmaIPhiIPhi[ipho]<0.002) return 127;

  if(phoHoverE[ipho]>0.05) passID |= 1<<HoverE_bit ;
  if(phohasPixelSeed[ipho]==1) passID |= 1<<Pixel_bit;

  Float_t rhofrac = 0.157/0.6;    
  Float_t combIso = phoEcalIsoDR04[ipho]+phoHcalIsoDR04[ipho]+phoTrkIsoHollowDR04[ipho];
  Float_t sieie = phoSigmaIEtaIEta[ipho];

  if( TMath::Abs(phoSCEta[ipho]) < 1.4442) {
    combIso -= rho25*rhofrac;
    if(isData!=1) {
      combIso -= 0.09;
      sieie -= 0.00003;
    }
    if(sieie>0.01) passID |= 1<<sieie_bit;
    if(combIso > 4.) passID |= 1<<combIso_bit;
  }else {
    rhofrac = 0.181/0.6;
    combIso -= rho25*rhofrac;
    if(isData!=1) {
      combIso -= 0.08;
      sieie -= 0.0003;
    }
    if(sieie>0.03) passID |= 1<<sieie_bit;
    if(combIso > 4.) passID |= 1<<combIso_bit;
  }

//   printf("photon eta %f, sieie %f, cutbit %d\n",  phoSCEta[ipho],  sieie, passID);

  return passID;
}

Int_t Ana_gj::PhotonID_3Iso(int ipho) {

  Int_t passID = 0;
  if(phoSigmaIEtaIEta[ipho]<0.002) return 127;
  if(phoSigmaIPhiIPhi[ipho]<0.002) return 127;

  if(phoHoverE[ipho]>0.05) passID |= 1<<HoverE_bit ;
  if(phohasPixelSeed[ipho]==1) passID |= 1<<Pixel_bit;

  Float_t rhofrac = 0.157/0.6;    
  Float_t sieie = phoSigmaIEtaIEta[ipho];

  if( TMath::Abs(phoSCEta[ipho]) < 1.4442) {
    if(isData!=1) {
      sieie -= 0.00003;
    }
    if(sieie>0.01) passID |= 1<<sieie_bit;
    if((phoEcalIsoDR04[ipho]-rho25*0.110) > 4.2*0.006*phoEt[ipho] ) passID |= 1<<ecalIso_bit;
    if((phoHcalIsoDR04[ipho]-rho25*0.037) > 2.2*0.0025*phoEt[ipho] ) passID |= 1<<hcalIso_bit;
    if((phoTrkIsoHollowDR04[ipho]-rho25*0.010) > 2.*0.001*phoEt[ipho] ) passID |= 1<<trkIso_bit;

  }else {
    rhofrac = 0.181/0.6;
    if(isData!=1) {
      sieie -= 0.0003;
    }
    if(sieie>0.03) passID |= 1<<sieie_bit;
    if((phoEcalIsoDR04[ipho]-rho25*0.054) > 4.2*0.006*phoEt[ipho] ) passID |= 1<<ecalIso_bit;
    if((phoHcalIsoDR04[ipho]-rho25*0.108) > 2.2*0.0025*phoEt[ipho] ) passID |= 1<<hcalIso_bit;
    if((phoTrkIsoHollowDR04[ipho]-rho25*0.019) > 2.*0.001*phoEt[ipho] ) passID |= 1<<trkIso_bit;

  }

//   printf("photon eta %f, sieie %f, cutbit %d\n",  phoSCEta[ipho],  sieie, passID);

  return passID;
}



Int_t Ana_gj::jetID(int ijet){
  
  Int_t passID=0;
  if(TMath::Abs(jetEta[ijet])>2.5) return 99;

  if(TMath::Abs(jetEta[ijet])<2.4){
    if(jetNHF[ijet]>=0.99) passID=1;
    if(jetNEF[ijet]>=0.99) passID=1;
    if(jetNConstituents[ijet]<=1) passID=1;
  }else {
    if(jetCHF[ijet]<=0.) passID=3; 
    if(jetNCH[ijet]<=0.) passID=3; 
    if(jetCEF[ijet]>=0.99) passID=3; 
  }
  return passID;
}

Int_t Ana_gj::eID(int iele) {
  Int_t passID=0;
  if(TMath::Abs(eleSCEta[iele])>2.5) return 100;
  if(TMath::Abs(eleSCEta[iele])>1.4442 && TMath::Abs(eleSCEta[iele])<1.556) return 100;     

  if(eleSigmaIEtaIEta[iele]<0.0002) return 100;
  if(eleSigmaIPhiIPhi[iele]<0.0002) return 100;

  int phoindex=-1;
  for(int ipho=0; ipho<nPho; ipho++){
    if(TMath::Abs(phoSCEta[ipho]-eleSCEta[iele])<0.1 &&
       TMath::Abs(phoSCPhi[ipho]-eleSCPhi[iele])<0.1) {
      phoindex = ipho;
      break;
    }
  }

  if(phoSigmaIEtaIEta[phoindex]<0.0002) return 100;
  if(phoSigmaIPhiIPhi[phoindex]<0.0002) return 100;


  float combIso = 100.;
  if(phoindex>=0) {
    combIso =  phoEcalIsoDR04[phoindex]+phoHcalIsoDR04[phoindex]+phoTrkIsoHollowDR04[phoindex];
  }

  float rhofrac = 0.157/0.6;//0.122; 
  if( TMath::Abs(eleSCEta[iele]) < 1.4442 ) {
    combIso -= rhofrac*rho25;
    //if(combIso/elePt[iele] >0.053) return 0;
    if(combIso>4.) passID |= 1<<0;
    if(eleSigmaIEtaIEta[iele]>0.01) passID |= 1<<1;
    if(TMath::Abs(eledEtaAtVtx[iele])>0.039) return 100;
    if(TMath::Abs(eledPhiAtVtx[iele])>0.005) return 100;
  }else if ( TMath::Abs(eleSCEta[iele]) >1.556 && 
	     TMath::Abs(eleSCEta[iele]) <2.5 ) {
    rhofrac = 0.181/0.6;
    combIso -= rhofrac*rho25;
    //if(combIso/elePt[iele] >0.042) return 0;
    if(combIso>4.) passID |= 1<<0;
    if(eleSigmaIEtaIEta[iele]>0.031) passID |= 1<<1;
    if(TMath::Abs(eledEtaAtVtx[iele])>0.028) return 100;
    if(TMath::Abs(eledPhiAtVtx[iele])>0.007) return 100;    
  }else{
    return 100;
  }
  return passID;
}
Int_t Ana_gj::phocat(int ipho) {

  Int_t r9 = (phoR9[ipho] > 0.94) ? 0 : 1;
  Int_t eta = (fabs(phoSCEta[ipho]) < 1.479) ? 0 : 1;

  return r9 + 2*eta;
}
