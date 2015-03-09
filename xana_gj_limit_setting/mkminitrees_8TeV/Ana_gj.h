//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 16 21:05:03 2011 by ROOT version 5.28/00
// from TTree EventTree/Event data
// found on file: /home/cmkuo/job_photon_2011a_Jul5rereco_diphoton_Jul6_JSON_specialskim_4283.root
//////////////////////////////////////////////////////////

#ifndef Ana_gj_h
#define Ana_gj_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
//======= LP ==================
float sysEres[4]={0.23e-2,0.42e-2,0.51e-2,0.38e-2};
float Eres[4] = {1.0e-2, 1.69e-2, 3.03e-2, 2.96e-2};

float EresPlus[4]={0,0,0,0};
float EresMinus[4]={0,0,0,0};
//
float Eoffset[4] = { 0.02e-2, 0.03e-2,  0.0e-2, 0.02e-2};

float sysEoffset[4] = { 0.09e-2, 0.44e-2,  0.29e-2, 0.41e-2};

float Eoffset_jul5[4] = {0.094, -0.06, 1.641, 0.150};
float Eoffset_aug5[4] = {-0.0846, -0.1263, 2.57, 0.974};
float Eoffset_pr[4] = {-0.225, -0.266, 3.54, 2.03};
//======= LP ==================


   const Int_t kMaxeleESEffSigmaRR = 1;
   const Int_t kMaxeleESE1 = 1;
   const Int_t kMaxeleESE3 = 1;
   const Int_t kMaxeleESE5 = 1;
   const Int_t kMaxeleESE7 = 1;
   const Int_t kMaxeleESE11 = 1;
   const Int_t kMaxeleESE21 = 1;
   const Int_t kMaxphoESEffSigmaRR = 1;
   const Int_t kMaxphoESE1 = 1;
   const Int_t kMaxphoESE3 = 1;
   const Int_t kMaxphoESE5 = 1;
   const Int_t kMaxphoESE7 = 1;
   const Int_t kMaxphoESE11 = 1;
   const Int_t kMaxphoESE21 = 1;
   const Int_t kMaxnPFPho = 627;
   const Int_t kMaxPFPhoE = 1;
   const Int_t kMaxPFPhoEt = 1;
   const Int_t kMaxPFPhoEta = 1;
   const Int_t kMaxPFPhoPhi = 1;
   const Int_t kMaxPFPhoIsConv = 1;
   const Int_t kMaxPFPhoTrkIsoHollowDR04 = 1;
   const Int_t kMaxPFPhoEcalIsoDR04 = 1;
   const Int_t kMaxPFPhoHcalIsoDR04 = 1;
   const Int_t kMaxPFPhoHoverE = 1;
   const Int_t kMaxnPFEle = 5;
   const Int_t kMaxPFElePt = 1;
   const Int_t kMaxPFEleEta = 1;
   const Int_t kMaxPFElePhi = 1;
   const Int_t kMaxPFEleEn = 1;
   const Int_t kMaxPFEleCharge = 1;
   const Int_t kMaxnPFMu = 6;
   const Int_t kMaxPFmuCharge = 1;
   const Int_t kMaxPFmuPhi = 1;
   const Int_t kMaxPFmuEta = 1;
   const Int_t kMaxPFmuPt = 1;
   const Int_t kMaxPFmuPz = 1;
   const Int_t kMaxPFmuChambers = 1;
   const Int_t kMaxPFmuD0 = 1;
   const Int_t kMaxPFmuDz = 1;
   const Int_t kMaxPFmuChi2NDF = 1;
   const Int_t kMaxPFmuNumberOfValidTrkHits = 1;
   const Int_t kMaxPFmuNumberOfValidPixelHits = 1;
   const Int_t kMaxPFmuNumberOfValidMuonHits = 1;

class Ana_gj {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           orbit;
   Int_t           bx;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nHLT;
   Int_t           HLT[385];   //[nHLT]
   Int_t           HLTIndex[50];
   Float_t         bspotPos[3];
   Int_t           nVtx;
   Float_t         vtx[31][3];   //[nVtx]
   Int_t           vtxNTrk[31];   //[nVtx]
   Int_t           vtxNDF[31];   //[nVtx]
   Float_t         vtxD0[31];   //[nVtx]
   Int_t           IsVtxGood;
   Int_t           nVtxBS;
   Float_t         vtxbs[31][3];   //[nVtxBS]
   Float_t         vtxbsPtMod[31];   //[nVtxBS]
   Float_t         vtxbsSumPt2[31];   //[nVtxBS]
   vector<vector<int> > *vtxbsTkIndex;
   vector<vector<float> > *vtxbsTkWeight;
   Int_t           nTrk;
   Int_t           nGoodTrk;
   Int_t           IsTracksGood;
   Float_t         MET;
   Float_t         METx;
   Float_t         METy;
   Float_t         METPhi;
   Float_t         METsumEt;
   Float_t         tcMET;
   Float_t         tcMETx;
   Float_t         tcMETy;
   Float_t         tcMETPhi;
   Float_t         tcMETsumEt;
   Float_t         tcMETmEtSig;
   Float_t         tcMETSig;
   Float_t         pfMET;
   Float_t         pfMETx;
   Float_t         pfMETy;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Int_t           nEle;
   Int_t           eleTrg[9][13];   //[nEle]
   Int_t           eleID[9][12];   //[nEle]
   Int_t           eleClass[9];   //[nEle]
   Int_t           eleCharge[9];   //[nEle]
   Float_t         eleEn[9];   //[nEle]
   Float_t         eleSCRawEn[9];   //[nEle]
   Float_t         eleSCEn[9];   //[nEle]
   Float_t         elePt[9];   //[nEle]
   Float_t         eleEta[9];   //[nEle]
   Float_t         elePhi[9];   //[nEle]
   Float_t         eleSCEta[9];   //[nEle]
   Float_t         eleSCPhi[9];   //[nEle]
   Float_t         eleSCEtaWidth[9];   //[nEle]
   Float_t         eleSCPhiWidth[9];   //[nEle]
   Float_t         eleVtx[9][3];   //[nEle]
   Float_t         eleD0[9];   //[nEle]
   Float_t         eleDz[9];   //[nEle]
   Float_t         eleHoverE[9];   //[nEle]
   Float_t         eleEoverP[9];   //[nEle]
   Float_t         elePin[9];   //[nEle]
   Float_t         elePout[9];   //[nEle]
   Float_t         eleBrem[9];   //[nEle]
   Float_t         eledEtaAtVtx[9];   //[nEle]
   Float_t         eledPhiAtVtx[9];   //[nEle]
   Float_t         eleSigmaIEtaIEta[9];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[9];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[9];   //[nEle]
   Float_t         eleEmax[9];   //[nEle]
   Float_t         eleE3x3[9];   //[nEle]
   Float_t         eleE5x5[9];   //[nEle]
   Float_t         eleE2x5Right[9];   //[nEle]
   Float_t         eleE2x5Left[9];   //[nEle]
   Float_t         eleE2x5Top[9];   //[nEle]
   Float_t         eleE2x5Bottom[9];   //[nEle]
   Float_t         eleSeedTime[9];   //[nEle]
   Int_t           eleRecoFlag[9];   //[nEle]
   Float_t         eleIsoTrkDR03[9];   //[nEle]
   Float_t         eleIsoEcalDR03[9];   //[nEle]
   Float_t         eleIsoHcalDR03[9];   //[nEle]
   Float_t         eleIsoHcalSolidDR03[9];   //[nEle]
   Float_t         eleIsoTrkDR04[9];   //[nEle]
   Float_t         eleIsoEcalDR04[9];   //[nEle]
   Float_t         eleIsoHcalDR04[9];   //[nEle]
   Float_t         eleIsoHcalSolidDR04[9];   //[nEle]
   Int_t           eleMissHits[9];   //[nEle]
   Float_t         eleConvDist[9];   //[nEle]
   Float_t         eleConvDcot[9];   //[nEle]
   Float_t         eleESEffSigmaRR[9][3];   //[nEle]
   Float_t         eleESE1[9][2];   //[nEle]
   Float_t         eleESE3[9][2];   //[nEle]
   Float_t         eleESE5[9][2];   //[nEle]
   Float_t         eleESE7[9][2];   //[nEle]
   Float_t         eleESE11[9][2];   //[nEle]
   Float_t         eleESE21[9][2];   //[nEle]
   Float_t         eleBrLinear[9];   //[nEle]
   Float_t         eleCetaCorrE[9];   //[nEle]
   Float_t         eleCetaCorrEt[9];   //[nEle]
   Float_t         eleBremCorrE[9];   //[nEle]
   Float_t         eleBremCorrEt[9];   //[nEle]
   Float_t         eleFullCorrE[9];   //[nEle]
   Float_t         eleFullCorrEt[9];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[38][8];   //[nPho]
   Bool_t          phoIsPhoton[38];   //[nPho]
   Float_t         phoCaloPos[38][3];   //[nPho]
   Float_t         phoE[38];   //[nPho]
   Float_t         phoEt[38];   //[nPho]
   Float_t         phoEta[38];   //[nPho]
   Float_t         phoVtx[38][3];   //[nPho]
   Float_t         phoPhi[38];   //[nPho]
   Float_t         phoEtVtx[38][50];   //[nPho]
   Float_t         phoEtaVtx[38][50];   //[nPho]
   Float_t         phoPhiVtx[38][50];   //[nPho]
   Float_t         phoR9[38];   //[nPho]
   Float_t         phoCetaCorrE[38];   //[nPho]
   Float_t         phoCetaCorrEt[38];   //[nPho]
   Float_t         phoBremCorrE[38];   //[nPho]
   Float_t         phoBremCorrEt[38];   //[nPho]
   Float_t         phoFullCorrE[38];   //[nPho]
   Float_t         phoFullCorrEt[38];   //[nPho]
   Float_t         phoTrkIsoSolidDR03[38];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[38];   //[nPho]
   Float_t         phoEcalIsoDR03[38];   //[nPho]
   Float_t         phoHcalIsoDR03[38];   //[nPho]
   Float_t         phoHcalIsoSolidDR03[38];   //[nPho]
   Float_t         phoTrkIsoSolidDR04[38];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[38];   //[nPho]
   Float_t         phoTrkIsoHollowNoDzDR03[38];   //[nPho]
   Float_t         phoTrkIsoHollowNoDzDR04[38];   //[nPho]
   Float_t         phoCiCTrkIsoDR03[38][50];   //[nPho]
   Float_t         phoCiCTrkIsoDR04[38][50];   //[nPho]
   Float_t         phoCiCdRtoTrk[38];   //[nPho]
   Float_t         phoEcalIsoDR04[38];   //[nPho]
   Float_t         phoHcalIsoDR04[38];   //[nPho]
   Float_t         phoHcalIsoSolidDR04[38];   //[nPho]
   Float_t         phoHoverE[38];   //[nPho]
   Float_t         phoSigmaIEtaIEta[38];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[38];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[38];   //[nPho]
   Float_t         phoEmax[38];   //[nPho]
   Float_t         phoE3x3[38];   //[nPho]
   Float_t         phoE5x5[38];   //[nPho]
   Float_t         phoE2x5Right[38];   //[nPho]
   Float_t         phoE2x5Left[38];   //[nPho]
   Float_t         phoE2x5Top[38];   //[nPho]
   Float_t         phoE2x5Bottom[38];   //[nPho]
   Float_t         phoSeedTime[38];   //[nPho]
   Int_t           phoSeedDetId1[38];   //[nPho]
   Int_t           phoSeedDetId2[38];   //[nPho]
   Int_t           phoRecoFlag[38];   //[nPho]
   Int_t           phoPos[38];   //[nPho]
   Float_t         phoSCE[38];   //[nPho]
   Float_t         phoSCRawE[38];   //[nPho]
   Float_t         phoESEn[38];   //[nPho]
   Float_t         phoSCEt[38];   //[nPho]
   Float_t         phoSCEta[38];   //[nPho]
   Float_t         phoSCPhi[38];   //[nPho]
   Float_t         phoSCEtaWidth[38];   //[nPho]
   Float_t         phoSCPhiWidth[38];   //[nPho]
   Float_t         phoSCBrem[38];   //[nPho]
   Int_t           phoOverlap[38];   //[nPho]
   Int_t           phohasPixelSeed[38];   //[nPho]
   Int_t           phoIsConv[38];   //[nPho]
   Int_t           phoNConv[38];   //[nPho]
   Float_t         phoConvInvMass[38];   //[nPho]
   Float_t         phoConvCotTheta[38];   //[nPho]
   Float_t         phoConvEoverP[38];   //[nPho]
   Float_t         phoConvZofPVfromTrks[38];   //[nPho]
   Float_t         phoConvMinDist[38];   //[nPho]
   Float_t         phoConvdPhiAtVtx[38];   //[nPho]
   Float_t         phoConvdPhiAtCalo[38];   //[nPho]
   Float_t         phoConvdEtaAtCalo[38];   //[nPho]
   Float_t         phoConvTrkd0[38][2];   //[nPho]
   Float_t         phoConvTrkPin[38][2];   //[nPho]
   Float_t         phoConvTrkPout[38][2];   //[nPho]
   Float_t         phoConvTrkdz[38][2];   //[nPho]
   Float_t         phoConvTrkdzErr[38][2];   //[nPho]
   Float_t         phoConvChi2[38];   //[nPho]
   Float_t         phoConvChi2Prob[38];   //[nPho]
   Float_t         phoConvCharge[38][2];   //[nPho]
   Float_t         phoConvValidVtx[38];   //[nPho]
   Float_t         phoConvLikeLihood[38];   //[nPho]
   Float_t         phoConvP4[38][4];   //[nPho]
   Float_t         phoConvVtx[38][3];   //[nPho]
   Float_t         phoConvVtxErr[38][3];   //[nPho]
   Float_t         phoConvPairMomentum[38][3];   //[nPho]
   Float_t         phoConvRefittedMomentum[38][3];   //[nPho]
   Float_t         phoESEffSigmaRR[38][3];   //[nPho]
   Float_t         phoESE1[38][2];   //[nPho]
   Float_t         phoESE3[38][2];   //[nPho]
   Float_t         phoESE5[38][2];   //[nPho]
   Float_t         phoESE7[38][2];   //[nPho]
   Float_t         phoESE11[38][2];   //[nPho]
   Float_t         phoESE21[38][2];   //[nPho]
/*    Int_t           nMu; */
/*    Int_t           muTrg[46][6];   //[nMu] */
/*    Float_t         muEta[46];   //[nMu] */
/*    Float_t         muPhi[46];   //[nMu] */
/*    Int_t           muCharge[46];   //[nMu] */
/*    Float_t         muPt[46];   //[nMu] */
/*    Float_t         muPz[46];   //[nMu] */
/*    Float_t         muIsoTrk[46];   //[nMu] */
/*    Float_t         muIsoCalo[46];   //[nMu] */
/*    Float_t         muIsoEcal[46];   //[nMu] */
/*    Float_t         muIsoHcal[46];   //[nMu] */
/*    Float_t         muChi2NDF[46];   //[nMu] */
/*    Float_t         muEmVeto[46];   //[nMu] */
/*    Float_t         muHadVeto[46];   //[nMu] */
/*    Int_t           muType[46];   //[nMu] */
/*    Bool_t          muID[46][6];   //[nMu] */
/*    Float_t         muD0[46];   //[nMu] */
/*    Float_t         muDz[46];   //[nMu] */
/*    Int_t           muNumberOfValidTrkHits[46];   //[nMu] */
/*    Int_t           muNumberOfValidPixelHits[46];   //[nMu] */
/*    Int_t           muNumberOfValidMuonHits[46];   //[nMu] */
/*    Int_t           muStations[46];   //[nMu] */
/*    Int_t           muChambers[46];   //[nMu] */
/*    Int_t           nPFPho_; */
/*    Float_t         PFPhoE_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoEt_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoEta_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoPhi_[kMaxnPFPho];   //[nPFPho_] */
/*    Int_t           PFPhoIsConv_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoTrkIsoHollowDR04_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoEcalIsoDR04_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoHcalIsoDR04_[kMaxnPFPho];   //[nPFPho_] */
/*    Float_t         PFPhoHoverE_[kMaxnPFPho];   //[nPFPho_] */
/*    Int_t           nPFEle_; */
/*    Float_t         PFElePt_[kMaxnPFEle];   //[nPFEle_] */
/*    Float_t         PFEleEta_[kMaxnPFEle];   //[nPFEle_] */
/*    Float_t         PFElePhi_[kMaxnPFEle];   //[nPFEle_] */
/*    Float_t         PFEleEn_[kMaxnPFEle];   //[nPFEle_] */
/*    Int_t           PFEleCharge[kMaxnPFEle];   //[nPFEle_] */
/*    Int_t           nPFMu_; */
/*    Int_t           PFmuCharge_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuPhi_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuEta_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuPt_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuPz_[kMaxnPFMu];   //[nPFMu_] */
/*    Int_t           PFmuChambers_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuD0_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuDz_[kMaxnPFMu];   //[nPFMu_] */
/*    Float_t         PFmuChi2NDF_[kMaxnPFMu];   //[nPFMu_] */
/*    Int_t           PFmuNumberOfValidTrkHits_[kMaxnPFMu];   //[nPFMu_] */
/*    Int_t           PFmuNumberOfValidPixelHits_[kMaxnPFMu];   //[nPFMu_] */
/*    Int_t           PFmuNumberOfValidMuonHits_[kMaxnPFMu];   //[nPFMu_] */
   Float_t         rho25;
   Float_t         rho44;
   Int_t           nJet;
   Int_t           jetTrg[34][14];   //[nJet]
   Float_t         jetEn[34];   //[nJet]
   Float_t         jetPt[34];   //[nJet]
   Float_t         jetEta[34];   //[nJet]
   Float_t         jetPhi[34];   //[nJet]
   Float_t         jetEt[34];   //[nJet]
   Float_t         jetRawPt[34];   //[nJet]
   Float_t         jetRawEn[34];   //[nJet]
   Float_t         jetCHF[34];   //[nJet]
   Float_t         jetNHF[34];   //[nJet]
   Float_t         jetCEF[34];   //[nJet]
   Float_t         jetNEF[34];   //[nJet]
   Int_t           jetNCH[34];   //[nJet]
   Float_t         jetHFHAE[34];   //[nJet]
   Float_t         jetHFEME[34];   //[nJet]
   Int_t           jetNConstituents[34];   //[nJet]
   Float_t         jetTrackCountHiEffBJetTags[34];   //[nJet]
   Float_t         jetTrackCountHiPurBJetTags[34];   //[nJet]
   Float_t         jetSimpleSVHiEffBJetTags[34];   //[nJet]
   Float_t         jetSimpleSVHiPurBJetTags[34];   //[nJet]
   Int_t           nZee;
   Float_t         ZeeMass[36];   //[nZee]
   Float_t         ZeePt[36];   //[nZee]
   Float_t         ZeeEta[36];   //[nZee]
   Float_t         ZeePhi[36];   //[nZee]
   Int_t           ZeeLeg1Index[36];   //[nZee]
   Int_t           ZeeLeg2Index[36];   //[nZee]
   Int_t           nZmumu;
   Float_t         ZmumuMass[6];   //[nZmumu]
   Float_t         ZmumuPt[6];   //[nZmumu]
   Float_t         ZmumuEta[6];   //[nZmumu]
   Float_t         ZmumuPhi[6];   //[nZmumu]
   Int_t           ZmumuLeg1Index[6];   //[nZmumu]
   Int_t           ZmumuLeg2Index[6];   //[nZmumu]
/*    Int_t           nWenu; */
/*    Float_t         WenuMassTPfMET[9];   //[nWenu] */
/*    Float_t         WenuEtPfMET[9];   //[nWenu] */
/*    Float_t         WenuACopPfMET[9];   //[nWenu] */
/*    Int_t           WenuEleIndex[9];   //[nWenu] */
/*    Int_t           nWmunu; */
/*    Float_t         WmunuMassTPfMET[6];   //[nWmunu] */
/*    Float_t         WmunuEtPfMET[6];   //[nWmunu] */
/*    Float_t         WmunuACopPfMET[6];   //[nWmunu] */
/*    Int_t           WmunuMuIndex[6];   //[nWmunu] */
/*    Int_t           nConv; */
/*    Float_t         convP4[305][4];   //[nConv] */
/*    Float_t         convVtx[305][3];   //[nConv] */
/*    Float_t         convVtxErr[305][3];   //[nConv] */
/*    Float_t         convPairMomentum[305][3];   //[nConv] */
/*    Float_t         convRefittedMomentum[305][3];   //[nConv] */
/*    Int_t           convNTracks[305];   //[nConv] */
/*    Float_t         convPairInvMass[305];   //[nConv] */
/*    Float_t         convPairCotThetaSep[305];   //[nConv] */
/*    Float_t         convEoverP[305];   //[nConv] */
/*    Float_t         convDistOfMinApproach[305];   //[nConv] */
/*    Float_t         convDPhiTrksAtVtx[305];   //[nConv] */
/*    Float_t         convDPhiTrksAtEcal[305];   //[nConv] */
/*    Float_t         convDEtaTrksAtEcal[305];   //[nConv] */
/*    Float_t         convDxy[305];   //[nConv] */
/*    Float_t         convDz[305];   //[nConv] */
/*    Float_t         convLxy[305];   //[nConv] */
/*    Float_t         convLz[305];   //[nConv] */
/*    Float_t         convZofPrimVtxFromTrks[305];   //[nConv] */
/*    Int_t           convNHitsBeforeVtx[305][2];   //[nConv] */
/*    Int_t           convNSharedHits[305];   //[nConv] */
/*    Int_t           convValidVtx[305];   //[nConv] */
/*    Float_t         convMVALikelihood[305];   //[nConv] */
/*    Float_t         convChi2[305];   //[nConv] */
/*    Float_t         convChi2Probability[305];   //[nConv] */
/*    Float_t         convTk1Dz[305];   //[nConv] */
/*    Float_t         convTk2Dz[305];   //[nConv] */
/*    Float_t         convTk1DzErr[305];   //[nConv] */
/*    Float_t         convTk2DzErr[305];   //[nConv] */
/*    Int_t           convCh1Ch2[305];   //[nConv] */
/*    Float_t         convTk1D0[305];   //[nConv] */
/*    Float_t         convTk1Pout[305];   //[nConv] */
/*    Float_t         convTk1Pin[305];   //[nConv] */
/*    Float_t         convTk2D0[305];   //[nConv] */
/*    Float_t         convTk2Pout[305];   //[nConv] */
/*    Float_t         convTk2Pin[305];   //[nConv] */

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_bspotPos;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vtxNTrk;   //!
   TBranch        *b_vtxNDF;   //!
   TBranch        *b_vtxD0;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nVtxBS;   //!
   TBranch        *b_vtxbs;   //!
   TBranch        *b_vtxbsPtMod;   //!
   TBranch        *b_vtxbsSumPt2;   //!
   TBranch        *b_vtxbsTkIndex;   //!
   TBranch        *b_vtxbsTkWeight;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_nGoodTrk;   //!
   TBranch        *b_IsTracksGood;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METsumEt;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_tcMETx;   //!
   TBranch        *b_tcMETy;   //!
   TBranch        *b_tcMETPhi;   //!
   TBranch        *b_tcMETsumEt;   //!
   TBranch        *b_tcMETmEtSig;   //!
   TBranch        *b_tcMETSig;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETx;   //!
   TBranch        *b_pfMETy;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleID;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleEmax;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleE2x5Right;   //!
   TBranch        *b_eleE2x5Left;   //!
   TBranch        *b_eleE2x5Top;   //!
   TBranch        *b_eleE2x5Bottom;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalSolidDR03;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalSolidDR04;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_eleESE1;   //!
   TBranch        *b_eleESE3;   //!
   TBranch        *b_eleESE5;   //!
   TBranch        *b_eleESE7;   //!
   TBranch        *b_eleESE11;   //!
   TBranch        *b_eleESE21;   //!
   TBranch        *b_eleBrLinear;   //!
   TBranch        *b_eleCetaCorrE;   //!
   TBranch        *b_eleCetaCorrEt;   //!
   TBranch        *b_eleBremCorrE;   //!
   TBranch        *b_eleBremCorrEt;   //!
   TBranch        *b_eleFullCorrE;   //!
   TBranch        *b_eleFullCorrEt;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoCaloPos;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoVtx;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoCetaCorrE;   //!
   TBranch        *b_phoCetaCorrEt;   //!
   TBranch        *b_phoBremCorrE;   //!
   TBranch        *b_phoBremCorrEt;   //!
   TBranch        *b_phoFullCorrE;   //!
   TBranch        *b_phoFullCorrEt;   //!
   TBranch        *b_phoTrkIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoSolidDR04;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoTrkIsoHollowNoDzDR03;   //!
   TBranch        *b_phoTrkIsoHollowNoDzDR04;   //!
   TBranch        *b_phoCiCTrkIsoDR03;   //!
   TBranch        *b_phoCiCTrkIsoDR04;   //!
   TBranch        *b_phoCiCdRtoTrk;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoSolidDR04;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoEmax;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedDetId1;   //!
   TBranch        *b_phoSeedDetId2;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoNConv;   //!
   TBranch        *b_phoConvInvMass;   //!
   TBranch        *b_phoConvCotTheta;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoConvZofPVfromTrks;   //!
   TBranch        *b_phoConvMinDist;   //!
   TBranch        *b_phoConvdPhiAtVtx;   //!
   TBranch        *b_phoConvdPhiAtCalo;   //!
   TBranch        *b_phoConvdEtaAtCalo;   //!
   TBranch        *b_phoConvTrkd0;   //!
   TBranch        *b_phoConvTrkPin;   //!
   TBranch        *b_phoConvTrkPout;   //!
   TBranch        *b_phoConvTrkdz;   //!
   TBranch        *b_phoConvTrkdzErr;   //!
   TBranch        *b_phoConvChi2;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvCharge;   //!
   TBranch        *b_phoConvValidVtx;   //!
   TBranch        *b_phoConvLikeLihood;   //!
   TBranch        *b_phoConvP4;   //!
   TBranch        *b_phoConvVtx;   //!
   TBranch        *b_phoConvVtxErr;   //!
   TBranch        *b_phoConvPairMomentum;   //!
   TBranch        *b_phoConvRefittedMomentum;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoESE1;   //!
   TBranch        *b_phoESE3;   //!
   TBranch        *b_phoESE5;   //!
   TBranch        *b_phoESE7;   //!
   TBranch        *b_phoESE11;   //!
   TBranch        *b_phoESE21;   //!
/*    TBranch        *b_nMu;   //! */
/*    TBranch        *b_muTrg;   //! */
/*    TBranch        *b_muEta;   //! */
/*    TBranch        *b_muPhi;   //! */
/*    TBranch        *b_muCharge;   //! */
/*    TBranch        *b_muPt;   //! */
/*    TBranch        *b_muPz;   //! */
/*    TBranch        *b_muIsoTrk;   //! */
/*    TBranch        *b_muIsoCalo;   //! */
/*    TBranch        *b_muIsoEcal;   //! */
/*    TBranch        *b_muIsoHcal;   //! */
/*    TBranch        *b_muChi2NDF;   //! */
/*    TBranch        *b_muEmVeto;   //! */
/*    TBranch        *b_muHadVeto;   //! */
/*    TBranch        *b_muType;   //! */
/*    TBranch        *b_muID;   //! */
/*    TBranch        *b_muD0;   //! */
/*    TBranch        *b_muDz;   //! */
/*    TBranch        *b_muNumberOfValidTrkHits;   //! */
/*    TBranch        *b_muNumberOfValidPixelHits;   //! */
/*    TBranch        *b_muNumberOfValidMuonHits;   //! */
/*    TBranch        *b_muStations;   //! */
/*    TBranch        *b_muChambers;   //! */
/*    TBranch        *b_nPFPho_;   //! */
/*    TBranch        *b_PFPhoE_;   //! */
/*    TBranch        *b_PFPhoEt_;   //! */
/*    TBranch        *b_PFPhoEta_;   //! */
/*    TBranch        *b_PFPhoPhi_;   //! */
/*    TBranch        *b_PFPhoIsConv_;   //! */
/*    TBranch        *b_PFPhoTrkIsoHollowDR04_;   //! */
/*    TBranch        *b_PFPhoEcalIsoDR04_;   //! */
/*    TBranch        *b_PFPhoHcalIsoDR04_;   //! */
/*    TBranch        *b_PFPhoHoverE_;   //! */
/*    TBranch        *b_nPFEle_;   //! */
/*    TBranch        *b_PFElePt_;   //! */
/*    TBranch        *b_PFEleEta_;   //! */
/*    TBranch        *b_PFElePhi_;   //! */
/*    TBranch        *b_PFEleEn_;   //! */
/*    TBranch        *b_PFEleCharge;   //! */
/*    TBranch        *b_nPFMu_;   //! */
/*    TBranch        *b_PFmuCharge_;   //! */
/*    TBranch        *b_PFmuPhi_;   //! */
/*    TBranch        *b_PFmuEta_;   //! */
/*    TBranch        *b_PFmuPt_;   //! */
/*    TBranch        *b_PFmuPz_;   //! */
/*    TBranch        *b_PFmuChambers_;   //! */
/*    TBranch        *b_PFmuD0_;   //! */
/*    TBranch        *b_PFmuDz_;   //! */
/*    TBranch        *b_PFmuChi2NDF_;   //! */
/*    TBranch        *b_PFmuNumberOfValidTrkHits_;   //! */
/*    TBranch        *b_PFmuNumberOfValidPixelHits_;   //! */
/*    TBranch        *b_PFmuNumberOfValidMuonHits_;   //! */
   TBranch        *b_rho25;   //!
   TBranch        *b_rho44;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetHFHAE;   //!
   TBranch        *b_jetHFEME;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetTrackCountHiEffBJetTags;   //!
   TBranch        *b_jetTrackCountHiPurBJetTags;   //!
   TBranch        *b_jetSimpleSVHiEffBJetTags;   //!
   TBranch        *b_jetSimpleSVHiPurBJetTags;   //!
   TBranch        *b_nZee;   //!
   TBranch        *b_ZeeMass;   //!
   TBranch        *b_ZeePt;   //!
   TBranch        *b_ZeeEta;   //!
   TBranch        *b_ZeePhi;   //!
   TBranch        *b_ZeeLeg1Index;   //!
   TBranch        *b_ZeeLeg2Index;   //!
   TBranch        *b_nZmumu;   //!
   TBranch        *b_ZmumuMass;   //!
   TBranch        *b_ZmumuPt;   //!
   TBranch        *b_ZmumuEta;   //!
   TBranch        *b_ZmumuPhi;   //!
   TBranch        *b_ZmumuLeg1Index;   //!
   TBranch        *b_ZmumuLeg2Index;   //!
/*    TBranch        *b_nWenu;   //! */
/*    TBranch        *b_WenuMassTPfMET;   //! */
/*    TBranch        *b_WenuEtPfMET;   //! */
/*    TBranch        *b_WenuACopPfMET;   //! */
/*    TBranch        *b_WenuEleIndex;   //! */
/*    TBranch        *b_nWmunu;   //! */
/*    TBranch        *b_WmunuMassTPfMET;   //! */
/*    TBranch        *b_WmunuEtPfMET;   //! */
/*    TBranch        *b_WmunuACopPfMET;   //! */
/*    TBranch        *b_WmunuMuIndex;   //! */
/*    TBranch        *b_nConv;   //! */
/*    TBranch        *b_convP4;   //! */
/*    TBranch        *b_convVtx;   //! */
/*    TBranch        *b_convVtxErr;   //! */
/*    TBranch        *b_convPairMomentum;   //! */
/*    TBranch        *b_convRefittedMomentum;   //! */
/*    TBranch        *b_convNTracks;   //! */
/*    TBranch        *b_convPairInvMass;   //! */
/*    TBranch        *b_convPairCotThetaSep;   //! */
/*    TBranch        *b_convEoverP;   //! */
/*    TBranch        *b_convDistOfMinApproach;   //! */
/*    TBranch        *b_convDPhiTrksAtVtx;   //! */
/*    TBranch        *b_convDPhiTrksAtEcal;   //! */
/*    TBranch        *b_convDEtaTrksAtEcal;   //! */
/*    TBranch        *b_convDxy;   //! */
/*    TBranch        *b_convDz;   //! */
/*    TBranch        *b_convLxy;   //! */
/*    TBranch        *b_convLz;   //! */
/*    TBranch        *b_convZofPrimVtxFromTrks;   //! */
/*    TBranch        *b_convNHitsBeforeVtx;   //! */
/*    TBranch        *b_convNSharedHits;   //! */
/*    TBranch        *b_convValidVtx;   //! */
/*    TBranch        *b_convMVALikelihood;   //! */
/*    TBranch        *b_convChi2;   //! */
/*    TBranch        *b_convChi2Probability;   //! */
/*    TBranch        *b_convTk1Dz;   //! */
/*    TBranch        *b_convTk2Dz;   //! */
/*    TBranch        *b_convTk1DzErr;   //! */
/*    TBranch        *b_convTk2DzErr;   //! */
/*    TBranch        *b_convCh1Ch2;   //! */
/*    TBranch        *b_convTk1D0;   //! */
/*    TBranch        *b_convTk1Pout;   //! */
/*    TBranch        *b_convTk1Pin;   //! */
/*    TBranch        *b_convTk2D0;   //! */
/*    TBranch        *b_convTk2Pout;   //! */
/*    TBranch        *b_convTk2Pin;   //! */

   Ana_gj(TTree *tree=0);
   virtual ~Ana_gj();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1); 
   virtual Int_t    PhotonID(int ipho);
   virtual Int_t    PhotonID_3Iso(int ipho);
   virtual Int_t    eID(int iele);
   virtual Int_t    jetID(int ijet);   
   virtual Int_t    phocat(int i);

};

#endif

#ifdef Ana_gj_cxx
Ana_gj::Ana_gj(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*    if (tree == 0) { */
/*       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/cmkuo/job_photon_2011a_Jul5rereco_diphoton_Jul6_JSON_specialskim_4283.root"); */
/*       if (!f) { */
/*          f = new TFile("/home/cmkuo/job_photon_2011a_Jul5rereco_diphoton_Jul6_JSON_specialskim_4283.root"); */
/*          f->cd("/home/cmkuo/job_photon_2011a_Jul5rereco_diphoton_Jul6_JSON_specialskim_4283.root:/ggNtuplizer"); */
/*       } */
/*       tree = (TTree*)gDirectory->Get("EventTree"); */

/*    } */
/*    Init(tree); */

  
  TChain *tc = new TChain("ggNtuplizer/EventTree","ana");
  tc->Add("/data08a/cmkuo/job_photon_2011a_Aug5rereco_photon_Aug5rerecoJSONv2_specialskim_4284p1.root");
  tc->Add("/data08a/cmkuo/job_photon_2011a_Jul5rereco_diphoton_Aug25_JSON_specialskim_4284p1.root");
  tc->Add("/data08a/cmkuo/job_photon_2011a_PRv6_photon_Aug26_JSON_specialskim_4284p1.root");
  Init(tc);



}

Ana_gj::~Ana_gj()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Ana_gj::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Ana_gj::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Ana_gj::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vtxbsTkIndex = 0;
   vtxbsTkWeight = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT", HLT, &b_HLT);
   fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
   fChain->SetBranchAddress("bspotPos", bspotPos, &b_bspotPos);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("vtxNTrk", vtxNTrk, &b_vtxNTrk);
   fChain->SetBranchAddress("vtxNDF", vtxNDF, &b_vtxNDF);
   fChain->SetBranchAddress("vtxD0", vtxD0, &b_vtxD0);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nVtxBS", &nVtxBS, &b_nVtxBS);
   fChain->SetBranchAddress("vtxbs", vtxbs, &b_vtxbs);
   fChain->SetBranchAddress("vtxbsPtMod", vtxbsPtMod, &b_vtxbsPtMod);
   fChain->SetBranchAddress("vtxbsSumPt2", vtxbsSumPt2, &b_vtxbsSumPt2);
   fChain->SetBranchAddress("vtxbsTkIndex", &vtxbsTkIndex, &b_vtxbsTkIndex);
   fChain->SetBranchAddress("vtxbsTkWeight", &vtxbsTkWeight, &b_vtxbsTkWeight);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("nGoodTrk", &nGoodTrk, &b_nGoodTrk);
   fChain->SetBranchAddress("IsTracksGood", &IsTracksGood, &b_IsTracksGood);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("METx", &METx, &b_METx);
   fChain->SetBranchAddress("METy", &METy, &b_METy);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METsumEt", &METsumEt, &b_METsumEt);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
   fChain->SetBranchAddress("tcMETx", &tcMETx, &b_tcMETx);
   fChain->SetBranchAddress("tcMETy", &tcMETy, &b_tcMETy);
   fChain->SetBranchAddress("tcMETPhi", &tcMETPhi, &b_tcMETPhi);
   fChain->SetBranchAddress("tcMETsumEt", &tcMETsumEt, &b_tcMETsumEt);
   fChain->SetBranchAddress("tcMETmEtSig", &tcMETmEtSig, &b_tcMETmEtSig);
   fChain->SetBranchAddress("tcMETSig", &tcMETSig, &b_tcMETSig);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETx", &pfMETx, &b_pfMETx);
   fChain->SetBranchAddress("pfMETy", &pfMETy, &b_pfMETy);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleID", eleID, &b_eleID);
   fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx", eleVtx, &b_eleVtx);
   fChain->SetBranchAddress("eleD0", eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", elePout, &b_elePout);
   fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleEmax", eleEmax, &b_eleEmax);
   fChain->SetBranchAddress("eleE3x3", eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleE2x5Right", eleE2x5Right, &b_eleE2x5Right);
   fChain->SetBranchAddress("eleE2x5Left", eleE2x5Left, &b_eleE2x5Left);
   fChain->SetBranchAddress("eleE2x5Top", eleE2x5Top, &b_eleE2x5Top);
   fChain->SetBranchAddress("eleE2x5Bottom", eleE2x5Bottom, &b_eleE2x5Bottom);
   fChain->SetBranchAddress("eleSeedTime", eleSeedTime, &b_eleSeedTime);
   fChain->SetBranchAddress("eleRecoFlag", eleRecoFlag, &b_eleRecoFlag);
   fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalSolidDR03", eleIsoHcalSolidDR03, &b_eleIsoHcalSolidDR03);
   fChain->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04, &b_eleIsoTrkDR04);
   fChain->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04, &b_eleIsoEcalDR04);
   fChain->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04, &b_eleIsoHcalDR04);
   fChain->SetBranchAddress("eleIsoHcalSolidDR04", eleIsoHcalSolidDR04, &b_eleIsoHcalSolidDR04);
   fChain->SetBranchAddress("eleMissHits", eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleConvDist", eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("eleESEffSigmaRR", eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("eleESE1", eleESE1, &b_eleESE1);
   fChain->SetBranchAddress("eleESE3", eleESE3, &b_eleESE3);
   fChain->SetBranchAddress("eleESE5", eleESE5, &b_eleESE5);
   fChain->SetBranchAddress("eleESE7", eleESE7, &b_eleESE7);
   fChain->SetBranchAddress("eleESE11", eleESE11, &b_eleESE11);
   fChain->SetBranchAddress("eleESE21", eleESE21, &b_eleESE21);
   fChain->SetBranchAddress("eleBrLinear", eleBrLinear, &b_eleBrLinear);
   fChain->SetBranchAddress("eleCetaCorrE", eleCetaCorrE, &b_eleCetaCorrE);
   fChain->SetBranchAddress("eleCetaCorrEt", eleCetaCorrEt, &b_eleCetaCorrEt);
   fChain->SetBranchAddress("eleBremCorrE", eleBremCorrE, &b_eleBremCorrE);
   fChain->SetBranchAddress("eleBremCorrEt", eleBremCorrEt, &b_eleBremCorrEt);
   fChain->SetBranchAddress("eleFullCorrE", eleFullCorrE, &b_eleFullCorrE);
   fChain->SetBranchAddress("eleFullCorrEt", eleFullCorrEt, &b_eleFullCorrEt);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoIsPhoton", phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoCaloPos", phoCaloPos, &b_phoCaloPos);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoVtx", phoVtx, &b_phoVtx);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEtVtx", phoEtVtx, &b_phoEtVtx);
   fChain->SetBranchAddress("phoEtaVtx", phoEtaVtx, &b_phoEtaVtx);
   fChain->SetBranchAddress("phoPhiVtx", phoPhiVtx, &b_phoPhiVtx);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoCetaCorrE", phoCetaCorrE, &b_phoCetaCorrE);
   fChain->SetBranchAddress("phoCetaCorrEt", phoCetaCorrEt, &b_phoCetaCorrEt);
   fChain->SetBranchAddress("phoBremCorrE", phoBremCorrE, &b_phoBremCorrE);
   fChain->SetBranchAddress("phoBremCorrEt", phoBremCorrEt, &b_phoBremCorrEt);
   fChain->SetBranchAddress("phoFullCorrE", phoFullCorrE, &b_phoFullCorrE);
   fChain->SetBranchAddress("phoFullCorrEt", phoFullCorrEt, &b_phoFullCorrEt);
   fChain->SetBranchAddress("phoTrkIsoSolidDR03", phoTrkIsoSolidDR03, &b_phoTrkIsoSolidDR03);
   fChain->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
   fChain->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03, &b_phoEcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03, &b_phoHcalIsoDR03);
   fChain->SetBranchAddress("phoHcalIsoSolidDR03", phoHcalIsoSolidDR03, &b_phoHcalIsoSolidDR03);
   fChain->SetBranchAddress("phoTrkIsoSolidDR04", phoTrkIsoSolidDR04, &b_phoTrkIsoSolidDR04);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoTrkIsoHollowNoDzDR03", phoTrkIsoHollowNoDzDR03, &b_phoTrkIsoHollowNoDzDR03);
   fChain->SetBranchAddress("phoTrkIsoHollowNoDzDR04", phoTrkIsoHollowNoDzDR04, &b_phoTrkIsoHollowNoDzDR04);
   fChain->SetBranchAddress("phoCiCTrkIsoDR03", phoCiCTrkIsoDR03, &b_phoCiCTrkIsoDR03);
   fChain->SetBranchAddress("phoCiCTrkIsoDR04", phoCiCTrkIsoDR04, &b_phoCiCTrkIsoDR04);
   fChain->SetBranchAddress("phoCiCdRtoTrk", phoCiCdRtoTrk, &b_phoCiCdRtoTrk);
   fChain->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoSolidDR04", phoHcalIsoSolidDR04, &b_phoHcalIsoSolidDR04);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoEmax", phoEmax, &b_phoEmax);
   fChain->SetBranchAddress("phoE3x3", phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE5x5", phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoE2x5Right", phoE2x5Right, &b_phoE2x5Right);
   fChain->SetBranchAddress("phoE2x5Left", phoE2x5Left, &b_phoE2x5Left);
   fChain->SetBranchAddress("phoE2x5Top", phoE2x5Top, &b_phoE2x5Top);
   fChain->SetBranchAddress("phoE2x5Bottom", phoE2x5Bottom, &b_phoE2x5Bottom);
   fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedDetId1", phoSeedDetId1, &b_phoSeedDetId1);
   fChain->SetBranchAddress("phoSeedDetId2", phoSeedDetId2, &b_phoSeedDetId2);
   fChain->SetBranchAddress("phoRecoFlag", phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoPos", phoPos, &b_phoPos);
   fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEn", phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoSCEt", phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoOverlap", phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
   fChain->SetBranchAddress("phoNConv", phoNConv, &b_phoNConv);
   fChain->SetBranchAddress("phoConvInvMass", phoConvInvMass, &b_phoConvInvMass);
   fChain->SetBranchAddress("phoConvCotTheta", phoConvCotTheta, &b_phoConvCotTheta);
   fChain->SetBranchAddress("phoConvEoverP", phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoConvZofPVfromTrks", phoConvZofPVfromTrks, &b_phoConvZofPVfromTrks);
   fChain->SetBranchAddress("phoConvMinDist", phoConvMinDist, &b_phoConvMinDist);
   fChain->SetBranchAddress("phoConvdPhiAtVtx", phoConvdPhiAtVtx, &b_phoConvdPhiAtVtx);
   fChain->SetBranchAddress("phoConvdPhiAtCalo", phoConvdPhiAtCalo, &b_phoConvdPhiAtCalo);
   fChain->SetBranchAddress("phoConvdEtaAtCalo", phoConvdEtaAtCalo, &b_phoConvdEtaAtCalo);
   fChain->SetBranchAddress("phoConvTrkd0", phoConvTrkd0, &b_phoConvTrkd0);
   fChain->SetBranchAddress("phoConvTrkPin", phoConvTrkPin, &b_phoConvTrkPin);
   fChain->SetBranchAddress("phoConvTrkPout", phoConvTrkPout, &b_phoConvTrkPout);
   fChain->SetBranchAddress("phoConvTrkdz", phoConvTrkdz, &b_phoConvTrkdz);
   fChain->SetBranchAddress("phoConvTrkdzErr", phoConvTrkdzErr, &b_phoConvTrkdzErr);
   fChain->SetBranchAddress("phoConvChi2", phoConvChi2, &b_phoConvChi2);
   fChain->SetBranchAddress("phoConvChi2Prob", phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvCharge", phoConvCharge, &b_phoConvCharge);
   fChain->SetBranchAddress("phoConvValidVtx", phoConvValidVtx, &b_phoConvValidVtx);
   fChain->SetBranchAddress("phoConvLikeLihood", phoConvLikeLihood, &b_phoConvLikeLihood);
   fChain->SetBranchAddress("phoConvP4", phoConvP4, &b_phoConvP4);
   fChain->SetBranchAddress("phoConvVtx", phoConvVtx, &b_phoConvVtx);
   fChain->SetBranchAddress("phoConvVtxErr", phoConvVtxErr, &b_phoConvVtxErr);
   fChain->SetBranchAddress("phoConvPairMomentum", phoConvPairMomentum, &b_phoConvPairMomentum);
   fChain->SetBranchAddress("phoConvRefittedMomentum", phoConvRefittedMomentum, &b_phoConvRefittedMomentum);
   fChain->SetBranchAddress("phoESEffSigmaRR", phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoESE1", phoESE1, &b_phoESE1);
   fChain->SetBranchAddress("phoESE3", phoESE3, &b_phoESE3);
   fChain->SetBranchAddress("phoESE5", phoESE5, &b_phoESE5);
   fChain->SetBranchAddress("phoESE7", phoESE7, &b_phoESE7);
   fChain->SetBranchAddress("phoESE11", phoESE11, &b_phoESE11);
   fChain->SetBranchAddress("phoESE21", phoESE21, &b_phoESE21);
/*    fChain->SetBranchAddress("nMu", &nMu, &b_nMu); */
/*    fChain->SetBranchAddress("muTrg", muTrg, &b_muTrg); */
/*    fChain->SetBranchAddress("muEta", muEta, &b_muEta); */
/*    fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi); */
/*    fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge); */
/*    fChain->SetBranchAddress("muPt", muPt, &b_muPt); */
/*    fChain->SetBranchAddress("muPz", muPz, &b_muPz); */
/*    fChain->SetBranchAddress("muIsoTrk", muIsoTrk, &b_muIsoTrk); */
/*    fChain->SetBranchAddress("muIsoCalo", muIsoCalo, &b_muIsoCalo); */
/*    fChain->SetBranchAddress("muIsoEcal", muIsoEcal, &b_muIsoEcal); */
/*    fChain->SetBranchAddress("muIsoHcal", muIsoHcal, &b_muIsoHcal); */
/*    fChain->SetBranchAddress("muChi2NDF", muChi2NDF, &b_muChi2NDF); */
/*    fChain->SetBranchAddress("muEmVeto", muEmVeto, &b_muEmVeto); */
/*    fChain->SetBranchAddress("muHadVeto", muHadVeto, &b_muHadVeto); */
/*    fChain->SetBranchAddress("muType", muType, &b_muType); */
/*    fChain->SetBranchAddress("muID", muID, &b_muID); */
/*    fChain->SetBranchAddress("muD0", muD0, &b_muD0); */
/*    fChain->SetBranchAddress("muDz", muDz, &b_muDz); */
/*    fChain->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits); */
/*    fChain->SetBranchAddress("muNumberOfValidPixelHits", muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits); */
/*    fChain->SetBranchAddress("muNumberOfValidMuonHits", muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits); */
/*    fChain->SetBranchAddress("muStations", muStations, &b_muStations); */
/*    fChain->SetBranchAddress("muChambers", muChambers, &b_muChambers); */
/*    fChain->SetBranchAddress("nPFPho_", &nPFPho_, &b_nPFPho_); */
/*    fChain->SetBranchAddress("PFPhoE_", PFPhoE_, &b_PFPhoE_); */
/*    fChain->SetBranchAddress("PFPhoEt_", PFPhoEt_, &b_PFPhoEt_); */
/*    fChain->SetBranchAddress("PFPhoEta_", PFPhoEta_, &b_PFPhoEta_); */
/*    fChain->SetBranchAddress("PFPhoPhi_", PFPhoPhi_, &b_PFPhoPhi_); */
/*    fChain->SetBranchAddress("PFPhoIsConv_", PFPhoIsConv_, &b_PFPhoIsConv_); */
/*    fChain->SetBranchAddress("PFPhoTrkIsoHollowDR04_", PFPhoTrkIsoHollowDR04_, &b_PFPhoTrkIsoHollowDR04_); */
/*    fChain->SetBranchAddress("PFPhoEcalIsoDR04_", PFPhoEcalIsoDR04_, &b_PFPhoEcalIsoDR04_); */
/*    fChain->SetBranchAddress("PFPhoHcalIsoDR04_", PFPhoHcalIsoDR04_, &b_PFPhoHcalIsoDR04_); */
/*    fChain->SetBranchAddress("PFPhoHoverE_", PFPhoHoverE_, &b_PFPhoHoverE_); */
/*    fChain->SetBranchAddress("nPFEle_", &nPFEle_, &b_nPFEle_); */
/*    fChain->SetBranchAddress("PFElePt_", PFElePt_, &b_PFElePt_); */
/*    fChain->SetBranchAddress("PFEleEta_", PFEleEta_, &b_PFEleEta_); */
/*    fChain->SetBranchAddress("PFElePhi_", PFElePhi_, &b_PFElePhi_); */
/*    fChain->SetBranchAddress("PFEleEn_", PFEleEn_, &b_PFEleEn_); */
/*    fChain->SetBranchAddress("PFEleCharge", PFEleCharge, &b_PFEleCharge); */
/*    fChain->SetBranchAddress("nPFMu_", &nPFMu_, &b_nPFMu_); */
/*    fChain->SetBranchAddress("PFmuCharge_", PFmuCharge_, &b_PFmuCharge_); */
/*    fChain->SetBranchAddress("PFmuPhi_", PFmuPhi_, &b_PFmuPhi_); */
/*    fChain->SetBranchAddress("PFmuEta_", PFmuEta_, &b_PFmuEta_); */
/*    fChain->SetBranchAddress("PFmuPt_", PFmuPt_, &b_PFmuPt_); */
/*    fChain->SetBranchAddress("PFmuPz_", PFmuPz_, &b_PFmuPz_); */
/*    fChain->SetBranchAddress("PFmuChambers_", PFmuChambers_, &b_PFmuChambers_); */
/*    fChain->SetBranchAddress("PFmuD0_", PFmuD0_, &b_PFmuD0_); */
/*    fChain->SetBranchAddress("PFmuDz_", PFmuDz_, &b_PFmuDz_); */
/*    fChain->SetBranchAddress("PFmuChi2NDF_", PFmuChi2NDF_, &b_PFmuChi2NDF_); */
/*    fChain->SetBranchAddress("PFmuNumberOfValidTrkHits_", PFmuNumberOfValidTrkHits_, &b_PFmuNumberOfValidTrkHits_); */
/*    fChain->SetBranchAddress("PFmuNumberOfValidPixelHits_", PFmuNumberOfValidPixelHits_, &b_PFmuNumberOfValidPixelHits_); */
/*    fChain->SetBranchAddress("PFmuNumberOfValidMuonHits_", PFmuNumberOfValidMuonHits_, &b_PFmuNumberOfValidMuonHits_); */
   fChain->SetBranchAddress("rho25", &rho25, &b_rho25);
   fChain->SetBranchAddress("rho44", &rho44, &b_rho44);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   fChain->SetBranchAddress("jetRawPt", jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetCHF", jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetHFHAE", jetHFHAE, &b_jetHFHAE);
   fChain->SetBranchAddress("jetHFEME", jetHFEME, &b_jetHFEME);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetTrackCountHiEffBJetTags", jetTrackCountHiEffBJetTags, &b_jetTrackCountHiEffBJetTags);
   fChain->SetBranchAddress("jetTrackCountHiPurBJetTags", jetTrackCountHiPurBJetTags, &b_jetTrackCountHiPurBJetTags);
   fChain->SetBranchAddress("jetSimpleSVHiEffBJetTags", jetSimpleSVHiEffBJetTags, &b_jetSimpleSVHiEffBJetTags);
   fChain->SetBranchAddress("jetSimpleSVHiPurBJetTags", jetSimpleSVHiPurBJetTags, &b_jetSimpleSVHiPurBJetTags);
   fChain->SetBranchAddress("nZee", &nZee, &b_nZee);
   fChain->SetBranchAddress("ZeeMass", ZeeMass, &b_ZeeMass);
   fChain->SetBranchAddress("ZeePt", ZeePt, &b_ZeePt);
   fChain->SetBranchAddress("ZeeEta", ZeeEta, &b_ZeeEta);
   fChain->SetBranchAddress("ZeePhi", ZeePhi, &b_ZeePhi);
   fChain->SetBranchAddress("ZeeLeg1Index", ZeeLeg1Index, &b_ZeeLeg1Index);
   fChain->SetBranchAddress("ZeeLeg2Index", ZeeLeg2Index, &b_ZeeLeg2Index);
   fChain->SetBranchAddress("nZmumu", &nZmumu, &b_nZmumu);
   fChain->SetBranchAddress("ZmumuMass", ZmumuMass, &b_ZmumuMass);
   fChain->SetBranchAddress("ZmumuPt", ZmumuPt, &b_ZmumuPt);
   fChain->SetBranchAddress("ZmumuEta", ZmumuEta, &b_ZmumuEta);
   fChain->SetBranchAddress("ZmumuPhi", ZmumuPhi, &b_ZmumuPhi);
   fChain->SetBranchAddress("ZmumuLeg1Index", ZmumuLeg1Index, &b_ZmumuLeg1Index);
   fChain->SetBranchAddress("ZmumuLeg2Index", ZmumuLeg2Index, &b_ZmumuLeg2Index);
/*    fChain->SetBranchAddress("nWenu", &nWenu, &b_nWenu); */
/*    fChain->SetBranchAddress("WenuMassTPfMET", WenuMassTPfMET, &b_WenuMassTPfMET); */
/*    fChain->SetBranchAddress("WenuEtPfMET", WenuEtPfMET, &b_WenuEtPfMET); */
/*    fChain->SetBranchAddress("WenuACopPfMET", WenuACopPfMET, &b_WenuACopPfMET); */
/*    fChain->SetBranchAddress("WenuEleIndex", WenuEleIndex, &b_WenuEleIndex); */
/*    fChain->SetBranchAddress("nWmunu", &nWmunu, &b_nWmunu); */
/*    fChain->SetBranchAddress("WmunuMassTPfMET", WmunuMassTPfMET, &b_WmunuMassTPfMET); */
/*    fChain->SetBranchAddress("WmunuEtPfMET", WmunuEtPfMET, &b_WmunuEtPfMET); */
/*    fChain->SetBranchAddress("WmunuACopPfMET", WmunuACopPfMET, &b_WmunuACopPfMET); */
/*    fChain->SetBranchAddress("WmunuMuIndex", WmunuMuIndex, &b_WmunuMuIndex); */
/*    fChain->SetBranchAddress("nConv", &nConv, &b_nConv); */
/*    fChain->SetBranchAddress("convP4", convP4, &b_convP4); */
/*    fChain->SetBranchAddress("convVtx", convVtx, &b_convVtx); */
/*    fChain->SetBranchAddress("convVtxErr", convVtxErr, &b_convVtxErr); */
/*    fChain->SetBranchAddress("convPairMomentum", convPairMomentum, &b_convPairMomentum); */
/*    fChain->SetBranchAddress("convRefittedMomentum", convRefittedMomentum, &b_convRefittedMomentum); */
/*    fChain->SetBranchAddress("convNTracks", convNTracks, &b_convNTracks); */
/*    fChain->SetBranchAddress("convPairInvMass", convPairInvMass, &b_convPairInvMass); */
/*    fChain->SetBranchAddress("convPairCotThetaSep", convPairCotThetaSep, &b_convPairCotThetaSep); */
/*    fChain->SetBranchAddress("convEoverP", convEoverP, &b_convEoverP); */
/*    fChain->SetBranchAddress("convDistOfMinApproach", convDistOfMinApproach, &b_convDistOfMinApproach); */
/*    fChain->SetBranchAddress("convDPhiTrksAtVtx", convDPhiTrksAtVtx, &b_convDPhiTrksAtVtx); */
/*    fChain->SetBranchAddress("convDPhiTrksAtEcal", convDPhiTrksAtEcal, &b_convDPhiTrksAtEcal); */
/*    fChain->SetBranchAddress("convDEtaTrksAtEcal", convDEtaTrksAtEcal, &b_convDEtaTrksAtEcal); */
/*    fChain->SetBranchAddress("convDxy", convDxy, &b_convDxy); */
/*    fChain->SetBranchAddress("convDz", convDz, &b_convDz); */
/*    fChain->SetBranchAddress("convLxy", convLxy, &b_convLxy); */
/*    fChain->SetBranchAddress("convLz", convLz, &b_convLz); */
/*    fChain->SetBranchAddress("convZofPrimVtxFromTrks", convZofPrimVtxFromTrks, &b_convZofPrimVtxFromTrks); */
/*    fChain->SetBranchAddress("convNHitsBeforeVtx", convNHitsBeforeVtx, &b_convNHitsBeforeVtx); */
/*    fChain->SetBranchAddress("convNSharedHits", convNSharedHits, &b_convNSharedHits); */
/*    fChain->SetBranchAddress("convValidVtx", convValidVtx, &b_convValidVtx); */
/*    fChain->SetBranchAddress("convMVALikelihood", convMVALikelihood, &b_convMVALikelihood); */
/*    fChain->SetBranchAddress("convChi2", convChi2, &b_convChi2); */
/*    fChain->SetBranchAddress("convChi2Probability", convChi2Probability, &b_convChi2Probability); */
/*    fChain->SetBranchAddress("convTk1Dz", convTk1Dz, &b_convTk1Dz); */
/*    fChain->SetBranchAddress("convTk2Dz", convTk2Dz, &b_convTk2Dz); */
/*    fChain->SetBranchAddress("convTk1DzErr", convTk1DzErr, &b_convTk1DzErr); */
/*    fChain->SetBranchAddress("convTk2DzErr", convTk2DzErr, &b_convTk2DzErr); */
/*    fChain->SetBranchAddress("convCh1Ch2", convCh1Ch2, &b_convCh1Ch2); */
/*    fChain->SetBranchAddress("convTk1D0", convTk1D0, &b_convTk1D0); */
/*    fChain->SetBranchAddress("convTk1Pout", convTk1Pout, &b_convTk1Pout); */
/*    fChain->SetBranchAddress("convTk1Pin", convTk1Pin, &b_convTk1Pin); */
/*    fChain->SetBranchAddress("convTk2D0", convTk2D0, &b_convTk2D0); */
/*    fChain->SetBranchAddress("convTk2Pout", convTk2Pout, &b_convTk2Pout); */
/*    fChain->SetBranchAddress("convTk2Pin", convTk2Pin, &b_convTk2Pin); */
   Notify();
}

Bool_t Ana_gj::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Ana_gj::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Ana_gj::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Ana_gj_cxx
