#include <TMath.h>
#include <TMVA/Reader.h>

using namespace std;

Int_t HggPreselection(TreeReader &data, Int_t ipho, Bool_t eleVeto=kTRUE) {

  Int_t phoID=1;

  // load relevant branches from TTree/TChain
  Int_t    nPho                = data.GetInt("nPho");
  Float_t* phoEt               = data.GetPtrFloat("phoEt");
  Float_t* phoSCEta            = data.GetPtrFloat("phoSCEta");
  Int_t*   phoEleVeto          = data.GetPtrInt("phoEleVeto");
  Float_t* phoSigmaIEtaIEta    = data.GetPtrFloat("phoSigmaIEtaIEta_2012");
  Float_t* phoR9               = data.GetPtrFloat("phoR9");
  Float_t* phoHoverE           = data.GetPtrFloat("phoHoverE");
  Float_t* phoHcalIsoDR03      = data.GetPtrFloat("phoHcalIsoDR03");
  Float_t* phoTrkIsoHollowDR03 = data.GetPtrFloat("phoTrkIsoHollowDR03");
  vector<float>* phoCiCPF4chgpfIso02 = data.GetPtrVectorFloat("phoCiCPF4chgpfIso02", nPho);
                              
  if (phoEt[ipho] < 15.) phoID = 0;
  if (phoEt[ipho] > 200.) phoID = 0;
  if (TMath::Abs(phoSCEta[ipho]) > 1.4442 && TMath::Abs(phoSCEta[ipho]) < 1.566) phoID = 0;
  if (TMath::Abs(phoSCEta[ipho]) > 2.5) phoID = 0;
  if (eleVeto && phoEleVeto[ipho] == 1) phoID = 0;  

  // Hgg photon Preselection-------------
  if (TMath::Abs(phoSCEta[ipho]) < 1.5) {
    if (phoSigmaIEtaIEta[ipho] > 0.014) phoID = 0;
    if (phoR9[ipho] > 0.9) {
      if (phoHoverE[ipho] >= 0.082) phoID = 0;
    } else {
      if (phoHoverE[ipho] >= 0.075) phoID = 0;
    }
  } else {
    if (phoSigmaIEtaIEta[ipho] > 0.034) phoID = 0;
    if (phoHoverE[ipho] >= 0.075) phoID = 0;
  }
  
  if (phoR9[ipho] > 0.9) {
    if (phoHcalIsoDR03[ipho]-0.005*phoEt[ipho] > 50.) phoID = 0;
    if (phoTrkIsoHollowDR03[ipho]-0.002*phoEt[ipho] > 50.) phoID = 0;
    if (phoCiCPF4chgpfIso02[ipho][0] >4.) phoID = 0;
  } else {
    if (phoHcalIsoDR03[ipho]-0.005*phoEt[ipho] > 4.) phoID = 0;
    if (phoTrkIsoHollowDR03[ipho]-0.002*phoEt[ipho] > 4.) phoID = 0;
    if (phoCiCPF4chgpfIso02[ipho][0] >4.) phoID = 0;
  }
  //end of Hgg preselection----------
  
  return phoID;
}

void select_photon(Int_t iWP, TreeReader &data, vector<int>& accepted) {

  accepted.clear();
  // load relevant branches from TTree/TChain
  Int_t    nPho                = data.GetInt("nPho");
  Float_t* phoEt               = data.GetPtrFloat("phoEt");
  Float_t* phoSCEta            = data.GetPtrFloat("phoSCEta");
  Int_t*   phoEleVeto          = data.GetPtrInt("phoEleVeto");
  Float_t* phoSigmaIEtaIEta    = data.GetPtrFloat("phoSigmaIEtaIEta_2012");
  Float_t* phoHoverE           = data.GetPtrFloat("phoHoverE");
  Float_t* phoPFChIso          = data.GetPtrFloat("phoPFChIso");
  Float_t* phoPFNeuIso         = data.GetPtrFloat("phoPFNeuIso");
  Float_t* phoPFPhoIso         = data.GetPtrFloat("phoPFPhoIso");
  Float_t  rho                = data.GetFloat("rho");

  for (int iPho = 0; iPho < nPho; ++iPho) {

    if (phoEt[iPho] < 10. ) continue ;
    
    Int_t phoEB = 0 ;
    if (iWP == -1 ) continue ;
    if ( fabs(phoSCEta[iPho]) >= 1.566 ) phoEB = 1 ;
    if ( phoEleVeto[iPho] != 0 ) continue ;
    if ( phoHoverE[iPho] > 0.05 ) continue ;
    Float_t sIEIECut[3][2] = { { 0.012 , 0.034 } , { 0.011 , 0.033 } , { 0.011 , 0.031 } } ;
    if ( phoSigmaIEtaIEta[iPho] >= sIEIECut[iWP][phoEB] ) continue ;
    Float_t effArea[3][7] = { //[Ch,Nu,Pho][iPhof_eta]
      { 0.012 , 0.010 , 0.014 , 0.012 , 0.016 , 0.020 , 0.012 } ,
      { 0.030 , 0.057 , 0.039 , 0.015 , 0.024 , 0.039 , 0.072 } ,
      { 0.148 , 0.130 , 0.112 , 0.216 , 0.262 , 0.260 , 0.266 } 
    } ;
    Int_t i_effArea = 0 ; // effective area for pile up correction for DR04 combine rel. Iso
    if      ( fabs(phoSCEta[iPho]) < 1.0                                        ) i_effArea = 0 ;
    else if ( fabs(phoSCEta[iPho]) >= 1.0   && fabs(phoSCEta[iPho]) < 1.479  ) i_effArea = 1 ;
    else if ( fabs(phoSCEta[iPho]) >= 1.479 && fabs(phoSCEta[iPho]) < 2.0    ) i_effArea = 2 ;
    else if ( fabs(phoSCEta[iPho]) >= 2.0   && fabs(phoSCEta[iPho]) < 2.2    ) i_effArea = 3 ;
    else if ( fabs(phoSCEta[iPho]) >= 2.2   && fabs(phoSCEta[iPho]) < 2.3    ) i_effArea = 4 ;
    else if ( fabs(phoSCEta[iPho]) >= 2.3   && fabs(phoSCEta[iPho]) < 2.4    ) i_effArea = 5 ;
    else if ( fabs(phoSCEta[iPho]) >= 2.4                                       ) i_effArea = 6 ;
    Float_t chIsoCut[3][2] = { {2.6,2.3} , {1.5,1.2} , {0.7,0.5} } ; //[Loose,Medium,Tight][EB,EE]
    float corrIso = (float) TMath::Max( float(0.) , phoPFChIso[iPho] - ( effArea[0][i_effArea] * rho ) ) ;
    if ( corrIso >= chIsoCut[iWP][phoEB] ) continue ;
    Float_t neuIsoCut[3][2] = { { 3.5, 2.9 } , { 1.0, 1.5 } , { 0.4, 1.5 } } ; //[Loose,Medium,Tight][EB,EE]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 2; j++) neuIsoCut[i][j] += ( 0.04 * phoEt[iPho] );
    corrIso = TMath::Max( float(0.) , phoPFNeuIso[iPho] - ( effArea[1][i_effArea] * rho ) ) ;
    if ( corrIso >= neuIsoCut[iWP][phoEB] ) continue ;
    Float_t gammaIsoCut[3][2] = { {1.3,999.} , {0.7,1.0} , {0.5,1.0} } ; //[Loose,Medium,Tight][EB,EE]
    for (int i = 0; i < 3; i++) for (int j = 0; j < 2; j++) gammaIsoCut[i][j] += ( 0.005 * phoEt[iPho] ) ;
    corrIso = TMath::Max( float(0.) , phoPFPhoIso[iPho] - ( effArea[2][i_effArea] * rho ) ) ;
    if ( corrIso >= gammaIsoCut[iWP][phoEB] ) continue ;
    accepted.push_back(iPho);
  }

}

float select_photon_mva(TreeReader &data, Int_t i) {

  /* Photon identification with the Zgamma MVA. Returns the MVA evaluated value.
   *
   * Documentation:
   * https://indico.cern.ch/getFile.py/access?contribId=3&resId=0&materialId=slides&confId=298231
   *
   * data = handle providing access to an input event;
   * i = index of a photon candidate to consider.
   */

  // load necessary tree branches
  float* phoEt             = data.GetPtrFloat("phoEt");
  float* phoEta            = data.GetPtrFloat("phoEta");
  float* phoPhi            = data.GetPtrFloat("phoPhi");
  float* phoR9             = data.GetPtrFloat("phoR9");
  float* phoSigmaIEtaIEta  = data.GetPtrFloat("phoSigmaIEtaIEta_2012");
  float* phoSigmaIEtaIPhi  = data.GetPtrFloat("phoSigmaIEtaIPhi_2012");
  float* phoE1x3           = data.GetPtrFloat("phoE1x3");
  float* phoE2x2           = data.GetPtrFloat("phoE2x2");
  float* phoE5x5           = data.GetPtrFloat("phoE5x5");
  float* phoE2x5Max        = data.GetPtrFloat("phoE2x5Max");
  float* phoSCEta          = data.GetPtrFloat("phoSCEta");
  float* phoSCRawE         = data.GetPtrFloat("phoSCRawE");
  float* phoSCEtaWidth     = data.GetPtrFloat("phoSCEtaWidth");
  float* phoSCPhiWidth     = data.GetPtrFloat("phoSCPhiWidth");
  float  rho               = data.GetFloat("rho");
  float* phoPFPhoIso       = data.GetPtrFloat("phoPFPhoIso");
  float* phoPFChIso        = data.GetPtrFloat("phoPFChIso");
  float* phoESEn           = data.GetPtrFloat("phoESEn");
  float* phoESEffSigmaRR   = data.GetPtrFloat("phoESEffSigmaRR");
  float* phoPFChWorstIso   = data.GetPtrFloat("phoPFChWorstIso");

  Float_t* phoSigmaIEtaIEta_2012  = data.GetPtrFloat("phoSigmaIEtaIEta_2012");
  Float_t* phoSigmaIEtaIPhi_2012  = data.GetPtrFloat("phoSigmaIEtaIPhi_2012");
  Float_t* phoSigmaIPhiIPhi_2012  = data.GetPtrFloat("phoSigmaIPhiIPhi_2012");  
  Float_t* phoE1x3_2012           = data.GetPtrFloat("phoE1x3_2012");
  Float_t* phoE2x2_2012            = data.GetPtrFloat("phoE2x2_2012");
  Float_t* phoE5x5_2012            = data.GetPtrFloat("phoE5x5_2012");
  Float_t* phoE2x5Max_2012        = data.GetPtrFloat("phoE2x5Max_2012");





  // classification variables
  static float phoEt_, phoEta_, phoPhi_, phoR9_;
  static float phoSigmaIEtaIEta_, phoSigmaIEtaIPhi_;
  static float phoS13_, phoS4_, phoS25_, phoSCEta_, phoSCRawE_;
  static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
  static float phoPFPhoIso_, phoPFChIso_, phoPFChIsoWorst_;
  static float phoESEnToRawE_, phoESEffSigmaRR_;

  static float sieie_2012, sieip_2012, s13_2012, s4_2012, s25_2012;

  // MVA classifiers for 0=ECAL barrel and 1=ECAL endcaps
  static TMVA::Reader* tmvaReader[2] = {NULL, NULL};

  // 0=ECAL barrel or 1=ECAL endcaps
  int iBE = (fabs(phoSCEta[i]) < 1.479) ? 0 : 1;

  // one-time MVA initialization
  if (!tmvaReader[iBE]) {
    tmvaReader[iBE] = new TMVA::Reader("!Color:Silent");

    // add classification variables
/*     tmvaReader[iBE]->AddVariable("phoPhi", &phoPhi_); */
/*     tmvaReader[iBE]->AddVariable("phoR9", &phoR9_); */
/*     tmvaReader[iBE]->AddVariable("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);     */
/*     tmvaReader[iBE]->AddVariable("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi_); */
/*     tmvaReader[iBE]->AddVariable("s13", &phoS13_); */
/*     tmvaReader[iBE]->AddVariable("s4ratio", &phoS4_); */
/*     tmvaReader[iBE]->AddVariable("s25", &phoS25_); */

/*     tmvaReader[iBE]->AddVariable("phoSCEta", &phoSCEta_); */
/*     tmvaReader[iBE]->AddVariable("phoSCRawE", &phoSCRawE_); */
/*     tmvaReader[iBE]->AddVariable("phoSCEtaWidth", &phoSCEtaWidth_); */
/*     tmvaReader[iBE]->AddVariable("phoSCPhiWidth", &phoSCPhiWidth_); */

/*     if (iBE == 1) { */
/*       tmvaReader[iBE]->AddVariable("phoESEn/phoSCRawE", &phoESEnToRawE_); */
/*       tmvaReader[iBE]->AddVariable("phoESEffSigmaRR", &phoESEffSigmaRR_); */
/*     } */

/*     tmvaReader[iBE]->AddVariable("rho_2012", &rho_); */
/*     tmvaReader[iBE]->AddVariable("phoPFPhoIso", &phoPFPhoIso_); */
/*     tmvaReader[iBE]->AddVariable("phoPFChIso", &phoPFChIso_); */
/*     tmvaReader[iBE]->AddVariable("phoPFChIsoWorst", &phoPFChIsoWorst_); */

/*     // FIXME: why do we need this? */
/*     tmvaReader[iBE]->AddSpectator("phoEt", &phoEt_); */
/*     tmvaReader[iBE]->AddSpectator("phoEta", &phoEta_); */


    tmvaReader[iBE]->AddVariable("recoPhi", &phoPhi_);
    tmvaReader[iBE]->AddVariable("r9", &phoR9_);
    tmvaReader[iBE]->AddVariable( "sieie_2012",       	      &sieie_2012 );     
    tmvaReader[iBE]->AddVariable( "sieip_2012",       	      &sieip_2012 );     
    tmvaReader[iBE]->AddVariable( "s13 := e1x3_2012/e5x5_2012",   &s13_2012 );	        
    tmvaReader[iBE]->AddVariable( "s4 := e2x2_2012/e5x5_2012",    &s4_2012 );	       
    tmvaReader[iBE]->AddVariable( "s25 := e2x5_2012/e5x5_2012",   &s25_2012 );	        
    tmvaReader[iBE]->AddVariable("recoSCEta", &phoSCEta_);
    tmvaReader[iBE]->AddVariable("rawE", &phoSCRawE_);
    tmvaReader[iBE]->AddVariable("scEtaWidth", &phoSCEtaWidth_);
    tmvaReader[iBE]->AddVariable("scPhiWidth", &phoSCPhiWidth_);
    if (iBE == 1) {
      tmvaReader[iBE]->AddVariable("ESEn := esEn/rawE", &phoESEnToRawE_);
      tmvaReader[iBE]->AddVariable("esRR", &phoESEffSigmaRR_);
    }
    tmvaReader[iBE]->AddVariable("rho", &rho_);
    tmvaReader[iBE]->AddVariable("phoIsoRaw", &phoPFPhoIso_);
    tmvaReader[iBE]->AddVariable("chIsoRaw", &phoPFChIso_);
    tmvaReader[iBE]->AddVariable("chWorstRaw", &phoPFChIsoWorst_);

    // FIXME: why do we need this?
    tmvaReader[iBE]->AddSpectator("recoPt", &phoEt_);
    tmvaReader[iBE]->AddSpectator("recoEta", &phoEta_);

    // read weight files
    if (iBE == 0){
      //tmvaReader[0]->BookMVA("BDT", "external/photonID_EB_BDT_20140202.weights.xml");
      tmvaReader[0]->BookMVA("BDT", "external/photonID_EB_BDT_phys14.weights.xml");
    }else{
      //tmvaReader[1]->BookMVA("BDT", "external/photonID_EE_BDT_20140202.weights.xml");
      tmvaReader[1]->BookMVA("BDT", "external/photonID_EE_BDT_phys14.weights.xml");
    }
  } // one-time initialization

  // set MVA variables
  phoPhi_ = phoPhi[i];
  phoR9_ = phoR9[i];
  phoSigmaIEtaIEta_ = phoSigmaIEtaIEta[i];
  phoSigmaIEtaIPhi_ = phoSigmaIEtaIPhi[i];
  phoS4_ = phoE2x2[i]/phoE5x5[i];
  phoS13_ = phoE1x3[i]/phoE5x5[i];
  phoS25_ = phoE2x5Max[i]/phoE5x5[i];
  phoSCEta_ = phoSCEta[i];
  phoSCRawE_ = phoSCRawE[i];
  phoSCEtaWidth_ = phoSCEtaWidth[i];
  phoSCPhiWidth_ = phoSCPhiWidth[i];
  rho_ = rho;
  phoPFPhoIso_ = phoPFPhoIso[i];
  phoPFChIso_ = phoPFChIso[i];
  phoESEnToRawE_ = phoESEn[i]/phoSCRawE[i];
  phoESEffSigmaRR_= phoESEffSigmaRR[i];
  phoEt_ = phoEt[i];
  phoEta_ = phoEta[i];
  phoPFChIsoWorst_ = phoPFChWorstIso[i];
  

  sieie_2012 = phoSigmaIEtaIEta_2012[i];
  sieip_2012 = phoSigmaIEtaIPhi_2012[i];

  s13_2012 = phoE1x3_2012[i]/phoE5x5_2012[i];
  s4_2012 = phoE2x2_2012[i]/phoE5x5_2012[i];
  s25_2012 = phoE2x5Max_2012[i]/phoE5x5_2012[i];

  return tmvaReader[iBE]->EvaluateMVA("BDT");

}
