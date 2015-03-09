/* Takes pre-minitrees from the xAnaZg-ng analysis and produces corresponding
 * minitrees, writing output into one and the same file output/minitrees.root.
 *
 * Actual configuration is given in the mkminitrees() steering function at the
 * end of this file.
 *
 * The directory "./output" is removed every time mkminitrees() is executed.
 * 
 * Usage: root -b -l -q mkminitrees.cc+
 */

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TLorentzVector.h>

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

// types of photon, electron and muon identification
enum {
   PHOIDLoose, PHOIDMedium, PHOIDTight, PHOIDMVA,
   ELEIDLoose, ELEIDMedium, ELEIDTight, ELEIDMVATrig,
   MUIDLoose, MUIDTight
};

//______________________________________________________________________________
void SetBranchAddress(TTree* tree, const char* bname, void* ptr)
{
   /* Associates tree branch with address of a variable.
    */

   // verify branch existence
   if (!tree->GetBranch(bname))
      FATAL(Form("tree branch \"%s\" does not exist", bname));

   // activate this branch
   tree->SetBranchStatus(bname, 1);

   // associate tree branch with address and check return value
   Int_t ret = tree->SetBranchAddress(bname, ptr);
   if (ret != 0 && ret != 4)
      FATAL("TTree::SetBranchAddress() returned bad code");
}

//______________________________________________________________________________
Int_t event_category(bool takeEle, float absEta1, float absEta2, float absEta, float phoR9)
{
   /* Returns event category number: any of 1, 2, 3 or 4.
    *
    * takeEle = true/electrons, false/muons;
    * absEta1, absEta2 = |eta| values for the first and the second lepton;
    * absEta = |eta| value for the photon;
    * phoR9 = photon's R9 value.
    */

   // electrons vs muons
   if (takeEle) {
      if (absEta1 < 1.479 && absEta2 < 1.479 && absEta < 1.479 && phoR9 >= 0.94)
         return 1;
      else if (absEta1 < 1.479 && absEta2 < 1.479 && absEta < 1.479 && phoR9 < 0.94)
         return 2;
      else if ((absEta1 >= 1.479 || absEta2 >= 1.479) && absEta < 1.479)
         return 3;
      else if (absEta >= 1.479)
         return 4;
   }
   else {
      if (absEta1 < 2.1 && absEta2 < 2.1 && (absEta1 < 0.9 || absEta2 < 0.9) && absEta < 1.479 && phoR9 >= 0.94)
         return 1;
      else if (absEta1 < 2.1 && absEta2 < 2.1 && (absEta1 < 0.9 || absEta2 < 0.9) && absEta < 1.479 && phoR9 < 0.94)
         return 2;
      else if (((absEta1 >= 0.9 && absEta2 >= 0.9) || (absEta1 >= 2.1 || absEta2 >= 2.1)) && absEta < 1.479)
         return 3;
      else if (absEta >= 1.479)
         return 4;
   }

   // never reached
   FATAL("inconsistency in event_category()");
   return -1;
}

//______________________________________________________________________________
void mkminitree(const char* inpath = "preminitree.root", const char* outname = "minitree",
                int lepIDtype = ELEIDMVATrig, int phoIDtype = PHOIDMVA,
                bool rmISRFSRPhotons = false)
{
   /* Processes given pre-minitree and produces minitree.
    *
    * Each event in the minitree contains:
    *   - run number and event number;
    *   - (overall) event weight (MC only);
    *   - event category;
    *   - 3-particle invariant mass;
    *   - pt, eta and phi of the two selected leptons and the photon.
    *
    * NOTE: MC events with pileup weight bigger than 20 are skipped
    * (corresponding warning is printed out).
    *
    * inpath = path to root file with pre-minitree to process;
    * outname = name of output minitree to be written into minitrees.root;
    * lepIDtype, phoIDtype = types of lepton and photon identification to take;
    *
    * rmISRFSRPhotons = (for MC) if true, ISR (initial state radiation) and FSR
    *     (final state radiation) photons will be removed from all events.
    */

   Printf("Processing %s: lepIDtype=%i, phoIDtype=%i, rmISRFSRPhotons=%i ...",
          inpath, lepIDtype, phoIDtype, rmISRFSRPhotons);

   // determine whether we are dealing with electrons or muons
   bool takeEle = (lepIDtype == MUIDLoose || lepIDtype == MUIDTight) ? false : true;

   // open root file, get input TTree
   TFile* fi = TFile::Open(inpath);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   TTree* intree = dynamic_cast<TTree*>(fi->Get("ggNtuplizer/EventTree"));
   if (!intree) FATAL("TFile::Get() failed");

   // determine whether MC branches are present or not
   bool hasMC = intree->GetBranch("phoGenMomPID") ? true : false;

   // variables to be associated with the input tree branches
   Int_t run;
   Long64_t event;
   float puwei;
   Int_t nLep, nPho, nJet;
   Int_t lepCharge[50];
   Int_t phoGenMomPID[50];
   float lepPt[50], lepEta[50], lepPhi[50], lepID[50], lepSF[50];
   float phoPt[50], phoEta[50], phoPhi[50], phoR9[50], phoID[50], phoSF[50];
   float jetPt[50], jetEta[50], jetPhi[50], jetEn[50];

   // disable all branches by default
   intree->SetBranchStatus("*", 0);

   // associate tree branches with variables
   SetBranchAddress(intree, "run", &run);
   SetBranchAddress(intree, "event", &event);

   // per-event weight due to pileup difference between data and MC
   if (hasMC)
      SetBranchAddress(intree, takeEle ? "puweiEle" : "puweiMu", &puwei);

   // electrons vs muons
   if (takeEle) {
      SetBranchAddress(intree, "nEle", &nLep);
      SetBranchAddress(intree, "elePt", lepPt);
      SetBranchAddress(intree, "eleEta", lepEta);
      SetBranchAddress(intree, "elePhi", lepPhi);
      SetBranchAddress(intree, "eleCharge", lepCharge);

      // ID
      if (lepIDtype == ELEIDLoose) {
         SetBranchAddress(intree, "eleIDLoose", lepID);
         if (hasMC) SetBranchAddress(intree, "eleSFLoose", lepSF);
      } else if (lepIDtype == ELEIDMedium) {
         SetBranchAddress(intree, "eleIDMedium", lepID);
         if (hasMC) SetBranchAddress(intree, "eleSFMedium", lepSF);
      } else if (lepIDtype == ELEIDTight) {
         SetBranchAddress(intree, "eleIDTight", lepID);
         if (hasMC) SetBranchAddress(intree, "eleSFTight", lepSF);
      } else if (lepIDtype == ELEIDMVATrig) {
         SetBranchAddress(intree, "eleIDMVATrig", lepID);
         if (hasMC) SetBranchAddress(intree, "eleSFMVATrig", lepSF);
      } else
         FATAL("wrong electron ID type provided");

   } else {
      SetBranchAddress(intree, "nMu", &nLep);
      SetBranchAddress(intree, "muPt", lepPt);
      SetBranchAddress(intree, "muEta", lepEta);
      SetBranchAddress(intree, "muPhi", lepPhi);
      SetBranchAddress(intree, "muCharge", lepCharge);

      // ID
      if (lepIDtype == MUIDLoose) {
         SetBranchAddress(intree, "muIDLoose", lepID);
         if (hasMC) SetBranchAddress(intree, "muSFLoose", lepSF);
      } else if (lepIDtype == MUIDTight) {
         SetBranchAddress(intree, "muIDTight", lepID);
         if (hasMC) SetBranchAddress(intree, "muSFTight", lepSF);
      } else
         FATAL("wrong muon ID type provided");
   }

   // photons
   SetBranchAddress(intree, "nPho", &nPho);
   SetBranchAddress(intree, "phoPt", phoPt);
   SetBranchAddress(intree, "phoEta", phoEta);
   SetBranchAddress(intree, "phoPhi", phoPhi);
   SetBranchAddress(intree, "phoR9", phoR9);

   if (rmISRFSRPhotons)
      SetBranchAddress(intree, "phoGenMomPID", phoGenMomPID);

   // photon ID
   if (phoIDtype == PHOIDLoose) {
      SetBranchAddress(intree, "phoIDLoose", phoID);
      if (hasMC) SetBranchAddress(intree, "phoSFLoose", phoSF);
   } else if (phoIDtype == PHOIDMedium) {
      SetBranchAddress(intree, "phoIDMedium", phoID);
      if (hasMC) SetBranchAddress(intree, "phoSFMedium", phoSF);
   } else if (phoIDtype == PHOIDTight) {
      SetBranchAddress(intree, "phoIDTight", phoID);
      if (hasMC) SetBranchAddress(intree, "phoSFTight", phoSF);
   } else if (phoIDtype == PHOIDMVA) {
      SetBranchAddress(intree, "phoIDMVA", phoID);
      if (hasMC) SetBranchAddress(intree, "phoSFMVA", phoSF);
   } else
      FATAL("wrong photon ID type provided");

   // jets
   SetBranchAddress(intree, "nJet", &nJet);
   SetBranchAddress(intree, "jetPt", jetPt);
   SetBranchAddress(intree, "jetEta", jetEta);
   SetBranchAddress(intree, "jetPhi", jetPhi);
   SetBranchAddress(intree, "jetEn", jetEn);

   // prepare output tree
   TFile* fo = TFile::Open("output/minitrees.root", "UPDATE");
   if (!fo || fo->IsZombie())
      FATAL("TFile::Open() failed");

   TTree* outtree = new TTree(outname, "Selected Zg events");

   // variables to be associated with the output tree branches
   Int_t category;
   float mcwei, hzg_mass;
   float lep1Pt, lep1Eta, lep1Phi;
   float lep2Pt, lep2Eta, lep2Phi;
   float phoPt_, phoEta_, phoPhi_;

   // associate variables with the output tree branches
   outtree->Branch("run", &run);
   outtree->Branch("event", &event);

   if (hasMC)
      outtree->Branch("mcwei", &mcwei);

   outtree->Branch("category", &category);
   outtree->Branch("hzg_mass", &hzg_mass);
   outtree->Branch("lep1Pt",  &lep1Pt);
   outtree->Branch("lep1Eta", &lep1Eta);
   outtree->Branch("lep1Phi", &lep1Phi);
   outtree->Branch("lep2Pt",  &lep2Pt);
   outtree->Branch("lep2Eta", &lep2Eta);
   outtree->Branch("lep2Phi", &lep2Phi);
   outtree->Branch("phoPt",   &phoPt_);
   outtree->Branch("phoEta",  &phoEta_);
   outtree->Branch("phoPhi",  &phoPhi_);

   // electron mass vs muon mass
   double lepMass = takeEle ? 0.000511 : 0.105658;

   // loop over events
   for (Long64_t ev = 0; ev < intree->GetEntriesFast(); ev++) {
      if (intree->GetEntry(ev) <= 0)
         FATAL("TTree::GetEntry() failed");

      TLorentzVector lep1, lep2, pho, jet1, jet2;

      int i1 = -1, i2 = -1;  // indices of best pair of leptons
      double dMZ = 1e+10;    // |two-lepton mass - true Z mass|

      // select best pair of leptons with invariant mass closest to the Z mass
      for (int i = 0; i < nLep; i++) {
         if (lepID[i] < 0.5) continue;
         if (lepPt[i] < 10) continue;

         for (int j = i + 1; j < nLep; j++) {
            if (lepID[j] < 0.5) continue;
            if (lepPt[j] < 10) continue;

            if (lepCharge[i] == lepCharge[j]) continue;
            if (lepPt[i] < 20 && lepPt[j] < 20) continue;

            lep1.SetPtEtaPhiM(lepPt[i], lepEta[i], lepPhi[i], lepMass);
            lep2.SetPtEtaPhiM(lepPt[j], lepEta[j], lepPhi[j], lepMass);

            // NOTE: loose muon ID should be complemented with DeltaR > 0.02 cut
            // in order to suppress contribution from split tracks, see
            // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId?rev=49#Loose_Muon
            if (lepIDtype == MUIDLoose && lep1.DeltaR(lep2) < 0.02)
               continue;

            TLorentzVector Z = lep1 + lep2;

            double delta = fabs(Z.M() - 91.188);
            if (delta < dMZ) {
               dMZ = delta;
               i1 = i;
               i2 = j;
            }
         }
      } // pair of leptons

      // skip event if no pair of leptons is accepted
      if (i1 < 0 || i2 < 0) continue;

      lep1.SetPtEtaPhiM(lepPt[i1], lepEta[i1], lepPhi[i1], lepMass);
      lep2.SetPtEtaPhiM(lepPt[i2], lepEta[i2], lepPhi[i2], lepMass);

      TLorentzVector Z = lep1 + lep2;
      if (Z.M() < 50) continue;

      // best photon selection;
      // NOTE: photons are already ordered in descending order of their pt
      int k;
      for (k = 0; k < nPho; k++) {
         if (phoPt[k] < 15) continue;

         // remove all ISR and FSR photons, if requested
         if (rmISRFSRPhotons && fabs(phoGenMomPID[k]) <= 22) continue;

         category = event_category(takeEle, fabs(lepEta[i1]), fabs(lepEta[i2]),
                                   fabs(phoEta[k]), phoR9[k]);

         // works both for cut-based photon ID as well as for MVA photon ID
         if (category == 1 && phoID[k] < 0.126) continue;
         if (category == 2 && phoID[k] < 0.107) continue;
         if (category == 3 && phoID[k] < 0.126) continue;
         if (category == 4 && phoID[k] < 0.135) continue;

         pho.SetPtEtaPhiM(phoPt[k], phoEta[k], phoPhi[k], 0);

         // deltaR cuts: ISR and FSR photon rejection
         if (lep1.DeltaR(pho) < 0.4 || lep2.DeltaR(pho) < 0.4)
            continue;

         // graceful photon pt cut
         if (phoPt[k] < (Z + pho).M() * 15./110) continue;

         // select photon with highest pt
         break;
      } // photon selection

      // skip event if no photon is accepted
      if (k == nPho) continue;

      TLorentzVector H = Z + pho;
      hzg_mass = H.M();

      if (hzg_mass < 100 || hzg_mass > 190) continue;

      // FSR photon rejection
      if (Z.M() + hzg_mass < 185) continue;

      // try to take two highest-pt jets (category 5);
      // NOTE: jets are already ordered in descending order of their pt
      if (nJet >= 2) {
         jet1.SetPtEtaPhiE(jetPt[0], jetEta[0], jetPhi[0], jetEn[0]);
         jet2.SetPtEtaPhiE(jetPt[1], jetEta[1], jetPhi[1], jetEn[1]);

         TLorentzVector dijet = jet1 + jet2;

         double zeppenfeld = H.Eta() - 0.5 * (jetEta[0] + jetEta[1]);

         if (jetPt[0] >= 30 && jetPt[1] >= 30 &&
             lep1.DeltaR(jet1) >= 0.5 && lep2.DeltaR(jet1) >= 0.5 && pho.DeltaR(jet1) >= 0.5 &&
             lep1.DeltaR(jet2) >= 0.5 && lep2.DeltaR(jet2) >= 0.5 && pho.DeltaR(jet2) >= 0.5 &&
             fabs(jetEta[0] - jetEta[1]) > 3.5 && fabs(zeppenfeld) < 2.5 &&
             dijet.M() > 500 && fabs(dijet.DeltaPhi(H)) > 2.4)
            category = 5;
      } // jets

      // skip MC events which lead to "statistical spikes"
      if (hasMC && puwei > 20) {
         fprintf(stderr, "WARNING: puwei = %.1f, event skipped\n", puwei);
         continue;
      }

      // fill output TTree
      lep1Pt  = lepPt[i1];
      lep2Pt  = lepPt[i2];
      lep1Eta = lepEta[i1];
      lep2Eta = lepEta[i2];
      lep1Phi = lepPhi[i1];
      lep2Phi = lepPhi[i2];

      phoPt_  = phoPt[k];
      phoEta_ = phoEta[k];
      phoPhi_ = phoPhi[k];

      if (hasMC)
         mcwei = puwei * lepSF[i1] * lepSF[i2] * phoSF[k];

      outtree->Fill();
   } // event loop

   // flush caches
   outtree->Write("", TObject::kOverwrite);

   // cleanup
   delete intree;
   delete outtree;
   delete fi;
   delete fo;
}

//______________________________________________________________________________
void mkminitrees()
{
   /* Steering function.
    */

   // - photon identification to take;
   // - common part in paths to pre-minitrees
   int phoIDtype = PHOIDMVA;
   TString basepath = "../xanaZg-ng/preminitrees/preminitree_";

   // delete and then recreate empty output directory
    if (gSystem->Exec("rm -fr output") != 0 || gSystem->mkdir("output") != 0)
        FATAL("failed to remove or create directory \"./output\"");

   // data 2012
   mkminitree(basepath + "2ele_data2012ABCD.root", "minitree_eeg_data", ELEIDMVATrig, phoIDtype);
   mkminitree(basepath + "2muo_data2012ABCD.root", "minitree_mmg_data", MUIDTight,    phoIDtype);

   // MC samples for description of background:
   //  - Z+gamma MC
   //  - DY+jets MC
   //  - DY+jets MC with ISR and FSR photons removed
   mkminitree(basepath + "Zg_MC.root", "minitree_eeg_Zg_MC", ELEIDMVATrig, phoIDtype);
   mkminitree(basepath + "Zg_MC.root", "minitree_mmg_Zg_MC", MUIDTight,    phoIDtype);
   mkminitree(basepath + "DYJets_MC.root", "minitree_eeg_DYJets_MC", ELEIDMVATrig, phoIDtype);
   mkminitree(basepath + "DYJets_MC.root", "minitree_mmg_DYJets_MC", MUIDTight,    phoIDtype);
   mkminitree(basepath + "DYJets_MC.root", "minitree_eeg_DYJets_MC_withoutISRFSRPhotons", ELEIDMVATrig, phoIDtype, true);
   mkminitree(basepath + "DYJets_MC.root", "minitree_mmg_DYJets_MC_withoutISRFSRPhotons", MUIDTight,    phoIDtype, true);

   const char* proc[] = {"ggH", "qqH",  "WH", "ZH", "ttH"};

   // HZg signal MCs with M_H = 120-160 GeV
   for (int p = 0; p < 5; p++)
      for (int m = 0; m < 9; m++) {
         TString str = TString::Format("HZg_%s_%i", proc[p], 120 + 5 * m);
         mkminitree(basepath + str + ".root", "minitree_eeg_" + str, ELEIDMVATrig, phoIDtype);
         mkminitree(basepath + str + ".root", "minitree_mmg_" + str, MUIDTight,    phoIDtype);
      }
}
