/* Performs unbinned fitting of invariant mass shapes and creates text datacards
 * for the "combine" tool from the HiggsAnalysis/CombinedLimit package.
 *
 * Minitrees are taken from output/minitrees.root. The following branches are
 * assumed to exist in each minitree: category, gj_mass, mcwei (for MC), see
 * fill_events().
 *
 * Signal PDFs are evaluated from minitrees corresponding to existing MC
 * productions (i.e. for certain Higgs mass values). At the same time, signal
 * PDFs for intermediate Higgs masses (for which no MC productions exist) are
 * extrapolated from two PDFs corresponding to closest existing MC productions.
 *
 * The PDF for background description as well as the expected background yield
 * are both obtained from real data.
 *
 * TH1 versions of invariant mass shapes, TF1 versions of produced PDFs as well
 * as normalization factors (RooDouble) which adapt MC scales to real data are
 * written into output/for_visualizations.root.
 *
 * Observed mass points (from real data) together with all obtained PDFs are
 * written into per-channel (eeg or mmg) and per-category (cat1, ..., cat4)
 * files output/datacards/for_datacards_*.root. Datacards themselves are
 * produced as plain text files output/datacards/datacard_*.txt.
 *
 * The directory "./output/datacards" is removed every time this macro is
 * executed.
 *
 * Actual configuration is encoded directly into functions mkdatacard(),
 * normfactor_SMHiggs_8TeV() and mkdatacards().
 *
 * Documentation:
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombination
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSDataAnalysisSchool2014HiggsCombPropertiesExercise?rev=11
 * https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit?rev=111
 * https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationConventions?rev=18
 * https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationInput?rev=11
 * https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/HiggsCombinationInputUnbinnedDraft?rev=2
 *
 * Usage: root -b -l -q mkdatacards.C
 */

// prints a message and exits gracefully
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)

//______________________________________________________________________________
//RooAbsData* fill_events(TH1* hMass, const char* treename, int cat = -1, bool usewei = true)
RooAbsData* fill_events(TH1* hMass, const char* fname, int cat = -1, bool usewei = true)
{
   /* Fills hMass + a RooDataSet with invariant mass values from given minitree.
    * Returns the filled RooDataSet object.
    *
    * treename = name of minitree in minitrees.root to take;
    * cat = if > 0, take only events from this particular category;
    * usewei = (for MC) if true, attribute weights to all events.
    */

   // open file and get requested minitree
  //TFile* fi = TFile::Open("output/minitrees.root");
   TFile* fi = TFile::Open(fname);
   if (!fi || fi->IsZombie())
      FATAL("TFile::Open() failed");

   TTree* tree = (TTree*)fi->Get("tt");//dynamic_cast<TTree*> (fi->Get(treename));
   if (!tree) FATAL("TFile::Get() failed");

   // variables to be associated with minitree branches
   Int_t category=0;
   float gjmass;
   float mcwei;

   // associate tree branches with variables
//    if (tree->SetBranchAddress("category", &category) != 0)
//       FATAL("TTree::SetBranchAddress() failed");

   if (tree->SetBranchAddress("gjmass", &gjmass) != 0)
      FATAL("TTree::SetBranchAddress() failed");

   if (usewei && tree->SetBranchAddress("mcwei", &mcwei) != 0)
      FATAL("TTree::SetBranchAddress() failed");

   // NOTE: the X axis variable below (invariant mass) is named exactly as "x"
   // since otherwise unbinned fitting will not work. This is because
   // TF1Wrapper::fX is named as "x"
   RooRealVar x("x", "", 400, 500, 4500);
   RooRealVar w("w", "", 1);  // NOTE: weights are ignored without this variable
   RooDataSet* rooMass;
   if (usewei)
      rooMass = new RooDataSet("rooMass", "", RooArgSet(x, w), RooFit::WeightVar(w));
   else
      rooMass = new RooDataSet("rooMass", "", RooArgSet(x));

   // fill hMass and rooMass
   for (Long64_t ev = 0; ev < tree->GetEntriesFast(); ev++) {
      if (tree->GetEntry(ev) <= 0)
         FATAL("TTree::GetEntry() failed");

      if (cat >= 0 && category != cat)
         continue;
      if (gjmass<500) continue;
      x = gjmass;

      if (usewei) {
         hMass->Fill(gjmass, mcwei);
         rooMass->add(RooArgSet(x), mcwei);
      } else {
         hMass->Fill(gjmass);
         rooMass->add(RooArgSet(x));
      }
   }

   delete tree;
   delete fi;

   return rooMass;
}

//______________________________________________________________________________
double normfactor_SMHiggs_8TeV(const char* channel, int p, int m)
{
   /* Values were calculated by normfactors_and_theoretical_uncertainties.py.
    *
    * channel = either of "eeg" (for electrons) or "mmg" (for muons);
    * p = production process number (0=ggH, 1=qqH, 2=WH, 3=ZH, 4=ttH);
    * mass = value of SM Higgs mass (from the region 120-160GeV with step 1GeV).
    */


//    int index = mass - 120;
   int index = m;
   if (index < 0 || index > 9)
      FATAL("wrong mass given");

//    if (TString(channel).Contains("ee"))
//       return norm_eeg[p][index];
//    return norm_mmg[p][index];

   //int mass_exist[9] = {700, 1000, 1200, 1500, 1700, 2000, 2500, 3000, 3500 };
//    double norm_gj[10]={18.497023, 2.753563, 0.931327, 0.224205, 0.092183, 
//  		       0.012509, 0.001691, 0.000238, 0.001206};
    double norm_gj[10]={0.98701, 0.16431, 0.06038, 0.01638, 0.007332,
 		       0.001122, 0.0001881, 0.000033791, 0.000230180};
   
   
   return norm_gj[m];

}

//______________________________________________________________________________
void mkdatacard(const char* channel = "gj", int cat = 0)
{
   /* Does everything for particular channel and particular event category.
    *
    * channel = either of "eeg" (for electrons) or "mmg" (for muons);
    * cat = 1, 2, 3 or 4.
    *
    * NOTE: we follow naming conventions of
    * https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
    *
    * NOTE: in this function, memory is leaked (via "new") for coding brevity.
    */

   // cat1, cat2, cat3 or cat4
   TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("gj_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   if (!fo || fo->IsZombie())
      FATAL("TFile::Open() failed");

   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_8TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");

   //
   // evaluate PDF of background from real data
   //

   TH1D hDataObs("hDataObs", "", 100, 500., 4500.);
   hDataObs.SetXTitle("M_{#gamma jet} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");

   // false = do not add per-event weights
   TString treename = TString::Format("tt", channel);
   //RooAbsData* dataObs = fill_events(&hDataObs, treename, cat, false);
   RooAbsData* dataObs = fill_events(&hDataObs, "8TeV_2fb/ana_gj_out_data.root", cat, false);

   TF1RooDsigmadm bgrfit(500, 4500); // 5 = polynomial degree
   bgrfit.FitTo(&hDataObs, dataObs);

   // save TH1 and TF1 objects into root file for later visualization
   wd->cd();
   hDataObs.Write("", TObject::kOverwrite);
   bgrfit.GetTF1()->Write("fit_hDataObs", TObject::kOverwrite);

   // rename the X axis variable and add observed mass points into the workspace
   dataObs->SetName("data_obs");
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_gj_mass")))
      FATAL("RooWorkspace::import() failed");

   // - rename parameters of the background PDF;
   // - declare parameters of the background PDF to be constants (otherwise the
   //   parameters will be considered as freely floating even without using the
   //   "flatParam" directives below)
   for (int i = 0; i < bgrfit.GetNPar(); i++) {
     bgrfit.GetPar(i)->Print();
      const char* name = bgrfit.GetPar(i)->GetName();
      bgrfit.GetPar(i)->SetName(TString::Format("CMS_gj_%s_%s_8TeV_%s",
                                                channel, cats, name));
      bgrfit.GetPar(i)->setConstant(true);
   }

   // - make normalization utilized by the combine tool to be equal to expected
   //   number of events written into datacards;
   // - add (extended version of) background PDF into the workspace under the
   //   name "pdf_bgr";
   // - add suffix "_bgr" to all PDF's parameters (and to subPDFs, if any);
   // - connect the PDF's X axis variable to the X axis variable of data_obs
   bgrfit.GetPar(0)->setVal(1);
   bgrfit.GetExtPdf()->SetName("pdf");
   if (wspace.import(*bgrfit.GetExtPdf(), RooFit::RenameAllNodes("bgr"),
                                          RooFit::RenameAllVariablesExcept("bgr", "x"),
                                          RooFit::RenameVariable("x", "CMS_gj_mass")))
      FATAL("RooWorkspace::import() failed");

   // total number of events observed/expected in the region [100, 190]
   int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
   double expectation = bgrfit.GetTF1()->GetParameter(0);

   // important consistency check (inconsistency may happen e.g. if by a mistake
   // the mass region in mkminitrees.cc is different from [100, 190]
   printf("observation %d, expectation %.2f, dataobs %.2f \n", observation, expectation,
	  dataObs->sumEntries());
   
   //if (observation != dataObs->sumEntries())
   if (TMath::Abs(observation - dataObs->sumEntries())/observation > 0.01)
      FATAL("observation != dataObs->sumEntries()");

   //
   printf(" evaluate signal PDFs for existing MC productions \n");
   //

   // - references to signal PDFs for existing MC productions
   // - expected signal yields for all mass points
   TF1Wrapper* sigfit[5][9]; // [proc][mass]
   double    expected[5][41];

   // Higgs production processes & Higgs masses for existing MC productions
   const char* proc[5] = {"gj", "qqH", "WH", "ZH", "ttH"};
   int mass_exist[9] = {700, 1000, 1200, 1500, 1700, 2000, 2500, 3000, 3500 };

   for (int p = 0; p < 1; p++)
       for (int m = 0; m < 9; m++) {
         int mass = mass_exist[m];

         TString hname = TString::Format("hMass_%s_%i", proc[p], mass);

         TH1D hMass(hname, "", 400, 500, 4500);
         hMass.SetXTitle("M_{#gamma jet} (GeV/c^{2})");
         hMass.SetYTitle("Counts");
         hMass.Sumw2();

//          TString treename = TString::Format("minitree_%s_gj_%s_%i",
//                                             channel, proc[p], mass);
         TString fname = TString::Format("8TeV_2fb/ana_gj_out_job_qstar_%i.root",
                                            mass_exist[m]);
         RooAbsData* rooMass = fill_events(&hMass, fname, cat);

         // true = use nuisance parameters for the energy scale and resolution
         TF1Wrapper* sfit = new TF1RooCBGaussian(500, 4500, true);

         // set initial values and limits
         sfit->GetTF1()->SetParameter(0, hMass.Integral(1, hMass.GetNbinsX()));
         sfit->GetTF1()->SetParameter(1, mass);
         //sfit->GetTF1()->SetParLimits(1, 0.5 * mass, 1.5 * mass);

         sfit->FitTo(&hMass, rooMass);

         // save TH1 and TF1 objects into root file for later visualization
         wd->cd();
         hMass.Write("", TObject::kOverwrite);
         sfit->GetTF1()->Write("fit_" + hname, TObject::kOverwrite);

         // save the normalization factor which adapts the MC scale to real
         // data; evaluate expected signal yield
         double norm = normfactor_SMHiggs_8TeV(channel, p, m);
         RooDouble(norm).Write("norm_" + hname, TObject::kOverwrite);
         expected[p][m] = norm * sfit->GetTF1()->GetParameter(0);

         // unique suffix
         TString sfx = TString::Format("sig_%s_%i", proc[p], mass[m]);

         // - rename parameters of the PDF;
         // - declare parameters of the PDF to be constants (otherwise the
         //   parameters will be considered as freely floating, and the combine
         //   tool will produce weird results)
         for (int i = 0; i < sfit->GetNPar(); i++) {
            const char* name = sfit->GetPar(i)->GetName();
            sfit->GetPar(i)->SetName(TString::Format("CMS_gj_%s_%s_8TeV_%s_%s",
                                                     channel, cats, name, sfx.Data()));
            sfit->GetPar(i)->setConstant(true);
         }

         // set names for the nuisance parameters of energy scale and resolution
         sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
         sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));

         // - add signal PDF into the workspace under the name "pdf_sig_...";
         // - add suffix "_sig_..." to all subPDFs;
         // - connect the PDF's X axis variable to the X axis variable of data_obs
         // - connect the nuisance parameters of different signal PDFs together
         sfit->GetPdf()->SetName("pdf");
         if (wspace.import(*sfit->GetPdf(), RooFit::RenameAllNodes(sfx),
                                            RooFit::RenameVariable("x", "CMS_gj_mass")))
            FATAL("RooWorkspace::import() failed");

         sigfit[p][m] = sfit;
      } // process and mass loops

   //
   printf("extrapolate signal PDFs for intermediate mass points \n");
   //

   // use 1 GeV steps
//    for (int p = 0; p < 5; p++)
//       for (int m = 0; m < 8; m++)
//          for (int k = 1; k <= 4; k++) {
   for (int p = 0; p < 1; p++)
      for (int m = 0; m < 9; m++)
         for (int k = 1; k <= 1; k++) {
	   int mass = mass_exist[m]; //120 + m * 5 + k;

            // neighbours
//             double m1 = mass_exist[m];
//             double m2 = mass_exist[m + 1];

            // proportions of first/second neighbour
//             double a = (m2 - mass)/(m2 - m1);
//             double b = 1 - a;

            // true = use nuisance parameters for the energy scale and resolution
            TF1Wrapper* sfit = new TF1RooCBGaussian(500, 4500, true);

            // set values of parameters from linear extrapolation
             for (int i = 0; i < sfit->GetNPar() - 2; i++) {
//                double val1 = sigfit[p][m]    ->GetTF1()->GetParameter(i);
//                double val2 = sigfit[p][m + 1]->GetTF1()->GetParameter(i);

	       sfit->GetPar(i)->setVal(sigfit[p][m]    ->GetTF1()->GetParameter(i));
             }
	     
            // evaluate expected signal yield
            double norm = normfactor_SMHiggs_8TeV(channel, p, m);
            //expected[p][m * 5 + k] = norm * sfit->GetPar(0)->getVal();
            expected[p][m] = norm * sfit->GetPar(0)->getVal();

            // unique suffix
            TString sfx = TString::Format("sig_%s_%i", proc[p], mass);

            // - rename parameters of the PDF;
            // - declare parameters of the PDF to be constants (otherwise the
            //   parameters will be considered as freely floating, and the combine
            //   tool will produce weird results)
            for (int i = 0; i < sfit->GetNPar(); i++) {
               const char* name = sfit->GetPar(i)->GetName();
               sfit->GetPar(i)->SetName(TString::Format("CMS_gj_%s_%s_8TeV_%s_%s",
                                                      channel, cats, name, sfx.Data()));
               sfit->GetPar(i)->setConstant(true);
            }

            // set names for the nuisance parameters of energy scale and resolution
            sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
            sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));

            // - add signal PDF into the workspace under the name "pdf_sig_...";
            // - add suffix "_sig_..." to all subPDFs;
            // - connect the PDF's X axis variable to the X axis variable of data_obs
            // - connect the nuisance parameters of different signal PDFs together
            sfit->GetPdf()->SetName("pdf");
            if (wspace.import(*sfit->GetPdf(), RooFit::RenameAllNodes(sfx),
                                               RooFit::RenameVariable("x", "CMS_gj_mass")))
               FATAL("RooWorkspace::import() failed");
         } // extrapolation

   wspace.writeToFile(TString::Format("output/datacards/for_datacards_gj_%s_%s_8TeV.root",
                                      channel, cats));

   // close output/for_visualizations.root
   delete fo;

   //
   printf(" produce datacards \n");
   //

   // cuttent UTC time
   TString timestamp(TTimeStamp().AsString());
   timestamp.Remove(timestamp.Last(':'));

   // name of final state
   TString binString = TString::Format("%s_%s_8TeV", channel, cats);
   const char* bin = binString.Data();

   //for (int m = 0; m < 41; m++) {
   for (int m = 0; m < 9; m++) {
      int mass = mass_exist[m];

      TString datacard;
      datacard += "# Datacard for the gj analysis for limit setting\n";
      datacard += "# National Central University, Taiwan\n";
      datacard += "# " + timestamp + " (UTC)\n";
      datacard += TString::Format("# Usage: combine -U -M Asymptotic -m %i datacard.txt\n", mass);
      datacard += "#\n";
      datacard += "imax 1  # number of final states\n";
      datacard += "jmax *  # number of yields given below minus one\n";
      datacard += "kmax *  # number of sources of systematical uncertainties (nuisance parameters)\n";
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      datacard += TString::Format("shapes  *         * for_datacards_gj_%s_%s_8TeV.root gj_workspace:pdf_sig_$PROCESS_%i\n", channel, cats, mass);
      datacard += TString::Format("shapes  bgr       * for_datacards_gj_%s_%s_8TeV.root gj_workspace:pdf_bgr\n", channel, cats);
      datacard += TString::Format("shapes  data_obs  * for_datacards_gj_%s_%s_8TeV.root gj_workspace:data_obs\n", channel, cats);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
      datacard += TString::Format("bin            %s\n", bin);
      datacard += TString::Format("observation    %i\n", observation);
      datacard += "------------------------------------------------------------------------------------------------------------\n";
//       datacard += TString::Format("bin            %-15s %-15s %-15s %-15s %-15s %-15s\n", bin, bin, bin, bin, bin, bin);
//       datacard +=                 "process        ttH             ZH              WH              qqH             ggH             bgr\n";
//       datacard +=                 "process        -4              -3              -2              -1              0               1\n";
//       datacard += TString::Format("rate           %-15f %-15f %-15f %-15f %-15f %-15f\n",
//                                   expected[4][m], expected[3][m], expected[2][m], expected[1][m], expected[0][m], expectation);
      datacard += TString::Format("bin            %-15s %-15s\n", bin, bin);
      datacard +=                 "process        gj             bgr\n";
      datacard +=                 "process        0               1\n";
      datacard += TString::Format("rate           %-15f %-15f\n", expected[0][m], expectation);
      datacard += "------------------------------------------------------------------------------------------------------------\n";

      // // theoretical uncertainties
//       ifstream in("theoretical_uncertainties_SM_Higgs.txt");
//       int count = 0; // simple error protection
//       char line[4096];

//       while (in.good()) {
//          in.getline(line, 4096);
//          if (!in.eof() && in.fail())
//             FATAL("ifstream::getline() failed");

//          TString s = line;

//          if (s.BeginsWith(TString::Format("%.1f ", (double)mass))) {
//             s.Remove(0, s.First(':') + 1); // remove e.g. "155.0   :"
//             datacard += s + "\n";
//             count++;
//          }
//       }
//       if (count != 7) FATAL("theoretical uncertainties: line count != 7");


      datacard += "pdf_gj          lnN   1.020            -\n";


      // CMS uncertainties.
      //
      // NOTE: naming conventions are from
      // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsCombinationConventions?rev=18
      // TODO: updated numbers?

      datacard += "lumi_8TeV       lnN   1.026            -\n";

      // electron vs muon efficiency uncertainties, including uncertainties of
      // triggering efficiencies (it is assumed that the two are not correlated)
      if (TString(channel).Contains("ee"))
         datacard += "CMS_eff_e_8TeV  lnN   1.022            -\n";
      else
         datacard += "CMS_eff_m_8TeV  lnN   1.038            -\n";

      // photon efficiency uncertainties: ECAL barrel vs ECAL endcaps
      if (cat == 1 || cat == 2 || cat == 3)
         datacard += "CMS_eff_g_8TeV  lnN  1.006            -\n";
      else
         datacard += "CMS_eff_g_8TeV  lnN  1.010            -\n";

      // uncertainty on photon's R9 cut
      if (cat == 1 || cat == 2)
         datacard += "CMS_eff_R9_EB_8TeV lnN 1.050            -\n";

      float UEPS[] = {1.026, 1.035, 1.018, 1.021};
      datacard += TString::Format("CMS_UEPS_8TeV   lnN   1.002            -\n", UEPS[cat]);

      float JEC[] = {1.028, 1.032, 1.022, 1.022};
      datacard += TString::Format("CMS_JEC_8TeV    lnN   1.001            -\n", JEC[cat]);

      float JER[] = {1.010, 1.011, 1.011, 0.014};
      datacard += TString::Format("CMS_JER_8TeV    lnN   1.001            -\n", JER[cat]);

      //if (TString(channel).Contains("ee"))
      //   datacard += "pdf_PU_mu_8TeV     1.008            -\n";
      //else
      //   datacard += "pdf_PU_mu_8TeV     1.004            -\n";

      datacard += "------------------------------------------------------------------------------------------------------------\n";

      // declare uncertainties for the energy scale and resolution
      // (as described by the nuisance parameters above)
      // TODO: errors are set by hand
      datacard += TString::Format("CMS_scale_%s                                param          1     0.05\n", channel);
      datacard += TString::Format("CMS_res_%s                                  param          1     0.01\n", channel);

      // declare parameters of the background PDF to be freely floating
      //for (int i = 0; i < bgrfit.GetNPar(); i++) {
      for (int i = 0; i < 2; i++) {
         TString nameStr = TString::Format("%s_bgr", bgrfit.GetPar(i)->GetName());
         datacard += TString::Format("%-44s flatParam\n", nameStr.Data());
      }

      // make datacard file
      FILE* out = fopen(TString::Format("output/datacards/datacard_gj_%s_%s_8TeV_%i.txt",
                                        channel, cats, mass[m]).Data(), "w");
      if (!out) FATAL("fopen() failed");
      if (fputs(datacard.Data(), out) < 0) FATAL("fputs() failed");
      if (fclose(out) != 0) FATAL("fclose() failed");
   } // mass loop
}

//______________________________________________________________________________
void mkdatacards()
{
   /* Steering function.
    */

   // FIXME: this line is required on lxplus and pcncu1X machines
   gSystem->SetIncludePath("-I$ROOFITSYS/include");

   gROOT->LoadMacro("fitting_functions/RooGaussStepBernstein.cxx+");
   gROOT->LoadMacro("fitting_functions/RooDsigmadm.cxx+");
   gROOT->LoadMacro("fitting_functions/fitting_functions.cc+");

   // suppress useless messages from RooFit
//    RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
   RooFit::RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

   // prevent segfaults due to various TH1 objects not accurately attributed to
   // correct TDirectory
   TH1::AddDirectory(0);

   // prepare empty directory for datacards
   if (gSystem->Exec("rm -fr output/datacards") != 0 || gSystem->mkdir("output/datacards") != 0)
      FATAL("failed to remove or create directory \"./output/datacards\"");

   // remove previous version of for_visualizations.root, if any
   if (!gSystem->AccessPathName("output/for_visualizations.root"))
      if (gSystem->Unlink("output/for_visualizations.root") != 0)
         FATAL("TSystem::Unlink() failed");

   // actual work
   for (int cat = 0; cat < 1; cat++) {
      mkdatacard("gj", cat);
   }
}
