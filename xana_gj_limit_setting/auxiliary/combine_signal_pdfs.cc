/* Visualization of signal PDFs.
 *
 * This code requires compilation due to a bug in CINT.
 *
 * Usage: root -l combine_signal_pdfs.cc+
 */

#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooDouble.h>

// linear combination of TF1 objects
class LinearCombination {
 public:
   LinearCombination(TF1* pdfs[], double norms[], int npdfs) {
      for (int i = 0; i < npdfs; i++) {
         fPdf[i] = (TF1*) pdfs[i]->Clone();
         fNorm[i] = norms[i];
      }
      fNPdfs = npdfs;

      // - associate "fun" with LinearCombination::operator();
      // - improve drawing quality (get rid of drawing artifacts)
      fun = new TF1("fun", this, 100, 190, 0);
      fun->SetNpx(2000);
   }

   double operator() (double* x, double* /* pars */) {
      double result = 0;

      for (int i = 0; i < fNPdfs; i++)
         result += fNorm[i] * fPdf[i]->Eval(x[0]);

      return result;
   }

   TF1* fun;

 private:
   int     fNPdfs;
   TF1*    fPdf[64];
   double  fNorm[64];
};

//______________________________________________________________________________
void combine_single_mass(const char* path, const char* channel, double mH)
{
   /* Combines PDFs for single Higgs mass point on one canvas.
    */

   const char* procs[] = {"ggH", "qqH", "WH", "ZH", "ttH"};

   TF1* pdfs[5];
   double norms[5];

   TFile f(path);

   // collect signal PDFs and their normalizations for given mass point
   for (int i = 0; i < 5; i++) {
      TString name;

      name.Form("%s/fit_hMass_%s_%.0f", channel, procs[i], mH);
      pdfs[i] = dynamic_cast<TF1*>(f.Get(name));

      name.Form("%s/norm_hMass_%s_%.0f", channel, procs[i], mH);
      norms[i] = (double) (*dynamic_cast<RooDouble*>(f.Get(name)));
   }

   TString cname = TString::Format("signals_%s_mH%.0f", channel, mH);
   new TCanvas(cname, cname, 700, 700);

   gPad->SetLeftMargin(0.11);
   gPad->SetRightMargin(0.04);
   gPad->SetTopMargin(0.06);
   gPad->SetBottomMargin(0.09);
   gPad->SetGridx();
   gPad->SetGridy();
   gPad->SetLogy();

   TF1* sum = (new LinearCombination(pdfs, norms, 5))->fun;

   sum->GetHistogram()->SetTitle(TString::Format("%s_mH%.0f", channel, mH));
   sum->GetHistogram()->SetXTitle("M_{ll#gamma}^{rec} (GeV/c^{2})");
   sum->GetHistogram()->SetYTitle("Number of events");
   sum->GetHistogram()->SetTitleOffset(1.2, "X");
   sum->GetHistogram()->SetTitleOffset(1.5, "Y");

   sum->SetLineColor(kBlack);
   sum->Draw();

   TLegend* leg = new TLegend(0.8, 0.7, 0.94, 0.92);

   leg->AddEntry(sum, "total sum", "l");

   int colors[] = {kBlue, 8, kPink, kOrange, kYellow};

   for (int i = 0; i < 5; i++) {
      TF1* one = (new LinearCombination(&pdfs[i], &norms[i], 1))->fun;

      one->SetLineColor(colors[i]);
      one->SetLineWidth(1);
      one->Draw("same");

      leg->AddEntry(one, TString::Format("%s", procs[i]), "l");
   }

   leg->SetFillColor(kWhite);
   leg->Draw("same");
}

//______________________________________________________________________________
void combine_pdfs(const char* path, const char* channel)
{
   /* Combines PDFs for several Higgs mass points on one canvas.
    */

   int nM = 9;
//    double mH[] = {120, 130, 140, 150, 160};
   double mH[] = {120, 125, 130, 135, 140, 145, 150, 155, 160};

   TF1* sums[64];

   TFile f(path);

   // collect sums of PDFs
   for (int k = 0; k < nM; k++) {
      const char* procs[] = {"ggH", "qqH", "WH", "ZH", "ttH"};

      TF1* pdfs[5];
      double norms[5];

      // collect signal PDFs and their normalizations for given mass point
      for (int i = 0; i < 5; i++) {
         TString name;

         name.Form("%s/fit_hMass_%s_%.0f", channel, procs[i], mH[k]);
         pdfs[i] = dynamic_cast<TF1*>(f.Get(name));

         name.Form("%s/norm_hMass_%s_%.0f", channel, procs[i], mH[k]);
         norms[i] = (double) (*dynamic_cast<RooDouble*>(f.Get(name)));
      }

      sums[k] = (new LinearCombination(pdfs, norms, 5))->fun;
   }

   // draw results
   TString cname = TString::Format("signals_%s", channel);
   new TCanvas(cname, cname, 700, 700);

   gPad->SetLeftMargin(0.11);
   gPad->SetRightMargin(0.04);
   gPad->SetTopMargin(0.06);
   gPad->SetBottomMargin(0.09);
   gPad->SetGridx();
   gPad->SetGridy();

   TLegend* leg = new TLegend(0.8, 0.7, 0.94, 0.92);

   for (int k = 0; k < nM; k++) {
      sums[k]->GetHistogram()->SetTitle(channel);
      sums[k]->GetHistogram()->SetXTitle("M_{ll#gamma}^{rec} (GeV/c^{2})");
      sums[k]->GetHistogram()->SetYTitle("Number of events");
      sums[k]->GetHistogram()->SetTitleOffset(1.2, "X");
      sums[k]->GetHistogram()->SetTitleOffset(1.5, "Y");
      sums[k]->GetHistogram()->SetAxisRange(0, 0.45, "Y");

      sums[k]->SetLineWidth(2);
      sums[k]->SetLineColor(kBlack);
      sums[k]->Draw(k == 0 ? "" : "same");

      leg->AddEntry(sums[k], TString::Format("mH = %.0f", mH[k]), "l");
   }

//    leg->SetFillColor(kWhite);
//    leg->Draw("same");
}

//______________________________________________________________________________
void combine_signal_pdfs()
{
   /* Steering function.
    */

   // configuration
   const char* path = "../output/for_visualizations.root";
   const char* channel = "mmg_cat4_8TeV";

   combine_single_mass(path, channel, 125);
   combine_single_mass(path, channel, 160);

   combine_pdfs(path, channel);
}
