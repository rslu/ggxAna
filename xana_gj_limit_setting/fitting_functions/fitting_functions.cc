/* Ugly but handy TF1 wrappers for several RooAbsPdf probability density
 * functions for doing unbinned fits of signal/background invariant mass shapes.
 *
 * This source code defines a common base class TF1Wrapper as well as its
 * sub-classes for the following PDFs:
 *   - Gaussian, double Gaussian, triple Gaussian, Crystal Ball + Gaussian with
 *     or without nuisance parameters for scale and resolution;
 *   - RooGaussStepBernstein of arbitrary (but reasonable) polynomial degree.
 *
 * The TF1Wrapper class adds the following functionality to RooAbsPdf:
 *   1. Allows one to treat RooAbsPdf objects as TF1 objects, thus bypassing
 *      RooFit's drawing limitations and releasing full drawing power of ROOT.
 *      In particular, the TF1 objects can be saved into root files for later
 *      drawing to be performed by separate macros. This greatly simplifies
 *      drawing-related programming.
 *   2. Enables pre-fitting using binned data and ROOT's chisquare minimization
 *      technique. The latter leads usually to much better convergence and thus
 *      circumvents RooFit's troubles with unstable/bad fits.
 *   3. Adds a consistency check of unbinned fitting results since fitting is
 *      performed to RooAbsData objects while drawing can be performed with TH1
 *      objects which were used during the pre-fitting procedure. This turns out
 *      to be useful, in particular, to verify that per-event weights were
 *      treated correctly indeed.
 *   4. Provides simple uniform interface for getting parameters' values and
 *      associated errors.
 *
 * This source code can be compiled in the following way:
 *   gROOT->LoadMacro("RooGaussStepBernstein.cxx+");
 *   gROOT->LoadMacro("fitting_functions.cc+");
 *
 * Usage example:
 *   TH1D hMass(...);
 *   RooDataSet rooMass(...);  // RooDataHist can also be used, if necessary
 *
 *   // fill hMass and rooMass with the very same data
 *   // ...
 *
 *   TF1RooDoubleGaussian fit(100, 190);  // [100, 190] = fitting region
 *
 *   // set allowed variation ranges of the fit parameters, if necessary
 *   // (they will be propagated to the unbinned fit)
 *   fit.GetTF1()->SetParLimits(1, 110, 140);
 *
 *   // set initial values of the fit parameters, if necessary
 *   // (below: total number of events and main peak position)
 *   fit.GetTF1()->SetParameter(0, hMass.Integral(1, hMass.GetNbinsX()));
 *   fit.GetTF1()->SetParameter(1, 125);
 *
 *   fit.FitTo(&hMass, &rooMass);
 *
 *   // access results, e.g.
 *   fit.GetTF1()->Write(...);
 *   fit.GetPdf()->Write(...);
 *   fit.GetExtPdf()->Write(...);
 *
 *   // in particular, access values of the fit parameters and their errors
 *   double mean = fit.GetTF1()->GetParameter(1);
 *   double err  = fit.GetTF1()->GetParError(1);
 *
 *   // parameters of the unbinned fit can be accessed directly
 *   int npar = fit.GetNPar();
 *   fit.GetPar(1)->SetName("CMS_hgz_mean");
 */

#include "fitting_functions.h"

//______________________________________________________________________________
TF1Wrapper::TF1Wrapper(double xmin, double xmax) :
   fX("x", "", xmin, xmax)
{
   /* Initialization.
    *
    * [xmin, xmax] = fitting region.
    *
    * NOTE: zeroth parameter of the fitting function (fFitUnbin and fFitBinned)
    * always corresponds to the total number of events in [xmin, xmax] region,
    * see SetPdf() for details.
    */

   fPar[0] = new RooRealVar("norm", "", 200, 0, 1e+6);

   fNPar = 1;
   fBinWidth = 1;

   fPdf = NULL;
   fFitUnbin = NULL;
   fFitBinned = NULL;
}

//______________________________________________________________________________
TF1Wrapper::~TF1Wrapper()
{
   // memory cleanup
   if (fFitBinned) delete fFitBinned;
   if (fFitUnbin) delete fFitUnbin;
   if (fPdf) delete fPdf;

   for (int i = 0; i < fNPar; i++)
      delete fPar[i];
}

//______________________________________________________________________________
void TF1Wrapper::FitTo(TH1* histo, RooAbsData* data)
{
   /* Performs unbinned (or binned in case of RooDataHist) fitting to "data"
    * after pre-fitting to "histo" with ROOT's chisquare and loglikelihood
    * minimization methods.
    *
    * NOTE: the (RooRealVar) X axis variable which was used for filling of
    * "data" must be named exactly as "x". Otherwise the fitting will not work.
    * This is because TF1Wrapper::fX is named as "x".
    *
    * This function must be called once per class instance.
    */

   // normalization for operator() return value
   fBinWidth = histo->GetBinWidth(1);

   double xmin = fX.getMin();
   double xmax = fX.getMax();

   // pre-fit binned data with chisquare and loglikelihood minimization methods
   histo->Fit(fFitBinned, "QEMN", "", xmin, xmax);
   histo->Fit(fFitBinned, "QEMNL", "", xmin, xmax);
   // propagate allowed variation ranges of parameters to fFitUnbin
   for (int i = 0; i < fNPar; i++) {
      double pmin, pmax;
      fFitBinned->GetParLimits(i, pmin, pmax);
      fPar[i]->setRange(pmin, pmax);
   }

   // perform unbinned fit
   fFitUnbin->fitTo(*data, RooFit::PrintLevel(-1),
		    RooFit::Range(xmin, xmax), RooFit::NumCPU(4),
		    RooFit::SumW2Error(kFALSE), RooFit::Hesse(kFALSE));

   // propagate back the unbinned fit parameters to fFitBinned parameters
   for (int i = 0; i < fNPar; i++) {
      fFitBinned->SetParameter(i, fPar[i]->getVal());
      fFitBinned->SetParError(i, fPar[i]->getError());
   }
}

//______________________________________________________________________________
double TF1Wrapper::operator()(double* x, double* par)
{
   /* Evaluates PDF value times normalization for given parameters and the X
    * axis value.
    *
    * NOTE: this function is not for public use. It is associated with
    * fFitBinned in SetPdf() and thus it must be declared as public for the TF1
    * wrapper to work.
    */

   fX = x[0];

   for (int i = 0; i < fNPar; i++)
      *fPar[i] = par[i];

   return fBinWidth * par[0] * fFitUnbin->getVal(RooArgSet(fX));
};

//______________________________________________________________________________
RooRealVar* TF1Wrapper::AddParameter(const char* name, double val, double vmin, double vmax)
{
   /* Adds new variable to fPar[].
    *
    * name, val = name and initial value of new variable;
    * [vmin, vmax] = its allowed range of variation during fitting.
    */

   // NOTE: check (fNpar < 64) is not done for simplicity

   fPar[fNPar] = new RooRealVar(name, "", val, vmin, vmax);
   fNPar++;

   return fPar[fNPar-1];
}

//______________________________________________________________________________
void TF1Wrapper::SetPdf(RooAbsPdf* pdf)
{
   /* Associates particular PDF with the wrapper.
    */

   fPdf = pdf;

   double xmin = fX.getMin();
   double xmax = fX.getMax();

   // define fPar[0] as the number of entries in [xmin, xmax] region
   fX.setRange("our_window", xmin, xmax);
   fFitUnbin = new RooExtendPdf("fit", "", *pdf, *fPar[0], "our_window");

   // NOTE: associate fFitBinned with TF1Wrapper::operator()
   fFitBinned = new TF1("fit", this, xmin, xmax, fNPar);

   // propagate initial values and limits to the TF1 function
   for (int i = 0; i < fNPar; i++) {
      fFitBinned->SetParameter(i, fPar[i]->getVal());
      fFitBinned->SetParLimits(i, fPar[i]->getMin(), fPar[i]->getMax());
   }

   // increase the number of sampling points used during drawings as well as in
   // all operations which require division of the X axis into small segments
   fFitBinned->SetNpx(2000);  // ROOT's default is 100
}
