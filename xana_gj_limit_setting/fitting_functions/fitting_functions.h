#ifndef FITTING_FUNCTIONS_H
#define FITTING_FUNCTIONS_H

#include <TH1.h>
#include <TF1.h>
#include <TString.h>

// RooFit
#include <RooAddPdf.h>
#include <RooRealVar.h>
#include <RooAbsData.h>
#include <RooCBShape.h>
#include <RooProduct.h>
#include <RooAddition.h>
#include <RooGaussian.h>
#include <RooExtendPdf.h>

#include "RooGaussStepBernstein.h"
#include "RooDsigmadm.h"

// common base wrapper class
class TF1Wrapper {
 public:
   TF1Wrapper(double xmin, double xmax);
   virtual ~TF1Wrapper();

   // main method
   virtual void FitTo(TH1* histo, RooAbsData* data);

   // getters
   int           GetNPar()     { return fNPar; }
   RooRealVar*   GetPar(int i) { return fPar[i]; }
   RooAbsPdf*    GetPdf()      { return fPdf; }
   RooExtendPdf* GetExtPdf()   { return fFitUnbin; }
   TF1*          GetTF1()      { return fFitBinned; }

   // not for public use: associated with fFitBinned
   double operator()(double* x, double* par);

 protected:

   RooRealVar*   AddParameter(const char* name, double val, double vmin, double vmax);
   void          SetPdf(RooAbsPdf* pdf);

   // variables
   RooRealVar    fX;           // variable associated with the X axis
   RooRealVar*   fPar[64];     // fFitUnbin's parameters
   int           fNPar;        // actual number of parameters in fPar[]
   double        fBinWidth;    // bin width of a TH1 histogram given to FitTo()
   RooAbsPdf*    fPdf;         // non-extended version of PDF
   RooExtendPdf* fFitUnbin;    // extended version of PDF
   TF1*          fFitBinned;   // TF1 version of fFitUnbin
};

// single Gaussian
class TF1RooGaussian : public TF1Wrapper {
 public:
   TF1RooGaussian(double xmin, double xmax, bool nuisance = false) :
      TF1Wrapper(xmin, xmax)
   {
      RooAbsReal* mean = AddParameter("mean", 125, xmin, xmax);
      RooAbsReal* sigma = AddParameter("sigma", 2, 0.5, 8);

      if (nuisance) {
         RooRealVar* scale = AddParameter("scale", 1, 0.5, 1.5);
         RooRealVar* resol = AddParameter("resolution", 1, 0.5, 1.5);

         scale->setConstant(true);
         resol->setConstant(true);

         // FIXME: memory is leaked below (not freed in destructor)
         mean = new RooProduct("mean_nui", "", RooArgList(*mean, *scale));
         sigma = new RooProduct("sigma_nui", "", RooArgList(*sigma, *resol));
      }

      SetPdf(new RooGaussian("Gaussian", "", fX, *mean, *sigma));

      if (nuisance) {
         fFitBinned->FixParameter(fNPar - 2, 1);
         fFitBinned->FixParameter(fNPar - 1, 1);
      }
   };
};

// double Gaussian
class TF1RooDoubleGaussian : public TF1Wrapper {
 public:
   TF1RooDoubleGaussian(double xmin, double xmax, bool nuisance = false) :
      TF1Wrapper(xmin, xmax)
   {
      RooAbsReal* mean1   = AddParameter("mean1", 125, xmin, xmax);
      RooAbsReal* sigma1  = AddParameter("sigma1", 2, 0.5, 8);
      RooRealVar* delta21 = AddParameter("delta21", 0, -50, 50); // mean2 = mean1 + delta21
      RooRealVar* s21     = AddParameter("s21", 3, 1, 30);       // sigma2/sigma1
      RooRealVar* frac    = AddParameter("frac", 0.9, 0, 1);

      // FIXME: memory is leaked below (not freed in destructor)

      RooAbsReal* mean2 = new RooAddition("mean2", "", RooArgList(*mean1, *delta21));
      RooAbsReal* sigma2 = new RooProduct("sigma2", "", RooArgList(*sigma1, *s21));

      if (nuisance) {
         RooRealVar* scale = AddParameter("scale", 1, 0.5, 1.5);
         RooRealVar* resol = AddParameter("resolution", 1, 0.5, 1.5);

         scale->setConstant(true);
         resol->setConstant(true);

         mean1 = new RooProduct("mean1_nui", "", RooArgList(*mean1, *scale));
         mean2 = new RooProduct("mean2_nui", "", RooArgList(*mean2, *scale));
         sigma1 = new RooProduct("sigma1_nui", "", RooArgList(*sigma1, *resol));
         sigma2 = new RooProduct("sigma2_nui", "", RooArgList(*sigma2, *resol));
      }

      RooGaussian* pdf1 = new RooGaussian("subpdf1", "", fX, *mean1, *sigma1);
      RooGaussian* pdf2 = new RooGaussian("subpdf2", "", fX, *mean2, *sigma2);

      SetPdf(new RooAddPdf("DoubleGaussian", "", *pdf1, *pdf2, *frac));

      if (nuisance) {
         fFitBinned->FixParameter(fNPar - 2, 1);
         fFitBinned->FixParameter(fNPar - 1, 1);
      }
   };
};

// triple Gaussian
class TF1RooTripleGaussian : public TF1Wrapper {
 public:
   TF1RooTripleGaussian(double xmin, double xmax, bool nuisance = false) :
      TF1Wrapper(xmin, xmax)
   {
      RooAbsReal* mean1   = AddParameter("mean1", 125, xmin, xmax);
      RooAbsReal* sigma1  = AddParameter("sigma1", 2, 0.5, 8);
      RooRealVar* delta21 = AddParameter("delta21", 0, -50, 50); // mean2 = mean1 + delta21
      RooRealVar* s21     = AddParameter("s21", 3, 1, 30);       // sigma2/sigma1
      RooRealVar* delta31 = AddParameter("delta31", 0, -50, 50); // mean3 = mean1 + delta31
      RooRealVar* s32     = AddParameter("s32", 3, 1, 30);       // sigma3/sigma2
      RooRealVar* frac23  = AddParameter("frac23", 0.9, 0, 1);
      RooRealVar* frac123 = AddParameter("frac123", 0.9, 0, 1);

      // FIXME: memory is leaked below (not freed in destructor)

      RooAbsReal* mean2  = new RooAddition("mean2", "", RooArgList(*mean1, *delta21));
      RooAbsReal* mean3  = new RooAddition("mean3", "", RooArgList(*mean1, *delta31));
      RooAbsReal* sigma2 = new RooProduct("sigma2", "", RooArgList(*sigma1, *s21));
      RooAbsReal* sigma3 = new RooProduct("sigma3", "", RooArgList(*sigma2, *s32));

      if (nuisance) {
         RooRealVar* scale = AddParameter("scale", 1, 0.5, 1.5);
         RooRealVar* resol = AddParameter("resolution", 1, 0.5, 1.5);

         scale->setConstant(true);
         resol->setConstant(true);

         mean1 = new RooProduct("mean1_nui", "", RooArgList(*mean1, *scale));
         mean2 = new RooProduct("mean2_nui", "", RooArgList(*mean2, *scale));
         mean3 = new RooProduct("mean3_nui", "", RooArgList(*mean3, *scale));
         sigma1 = new RooProduct("sigma1_nui", "", RooArgList(*sigma1, *resol));
         sigma2 = new RooProduct("sigma2_nui", "", RooArgList(*sigma2, *resol));
         sigma3 = new RooProduct("sigma3_nui", "", RooArgList(*sigma3, *resol));
      }

      RooGaussian* pdf1 = new RooGaussian("subpdf1", "", fX, *mean1, *sigma1);
      RooGaussian* pdf2 = new RooGaussian("subpdf2", "", fX, *mean2, *sigma2);
      RooGaussian* pdf3 = new RooGaussian("subpdf3", "", fX, *mean3, *sigma3);

      RooAbsPdf* pdf23 = new RooAddPdf("subpdf23", "", *pdf2, *pdf3, *frac23);

      SetPdf(new RooAddPdf("TripleGaussian", "", *pdf1, *pdf23, *frac123));

      if (nuisance) {
         fFitBinned->FixParameter(fNPar - 2, 1);
         fFitBinned->FixParameter(fNPar - 1, 1);
      }
   };
};

// Crystal Ball + Gaussian
class TF1RooCBGaussian : public TF1Wrapper {
 public:
   TF1RooCBGaussian(double xmin, double xmax, bool nuisance = false) :
      TF1Wrapper(xmin, xmax)
   {
      RooAbsReal* mean1   = AddParameter("mean1", 125, xmin, xmax);
      RooAbsReal* sigma1  = AddParameter("sigma1", 20, 2., 200);
      RooRealVar* alpha   = AddParameter("alpha", 2, 0.5, 10);
      RooRealVar* power   = AddParameter("power", 3, 0.5, 10);
      RooRealVar* delta21 = AddParameter("delta21", 0, -500, 500); // mean2 = mean1 + delta21
      RooRealVar* s21     = AddParameter("s21", 3, 1, 30);       // sigma2/sigma1
      RooRealVar* frac    = AddParameter("frac", 0.9, 0, 1);

      // FIXME: memory is leaked below (not freed in destructor)

      RooAbsReal* mean2  = new RooAddition("mean2", "", RooArgList(*mean1, *delta21));
      RooAbsReal* sigma2 = new RooProduct("sigma2", "", RooArgList(*sigma1, *s21));

      if (nuisance) {
         RooRealVar* scale = AddParameter("scale", 1, 0.5, 1.5);
         RooRealVar* resol = AddParameter("resolution", 1, 0.5, 1.5);

         scale->setConstant(true);
         resol->setConstant(true);

         mean1 = new RooProduct("mean1_nui", "", RooArgList(*mean1, *scale));
         mean2 = new RooProduct("mean2_nui", "", RooArgList(*mean2, *scale));
         sigma1 = new RooProduct("sigma1_nui", "", RooArgList(*sigma1, *resol));
         sigma2 = new RooProduct("sigma2_nui", "", RooArgList(*sigma2, *resol));
      }

      RooCBShape*  pdf1 = new RooCBShape ("subpdf1", "", fX, *mean1, *sigma1, *alpha, *power);
      RooGaussian* pdf2 = new RooGaussian("subpdf2", "", fX, *mean2, *sigma2);

      SetPdf(new RooAddPdf("CBGaussian", "", *pdf1, *pdf2, *frac));

      if (nuisance) {
         fFitBinned->FixParameter(fNPar - 2, 1);
         fFitBinned->FixParameter(fNPar - 1, 1);
      }
   };
};

// RooGaussStepBernstein
class TF1RooGaussStepBernstein : public TF1Wrapper {
 public:
   TF1RooGaussStepBernstein(double xmin, double xmax, int poldeg) :
      TF1Wrapper(xmin, xmax)
   {
      RooRealVar* mean = AddParameter("mean", 0, -50, 50);
      RooRealVar* sigma = AddParameter("sigma", 4, 2, 10);
      RooRealVar* stepval = AddParameter("stepval", 110, 50, 200);

      RooArgList coeflist;

      // collect polynomial's coefficients (the zeroth coefficient is always 1)
      for (int i = 1; i <= poldeg; i++)
         coeflist.add(*AddParameter(TString::Format("coef%i", i), 3, 0, 30));

      // may improve convergence
      if (poldeg > 0) *fPar[4] = 15;

      SetPdf(new RooGaussStepBernstein("GaussStepBernstein", "", fX, *mean, *sigma, *stepval, coeflist));
   };
};

// RooDsigmadm
class TF1RooDsigmadm : public TF1Wrapper {
 public:
   TF1RooDsigmadm(double xmin, double xmax) :
      TF1Wrapper(xmin, xmax)
   {
      RooRealVar* p1 = AddParameter("p1", 1., 0., 10.);
      RooRealVar* p2 = AddParameter("p2", 10, 0., 30);
      RooRealVar* p3 = AddParameter("p3", 1.2, 1.1, 1.3);


      SetPdf(new RooDsigmadm("Dsigmadm", "", fX, *p1, *p2, *p3));
   };
};

#endif
