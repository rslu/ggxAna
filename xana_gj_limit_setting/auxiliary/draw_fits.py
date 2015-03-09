#!/usr/bin/env python
"""Visualization of fits for background as well as for HZg signals.
"""

# python-2 compatibility
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

import ROOT

# for keeping drawed objects in memory
saves = []

def main():
    """Steering function.

    TODO: error bars are not quite correct since non-unit event weights of MCs
    are not propagated accurately.
    """
    # configuration
    channel = 'gj_cat0_8TeV'
    mH = 2000
    procs = ['gj']

    f = ROOT.TFile('../output/for_visualizations.root')

    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(0)

    # disable RooHist warnings 'non-integer bin entry with Poisson errors'
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

    # background
    bgr = f.Get(channel + '/hDataObs')
    bgrfit = f.Get(channel + '/fit_hDataObs')

    cname = 'background_{0}'.format(channel)
    draw(bgr, bgrfit, cname, channel, 'data 2012, 2.0 fb^{-1}')

    # signals
    for proc in procs:
        sig = f.Get('{0}/hMass_{1}_{2}'.format(channel, proc, mH))
        sigfit = f.Get('{0}/fit_hMass_{1}_{2}'.format(channel, proc, mH))

        cname = 'signal_{0}_{1}_{2}'.format(channel, proc, mH)
        legtxt = '{0} signal, mH={1} GeV/c^{{2}}'.format(proc, mH)
        draw(sig, sigfit, cname, channel, legtxt)

def draw(histo, fit, cname, title, legtxt):
    """Draws histogram with its fit as well as histogram/fit ratio.
    """
    c = ROOT.TCanvas(cname, cname, 1000, 700)
    saves.append(c)

    pad1 = ROOT.TPad('top', 'top', 0, 0.3, 1, 1)
    pad2 = ROOT.TPad('bottom', 'bottom', 0, 0, 1, 0.3)
    pad1.Draw()
    pad2.Draw()
    saves.append((pad1, pad2))

    # top pad
    pad1.cd()
    pad1.SetGridx()
    if 'Data' in legtxt:
        pad1.SetLogy()
    
    pad1.SetGridy()
    pad1.SetTopMargin(0.07)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.05)
    pad1.SetRightMargin(0.01)

    xtitle = 'M_{ee#gamma}^{rec} (GeV/c^{2})'
    if 'mmg' in title:
        xtitle = 'M_{#mu#mu#gamma}^{rec} (GeV/c^{2})'

    gr = ROOT.RooHist(histo) # make errors asymmetric
    gr.SetLineColor(ROOT.kBlack)
    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.7)
    gr.GetHistogram().SetTitle(title)
    gr.GetHistogram().SetXTitle(xtitle)
    gr.GetHistogram().SetYTitle('Entries')
    gr.GetHistogram().SetTitleSize(0.037, 'XY')
    gr.GetHistogram().SetLabelSize(0.037, 'XY')
    gr.GetHistogram().SetTitleOffset(0.7, 'Y')
    gr.GetHistogram().SetAxisRange(500, 3700, 'X')
    gr.GetHistogram().SetAxisRange(0, histo.GetMaximum() * 1.2, 'Y')
    if 'Data' in histo:
        gr.GetHistogram().SetAxisRange(0, histo.GetMaximum() * 1.5, 'Y')
    gr.Draw('APZ')

    fit.SetLineColor(ROOT.kRed)
    fit.SetLineWidth(1)
    fit.Draw('same')

    # legend
    leg = ROOT.TLegend(0.69, 0.73, 0.98, 0.91)
    leg.AddEntry(gr, legtxt, 'p')
    leg.AddEntry(fit, 'background model', 'l')
    leg.AddEntry(fit, '#chi^{{2}}/ndf of binned fit = {0:.3f}'.format(fit.GetChisquare()/fit.GetNDF()), 'l')
    leg.SetFillColor(0)
    leg.Draw('same')

    saves.append((gr, fit, leg))

    # bottom pad
    pad2.cd()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetLeftMargin(0.05)
    pad2.SetRightMargin(0.01)
    pad2.SetGridx()
    pad2.SetGridy()

    # gr/fit ratio with asymmetric errors
    grRat = ROOT.TGraphAsymmErrors()
    saves.append(grRat)

    # fill grRat
    for i in range(gr.GetN()):
        (x, exl, exh) = (gr.GetX()[i], gr.GetEXlow()[i], gr.GetEXhigh()[i])
        (y, eyl, eyh) = (gr.GetY()[i], gr.GetEYlow()[i], gr.GetEYhigh()[i])

        y0 = fit.Eval(x)

        # do not show points with no entries
        if y == 0 or y0 == 0:
            grRat.SetPoint(i, x, -1)
            grRat.SetPointError(i, exl, exh, 0, 0)
            continue

        grRat.SetPoint(i, x, y/y0)
        grRat.SetPointError(i, exl, exh, eyl/y0, eyh/y0)

    grRat.SetMarkerStyle(20)
    grRat.SetMarkerSize(0.7)
    grRat.SetLineColor(ROOT.kBlack)
    grRat.GetHistogram().SetXTitle(xtitle)
    grRat.GetHistogram().SetYTitle('Point/Fit')
    grRat.GetHistogram().SetTitleSize(0.085, 'XY')
    grRat.GetHistogram().SetLabelSize(0.09, 'XY')
    grRat.GetHistogram().SetTitleOffset(0.3, 'Y')
    grRat.GetHistogram().SetAxisRange(500, 3700, 'X')
    grRat.GetHistogram().SetAxisRange(0, 2, 'Y')
    grRat.GetYaxis().SetNdivisions(5, 5, 1)
    grRat.Draw('APZ')

    # fit with a constant
    fit = ROOT.TF1('fit', '[0]', 100, 3000)
    fit.SetParameter(0, 1)
    fit.SetLineWidth(1)
    grRat.Fit(fit, 'QEM', 'same', 100, 3000)

    # legend
    leg = ROOT.TLegend(0.48, 0.25, 0.97, 0.36)
    saves.append(leg)
    fmt = 'Average in [500, 3500]: {0:.3f} #pm {1:.3f}'
    leg.AddEntry(fit, fmt.format(fit.GetParameter(0), fit.GetParError(0)), 'l')
    leg.SetFillColor(0)
    leg.Draw('same')


if __name__ == '__main__':
    main()
