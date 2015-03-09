#!/usr/bin/env python
"""Visualization of parameters of background fitting function.
"""

# python-2 compatibility
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

import ROOT

# for keeping drawed objects in memory
saves = []

def main():
    """Steering function.
    """
    f = ROOT.TFile('../output/for_visualizations.root')

    # find out the number of parameters in the fitting function
    npars = f.Get('eeg_cat1_8TeV/fit_hDataObs').GetNpar()
    if npars > 11:
        raise Exception('The canvas is of size 4x3 - 1')

    canv = ROOT.TCanvas('params_background', 'params_background', 1280, 720)
    saves.append(canv)
    canv.SetLeftMargin(0)
    canv.SetRightMargin(0)
    canv.SetTopMargin(0)
    canv.SetBottomMargin(0)
    canv.Divide(4, 3)

    # names of parameters
    parnames = (['nent', 'mean', 'sigma', 'stepval'] +
                ['coef{0}'.format(i + 1) for i in range(16)])

    for par in range(npars):
        canv.cd(par + 1)

        ROOT.gPad.SetLeftMargin(0.14)
        ROOT.gPad.SetRightMargin(0.01)
        ROOT.gPad.SetTopMargin(0.08)
        ROOT.gPad.SetBottomMargin(0.09)
        ROOT.gPad.SetGridy()

        gr = [ROOT.TGraphErrors() for c in range(2)]
        saves.append(gr)

        clrs = [ROOT.kBlack, ROOT.kBlue]

        # fill graphs with per-category points
        for (c, channel) in enumerate(['eeg', 'mmg']):
            for cat in [1, 2, 3, 4]:
                fit = f.Get('{0}_cat{1}_8TeV/fit_hDataObs'.format(channel, cat))

                gr[c].SetPoint(cat - 1, cat, fit.GetParameter(par))
                gr[c].SetPointError(cat - 1, 0, fit.GetParError(par))

            # bin labels
            gr[c].GetHistogram().SetBins(4, 0.5, 4.5)
            for cat in [1, 2, 3, 4]:
                gr[c].GetXaxis().SetBinLabel(cat, 'cat{0}'.format(cat))

            gr[c].SetMarkerStyle(20)
            gr[c].SetMarkerSize(0.8)
            gr[c].SetMarkerColor(clrs[c])
            gr[c].SetLineColor(clrs[c])
            gr[c].GetXaxis().SetLabelSize(0.09)
            gr[c].GetYaxis().SetLabelSize(0.05)
            gr[c].GetYaxis().SetTitleSize(0.05)

            gr[c].GetHistogram().SetTitle(parnames[par])
            gr[c].GetHistogram().SetYTitle('Value')
            gr[c].GetHistogram().SetTitleOffset(1.6, 'Y')

        # set appropriate Y axis range
        ymin = min(gr[c].GetY()[i] - gr[c].GetEY()[i] for i in range(4) for c in range(2))
        ymax = max(gr[c].GetY()[i] + gr[c].GetEY()[i] for i in range(4) for c in range(2))
        dy = ymax - ymin
        gr[0].GetHistogram().SetAxisRange(ymin - 0.1 * dy, ymax + 0.1 * dy, 'Y')

        gr[0].Draw('APZC')
        gr[1].Draw('PZC,same')

    canv.cd(12)
    leg = ROOT.TLegend(0.3, 0.2, 0.8, 0.7)
    saves.append(leg)

    leg.AddEntry(gr[0], 'Z#gamma #rightarrow e^{+}e^{-}#gamma', 'p')
    leg.AddEntry(gr[1], 'Z#gamma #rightarrow #mu^{+}#mu^{-}#gamma', 'p')

    leg.SetFillColor(ROOT.kWhite)
    leg.Draw('same')


if __name__ == '__main__':
    main()
