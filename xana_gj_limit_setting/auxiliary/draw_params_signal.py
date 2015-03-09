#!/usr/bin/env python
"""Visualization of parameters of signal fitting functions.
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
    combine_params('../output/for_visualizations.root', 'eeg', cat=4)

def combine_params(path, channel, cat):
    """Combines parameters for given channel and category.
    """
    # names of parameters
    parnames = (['nent', 'mean1', 'sigma1', 'delta21', 's21', 'delta31', 's32',
                 'frac23', 'frac123', 'chi2ndf'])

    f = ROOT.TFile(path)

    # loop over parameters; the parameter number -1 is chi2/ndf
    for par in range(-1, 9):
        cname = '{0}_cat{1}_{2}'.format(channel, cat, parnames[par])
        canv = ROOT.TCanvas(cname, cname, 700, 700)
        saves.append(canv)

        canv.SetLeftMargin(0.10)
        canv.SetRightMargin(0.02)
        canv.SetTopMargin(0.06)
        canv.SetBottomMargin(0.08)
        canv.SetGridy()

        leg = ROOT.TLegend(0.85, 0.75, 0.96, 0.92)
        saves.append(leg)

        gr = [ROOT.TGraphErrors() for c in range(5)]
        saves.append(gr)

        clrs = [ROOT.kBlack, ROOT.kBlue, 8, ROOT.kRed, ROOT.kOrange, ROOT.kYellow]
        stls = [20, 20, 24, 24, 25, 25]

        for (c, proc) in enumerate(['ggH', 'qqH', 'WH', 'ZH', 'ttH']):
            i = 0 # counter of mass points

            # fill graphs
            for mass in range(120, 161, 1):
                fmt = '{0}_cat{1}_8TeV/fit_hMass_{2}_{3}'
                fit = f.Get(fmt.format(channel, cat, proc, mass))
                if not fit:
                    continue

                # chi2/ndf
                if par == -1:
                    gr[c].SetPoint(i, mass, fit.GetChisquare()/fit.GetNDF())
                    gr[c].SetPointError(i, 0, 0)
                    i += 1

                # fit parameters
                elif par < fit.GetNpar():
                    gr[c].SetPoint(i, mass, fit.GetParameter(par))
                    gr[c].SetPointError(i, 0, fit.GetParError(par))
                    i += 1

            gr[c].SetMarkerStyle(stls[c])
            gr[c].SetMarkerSize(0.8)
            gr[c].SetMarkerColor(clrs[c])
            gr[c].SetLineColor(clrs[c])

            title = '{0}_cat{1}, {2}'.format(channel, cat, parnames[par])
            gr[c].GetHistogram().SetTitle(title)
            gr[c].GetHistogram().SetXTitle('M_{ll#gamma}^{rec} (GeV/c^{2})')
            gr[c].GetHistogram().SetYTitle('Value')
            gr[c].GetHistogram().SetTitleOffset(1.4, 'Y')

            gr[c].Draw('APZC' if c == 0 else 'PZC,same')

            leg.AddEntry(gr[c], proc, 'p')

        # set appropriate Y axis range
        ymin = min(gr[c].GetY()[i] - gr[c].GetEY()[i] for i in range(gr[c].GetN()) for c in range(5))
        ymax = max(gr[c].GetY()[i] + gr[c].GetEY()[i] for i in range(gr[c].GetN()) for c in range(5))
        dy = ymax - ymin
        gr[0].GetHistogram().SetAxisRange(ymin - 0.1 * dy, ymax + 0.1 * dy, 'Y')

        leg.SetFillColor(ROOT.kWhite)
        leg.Draw('same')


if __name__ == '__main__':
    main()
