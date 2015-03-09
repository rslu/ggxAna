#!/usr/bin/env python
"""Visualization of two-lepton + photon invariant masses.
"""

# python-2 compatibility
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

import ROOT

# for keeping drawed objects in memory
saves = []

# normalization factors which adapt scales of the Zg and DYJets MCs to real data
normZg = 156.2 * 19693 / 6588161     # xs = 156.2 pb, lumi = 19.7 fb^{-1}
normDY = 3503.71 * 19693 / 30458871  # xs = 3503.71 pb

def main():
    """Steering function.

    TODO: in the ratio (bottom pad), error bars are not quite correct since
    non-unit event weights of MCs are not propagated accurately.
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(0)

    # disable RooHist warnings 'non-integer bin entry with Poisson errors'
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

    for cat in [1, 2, 3, 4]:
        combine('../output/minitrees.root', 'eeg', cat, 45)
        combine('../output/minitrees.root', 'mmg', cat, 45)

    combine('../output/minitrees.root', 'eeg', -1, 90)
    combine('../output/minitrees.root', 'mmg', -1, 90)

def combine(path, channel, cat=-1, nbins=90):
    """Draws the combination plot for particular channel and category.
    """
    hData = ROOT.TH1D('h', '', nbins, 100, 190)
    hZgMC = ROOT.TH1D('h', '', nbins, 100, 190)
    hDYZg = ROOT.TH1D('h', '', nbins, 100, 190)

    base = 'minitree_{0}_'.format(channel)

    # data 2012, Zg MC
    fill(hData, path, base + 'data', cat)
    fill(hZgMC, path, base + 'Zg_MC', cat, True, normZg)

    # Zg MC + DYJets MC with ISR/FSR photons removed
    fill(hDYZg, path, base + 'Zg_MC', cat, True, normZg)
    fill(hDYZg, path, base + 'DYJets_MC_withoutISRFSRPhotons', cat, True, normDY)

    cname = 'backgrounds_{0}_cat{1}'.format(channel, cat if cat > 0 else 'All')
    c = ROOT.TCanvas(cname, cname, 700, 700)
    saves.append(c)

    pad1 = ROOT.TPad('top', 'top', 0, 0.3, 1, 1)
    pad2 = ROOT.TPad('bottom', 'bottom', 0, 0, 1, 0.3)
    pad1.Draw()
    pad2.Draw()
    saves.append((pad1, pad2))

    # top pad
    pad1.cd()
    pad1.SetGridx()
    pad1.SetGridy()
    pad1.SetTopMargin(0.07)
    pad1.SetBottomMargin(0)
    pad1.SetLeftMargin(0.08)
    pad1.SetRightMargin(0.02)

    title = 'electrons, ' if channel == 'eeg' else 'muons, '
    title += 'cat{0}'.format(cat) if cat > 0 else 'no categorization'

    xtitle = 'M_{ee#gamma}^{rec} (GeV/c^{2})'
    if channel == 'mmg':
        xtitle = 'M_{#mu#mu#gamma}^{rec} (GeV/c^{2})'

    # draw the TGraphAsymmErrors object first in order to make its X axis range
    # exactly the same as used in the bottom pad below
    grData = ROOT.RooHist(hData) # make errors asymmetric
    grData.SetLineColor(ROOT.kBlack)
    grData.SetMarkerStyle(20)
    grData.SetMarkerSize(0.7)
    grData.GetHistogram().SetTitle(title)
    grData.GetHistogram().SetXTitle(xtitle)
    grData.GetHistogram().SetYTitle('Entries')
    grData.GetHistogram().SetTitleSize(0.037, 'XY')
    grData.GetHistogram().SetLabelSize(0.037, 'XY')
    grData.GetHistogram().SetTitleOffset(1.1, 'Y')
    grData.GetHistogram().SetAxisRange(90, 190, 'X')
    grData.GetHistogram().SetAxisRange(0, hData.GetMaximum() * 1.15, 'Y')
    grData.Draw('APZ')

    hDYZg.SetLineColor(ROOT.kGreen + 2)
    hDYZg.SetFillColor(ROOT.kGreen + 2)
    hDYZg.SetFillStyle(1001)
    hDYZg.Draw('same,hist')

    hZgMC.SetLineColor(ROOT.kOrange + 1)
    hZgMC.SetFillColor(ROOT.kOrange + 1)
    hZgMC.SetFillStyle(1001)
    hZgMC.Draw('same,hist')

    # redraw data points in order to put them on top of other objects
    grData.Draw('PZ')

    saves.append((grData, hDYZg, hZgMC))

    # legend
    leg = ROOT.TLegend(0.65, 0.77, 0.98, 0.93)
    saves.append(leg)
    leg.AddEntry(grData, 'data 2012, 19.7 fb^{-1}', 'p')
    leg.AddEntry(hDYZg, 'Z + jets MC', 'f')
    leg.AddEntry(hZgMC, 'Z#gamma MC', 'f')
    leg.SetFillColor(0)
    leg.Draw('same')

    # bottom pad
    pad2.cd()
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.SetLeftMargin(0.08)
    pad2.SetRightMargin(0.02)
    pad2.SetGridx()
    pad2.SetGridy()

    # grData/hDYZg ratio with asymmetric errors
    gr = ROOT.TGraphAsymmErrors()
    saves.append(gr)

    # fill gr
    grDnmr = ROOT.RooHist(hDYZg)
    for i in range(grData.GetN()):
        (x, exl, exh)    = (grData.GetX()[i], grData.GetEXlow()[i], grData.GetEXhigh()[i])
        (y1, eyl1, eyh1) = (grData.GetY()[i], grData.GetEYlow()[i], grData.GetEYhigh()[i])
        (y2, eyl2, eyh2) = (grDnmr.GetY()[i], grDnmr.GetEYlow()[i], grDnmr.GetEYhigh()[i])

        if y1 == 0 or y2 == 0:
            gr.SetPoint(i, x, 0)
            gr.SetPointError(i, exl, exh, 0, 0)
            continue

        eyl = y1/y2 * ((eyl1/y1)**2 + (eyl2/y2)**2)**0.5
        eyh = y1/y2 * ((eyh1/y1)**2 + (eyh2/y2)**2)**0.5

        gr.SetPoint(i, x, y1/y2)
        gr.SetPointError(i, exl, exh, eyl, eyh)

    gr.SetMarkerStyle(20)
    gr.SetMarkerSize(0.7)
    gr.SetLineColor(ROOT.kBlack)
    gr.GetHistogram().SetXTitle(xtitle)
    gr.GetHistogram().SetYTitle('Data/Z+jets')
    gr.GetHistogram().SetTitleSize(0.085, 'XY')
    gr.GetHistogram().SetLabelSize(0.09, 'XY')
    gr.GetHistogram().SetTitleOffset(0.45, 'Y')
    gr.GetHistogram().SetAxisRange(90, 190, 'X')
    gr.GetHistogram().SetAxisRange(0, 2, 'Y')
    gr.Draw('APZ')

    # fit with a constant
    fit = ROOT.TF1('fit', '[0]', 100, 190)
    fit.SetParameter(0, 1)
    fit.SetLineWidth(1)
    gr.Fit(fit, 'QEM', 'same', 100, 190)

    # legend
    leg = ROOT.TLegend(0.48, 0.25, 0.97, 0.36)
    saves.append(leg)
    fmt = 'Average in [100, 190]: {0:.3f} #pm {1:.3f}'
    leg.AddEntry(fit, fmt.format(fit.GetParameter(0), fit.GetParError(0)), 'l')
    leg.SetFillColor(0)
    leg.Draw('same')

    #c.Print(cname + '.png')

def fill(histo, path, treename, cat=-1, isMC=False, norm=1):
    """Fills histogram with invariant mass points from a minitree.
    """
    f = ROOT.TFile(path)
    tree = f.Get(treename)

    for ev in range(tree.GetEntriesFast()):
        if tree.GetEntry(ev) <= 0:
            raise Exception

        if cat > 0 and tree.category != cat:
            continue

        mcwei = tree.mcwei if isMC else 1
        histo.Fill(tree.hzg_mass, norm * mcwei)


if __name__ == '__main__':
    main()
