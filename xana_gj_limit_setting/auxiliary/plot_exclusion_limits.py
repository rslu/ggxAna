#!/usr/bin/env python
"""Draws final limit setting plot.
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
    # configuration
    path = '../output/higgsCombine_gj_cat0_8TeV.Asymptotic.merged.root'  #'../output/higgsCombine_eemmg_cat1234_8TeV.Asymptotic.merged.root'
    message = 'gj channels'

    limits = get_limits(path)

    grObs = ROOT.TGraph()
    grAv  = ROOT.TGraph()
    gr68  = ROOT.TGraphAsymmErrors()
    gr95  = ROOT.TGraphAsymmErrors()
    grXS = ROOT.TGraph()
    saves.append((grObs, grAv, gr68, gr95, grXS))

    for (i, (mH, obs, lo95, lo68, mean, hi68, hi95, XS)) in enumerate(limits):
        print('mH={0:.1f}, obs={1:6.2f}, median={2:6.2f}, SM XS={3:6.2f}'.format(mH, obs, mean, XS))

        # average & observed
        grObs.SetPoint(i, mH, obs)
        grAv.SetPoint(i, mH, mean)

        # 68
        mean68 = 0.5 * (lo68 + hi68)
        gr68.SetPoint(i, mH, mean68)
        gr68.SetPointError(i, 0, 0, mean68 - lo68, hi68 - mean68)

        # 95
        mean95 = 0.5 * (lo95 + hi95)
        gr95.SetPoint(i, mH, mean95)
        gr95.SetPointError(i, 0, 0, mean95 - lo95, hi95 - mean95)

        #XS
        grXS.SetPoint(i, mH, XS)
        

    c = ROOT.TCanvas('CL', 'CL', 800, 500)
    saves.append(c)

    c.SetLeftMargin(0.07)
    c.SetRightMargin(0.015)
    c.SetTopMargin(0.04)
    c.SetBottomMargin(0.08)
    c.SetGrid()
    c.SetLogy()
    
    # find out Y axis range
    #ymax = max(max(lim[1:]) for lim in limits)

    #frame = c.DrawFrame(500, 0.01, 4500, ymax * 1.1)
    frame = c.DrawFrame(500, 0.1, 4000, 30000.)    

    #frame.SetXTitle('M_{H} (GeV/c^{2})')
    #frame.SetYTitle('[#sigma#times BR]_{95%CL}/[#sigma#times BR]_{SM}')
    frame.SetXTitle('M_{#gammaj} (GeV/c^{2})')
    frame.SetYTitle('#sigma (fb)')
    frame.SetTitleOffset(0.9, 'Y')

    gr95.SetFillColor(ROOT.kGreen)
    gr95.Draw('e3')

    gr68.SetFillColor(ROOT.kYellow)
    gr68.Draw('e3')

    grAv.SetLineStyle(2)
    grAv.SetLineWidth(2)
    grAv.SetLineColor(ROOT.kRed)
    grAv.SetMarkerStyle(8)
    grAv.SetMarkerSize(1)
    grAv.Draw('PC')

    grObs.SetLineWidth(2)
    grObs.SetLineColor(ROOT.kBlack)
    grObs.Draw('PC,same')

    grXS.SetLineWidth(2)
    grXS.SetLineColor(ROOT.kBlue)
    grXS.Draw('PC,same')

    # put the grid on top of all graphs
    c.RedrawAxis('g')

    # legend
    leg = ROOT.TLegend(0.65, 0.73, 0.9, 0.93)
    saves.append(leg)

    leg.AddEntry(grObs, 'Observed', 'l')
    leg.AddEntry(grAv, 'Median expected', 'l')
    leg.AddEntry(gr68, 'Expected #pm 1#sigma', 'f')
    leg.AddEntry(gr95, 'Expected #pm 2#sigma', 'f')
    leg.AddEntry(grXS, 'SM q*', 'l')

    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(0)
    leg.Draw()

    # additional information
    txt = ROOT.TPaveText(0.18, 0.77, 0.48, 0.96, 'NDC')
    saves.append(txt)

    txt.AddText('q* #rightarrow #gamma jet')
    txt.AddText(message)
    txt.AddText('#sqrt{s} = 8 TeV, L = 2. fb^{-1}')

    txt.SetTextAlign(12)
    txt.SetTextFont(42)
    txt.SetTextSize(0.045)
    txt.SetFillStyle(0)
    txt.SetBorderSize(0)
    txt.Draw('same')

    c.Update()
    #c.Print('CL.png')

def get_limits(path):
    """Opens combine's output TTree, collects and returns confidence levels.
    """
    # {mass: value} dictionaries
    limObs  = {}
    lim95lo = {}
    lim68lo = {}
    limMean = {}
    lim68hi = {}
    lim95hi = {}
    XS = {}

    #XS[700] = 465.7
    #XS[1000] = 70.05
    #XS[1200] = 23.94
    #XS[1500] = 5.645
    #XS[1700] = 2.317
    #XS[2000] = 0.6531
    #XS[2500] = 0.08782
    #XS[3000] = 0.01234
    #XS[3500] = 0.001698

    XS[700] = 24850.
    XS[1000] = 4180.
    XS[1200] = 1552.
    XS[1500] = 412.4
    XS[1700] = 184.3
    XS[2000] = 58.581
    XS[2500] = 9.768
    XS[3000] = 1.755
    XS[3500] = 0.3241


    f = ROOT.TFile(path)
    tree = f.Get('limit')

    for entry in range(tree.GetEntriesFast()):
        if tree.GetEntry(entry) <= 0:
            raise Exception

        mH = tree.mh
        quant = tree.quantileExpected

        if -1.001 < quant < -0.999:
            limObs[mH] = tree.limit * XS[mH]
        elif 0.024 < quant < 0.026:
            lim95lo[mH] = tree.limit * XS[mH]
        elif 0.15 < quant < 0.17:
            lim68lo[mH] = tree.limit * XS[mH]
        elif 0.49 < quant < 0.51:
            limMean[mH] = tree.limit * XS[mH]
        elif 0.83 < quant < 0.85:
            lim68hi[mH] = tree.limit * XS[mH]
        elif 0.974 < quant < 0.976:
            lim95hi[mH] = tree.limit * XS[mH]
        else:
            raise Exception

    # simple test for duplicated entries
    if len(limObs.keys()) != tree.GetEntriesFast() // 6:
        raise Exception

    result = []

    # repack dictionaries into more robust representation
    for m in sorted(limObs.keys()):
        # TODO: is this possible?
        if lim95lo[m] > lim95hi[m] or lim68lo[m] > lim68hi[m]:
            raise Exception

        result.append((m, limObs[m],
                       min(lim95lo[m], lim95hi[m]),
                       min(lim68lo[m], lim68hi[m]),
                       limMean[m],
                       max(lim68lo[m], lim68hi[m]),
                       max(lim95lo[m], lim95hi[m]), 
                       XS[m]
                       ),                     
                      )

    return result


if __name__ == '__main__':
    main()
