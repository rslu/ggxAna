#!/usr/bin/env python
"""Manipulates with SM Higgs cross-sections, branching ratios and uncertainties.

Writes ready-to-be-used (in the HZg analysis datacards) lines with theoretical
uncertainties into text file for all Higgs mass values.

Extracts numbers of events from signal MCs, then evaluates and prints out to the
standard output normalization factors which adapt MC scales to real data.
"""

# python-2 support
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

import sys

def main():
    """Steering function.
    """
    # SM Higgs production processes
    procs = ['ggH', 'qqH', 'WH', 'ZH', 'ttH']

    # collect cross-sections, branching ratios, uncertainties into one place
    A = parse_files(procs)

    # add (with the aid of linear extrapolation) non-existent mass points in the
    # range 120-160GeV with step 1GeV
    masses_orig = A['masses'][:]
    for mass in range(120, 161):
       if '{0:.1f}'.format(mass) not in masses_orig:
         inject_mass(A, mass, masses_orig)

    # print out to the standard output numbers of expected H=>Zg=>eeg events
    expected_yields(A, lumi=19.7, masses=range(100, 195, 5))

    # write specially formatted lines with theoretical uncertainties into file
    write_uncertainties(A, 'theoretical_uncertainties_SM_Higgs.txt')

    # collected luminosities in 2012 (the sums over total recorded luminosities
    # extracted from /data3/ggNtuples/V05-03-12-03/job_2*_2012*.lumi files)
    lumi = {'eeg': 19712.2, 'mmg': 19674.2} # pb^(-1)

    # evaluate and print out normalization factors for signal MCs
    # NOTE: the code below requires access to the signal ntuples
    eval_normfactors(A, lumi, masses=range(120, 161))

def parse_files(procs):
    """Returns values of cross-sections, branching ratios, uncertainties.
    """
    xsections = {}   # {(proc, mass): cross-section}
    qcd_unc = {}     # {(proc, mass): (unc_down, unc_up)}
    pdf_unc = {}     # {(proc, mass): (unc_down, unc_up)}
    branchings = {}  # {mass: branching}
    br_unc = {}      # {mass: (unc_down, unc_up)}

    # parse text files with cross-sections and their uncertainties
    for proc in procs:
        with open('SM_Higgs_8TeV_{0}.txt'.format(proc), 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue

                (mass, xs, scup, scdn, pdfup, pdfdn) = line.split()
                ind = (proc, mass)

                xsections[ind] = float(xs)
                qcd_unc[ind] = (1 + 0.01 * float(scdn), 1 + 0.01 * float(scup))
                pdf_unc[ind] = (1 + 0.01 * float(pdfdn), 1 + 0.01 * float(pdfup))

    # parse text file with branching ratios and their uncertainties
    with open('SM_Higgs_branching_ratios.txt', 'r') as f:
        for line in f.readlines():
            if line.startswith('#'):
                continue

            # take branching ratios for the HZg analysis
            (mass, _, _, _, _, _, _, br, uncup, uncdn) = line.split()

            branchings[mass] = float(br)
            br_unc[mass] = (1 + 0.01 * float(uncdn), 1 + 0.01 * float(uncup))

    # sorted list with mass values
    masses = sorted(branchings.keys(), key=lambda x: float(x))

    return {'masses': masses, 'procs': procs, 'xs': xsections, 'br': branchings,
            'qcd_unc': qcd_unc, 'pdf_unc': pdf_unc, 'br_unc': br_unc}

def inject_mass(A, mass, masses_orig):
    """Inserts new mass point into A.

    Numbers are linearly extrapolated from neighbours corresponding to existing
    MC productions.

    mass = mass value to insert;
    masses_orig = A['masses'] just after parse_files() returned.
    """
    # find neighbours
    for (m1, m2) in zip(masses_orig, masses_orig[1:]):
        if float(m1) < mass < float(m2):
            break
    else:
        raise Exception

    # text form of mass value to insert
    m0 = '{0:.1f}'.format(mass)
    if m0 in A['masses']:
        raise Exception

    # proportions of first/second neighbour
    a = (float(m2) - mass)/(float(m2) - float(m1))
    b = 1 - a

    fmt = 'WARNING: injecting {0} as average of {1}({2:.3f}) and {3}({4:.3f})'
    print(fmt.format(m0, m1, a, m2, b), file=sys.stderr)

    for proc in A['procs']:
        i0 = (proc, m0)
        i1 = (proc, m1)
        i2 = (proc, m2)

        A['xs'][i0] = a * A['xs'][i1] + b * A['xs'][i2]
        A['br'][m0] = a * A['br'][m1] + b * A['br'][m2]

        A['qcd_unc'][i0] = (a * A['qcd_unc'][i1][0] + b * A['qcd_unc'][i2][0],
                            a * A['qcd_unc'][i1][1] + b * A['qcd_unc'][i2][1])

        A['pdf_unc'][i0] = (a * A['pdf_unc'][i1][0] + b * A['pdf_unc'][i2][0],
                            a * A['pdf_unc'][i1][1] + b * A['pdf_unc'][i2][1])

        A['br_unc'][m0] = (a * A['br_unc'][m1][0] + b * A['br_unc'][m2][0],
                           a * A['br_unc'][m1][1] + b * A['br_unc'][m2][1])

    A['masses'] = sorted(A['masses'] + [m0], key=lambda x: float(x))

def expected_yields(A, lumi, masses):
    """Prints out expected numbers of H=>Zg=>eeg events for given luminosity.
    """
    print('Expected numbers of H=>Zg=>eeg events for lumi={0}/fb:'.format(lumi))
    print('       ' + ' '.join('{0:6.0f}'.format(x) for x in masses))
    for proc in A['procs']:
        print('  {0:4s}'.format(proc + ':'), end='')

        for mass in masses:
            m = '{0:.1f}'.format(mass)
            nev_exp = A['xs'][(proc, m)] * lumi * 1000 * A['br'][m] * 0.033658

            print(' {0:6.1f}'.format(nev_exp), end='')

        print()

def write_uncertainties(A, outpath):
    """Writes specially formatted lines with theoretical uncertainties into file.
    """
    # the WH/ZH/ttH Higgs production channels are covered only up to 400 GeV
    masses = [x for x in A['masses'] if float(x) < 400.01]

    # make and collect lines with theoretical uncertainties
    res = '# Created by normfactors_and_theoretical_uncertainties.py\n'
    for mass in masses:
        fmt = ('{0:7} :'.format(mass) +
               '{0:<17}  lnN   {1:^15}{2:^15}{3:^15}{4:^15}{5:^15}{6:^15}\n')
        z = '-'

        # PDF+alphaS uncertainties
        pair = lambda x: '{0:.3f}/{1:.3f}'.format(*A['pdf_unc'][(x, mass)])

        res += fmt.format('pdf_gg', pair('ttH'), z, z, z, pair('ggH'), z)
        res += fmt.format('pdf_qqbar', z, pair('ZH'), pair('WH'), pair('qqH'), z, z)

        # QCD scale uncertainties
        pair = lambda x: '{0:.3f}/{1:.3f}'.format(*A['qcd_unc'][(x, mass)])

        res += fmt.format('QCDscale_ggH', z, z, z, z, pair('ggH'), z)
        res += fmt.format('QCDscale_qqH', z, z, z, pair('qqH'), z, z)
        res += fmt.format('QCDscale_VH',  z, pair('ZH'), pair('WH'), z, z, z)
        res += fmt.format('QCDscale_ttH', pair('ttH'), z, z, z, z, z)

        # branching ratios uncertainties
        unc = '{0:.4f}/{1:.4f}'.format(*A['br_unc'][mass])
        res += fmt.format('br_higgs_to_Zg', unc, unc, unc, unc, unc, z)

    print('Writing {0} ...'.format(outpath))
    with open(outpath, 'w') as f:
        f.write(res)

def eval_normfactors(A, lumi, masses):
    """Evaluates and prints out normalization factors for signal MCs.
    """
    # Higgs masses for existing MC productions
    mass_existing = [120, 125, 130, 135, 140, 145, 150, 155, 160]

    # template for paths to the signal ntuples
    pathfmt = '/data5/ggNtuples/V05-03-12-05/job_summer12_HZg_{0}_{1}.root'

    # production process <=> ntuple suffix map
    sfx = {'ggH':'ggH', 'qqH':'VBFH', 'WH':'WH_ZH', 'ZH':'WH_ZH', 'ttH':'TTH'}

    nevMC = {}  # {(proc, mass): number of events}

    print('Counting numbers of events in HZg signal MCs:')

    # collect numbers of events in MCs
    for proc in A['procs']:
        for mass in mass_existing:
            ntuple = pathfmt.format(sfx[proc], mass)

            if proc == 'WH':
                nev = get_entries(ntuple, takeWHZHttH=1)
            elif proc == 'ZH':
                nev = get_entries(ntuple, takeWHZHttH=2)
            elif proc == 'ttH':
                nev = get_entries(ntuple, takeWHZHttH=3)
            else:
                nev = get_entries(ntuple)

            print('  {0}: nevents_{1}={2}'.format(ntuple, proc, nev))

            nevMC[(proc, mass)] = nev

    # insert linearly extrapolated values for all intermediate mass points
    for mass in masses:
        if mass not in mass_existing:
            # find closest neighbours
            for (m1, m2) in zip(mass_existing, mass_existing[1:]):
                if m1 < mass < m2:
                    break
            else:
                raise Exception

            # proportions of first/second neighbour
            a = (m2 - mass)/(m2 - m1)
            b = 1 - a

            for proc in A['procs']:
                nevMC[(proc, mass)] = a * nevMC[(proc, m1)] + b * nevMC[(proc, m2)]

    print('Process order: ' + ' '.join(A['procs']))
    print('Mass order: ' + ' '.join(str(x) for x in masses))

    for channel in ['eeg', 'mmg']:
        # {proc: list with per-mass numbers of events}
        norm = dict((proc, []) for proc in A['procs'])

        for proc in A['procs']:
            for mass in masses:
                m = '{0:.1f}'.format(mass)

                # expected number of signal events in real data
                nev_exp = A['xs'][(proc, m)] * lumi[channel] * A['br'][m]

                # account for Z=>ee + Z=>mm + Z=>tt branching ratio
                #if proc not in ['WH', 'ZH']: # TODO: bug in WH_ZH samples
                nev_exp *= 0.033658 * 3

                norm[proc].append('{0:.5e}'.format(nev_exp/nevMC[(proc, mass)]))

        # print out the normalization factors
        arr = []
        for proc in A['procs']:
            arr.append('{' + ', '.join(norm[proc]) + '}')
        print('double norm_{0}[{1}][{2}] = {{\n   '.format(
                       channel, len(A['procs']), len(masses)) +
              ', \n   '.join(arr) + '};')

def get_entries(path, takeWHZHttH=-1):
    """Returns number of entries in ggNtuplizer's ntuple.

    path = path to root file with ggNtuplizer's TTree;
    takeWHZHttH =
        1 = count WH events only;
        2 = count ZH events only;
        3 = count events where Z decays into a pair of leptons (ee, mm or tt);
        any other value = count all events.
    """
    import ROOT

    f = ROOT.TFile(path)
    tree = f.Get('ggNtuplizer/EventTree')

    if takeWHZHttH not in [1, 2, 3]:
        return tree.GetEntriesFast()

    nent = 0

    # treatment of WH_ZH and ttH signal samples
    for ev in range(tree.GetEntriesFast()):
        if (tree.GetEntry(ev) <= 0):
            raise Exception('TTree::GetEntry() failed')

        mcPID = tree.mcPID

        # ttH
        if takeWHZHttH == 3:
            mcMomPID = tree.mcMomPID

            nlep = 0
            for k in range(tree.nMC):
                if mcMomPID[k] == 23 and abs(mcPID[k]) in [11, 13, 15]:
                    nlep += 1

            if (nlep >= 2):
                nent += 1

        # WH_ZH
        else:
            # count only WH events vs count only ZH events
            for k in range(tree.nMC):
                if abs(mcPID[k]) == 24:
                    if takeWHZHttH == 1:
                        nent += 1
                    break
            else:
                if takeWHZHttH == 2:
                    nent += 1

    return nent


if __name__ == '__main__':
    main()
