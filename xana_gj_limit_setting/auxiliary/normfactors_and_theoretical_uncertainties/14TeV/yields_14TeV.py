#!/usr/bin/env python
"""Prints out expected yields for the SM Higgs boson at 14TeV.
"""

# python-2 support
from __future__ import division        # 1/2 = 0.5, not 0
from __future__ import print_function  # print() syntax from python-3

def main():
    """Steering function.
    """
    # SM Higgs production processes
    procs = ['ggH', 'qqH', 'WH', 'ZH', 'ttH', 'bbH']

    # collect cross-sections, branching ratios, uncertainties into one place
    A = parse_files(procs)

    # print out to the standard output numbers of expected H=>Zg=>eeg events
    expected_yields(A, lumi=19.7, masses=[125, 125.5, 126])

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
        with open('SM_Higgs_14TeV_{0}.txt'.format(proc), 'r') as f:
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

def expected_yields(A, lumi, masses):
    """Prints out expected numbers of H=>Zg=>eeg events for given luminosity.
    """
    fmt = 'Expected numbers of H=>Zg=>eeg events for sqrt(s)=14TeV and lumi={0}/fb:'
    print(fmt.format(lumi))
    print('       ' + ' '.join('{0:6.1f}'.format(x) for x in masses))
    for proc in A['procs']:
        print('  {0:4s}'.format(proc + ':'), end='')

        for mass in masses:
            m = '{0:.1f}'.format(mass)
            nev_exp = A['xs'][(proc, m)] * lumi * 1000 * A['br'][m] * 0.033658

            print(' {0:6.1f}'.format(nev_exp), end='')

        print()


if __name__ == '__main__':
    main()
