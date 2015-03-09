#!/bin/bash
#
# This script obtains limit setting plots for the HZg analysis in one shot and
# serves as a source of documentation.
#
# NOTE: there is a bug in python's interface to ROOT due to which canvases are
# drawn with empty contents. Just force such a canvas to be redrawn by resizing
# it (e.g. by maximizing/restoring the canvas size).

# stop on first error
set -e

# Process pre-minitrees and create output/minitrees.root.
#
# Previous version (if any) of the directory ./output will be removed.
#
# Actual configuration for the code is given in the mkminitrees() steering
# function at the end of mkminitrees.cc file.
#--disabled by rslu Aug14.
#root -b -l -q mkminitrees.cc+

# Draw comparison of data/MC 3-particle invariant mass shapes.
#
# cd auxiliary
# python -i combine_backgrounds.py
# cd ..

# 1. Process text files with cross-sections, branching ratios and uncertainties
# for the SM Higgs and produce theoretical_uncertainties_SM_Higgs.txt with lines
# ready-to-be-used in the HZg analysis datacards.
#
# 2. Evaluate and print out a few numbers of expected H=>Zg=>eeg yields.
#
# 3. Extract numbers of events from ggNtuplizer's ntuples with HZg signal
# samples, evaluate and print out normalization factors which adapt MC scales to
# real data. For this step to succeed, the signal ntuples must be accessible.
# The printed arrays with normalization factors must be added to mkdatacards.C
# by hand.
#
# Those three steps were already done, this is just for documentation.
#
# cd auxiliary/normfactors_and_theoretical_uncertainties
# python normfactors_and_theoretical_uncertainties.py
# mv theoretical_uncertainties_SM_Higgs.txt ../..
# cd ../..

# Print out expected H=>Zg=>eeg yields for sqrt(s) = 14TeV.
#
# cd auxiliary/normfactors_and_theoretical_uncertainties/14TeV
# python yields_14TeV.py
# cd ../../..

# Process minitrees from output/minitrees.root: perform unbinned fitting of
# invariant mass shapes, create text datacards and corresponding root files with
# PDFs and observed mass points.
#
# Output: output/datacards/* and output/for_visualizations.root. The latter
# contains objects for easy drawing.
#
# Actual configuration is encoded in the macro itself.

#root -b -l -q mkdatacards.C

# Draw many interesting plots.
#
# Actual configuration on what to draw can be changed in a steering function of
# every python script or root macro.
#
# cd auxiliary
# python -i draw_params_background.py
# python -i draw_fits.py
# python -i draw_params_signal.py
# root -l combine_signal_pdfs.cc+
# cd ..

cd output/datacards
# Produce merged datacards.
#for ((mass=120;mass<161;mass++)); do
#    # merged categories
#    for ch in eeg mmg; do
#        combineCards.py cat1=datacard_hzg_${ch}_cat1_8TeV_${mass}.txt \
#                        cat2=datacard_hzg_${ch}_cat2_8TeV_${mass}.txt \
#                        cat3=datacard_hzg_${ch}_cat3_8TeV_${mass}.txt \
#                        cat4=datacard_hzg_${ch}_cat4_8TeV_${mass}.txt \
#                    > datacard_hzg_${ch}_cat1234_8TeV_${mass}.txt
#    done
#
#    # merged channels and categories
#    combineCards.py eeg_cat1=datacard_hzg_eeg_cat1_8TeV_${mass}.txt \
#                    eeg_cat2=datacard_hzg_eeg_cat2_8TeV_${mass}.txt \
#                    eeg_cat3=datacard_hzg_eeg_cat3_8TeV_${mass}.txt \
#                    eeg_cat4=datacard_hzg_eeg_cat4_8TeV_${mass}.txt \
#                    mmg_cat1=datacard_hzg_mmg_cat1_8TeV_${mass}.txt \
#                    mmg_cat2=datacard_hzg_mmg_cat2_8TeV_${mass}.txt \
#                    mmg_cat3=datacard_hzg_mmg_cat3_8TeV_${mass}.txt \
#                    mmg_cat4=datacard_hzg_mmg_cat4_8TeV_${mass}.txt \
#                > datacard_hzg_eemmg_cat1234_8TeV_${mass}.txt
#done

cd ..

#
# Run the combine tool from the HiggsAnalysis/CombinedLimit package.
#

# maximum number of jobs to run in parallel
nparallel=1

echo "Spawning per-datacard jobs, waiting for their termination:"

#for cat in 0 ; do
for cat in 0; do
    for ch in gj; do
        # do not iterate over categories if ch=eemmg
        [ "$ch" == "eemmg" ] && [ "$cat" != "1234" ] && continue

        outdir="results_gj_${ch}_cat${cat}_8TeV"
        mkdir $outdir

        #for ((mass=2000;mass<2001;mass++)); do
        for mass in 700 1000 1200 1500 1700 2000 2500 3000 3500; do
            card="datacard_gj_${ch}_cat${cat}_8TeV_${mass}.txt"

            echo "Spawning job for $card (`date`) ..."
            cd $outdir

            log="combine_gj_${ch}_cat${cat}_8TeV_${mass}.log"

            # NOTE: "tee" forbids the combine tool to delete its own output
#            nice combine -v1 -U -M Asymptotic -m $mass -n "_${ch}_cat${cat}_8TeV" \
#                ../datacards/$card 2>&1 | tee -a ${log} >/dev/null &
            combine -v1 -U -M Asymptotic -m $mass -n "_${ch}_cat${cat}_8TeV" \
                ../datacards/$card 2>&1 | tee -a ${log} >/dev/null 

            cd ..

            # do not execute more than $nparallel jobs in parallel
            while [ "`jobs -p | wc -l`" -ge "$nparallel" ]; do
                sleep 2
            done
        done
    done
done

echo "All jobs spawned. Waiting for all jobs to terminate ..."
#wait

# Merge per-mass point results together.
for cat in 0; do
    for ch in gj; do
        # do not iterate over categories if ch=eemmg
        [ "$ch" == "eemmg" ] && [ "$cat" != "1234" ] && continue

        outdir="results_gj_${ch}_cat${cat}_8TeV"

        hadd -v 1 higgsCombine_${ch}_cat${cat}_8TeV.Asymptotic.merged.root \
             $outdir/higgsCombine_${ch}_cat${cat}_8TeV.Asymptotic.mH*.root
    done
done

echo "Finished (`date`)."

# Draw final limit setting plot.
# cd ../auxiliary
# python -i plot_exclusion_limits.py
