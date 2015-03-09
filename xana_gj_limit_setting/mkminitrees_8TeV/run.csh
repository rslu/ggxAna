#!/bin/csh -f

source ~/setup_root.csh

echo "process data"
#root.exe -b -q run.C
#cp ana_gj_out_FNAME_MASK.root ana_gj_out_data.root

echo "process signal for pho"

./make_mc_code.pl --mc  job_qstar_700.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_1000_half.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_1000.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_1200.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_1500_half.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_1500.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_1700.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_2000_half.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_2000.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_2500_half.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_2500.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_3000.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qstar_3500.root
root.exe -b -q run_mc.C


echo "process pho+jet"
./make_mc_code.pl --mc  job_gjet_80to120.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_gjet_120to170.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_gjet_170to300.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_gjet_300to470.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_gjet_470to800.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_gjet_800to1000.root
root.exe -b -q run_mc.C

echo "process qcd "
./make_mc_code.pl --mc  job_qcd_80to120.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qcd_120to170.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qcd_170to300.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qcd_300to470.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qcd_470to600.root
root.exe -b -q run_mc.C
./make_mc_code.pl --mc  job_qcd_600to800.root
root.exe -b -q run_mc.C


hadd -f ana_gj_out_mc.root ana_gj_out_job_gjet*.root ana_gj_out_job_qcd_*.root 

