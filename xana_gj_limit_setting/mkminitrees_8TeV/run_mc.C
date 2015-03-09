void run_mc(){
  gROOT->LoadMacro("Ana_gj_mc.C++");
  Ana_gj_mc t;
  t.Loop();
}
