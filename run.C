void run(int option){
  gROOT->LoadMacro("xAna_C.so");
  xAna(option);
}
