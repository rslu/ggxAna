
void rootlogon()
{
   /* This function is automatically executed by root on start if no customized
    * rootlogon configuration is set in .rootrc.
    *
    * Useful for keeping the current working directory clean from .d and .so
    * files.
    */

   gSystem->SetBuildDir("tmpdir", kTRUE);
}
