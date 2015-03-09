#!/usr/bin/env perl

use strict;
use DBI;
use Getopt::Long;


open (INF, "<./Ana_gj.C");
open (OUTF, ">./Ana_gj_mc.C");
while (<INF>) {
  if ( $_ =~ s/if\(HLT/\/\/if(HLT/ ) {
    print OUTF $_;
  }elsif ( $_ =~ s/\/\/$ARGV[0]_mask// ) {
    print OUTF $_;
  }elsif ( $_ =~ s/FNAME_MASK.root/$ARGV[1]/ ) {
    print OUTF $_;
  }elsif ( $_ =~ s/Ana_gj/Ana_gj_mc/ ){
    print OUTF $_;
  }elsif ( $_ =~ s/\/\/MC_$ARGV[1]/ / ){
    print OUTF $_;
  }elsif ( $_ =~ s/MCTYPE 1/MCTYPE 2/ ){
    print OUTF $_;
  }else {
    print OUTF $_;
  }
}
    
close (INF);
close (OUTF);
  
open (INF, "<./Ana_gj_mc_template.h");
open (OUTF, ">./Ana_gj_mc.h");
while (<INF>) {
  if ( $_ =~ s/\/\/$ARGV[0]_mask// ) {
    print OUTF $_;
  }elsif ( $_ =~ s/FNAME_MASK/$ARGV[1]/ ) {
    print OUTF $_;
  }else {
    print OUTF $_;
  }
}

close (INF);
close (OUTF);
