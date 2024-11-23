#! \usr\bin\perl
## Copyright (C) 2020 Novelbio
## Running Platform: Ubuntu / Novelbrain
## Updator: LI, Yadang          Update time: 2020.03.16
## Author:  LI, Yadang          Create time: 2020.03.16
## Instruction: The function of this script is to check if input GeneSets files contain any illegal characters including quotes or apostrophes

use 5.014;
use strict;
use warnings;
use Cwd;

&main;
sub main {
    say "Running QuSAGE_GeneSetsFile_check.pl";
    &get_args;
    &check_Gset;
}
##  Func: Get arguments from outside of this program
sub get_args {
    my  %para = @ARGV;
    our $file_Gsets = $para{-t}; ## GeneSets Files
    our @files = split ",", $file_Gsets if exists $para{-t};
}

##  Func: Check illegal characters
sub check_Gset {
    our @files;
    my  ( $f, $filename, @info, @Gset, @illegal, $ErrInfo, $flag_err );
    $ErrInfo = "!!! Find illegal characters including quotes or apostrophes.\n";
    for $f ( @files ) {
      @info = &read_file ( $f );
      push @Gset, (split "\t", $_)[0] for @info[1..$#info];
      @Gset = &rm_dupl ( @Gset );
      (/['"]/) ? (push @illegal, "      $_") : () for @Gset;
      if ( @illegal > 0 ) {
        $flag_err++;
        $filename = (split '/', $f)[-1];
        $ErrInfo .= join "\n",
                         "    Please Check the following GeneSets Name:", @illegal,
                         "    in input GeneSetsFile: $filename\n";
        undef @illegal;
      }
    }
    die $ErrInfo if $flag_err > 0 ;
}


##  ---TOOLS------------------------------------------------------------TOOLS---

sub read_file {
    my ( $file, @read ) = $_[0];
    open FH, $file or die "Couldn't open $file for reading: $!";
      chomp for @read = <FH>;
    close FH;
    @read;
}

sub rm_dupl { ## remove duplicate of input array
    my ( @dupl, %count, @uniq ) = @_;
    @uniq = grep { ++$count{$_} < 2 } @dupl;
}

##  ---END----------------------------------------------------------------END---
