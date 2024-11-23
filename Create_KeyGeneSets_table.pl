#! usr\bin\perl

## Copyright (C) 2020 Novelbio
## Running Platform: Ubuntu / Novelbrain
## Updator: LI, Yadang          Update time: 2020.07.27
## Author:  LI, Yadang          Create time: 2020.07.27
## Instruction: The function of this script is to create a table of Top 3 GeneSets of each cluster

use 5.016;
use strict;
use warnings;
$| = 1;

##  Func: Global variable list
sub set_Global_Variable {
    our $dir_in;
    our $dir_out;
}

&main;
sub main {
    say "Run Create_KeyGeneSets_table.pl";
    our ( $dir_in, $dir_out );
    my  ( @files_summary, $f, $k );
    &set_Global_Variable;
    &get_CMD_line_ARGV;
    @files_summary = &get_input_file ( $dir_in );
    for $f( @files_summary ) { say $f;
        ($k = $f) =~ s/GeneSetInfoSummary/KeyGeneSet/;
        &create_KeyGeneSets_table ( $dir_in.$f, $dir_out.$k );
    }
}

##  Func: get arguments from outside of this program
sub get_CMD_line_ARGV {	#&get_CMD_line_ARGV;
    &help and exit if $ARGV[0] =~ /(-h)|(--h)|(-help)|(--help)/;
    my  %para = @ARGV;
    our $dir_in  = $para{-i};
    our $dir_out = $para{-o};
    say "Input Arguments:";
    say "@{[($_+1)/2]} $ARGV[$_-1] $para{$ARGV[$_-1]}" for grep { $_ & 1 } 0..$#ARGV;
}

sub get_input_file {
    our ( $dir_in );
    my  ( @files_summary );
    my  ( @files );
    @files = &read_dir ( $dir_in );
    (/.GeneSetInfoSummary.tsv/) ? (push @files_summary, $_) : () for @files;
    @files_summary;
}

##  Func: create table of top 3 genesets of each cluster
##
##  In:
##        file summary
##
##        GeneSet.name log.fold.change    p.Value              FDR                  Cluster    CI low            CI up
##        Ribosome     -0.993370528474239 4.21884749357559e-15 7.02438107680337e-13 Cluster 0 -1.22073936454752  -0.766558523732438
##        Proteasome   -0.29485343885929  2.89925281560244e-08 4.82725593797806e-07 Cluster 0 -0.398431190812074 -0.191782306712171
##        Ferroptosis  -0.184095657328657 4.33418345657799e-07 5.77313236416188e-06 Cluster 0 -0.254796924813815 -0.113567482403974
##        ...          ...                ...                  ...                  ...       ...                ...
##
##  Out:
##        file KeyGeneSet
##
##        ClusterID  GeneSet_1    GeneSet_2   GeneSet_3  logFC_1            logFC_2           logFC_3
##        Cluster 0  Ferroptosis  Proteasome  Ribosome   -0.184095657328657 -0.29485343885929 -0.993370528474239
##        Cluster 1  ...          ...         ...        ...                ...               ...
##        ...        ...          ...         ...        ...                ...               ...
sub create_KeyGeneSets_table {
    my ( $file_in, $file_out ) = @_[0,1];
    my ( @info, $Gset, $lgFC, $Clus, %all, $c, @key, @out );
    @info = &read_file_into_array ( $file_in, 'rm_header' );
    for ( @info ) {
        ($Gset, $lgFC, $Clus) = (split "\t", $_)[0,1,4];
        $all{$Clus}{$Gset} = $lgFC;
    }
    for $c( sort keys %all ) { #say $c;
        @key = sort { $all{$c}{$b} <=> $all{$c}{$a} } keys %{$all{$c}};
        push @out, (join "\t", $c, @key[0..2], $all{$c}{$key[0]}, $all{$c}{$key[1]}, $all{$c}{$key[2]});
    }
    unshift @out, "ClusterID\tGeneSet_1\tGeneSet_2\tGeneSet_3\tlogFC_1\tlogFC_2\tlogFC_3";
    &write_file ( $file_out, @out );
}


##  ---TOOLS------------------------------------------------------------TOOLS---

sub read_dir {
    my ( $dir, @f ) = $_[0];
    opendir DIR, $dir or die "Cannot open $dir: $!";
        (/^\./) ? () : (push @f, $_) for readdir DIR;	#say @f*1, "files";	#say for @f;
    closedir DIR;
    @f;
}

sub read_file_into_array { #&read_file_into_array (  );
    my ( $file, $rm_header, @info ) = @_[0,1]; #say "Reading file $file";
    open FH, $file or die "Cannot open $file for reading: $!";
        chomp for @info = <FH>;
        #s/\r?\n$// for @info = <FH>;
    close FH;
    shift @info if $rm_header;
    @info;
}

sub write_file { ## write info to file  #&write_file (  );
  	my ( $file, @info ) = @_;
  	open FH, '>', $file or die "Cannot open $file for writing: $!";
  	    say FH for @info;
  	close FH;  say "Finsh writing $file";
}

sub rm_dupl { ## remove duplicate of input array #&rm_dupl (  );
    my ( @dupl, %count, @uniq, @rm ) = @_;
  	#@rm = grep { ++$count{$_} > 1 } @dupl; say for "remove dupl element:", @rm;
    @uniq = grep { ++$count{$_} < 2 } @dupl;
}

#  ---HELP-------------------------------------------------------------------------------------------------------HELP---

sub help {
    print <<End;
-----------------------------
Instruction
-----------------------------
Usage:

Example:

perl 'ScriptPath/Create_KeyGeneSets_table.pl' -i FilePath/ -o FilePath/

End
}

## ---END---------------------------------------------------------------------------------------------------------END---
