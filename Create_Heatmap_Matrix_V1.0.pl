#! /usr/bin/perl
## Copyright (C) 2020 Novelbio
## Running Platform: Linux / Novelbrain
## Updator: LI, Yadang          Update time: 2020.05.20
## Author:  LI, Yadang          Create time: 2020.05.20
## Instruction: The function of this script is to create heatmap matrix table

use 5.018;
use strict;
use warnings;
use Cwd;
#say our $cwd = getcwd;

##  Func: Global variable list
sub set_Global_Variable {
    our $file_sig_result;   ## pvalues_0.05.txt
    our $file_cell_cluster; ## file of list of cells with their cluster
    our @cluster;           ## list of all cluster ordered in alphanumerical way
    our $file_CCI_matrix;   ## file of matrix of interaction numbers of each cluster pair
}

&main;
sub main {
    say "Run Create_Heatmap_Matrix_V1.0.pl";
    &set_Global_Variable;
    &read_ARGV;
    &get_cluster_info;
    &create_CCI_matrix;
}

##  Func: get arguments from outside of this program
sub read_ARGV {
    my %para = @ARGV;
    our ( $file_sig_result, $file_cell_cluster, $file_CCI_matrix ) = (
        $para{-i}, $para{-c}, $para{-o},
    );
    say "Input Arguments:";
    say "@{[($_+1)/2]} $ARGV[$_-1] $para{$ARGV[$_-1]}" for grep { $_ & 1 } 0..$#ARGV;
}

sub get_cluster_info {
    our ( $file_cell_cluster, @cluster );
    my  ( @info );
    @info = &read_file_into_array ( $file_cell_cluster );
    push @cluster, (split "\t", $_)[1] for @info[1..$#info];
    @cluster = sort { &cmp_alphanumerical ($a, $b) } &rm_dupl ( @cluster );
}

sub create_CCI_matrix {
    our ( @cluster, $file_sig_result, $file_CCI_matrix );
    my  ( $l, $r, @info, $ligand, $receptor, %CCI_num, @output );
    for $l ( @cluster ) {
        for $r ( @cluster ) {
            $CCI_num{$l}{$r} = 0; ## Initialization of %CCI_num;
        }
    }
    @info = &read_file_into_array ( $file_sig_result );
    for ( @info[1..$#info] ) {
        ( $ligand, $receptor ) = (split "\t", $_)[-5,-4];
        $CCI_num{$ligand}{$receptor}++;
    }
    for $r ( @cluster ) {
        push @output, $r;
        for $l ( @cluster ) {
            $output[-1] .= "\t$CCI_num{$l}{$r}";
        }
    }
    unshift @output, "ClusterID";
    $output[0] .= "\t$_" for @cluster;
    &write_file ( $file_CCI_matrix, @output );
}

#  ---TOOL-------------------------------------------------------------------------------------------------------TOOL---

sub read_file_into_array { #&read_file_into_array (  );
  	my ( $file, @info ) = $_[0]; #say "Reading file $file";
  	open FH, $file or die "Cannot open $file for reading: $!";
  		  #chomp for @info = <FH>;
        s/\r?\n$// for @info = <FH>;
  	close FH;
  	@info;
}

sub rm_dupl { ## remove duplicate of input array #= &rm_dupl (  );
    my ( @dupl, %count, @uniq, @rm ) = @_;
  	#@rm = grep { ++$count{$_} > 1 } @dupl; say for "remove dupl element:", @rm;
    @uniq = grep { ++$count{$_} < 2 } @dupl;
}

sub write_file { ## write info to file  #&write_file (  );
  	my ( $file, @info ) = @_;
  	open FH, '>', $file or die "Cannot open $file for writing: $!";
  	    say FH for @info;
  	close FH;  say "Finsh writing $file";
}

## Func: cmp 2 elements in alphanumerical way
##       element could be a number or non-digit string or non-digit string plus a numbers (e.g. "0", "A", "A0")
##       Example:
##           Input1 Input2 Output
##               1      1       0
##               1      2      -1
##               1      0       1
##              b      b        0
##              b      c       -1
##              b      a        1
##              b1     b1       0
##              b1     b2      -1
##              b1     b0       1
##              b1     a2      -1
##              b1     c0       1
##              b1     b        1
##              b1     a        1
##              b1     c       -1
##              b1      1       1
##              b       1       1
sub cmp_alphanumerical {
    my ( $e1, $e2, $c1, $c2, $d1, $d2, $result ) = @_[0,1];
    ($c1, $d1) = $e1 =~ /(\D+)?(\d+)?/;
    ($c2, $d2) = $e2 =~ /(\D+)?(\d+)?/;
    if ( defined $c1 and defined $c2 ) {
        $result = $c1 cmp $c2;
        if ( $c1 eq $c2 and defined $d1 and defined $d2 ) {
             $result = $d1 <=> $d2;
        }
        elsif ( !defined $d1 and !defined $d2 ) {}
        else {
            $result = $e1 cmp $e2;
        }
    }
    elsif ( !defined $c1 and !defined $c2 ) {
        $result = $d1 <=> $d2;
    }
    else {
        $result = $e1 cmp $e2;
    }
    $result;
}

#  ---END---------------------------------------------------------------------------------------------------------END---
