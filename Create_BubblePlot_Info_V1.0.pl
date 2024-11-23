#! \usr\bin\perl
## Copyright (C) 2019-2020 Novelbio
## Running Platform: Win / Novelbrain / Ubuntu
## Updator: LI, Yadang          Update time: 2020.05.18
## Author:  LI, Yadang          Create time: 2019.03.27
## Instruction: The function of this script is to create input files(including Means_for_Plot.txt; Pvalues_for_Plot.txt; Cluster_name_Ligand.gene.txt; Cluster_name_Ligand.cluster.txt; Cluster_name_Receptor.gene.txt; Cluster_name_Receptor.cluster.txt) for ScBubblePlot Task based on pvalues_0.05.txt, pvalues_revise.txt, significant_means_revise.txt, species, table file of comparison of gene naems of different species and Cluster_Type_file(optional)
## Function List:
##   main
##   get_cellphone_data
##   get_cluster_info
##   get_gene_name_switch_table
##   create_bubble_info
##     get_GG
##     get_CC_based_on_GG
##     output_Bubble_info
##     switch_gene_by_species

use 5.018;
use strict;
use warnings;

##  Func: Global variable list
sub set_Global_Variable {
    our $sep = ' >> ';        ## separator between cluster or gene
    our $dir_out;             ## directory of output
    our $MaxGenePair;         ## max gene pairs shown in each GenePair result file
    our $file_sig_result;     ## pvalues_0.05.txt
    our $file_pval;           ## pvalues_revise.txt
    our $file_mean;           ## significant_means_revise.txt
    our $species;             ## human / any other non human
    our $file_cluster_rename; ## file of list of clusters with their new names
    our $file_Gene_table;     ## file of list of gene names of the species being analyzed with corresponding Human ID
}

&main;
sub main {
    say "Run Create_BubblePlot_Info_V1.0.pl";
    &set_Global_Variable;
    our ( $dir_out, @cluster_type, @cluster_name, @cellphone_data, %human2other );
    our ( $MaxGenePair, $file_sig_result, $file_pval, $file_mean, $species, $file_Gene_table, $file_cluster_rename );
    &read_ARGV;
    my $dir_plot = $dir_out.'/BubblePlot_Files';
    mkdir $dir_plot unless -e $dir_plot;
    my ( $plotfile_means, $plotfile_pvalues, $plotfile_genepair, $plotfile_clusterpair ) = (
        $dir_plot.'/Means_for_Plot.txt',
        $dir_plot.'/Pvalues_for_Plot.txt',
        $dir_plot.'/GenePair.txt',
        $dir_plot.'/ClusterPair.txt'
    );
    my ( @type_class );
    &get_cellphone_data ( $file_sig_result );
    &get_cluster_info ( $file_sig_result, $file_cluster_rename );
    @type_class = &remove_duplicate ( @cluster_type );
    system "mkdir $dir_plot" unless -e $dir_plot;
    if ( $species eq "human" ) {
        system "cp pvalues_revise.txt $dir_plot/Pvalues_for_Plot.txt";
        system "cp significant_means_revise.txt $dir_plot/Means_for_Plot.txt";
        &create_bubble_info ( $MaxGenePair, $_, $plotfile_genepair, $plotfile_clusterpair ) for @type_class;
    } else {  #say $species;
        &get_gene_name_switch_table ( $species, $file_Gene_table ); say "finish getting gene name switch table";
        &switch_gene_name_of_CellPhoneDB_result ( $file_mean, $plotfile_means ); say "finish getting plotfile_means";
        &switch_gene_name_of_CellPhoneDB_result ( $file_pval, $plotfile_pvalues ); say "finish getting plotfile_pvalues";
        &create_bubble_info ( $MaxGenePair, $_, $plotfile_genepair, $plotfile_clusterpair, $species ) for @type_class;
    }
}

##  Func: get arguments from outside of this program
sub read_ARGV {	#&read_ARGV;
    my  %para = @ARGV;
    our ( $MaxGenePair, $file_sig_result, $file_pval, $file_mean, $species, $file_Gene_table, $file_cell_cluster, $dir_out ) = (
        $para{-n}, $para{-r}, $para{-p}, $para{-m}, $para{-s}, $para{-g}, $para{-c}, $para{-o},
    );
    say "Input Arguments:";
    say "@{[($_+1)/2]} $ARGV[$_-1] $para{$ARGV[$_-1]}" for grep { $_ & 1 } 0..$#ARGV;
}

# Function: to get the cellphone data from pvalues_0.05.txt and sort them from small to large by value of rank
# Input : $cellphone_file (pvalues_0.05.txt)
# Output: @cellphone_data
sub get_cellphone_data { say "start to get CellPhoneDB data from pvalues_0.05.txt";
    our ( @cellphone_data );
    my ( $cellphone_file ) = $_[0];
    @cellphone_data = &readfile ( $cellphone_file );
    shift @cellphone_data;
    @cellphone_data = sort {$a =~ /([^\t]+)\t[^\t]+$/; my $A = $1; $b =~ /([^\t]+)\t[^\t]+$/; my $B = $1; $A <=> $B} @cellphone_data;
}

# Function: to get the cluster classification info based on Cluster_rename_file(if exist) or pvalues_0.05.txt
# Input : $file_cell_cluster.txt(if exist)
# Output: @cluster_type and @cluster_name (these two arrays have same number of elements with one-to-one corresponding relaitons)
# Example:
#   output data from Cluster_Type_file:
#   @cluster_type @cluster_name
#           Tumor Tumor1
#           Tumor Tumor2
#           Bcell B1
#           Bcell B2
#             ... ...
#   output data from pvalues_0.05.txt
#   @cluster_type @cluster_name
#          Tumor1 Tumor1
#          Tumor2 Tumor2
#              B1 B1
#              B2 B2
#             ... ...
sub get_cluster_info {
    our ( @cellphone_data, @cluster_type, @cluster_name );
    my ( $file_sig_result, $cluster_type_file, @info, $t, $c, @cluster ) = @_[0,1];
    if ( defined $cluster_type_file ) { say "get info based on $cluster_type_file";
        @info = &readfile ( $cluster_type_file );
        ( $t, $c ) = split "\t" and push @cluster_name, $c and push @cluster_type, $t for @info[1..$#info];
    }
    else { say "get info based on $file_sig_result";
        for ( @cellphone_data ) {
            @info = split "\t";
            push @cluster, $info[10], $info[11];
        }
        @cluster = &remove_duplicate ( @cluster );
        @cluster_name = @cluster; @cluster_type = @cluster;
    }
}

# Function: get gene name comparison list between specified species and human
# Input:  $speices $file_Gene_table
# Output: %human2other
sub get_gene_name_switch_table {
    our ( %human2other );
    my ( $species, $file_Gene_table ) = @_[0,1];
    my ( @info, @list );
    @info = &readfile ( $file_Gene_table ); #say scalar @info;
    (/\t(.+)\t(.+)/) ? ($human2other{$2} .= ",$1") : () for @info[1..$#info];
    for ( keys %human2other ) {
        $human2other{$_} =~ s/^,//;
        @list = split ",", $human2other{$_};
        @list = &remove_duplicate ( @list );
    #    s/NEWGENE_/NEWGENE/ for @list;
    #    (/_/) ? (say "!!!Find illegal character \"_\" in gene name $_") : () for @list;
        $human2other{$_} = join ",", @list;
    } #say for keys %human2other; say for values %human2other;
    say "Size of gene name switch table: @{[scalar keys %human2other]}";
}

# Function: create the means/pvalues files for ScBubblePlot Task
# Input : $file_pval(pvalues_revise.txt)/$file_mean(significant_means_revise.txt)
sub switch_gene_name_of_CellPhoneDB_result {
    our ( @gene_pair );
    my ( $file_input, $file_output, @data, $pre, $suf ) = @_[0,1];
    @data = &readfile ( $file_input );
    ($pre, $gene_pair[0], $suf) = /^([^\t]+\t)([^\t]+)(\t.+)$/ and &switch_gene_by_species and $_ = $pre.$gene_pair[0].$suf for @data;
    undef @gene_pair;
    &write_file ( $file_output, @data );
}

# Function: create the gene pair/cluster pair files for ScBubblePlot Task
# Input : @cellphone_data
sub create_bubble_info {
    our ( @cellphone_data, @subdata, @gene_pair );
    my ( $MaxGenePair, $type_class, $gene_pair_file, $cluster_pair_file, $species, @info, $d, $i1, $i2 ) = @_;
    for $d( @cellphone_data ) { #&show_progress_bar (scalar @cellphone_data, ++$i1, 3);
        @info = split "\t", $d;
        &get_GG ( $type_class, $info[10], $d, $info[1] );
    } #say "finish getting gene pair"; #say "@{[$#gene_pair+1]} Gene Pair";
    @gene_pair = &remove_duplicate ( @gene_pair ); #say "@{[$#gene_pair+1]} Gene Pair";
    &get_CC_based_on_GG ( "L", $MaxGenePair );
    &switch_gene_by_species if $species;
    &output_Bubble_info ( $gene_pair_file, $cluster_pair_file, $type_class, "Ligand", $MaxGenePair );
    for $d( @cellphone_data ) { #&show_progress_bar (scalar @cellphone_data, ++$i2, 3);
        @info = split "\t", $d;
        &get_GG ( $type_class, $info[11], $d, $info[1] );
    } #say "finish getting gene pair"; #say "@{[$#gene_pair+1]} Gene Pair";
    @gene_pair = &remove_duplicate ( @gene_pair ); #say "@{[$#gene_pair+1]} Gene Pair";
    &get_CC_based_on_GG ( "R", $MaxGenePair );
    &switch_gene_by_species if $species;
    &output_Bubble_info ( $gene_pair_file, $cluster_pair_file, $type_class, "Receptor", $MaxGenePair );
}

# Function: to create gene_gene relation by specified cluster type and extract corresponding cellphone data
# Input : $type_class (specified cluster type) $cluster (actual cluster) $data (corresponding cellphone data) $GG(gene_gene info)
# Output: @gene_pair (gene_gene list) @subdata (corresponding cellphone data)
sub get_GG { #say "start to get Gene Pair";
    our ( @cluster_type, @cluster_name, @gene_pair, @subdata );
    my ( $type_class, $cluster, $data, $GG ) = @_[0..3];
      ($cluster eq $cluster_name[$_] and $cluster_type[$_] eq $type_class)
    ? (push @subdata, $data and push @gene_pair, $GG)
    : ()
    for 0 .. $#cluster_name; #say "@{[$#gene_pair+1]} Gene Pair";
}

# Function: to create cluster_pair based on gene_pair and corresponding cellphone data
# Input : @gene_pair (gene_pair list) @subdata (corresponding cellphone data) $LR_type (Ligand or Receptor)
# Output: @clus_pair (cluster_pair list)
# Example:
#                 for ligand:   for Receptor:
#   input         output        output
#   clus1 clus3   clus1 clus3   clus1 clus3
#   clus2 clus4   clus1 clus4   clus2 clus3
#                 clus2 clus3   clus1 clus4
#                 clus2 clus4   clus2 clus4
sub get_CC_based_on_GG { #say "start to get cluster Pair";
    our ( $sep, @subdata, @gene_pair, @clus_pair );
    my ( $LR_type, $MaxGenePair, @info, $line, $g, @clus_a, @clus_b, $l, $r ) = @_[0,1]; #say for @gene_pair;
    $MaxGenePair = $#gene_pair + 1 if $MaxGenePair > $#gene_pair + 1; #say for @gene_pair[0..$MaxGenePair];
    for $g( 0 .. $MaxGenePair-1 ) { #&show_progress_bar ($MaxGenePair, $g, 2);
        for $line( @subdata ) {
            @info = split "\t", $line;
            push @clus_a, $info[10] and push @clus_b, $info[11] if $info[1] eq $gene_pair[$g];
        }
    }
    @clus_a = &remove_duplicate ( @clus_a );
    @clus_b = &remove_duplicate ( @clus_b );
    if ( $LR_type eq "L" ) {
        for $l( @clus_a ) {
            push @clus_pair, $l.$sep.$_ for @clus_b;
        }
    }
    if ( $LR_type eq "R" ) {
        for $r( @clus_b ) {
            push @clus_pair, $_.$sep.$r for @clus_a;
        }
    }
    undef @subdata;
}

# Function: switch gene names from human to other Species
# Input: @gene_pair (gene_pair list)
sub switch_gene_by_species {
    our ( $sep, %human2other, @gene_pair );
    my ( @gene, $g, $p, $a, $b );
    for ( @gene_pair ) {
        s/^(.+)($sep.+)$/$human2other{$1}$2/ if /^(.+)$sep(.+)$/ and exists $human2other{$1};
        s/^(.+$sep)(.+)$/$1$human2other{$2}/ if /^(.+)$sep(.+)$/ and exists $human2other{$2};
    }
   1;
}

# Function: to create gene_pair list file and clus_pair list file
# Input : @gene_pair (gene_pair list) @clus_pair (cluster_pair list) $file_GR $file_CR (file name) $clus_type $LR_type
sub output_Bubble_info {
    our ( @gene_pair, @clus_pair );
    my ( $file_GR, $file_CR, $clus_type, $LR_type, $MaxGenePair ) = @_[0..4];
    ($#gene_pair + 1 < $MaxGenePair) ? ($MaxGenePair = $#gene_pair) : ($MaxGenePair--); say "Number of Gene Pair: @{[$MaxGenePair+1]}";
    $file_GR =~ s/GenePair.txt/$clus_type\_$LR_type\.GenePair.txt/;
    $file_CR =~ s/ClusterPair.txt/$clus_type\_$LR_type\.ClusterPair.txt/;
    &write_file ( $file_GR, 'GenePair', @gene_pair[0..$MaxGenePair] );
    &write_file ( $file_CR, 'ClusterPair', @clus_pair );
    undef @gene_pair; undef @clus_pair;
}

#  ---TOOL-------------------------------------------------------------------------------------------------------TOOL---

sub readfile {
    my ( $file, @read ) = $_[0];
    open FH, $file;
        chomp for @read = <FH>;
    close FH;
    @read;
}

sub write_file { ## write info to file  #&write_file (  );
  	my ( $file, @info ) = @_;
  	open FH, '>', $file or die "Cannot open $file for writing: $!";
  	    say FH for @info;
  	close FH;  say "Finsh writing $file";
}

sub remove_duplicate {
    my ( @array, %count, @uniq_itmes ) = @_;
    @uniq_itmes = grep { ++$count{ $_ } < 2 } @array; #say for @uniq_itmes; @uniq_itmes;
}

sub show_progress_bar{	# Usage: show the progress of a loop in a bar	#&show_progress_bar (scalar @array, ++$i, 10);
  	my ( $l, $i, $f, $d, $ii, $p ) = @_[0..2];	# input: list_length, cycle_index, display_frequency into $l $i $f
  	$f = 10 unless defined $f;
  	$d = length $l;
  	$ii = sprintf "%${d}d", $i;
  	$p = sprintf "%3d", int($i/$l*100);
  	say "[$ii/$l  - $p% ]" if int($l/$f) != 0 and $i % (int ($l / $f) +0) == 0;
}
#  ---END---------------------------------------------------------------------------------------------------------END---
