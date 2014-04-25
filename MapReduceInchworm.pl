#!/usr/bin/env perl

use strict;
use warnings;
use threads;
no strict qw(subs refs);

use FindBin;
use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/PerlLibAdaptors");
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case pass_through);

open (STDERR, ">&STDOUT"); 

my $CPU = 1;
my $IWORM_KMER_SIZE = 25;
my $MIN_IWORM_LEN = 25;

my $INCHWORM_CUSTOM_PARAMS;

# option list:
my ($seqType, @left_files, @right_files, @single_files, 
    $SS_lib_type, $min_contig_length, $output_directory
   );

my %allowed =
    ( seqType       => 'fa, or fq'
    );

my %allowed_check;
foreach my $all (keys %allowed) {
    my %h = map { (my $s = $_) =~ s/^or //; $s => 1 } split ', ', $allowed{$all};
    $allowed_check{$all} = \%h;
}

# default
#
$output_directory = &create_full_path("mrInchworm_out_dir");

$min_contig_length = 200;
my $path_reinforcement_distance;
my $PE_path_reinforcement_distance = 75;
my $SE_path_reinforcement_distance = 25;

my $min_kmer_cov = 1;
my $min_percent_read_iworm_kmers = -1;

my $usage = <<_EOUSAGE_;

####################################################################
#    
#   MapReduce-Inchworm
#
###################################################################
#
#  -seqType <string>      :type of reads: ( $allowed{seqType} ) 
#  
#
#
#

_EOUSAGE_

;

my $ROOTDIR = "$FindBin::RealBin";
my $MR_INCHWORM_DIR = "$ROOTDIR/MR_Inchworm";
my $FASTA_SPLITTER_DIR = "$ROOTDIR/Fasta_Splitter";

my $FASTOOL_DIR = "$ROOTDIR/trinity-plugins/fastool";

#unless (@ARGV) {
#    die "$usage\n";
#}

&GetOptions( 

    "seqType=s" => \$seqType,
    "left=s{,}" => \@left_files,
    "right=s{,}" => \@right_files,
    "single=s{,}" => \@single_files,
    
    "SS_lib_type=s" => \$SS_lib_type,

    "output=s" => \$output_directory,
    
    "min_contig_length=i" => \$min_contig_length,

    'CPU=i' => \$CPU,

    'min_kmer_cov=i'        => \$min_kmer_cov,
    'INCHWORM_CUSTOM_PARAMS=s' => \$INCHWORM_CUSTOM_PARAMS,

    'KMER_SIZE=i' => \$IWORM_KMER_SIZE,

);



sub check_option {
    my ($option, $name) = @_;
    $$option = lc $$option;
    if ($$option eq '') {
        die "Error, option '--$name' is required.\n";
    }
    if (!defined $allowed_check{$name}{$$option}) {
        die "Error, option '--$name' ($$option) not one of $allowed{$name}\n";
    }
}

check_option( \$seqType,     'seqType'     );

if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(R|F|RF|FR)$/) {
        die "Error, unrecognized SS_lib_type value of $SS_lib_type. Should be: F, R, RF, or FR\n";
    }
}

main: {



exit(0);

}



####
sub create_full_path {
    my ($file) = shift;
    if (ref($file) eq "ARRAY"){
       for (my $i=0;$i<scalar(@$file);$i++){
         $file->[$i] = &create_full_path($file->[$i]);
       }
       return @$file;
    }else{
      my $cwd = cwd();
      if ($file !~ m|^/|) { # must be a relative path
          $file = $cwd . "/$file";
      }
      return($file);
    }
}

