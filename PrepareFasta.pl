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

# option list:
my ($seqType, @left_files, @right_files, @single_files, 
    $SS_lib_type
   );

my %allowed =
    ( seqType       => 'fa, or fq'
    );

my %allowed_check;
foreach my $all (keys %allowed) {
    my %h = map { (my $s = $_) =~ s/^or //; $s => 1 } split ', ', $allowed{$all};
    $allowed_check{$all} = \%h;
}


## Performance monitoring options 
my $pm_left_fa_size=0;
my $pm_right_fa_size=0;
my $pm_single_fa_size=0;

my $ROOTDIR            = "$FindBin::RealBin";

&GetOptions( 

    "seqType=s" => \$seqType,
    "left=s{,}" => \@left_files,
    "right=s{,}" => \@right_files,
    "single=s{,}" => \@single_files,
    
    "SS_lib_type=s" => \$SS_lib_type,
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

unless ( (@left_files && @right_files) || @single_files ) {
    die "Error, need either options 'left' and 'right' or option 'single'\n";
}


my $curr_limit_settings = `/bin/sh -c 'ulimit -a' `; 
unless ($curr_limit_settings && $curr_limit_settings =~ /\w/) {
    $curr_limit_settings = `/bin/csh -c limit`; # backup, probably not needed.
}

print "Current settings:\n$curr_limit_settings\n\n";


##################################################################################
#

main: {

    @left_files = &create_full_path(\@left_files) if @left_files;
    @right_files = &create_full_path(\@right_files) if @right_files;
    @single_files = &create_full_path(\@single_files) if @single_files;

    my $trinity_target_fa = (@single_files) ? "single.fa" : "both.fa"; 

	## Prep data for Inchworm
	if (@left_files && @right_files) {

            unless (-s $trinity_target_fa && !-e "left.fa" && !-e "right.fa") {
                
                my ($left_SS_type, $right_SS_type);
                if ($SS_lib_type) {
                    ($left_SS_type, $right_SS_type) = split(//, $SS_lib_type);
                }
                print("Converting input files. (in parallel)");
                my $thr1;
                my $thr2;
                if (!(-s "left.fa")) {
                    $thr1 = threads->create('prep_seqs', \@left_files, $seqType, "left", $left_SS_type);
                } else {
                    $thr1 = threads->create(sub { print ("left file exists, nothing to do");});
                }
                if (!(-s "right.fa")) {
                    $thr2 = threads->create('prep_seqs', \@right_files, $seqType, "right", $right_SS_type);
                } else {
                    $thr2 = threads->create(sub { print ("right file exists, nothing to do");});
                }
                @left_files = @{$thr1->join()};
                @right_files =@{$thr2->join()};
                print("Done converting input files.");
		## Calculate input file sizes for performance monitoring
		## this should be set as the created fasta otherwise results will differ for same data passed as .fq and .fa?
		my $pm_temp = -s "left.fa";
                $pm_temp = $pm_temp / 1024 / 1024;
                $pm_left_fa_size = sprintf('%.0f', $pm_temp);
                $pm_temp = -s "right.fa";
                $pm_temp = $pm_temp / 1024 / 1024;
                $pm_right_fa_size = sprintf('%.0f', $pm_temp);
                
                &process_cmd("cat left.fa right.fa > $trinity_target_fa") unless (-s $trinity_target_fa && (-s $trinity_target_fa == ((-s "left.fa") + (-s "right.fa"))));
                unless (-s $trinity_target_fa == ((-s "left.fa") + (-s "right.fa"))){
                    die "$trinity_target_fa is smaller (".(-s $trinity_target_fa)." bytes) than the combined size of left.fa and right.fa (".((-s "left.fa") + (-s "right.fa"))." bytes)\n";
                }

		# we keep if we have jaccard; delete later
		unlink ("left.fa", "right.fa"); # unless $jaccard_clip; # no longer needed now that we have 'both.fa', which is needed by chryaslis
            }
        }
	elsif (@single_files) {
            
            @single_files = @{&prep_seqs(\@single_files, $seqType, "single", $SS_lib_type) unless (-s "single.fa")};
	    ## Calculate input file sizes for performance monitoring
	    my $pm_temp = -s "single.fa";
            $pm_temp = $pm_temp / 1024 / 1024;
            my $pm_single_fa_size = sprintf('%.0f', $pm_temp);
        }
        
        else {
            die "not sure what to do. "; # should never get here.
        }

exit(0);

}

#####################################################################################

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

###
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $start_time = time();
    my $ret = system($cmd);
    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

    return;
}

###
sub prep_seqs {
    my ($initial_files_ref, $seqType, $file_prefix, $SS_lib_type) = @_;
    my @initial_files = @$initial_files_ref;
    return if -e "$file_prefix.fa";

        for (my $i=0;$i<scalar(@initial_files);$i++){
         my $f = $initial_files[$i];
         if ($f=~/\.gz$/){
          my $new = $f;
          $new=~s/\.gz$//;
          unlink($new);
          &process_cmd("gunzip -c $f > $new");
          $initial_files[$i] = $new;
         }elsif ($f=~/\.bz2$/){
          my $new = $f;
          $new=~s/\.bz2$//;
          unlink($new);
          &process_cmd("bunzip2 -dkc $f > $new");
          $initial_files[$i] = $new;
         }
        }

        my $initial_file_str = join(" ",@initial_files);
    if ($seqType eq "fq") {
       # make fasta
       foreach my $f (@initial_files){
         my $perlcmd = "$UTILDIR/fastQ_to_fastA.pl -I $f ";
         my $fastool_cmd = "$FASTOOL_DIR/fastool";
         if ($SS_lib_type && $SS_lib_type eq "R") {
             $perlcmd .= " --rev ";
             $fastool_cmd .= " --rev ";
         }
         $fastool_cmd .= " --illumina-trinity --to-fasta $f >> $file_prefix.fa";
         $perlcmd .= " >> $file_prefix.fa";
         my $cmd = ($USE_FASTOOL) ? $fastool_cmd : $perlcmd;
         &process_cmd($cmd);
        }
    }
    elsif ($seqType eq "fa") {
        if (scalar(@initial_files) == 1 && (!$SS_lib_type || $SS_lib_type ne "R")) {
            ## just symlink it here:
            my $cmd = "ln -s $initial_file_str $file_prefix.fa";
            &process_cmd($cmd);
        }elsif(scalar(@initial_files) > 1 && (!$SS_lib_type || $SS_lib_type ne "R")){
                my $cmd = "cat $initial_file_str > $file_prefix.fa";
                &process_cmd($cmd);
        }else {
          #if ($SS_lib_type && $SS_lib_type eq "R") {
                  foreach my $f (@initial_files){
                my $cmd = "$UTILDIR/revcomp_fasta.pl $f >> $file_prefix.fa";
                &process_cmd($cmd);
                }
        }
    }
    elsif (($seqType eq "cfa") | ($seqType eq "cfq")) {
        # make double-encoded fasta
        foreach my $f (@initial_files){
         my $cmd = "$UTILDIR/csfastX_to_defastA.pl -I $f ";
         if ($SS_lib_type && $SS_lib_type eq "R") {
             $cmd .= " --rev ";
         }
         $cmd .= ">> $file_prefix.fa";
         &process_cmd($cmd);
        }
  }
  return \@initial_files;
}




