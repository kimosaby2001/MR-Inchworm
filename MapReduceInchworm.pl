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

my $VERSION = "MR_Inchworm_r2015";

my $CPU = 1;
my $MR_PAGE_SIZE = 1024;


my $IWORM_KMER_SIZE = 25;
my $MIN_IWORM_LEN = 25;

my $long_reads = "";

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
my $min_edge_cov = 1;
my $min_percent_read_iworm_kmers = -1;

## Performance monitoring options 
my $pm_logfile = "MR_Inchworm.timing";
my $pm_mrInchworm_start=0;
my $pm_mrInchworm_end=0;
my $pm_left_fa_size=0;
my $pm_right_fa_size=0;
my $pm_single_fa_size=0;
my $pm_trinity_fa_size=0;
my $pm_trinity_arguments="";
my $pm_inchworm_kmers=0;
my $pm_read_count=0;

my $run_with_collectl = 0;
# flush each second, record procs+rest every 5 secs, use only process subsystem
#my $collectl_param = "-F1 -i5:5 -sZ";
my $collectl_param = "-F1 -i5:5 -scfnmZ";
my $collectl_output_directory = "collectl";
my $collectl_pid = 0;
my $collectl_out = "";
my $collectl_titlename = "";
my $start_dir = cwd();

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

my $ROOTDIR            = "$FindBin::RealBin";
my $MR_INCHWORM_DIR    = "$ROOTDIR/MR_Inchworm";
my $FASTA_SPLITTER_DIR = "$ROOTDIR/Fasta_Splitter";
my $FASTOOL_DIR        = "$ROOTDIR/fastool";
my $UTILDIR            = "$ROOTDIR/util"; 
my $COLLECTL_DIR       = "$ROOTDIR/Collectl/bin";

unless (@ARGV) {
    die "$usage\n";
}

my $FULL_CLEANUP = 0;
my $NO_FASTOOL = 0;

&GetOptions( 

    "seqType=s" 	=> \$seqType,
    "left=s{,}" 	=> \@left_files,
    "right=s{,}" 	=> \@right_files,
    "single=s{,}" 	=> \@single_files,
    
    "SS_lib_type=s" 	=> \$SS_lib_type,

    "long_reads=s" 	=> \$long_reads,

    "output=s" 		=> \$output_directory,
    
    "min_contig_length=i" => \$min_contig_length,

    'CPU=i' 		=> \$CPU,

    'min_kmer_cov=i'        => \$min_kmer_cov,
    'min_edge_cov=i'        => \$min_edge_cov,
    'INCHWORM_CUSTOM_PARAMS=s' => \$INCHWORM_CUSTOM_PARAMS,

    'no_fastool' 	=> \$NO_FASTOOL,

    'KMER_SIZE=i' 	=> \$IWORM_KMER_SIZE,

    'page_size=i' 	=> \$MR_PAGE_SIZE,

    'monitoring' 	=> \$run_with_collectl,
    'collectl_dir=s' 	=> \$collectl_output_directory,

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

if ($run_with_collectl && $^O !~ /linux/i) {
    print STDERR "WARNING, --monitoring can only be used on linux. Turning it off.\n\n";
    $run_with_collectl = 0;
}


if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(R|F|RF|FR)$/) {
        die "Error, unrecognized SS_lib_type value of $SS_lib_type. Should be: F, R, RF, or FR\n";
    }
}

my $USE_FASTOOL = 1;
if ($NO_FASTOOL) {
    $USE_FASTOOL = 0;
}

unless ( (@left_files && @right_files) || @single_files ) {
    die "Error, need either options 'left' and 'right' or option 'single'\n";
}


my $curr_limit_settings = `/bin/sh -c 'ulimit -a' `; 
unless ($curr_limit_settings && $curr_limit_settings =~ /\w/) {
    $curr_limit_settings = `/bin/csh -c limit`; # backup, probably not needed.
}

print "Current settings:\n$curr_limit_settings\n\n";


sub collectl_start {
    # install signal handler to stop collectl on interrupt
    $SIG{INT} = sub { print "Trinity interrupted\n"; &collectl_stop(); exit(1); };

    if ($run_with_collectl){
        warn "STARTING COLLECTL\n";
        $collectl_output_directory = "$start_dir/collectl";
        `rm -rf $collectl_output_directory `;
        $collectl_output_directory = &create_full_path($collectl_output_directory);
        unless (-d $collectl_output_directory) {
            mkdir $collectl_output_directory or die "Error, cannot mkdir $collectl_output_directory";
        }
        my $collectl_userid = qx(id --user --real);
        chomp($collectl_userid);
        my $cmd = "cd $collectl_output_directory && exec ${COLLECTL_DIR}/collectl $collectl_param --procfilt u$collectl_userid -f $collectl_output_directory/y";
        ## fork a child to run collectl
        $collectl_pid = fork();
        if (not defined $collectl_pid) {
            warn "FORK FAILED - NO COLLECTL PROCESS STARTED\n";
        } elsif ($collectl_pid == 0) {
            warn "I'M THE CHILD RUNNING TRINITY\n";
            exec($cmd);
            warn "COLLECTL FINISHED BEVORE KILL WAS CALLED\n";
            exit(0);
        } else {
        warn "I'M THE PARENT, COLLECTL_PID=$collectl_pid\n";
        }
    }
}

sub collectl_stop {
    if ($run_with_collectl && $collectl_pid>0) {
        warn "TERMINATING COLLECTL, PID = $collectl_pid\n";
        # try to be nice here as a hard kill will result in broken/unusable raw.gz file
        system("sync");
        kill("INT", $collectl_pid);
        kill("TERM", $collectl_pid);
        waitpid($collectl_pid,0);
        chdir($collectl_output_directory) or return;
        system("$COLLECTL_DIR/make_data_files.sh");
        system("$COLLECTL_DIR/timetable.sh");
        $collectl_titlename = "${VERSION} ${CPU} @{left_files}@{single_files}";
        system("$COLLECTL_DIR/plot.sh \"$collectl_titlename\" ${CPU}");
    }
}


##################################################################################
#

my $MKDIR_OUTDIR_FLAG = 0;

main: {

    @left_files = &create_full_path(\@left_files) if @left_files;
    @right_files = &create_full_path(\@right_files) if @right_files;
    @single_files = &create_full_path(\@single_files) if @single_files;
    $output_directory = &create_full_path($output_directory);
    $long_reads = &create_full_path($long_reads) if $long_reads;

    unless (-d $output_directory) {
        &process_cmd("mkdir -p $output_directory");
        $MKDIR_OUTDIR_FLAG = 1;
    }

    chdir ($output_directory) or die "Error, cannot cd to $output_directory"; 

    ## create inchworm file name
    my $inchworm_file = "inchworm.K$IWORM_KMER_SIZE.L$MIN_IWORM_LEN";
    unless ($SS_lib_type) {
        $inchworm_file .= ".DS";
    }
    $inchworm_file .= ".fa";
    $inchworm_file = &create_full_path($inchworm_file);

    my $trinity_target_fa = (@single_files) ? "single.fa" : "both.fa"; 
    my $inchworm_target_fa = $trinity_target_fa; 

    ## Don't prep the inputs if Inchworm already exists.... Resuming earlier operations.
    my $inchworm_finished_checkpoint_file = "$inchworm_file.finished";
    if (-s $inchworm_file && -e $inchworm_finished_checkpoint_file) {
        print "\n\n#######################################################################\n"
            . "Inchworm file: $inchworm_file detected.\n"
            . "Skipping Inchworm Step, Using Previous Inchworm Assembly\n"
            . "#######################################################################\n\n";
        sleep(2);
    }
    else {

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

	my $count_of_reads = `wc -l < $inchworm_target_fa`;chomp($count_of_reads); #AP: grep is  expensive; one test took 2h...!
        $count_of_reads/=2;

        if ($long_reads) {
            $inchworm_target_fa .= ".wLongReads.fa";
            $count_of_reads += `grep -c '^>' $long_reads | wc -l`; #AP we don't know if these will be one single line
            &process_cmd("cat $long_reads $trinity_target_fa > $inchworm_target_fa");
        }
            
        open (my $ofh, ">$inchworm_target_fa.read_count") or die $!;
        print $ofh $count_of_reads."\n";
        close $ofh;
    }

    collectl_start() unless ($FULL_CLEANUP);

    ##############################################
    # MR-Inchworm
    #

    $pm_mrInchworm_start = `date +%s`;
    unless (-s $inchworm_file && -e $inchworm_finished_checkpoint_file) {
        &run_mrInchworm($inchworm_file, $inchworm_target_fa, $SS_lib_type);
        &process_cmd("touch $inchworm_finished_checkpoint_file");
    }
    $pm_mrInchworm_end = `date +%s`;

    my $pm_mrInchworm_time = $pm_mrInchworm_end - $pm_mrInchworm_start;
    print "\n  mrInchworm took  $pm_mrInchworm_time seconds\n";

exit(0);

}

#####################################################################################

sub run_mrInchworm {
    my ($inchworm_outfile, $reads, $strand_specific_flag) = @_;

    ## get count of number of reads to be assembled.
    my $read_count_file = "$reads.read_count";
    if (! -s $read_count_file) {
        my $count_of_reads = `wc -l < $reads`;chomp($count_of_reads); #AP: grep is  expensive; one test took 2h...! 
        $count_of_reads/=2;  # assume fasta; two lines per read
        $pm_read_count = $count_of_reads;
        open (my $ofh, ">$read_count_file") or die $!;
        print $ofh $count_of_reads."\n";
        close $ofh;
    }

    my $inchworm_cmd;
    
    my @tmp_files; # to be deleted after successful inchworm run.

    ######################################################
    #   C-FastaSplitter
    #################################################### 

    my $reads_header = "$reads.header";
    my $cmd_header = "grep -nr \">\" $reads | awk -F \":>\" '{print \$1}' > $reads_header"; 
    &process_cmd($cmd_header);

    my $sKmerDir = "sKmer_tmp";
    if (-d $sKmerDir) {
        &process_cmd("rm -r $sKmerDir");
    }
    &process_cmd("mkdir -p $sKmerDir");

    my $cmd_splitter = "mpirun -mca mpi_warn_on_fork 0 -np $CPU $FASTA_SPLITTER_DIR/Fasta_Splitter -r $reads";
    $cmd_splitter = $cmd_splitter . " -i $reads_header -o $sKmerDir";
    &process_cmd($cmd_splitter);

    ########################################################
    # MR-Inchworm
    #######################################################


    print $MR_PAGE_SIZE, "\n";


#    my $cmd_mrInchworm = "mpirun --mca btl self,openib -np $CPU $MR_INCHWORM_DIR/mr_inchworm -K $IWORM_KMER_SIZE -L $MIN_IWORM_LEN";
    my $cmd_mrInchworm = "mpirun --mca orte_base_help_aggregate 0 -np $CPU $MR_INCHWORM_DIR/mr_inchworm -K $IWORM_KMER_SIZE -L $MIN_IWORM_LEN";
    $cmd_mrInchworm .= " --PageSize $MR_PAGE_SIZE";

    if ($min_kmer_cov > 1) {
            $cmd_mrInchworm .= " --minKmerCount $min_kmer_cov";
    } 

    if ($min_edge_cov > 1) {
            $cmd_mrInchworm .= " --minEdgeCount $min_edge_cov";
    }
   
    unless ($strand_specific_flag) {
        $cmd_mrInchworm .= " --DS ";
    }

    if ($INCHWORM_CUSTOM_PARAMS) {
        $cmd_mrInchworm .= " $INCHWORM_CUSTOM_PARAMS";
    } 

    $cmd_mrInchworm .= " $sKmerDir/*"; 

    &process_cmd($cmd_mrInchworm);
    &process_cmd("$UTILDIR/combine_clustered_iworm.sh iworm $CPU $inchworm_outfile");

    ####
    &process_cmd("rm $reads_header");
    &process_cmd("rm -r $sKmerDir");

    return;
}

########
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

#    print "CMD: $cmd\n";

#    my $start_time = time();
    my $ret = system($cmd);
#    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
#    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

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

END {
    &collectl_stop();
}
