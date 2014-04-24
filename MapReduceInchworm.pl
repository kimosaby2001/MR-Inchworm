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



