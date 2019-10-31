#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

$| = 1;

my ($file_input, $file_output, $help, $verbose, $mverbose);

my $start_time = time();

######## options/help ########
GetOptions(
    'help|h'                    => \$help,
    'input|i=s'                 => \$file_input,
    'output|o=s'                => \$file_output,
) or die "\n Invalid commmand line options.\n";

my $USAGE =<<USAGE;
    Usage:
         $0 [-h|--help]
            [-v|--verbose]
            [-vv|--mverbose]
            [-i|--input]
            [-o|--output]

    Options:
          -h,  --help          show this help
          -i,  --input         input file
          -o,  --output        output file

USAGE

## check arguments
if($help) {
    print "$USAGE\n";
    exit 1;
}
if(!$file_input){
    print "  $0: Argument required. \n";
    print "$USAGE\n";
    exit 1
}
if(!$file_output){
    print "  $0: Argument required. \n";
    print "$USAGE\n";
    exit 1
}

open INPUT,  '<',  $file_input  or die "Can't read $file_input: $!";
open OUTPUT, '>', $file_output or die "Can't write to $file_output: $!";

my @sample = ();
my $header = '';
my $cell = '';
my $n = 0;
while( my $current_line = <INPUT>)  {
    if($current_line =~/^EstimateS/){
        ## process previous sample
        if($#sample + 1 > 0){
            ## get last line
            my $last_element = $sample[-1];

            if($n == 1){
                print OUTPUT $header;
            }
            print OUTPUT "$last_element";
        }
        @sample = ();
        $cell = '';
    }

    ## get cell id
    if($current_line =~/^Diversity Output from Input File:  Cell (\w+) \(/){
        $cell = $1;
    }

    ## save line with ^number
    if($current_line =~/^\d+/){
        push(@sample, $cell ."\t".$current_line);
    }

    ## header file
    if($current_line =~/^Samples/){
        $n = $n + 1;
        $header = $current_line;
    }
}

## process previous sample
if($#sample + 1 > 0){
    ## get last line
    my $last_element = $sample[-1];
    print OUTPUT "$last_element";
}

close INPUT;
close OUTPUT;
