#!/usr/bin/perl

#===================================================================#
#
# updateControlDict: update OpenFOAM controlDict with new start/endtimes
# Tim MJ Nijssen - Oktober 2019
#
# called as:
#    perl updateControlDict.pl oPath nPath startTime endTime
#
# oPath: path to original constroldict
# nPath: path where the updated controlDict will be written to
# nStartTime: new value of startTime. Leave empty (use ' ') to leave unchanged.
# nEndTime:   new value of endTime.   Leave empty (use ' ') to leave unchanged.
#
# explaination of variable names:
# o refers to the original file
# r refers to the replacement file
# n refers to the new file
#
#===================================================================#

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

# read input from commandline
my $oPath     = $ARGV[0];
my $nPath     = $ARGV[1];
my $nStartTime = $ARGV[2];
my $nEndTime   = $ARGV[3];

my $updateStartTime = 0;
my $updateEndTime   = 0;

# throw error if one of the files is missing
if ($oPath eq '') {
    die "No path to original file supplied";
}
if ($nPath eq '') {
    die "No path to new file supplied";
}

# read startTime input
if ($nStartTime eq ''|| $nStartTime eq ' ') {
    print("updateControlDict: Leaving startTime unchanged\n");
}
else {
    if (looks_like_number($nStartTime)) {
	$updateStartTime = 1;
    }
    else {
	die "startTime is non-numeric";
    }
}

# read endTime input
if ($nEndTime eq ''|| $nEndTime eq ' ') {
    print("updateControlDict: Leaving endTime unchanged\n");
}
else {
    if (looks_like_number($nEndTime)) {
	$updateEndTime = 1;
    }
    else {
	die "endTime is non-numeric";
    }
}

# read file into array
open(my $input1, "<", $oPath) or die "Can't open original file";
my @oFile = <$input1>;

open(my $output, ">", $nPath) or die "Can't open new file";
my @nFile = @oFile;

for (my $i=0; $i<@nFile; $i++) {
    my $line  = $nFile[$i]; # current line
    my @words = split(' ',$line);

    if (@words == 2) {

	if ($updateStartTime && $words[0]=~'startTime') {
	    my $oStartTime = $words[1];
	    $oStartTime =~ s/;//; # remove ; from number

	    # throw error if original startTime is non-numeric
	    if (!looks_like_number($oStartTime)) {
		die "startTime in original file is non-numeric";
	    }

	    # replace number
	    print("updateControlDict: Updating startTime from $oStartTime to $nStartTime\n");
	    $line =~ s/$oStartTime/$nStartTime/;

	    # write new line to array
	    $nFile[$i] = $line;
	}
	if ($updateEndTime && $words[0]=~'endTime') {
	    my $oEndTime = $words[1];
	    $oEndTime =~ s/;//; # remove ; from number

	    # throw error if original endTime is non-numeric
	    if (!looks_like_number($oEndTime)) {
		die "startTime in original file is non-numeric";
	    }

	    # replace number
	    print("updateControlDict: Updating endTime from $oEndTime to $nEndTime\n");
	    $line =~ s/$oEndTime/$nEndTime/;

	    # write new line to array
	    $nFile[$i] = $line;
	}
    }
}

# print new file
print $output @nFile;
