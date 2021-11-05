#!/usr/bin/perl

#===================================================================#
#
# replaceBC: replace boundary conditions in OpenFOAM files
# Tim MJ Nijssen - Oktober 2019
#
# called as:
#    perl replaceBC.pl oPath rPath nPath
#
# oPath: path to original file containing BCs
# rPath: path to file containing BCs to be substituted in oFile
# nPath: path where the updated file will be written to
#
# explaination of variable names:
# o refers to the original file
# r refers to the replacement file
# n refers to the new file
#
#===================================================================#

use strict;
use warnings;

# read input from commandline
my $oPath = $ARGV[0];
my $rPath = $ARGV[1];
my $nPath = $ARGV[2];

# throw error if one of the files is missing
if ($oPath eq '') {
    die "No path to original file supplied";
}
if ($rPath eq '') {
    die "No path to replacement file supplied";
}
if ($nPath eq '') {
    die "No path to new file supplied";
}

# read file into array
open(my $input1, "<", $oPath) or die "Can't open original file";
my @oFile = <$input1>;

open(my $input2, "<", $rPath) or die "Can't open replacement file";
my @rFile = <$input2>;

open(my $output, ">", $nPath) or die "Can't open new file";
my @nFile = @oFile;

my $rStartLine = -1; # line index patch starts at (-1 is not found)
my $rEndLine   = -1; # line index patch ends  at (-1 is not found)
my $rPatchName = ""; # name of patch
my $rMultiComment  = 0; # inside multi-line comment
my $rCodeBlock     = 0; # inside code block

for (my $i=0; $i<@rFile; $i++) {
    my $rLine  = $rFile[$i]; # current line
    my $rLineCompress = join('',split(" ",$rLine)); # remove all whitespace

    # remove single-line comments
    if ($rLineCompress=~'//') {
	$rLineCompress = substr($rLineCompress, 0, index($rLineCompress, '//'));
    }
    
    # remove multi-line comments
    if (!$rMultiComment && $rLineCompress=~'/[*]') {
	# get part before /*
	$rLineCompress = substr($rLineCompress, 0, index($rLineCompress, '/[*]')+1); 
	$rMultiComment = 1;
    }
    if ($rMultiComment && $rLineCompress=~'[*]/') {
	 # get part after */
        (undef, $rLineCompress) = split('[*]/',$rLineCompress);
	$rMultiComment = 0;
    }

    # remove code blocks
    if (!$rCodeBlock && $rLineCompress=~'#\{') {
	# get part before #{
	$rLineCompress = substr($rLineCompress, 0, index($rLineCompress, '#\{')+1); 
	$rCodeBlock = 1;
    }
    if ($rCodeBlock && $rLineCompress=~'#\};') {
	 # get part after #};
        (undef, $rLineCompress) = split('#\};',$rLineCompress);
	$rCodeBlock = 0;
    }

    if (!($rMultiComment|$rCodeBlock)) {

	# if no patch start found, look for one
	if ($rStartLine==-1) {
	    if (!length($rLineCompress)==0) { # if line is not empty
		$rStartLine = $i;
		$rEndLine   = -1;
		$rPatchName = $rLineCompress;
		$rPatchName =~ s/{//; # remove { from patch name
		print("Patch $rPatchName starts at line $rStartLine of replacement file\n");
	    }
	}
	
	# if patch start is found, look for end
	else {
	    if ($rLineCompress=~/}/) { # if line contains }	    
		$rEndLine   = $i;
		print("Patch $rPatchName ends at line $rEndLine of replacement file\n");

		# find patch in original file

		my $oStartLine = -1; # line index patch starts at (-1 is not found)
		my $oEndLine   = -1; # line index patch ends  at (-1 is not found)

		my $oMultiComment  = 0; # inside multi-line comment
		my $oCodeBlock     = 0; # inside code block

		for (my $j=0; $j<@nFile; $j++) { # nFile is the latest version of oFile
		    my $oLine  = $nFile[$j]; # current line
		    my $oLineCompress = join('',split(" ",$oLine)); # remove all whitespace

		    # remove single-line comments
		    if ($oLineCompress=~'//') {
			$oLineCompress = substr($rLineCompress, 0, index($rLineCompress, '//'));
		    } 

		    # remove multi-line comments
		    if (!$oMultiComment && $oLineCompress=~'/[*]') {
			# get part before /*
			$oLineCompress = substr($oLineCompress, 0, index($oLineCompress, '/[*]')+1); 
			$oMultiComment = 1;
		    }
		    if ($oMultiComment && $oLineCompress=~'[*]/') {
			# get part after */
			(undef, $oLineCompress) = split('[*]/',$oLineCompress);
			$oMultiComment = 0;
		    }

		    # remove code blocks
		    if (!$oCodeBlock && $oLineCompress=~'#\{') {
			# get part before #{
			$oLineCompress = substr($oLineCompress, 0, index($oLineCompress, '#\{')+1); 
			$oCodeBlock = 1;
		    }
		    if ($oCodeBlock && $oLineCompress=~'#\};') {
			# get part after #};
			(undef, $oLineCompress) = split('#\};',$oLineCompress);
			$oCodeBlock = 0;
		    }

		    if (!($oMultiComment|$oCodeBlock)) {

			# if no matching patch is found, look for one
			if ($oStartLine==-1) {
			    if ($oLine=~/$rPatchName/) {
				$oStartLine = $j;
				print("Patch $rPatchName starts at line $oStartLine of original file\n");
			    }
			    
			    # if end of file is reached without finding matching patch, throw warning
			    if ($j==(@nFile-1)) {
				print("WARNING: Patch $rPatchName not found in original file\n");
			    }
			}

			# if matching patch is found, look for end
			else {
			    if ($oLineCompress=~/}/) {
				$oEndLine = $j;
				print("Patch $rPatchName ends at line $oEndLine of original file\n");
				print("replaceBC: replacing patch $rPatchName\n");
				
				# replace patch in original file with patch in replacement file
				
				# start with original part before patch
				my @tempFile = @nFile[0..($oStartLine-1)];
				# add replacement patch
				push(@tempFile,@rFile[($rStartLine)..($rEndLine)]);
				# add original part after patch
				push(@tempFile,@nFile[($oEndLine+1)..(@nFile-1)]);
				# update nFile
				@nFile = @tempFile;

				last; # break loop, return to scanning replacement file for next patch
			    }
			    else {
				# throw error if end of patch is not found
				if ($j==(@nFile-1)) {
				    die "End of patch $rPatchName not found in original file"
				}
			    }
			}
		    }
		}
		# reset for next replacement
		$rPatchName = "";
		$rStartLine = -1;
	    }
	    else {
		# throw error if end of patch is not found
		if ($i==(@rFile-1)) {
		    die "End of patch $rPatchName not found in replacement file"
		}
	    }
	}
    }
}

# print new file
print $output @nFile;

