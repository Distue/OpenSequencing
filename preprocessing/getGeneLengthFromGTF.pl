#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas.schwarzl@ucd.ie
# Info: add chr to the beginning of each line in a  GTF file
# DevInfo: 
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

use warnings;
use strict;             
use Getopt::Long; # options

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------

my $usage = << "USE";

  Usage: 
            getGeneLengthFromGTF.pl <inputin>
           
  Options:                      
                                      
            --help                   show help message

  Examples:
     Output to console:
           cat genes.gtf | convertGTF.pl 
     Output to file (*nix)
           cat genes.gtf | convertGTF.pl > file2.txt
     Append to file (*nix)
           cat genes.gtf | convertGTF.pl >> file2.txt     
USE


# -----------------------------------------------------------------------------
# Arguments and options
# -----------------------------------------------------------------------------

# Options:
my $help = 0;
$help == 0 or die $usage;

# -----------------------------------------------------------------------------
# Assertions
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

my %haGenes;

while(my $line = <>) {  # pipe in
    # split the line
    my @elements = split("\t", $line);
    # i
    my $annotationField = $elements[8];
    my @annotation = split( ";" , $annotationField);
    my $id;
    if( $annotation[0] =~ /(ENSG[0-9]+)/ ) {
        $id = $1;
    }
    else {
        die( "could not find the id for $annotation[0]" );        
    }
    my $length = $elements[4] - $elements[3];

    if(!defined($haGenes{ $id }) ) {
        $haGenes{$id} = $length;
    }
    else {
        $haGenes{ $id } += $length;
    }
}

# print 
while ( my ( $id, $length) = each( %haGenes)){
    print $id . "\t" . $length . "\n"; 
}

# -----------------------------------------------------------------------------
# Close 
# -----------------------------------------------------------------------------

exit;
