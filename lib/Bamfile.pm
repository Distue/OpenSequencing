#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package Bamfile;
use strict;
use warnings;
use Carp qw(croak);

# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    
    # command 
    my $refVars = shift;
    
    checkFile($refVars->{'fullFileName'});

    my $self = {  # this variable stores the varibles from the object
        'fullFileName' => $refVars->{'fullFileName'},
        'samtools'     => $refVars->{'samtools'}
    };
    
    bless ($self, $class);
    return $self;
}


# destructor - is called when object is destroyed

sub DESTROY {
    my $self = shift;
    
    # check for an overridden destructor
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}

# -----------------------------------------------------------------------------
# Static Methods & Variables
# -----------------------------------------------------------------------------

sub checkFile {
    my $fullFileName = shift;
    
# TODO: check if not debug! 
#    # check if file is accessable
#    (-e $fullFileName) or croak("Cannot access Bamfile, or file does not exist '" . $fullFileName . "'");
#    
#    # if file name has 0 bytes, delete it
#    if (-s $fullFileName == 0) {
#        croak("Bam file is empty: '" . $fullFileName . "'");   
#    }
}


# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------

# setter
sub set {
    my $self = shift;
    
    my $fullFileName = shift;
    checkSamFile($fullFileName);
}


sub get {
    my $self = shift;
    
    return ($self->{'fullFileName'});
}


1;
