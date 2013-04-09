#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package Samfile;
use strict;
use warnings;
use Carp qw(croak);
use Bamfile;
use Command;
use CommandStack;
use OpenSeqConfig;

# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    
    # command 
    my $refVars = shift;
    
    checkDefined($refVars->{'config'});
    
    if ($refVars->{'config'}->isDebugOff()) {
        checkFile($refVars->{'fullFileName'});
    }

    #TODO check variables
    
    my $self = {  # this variable stores the varibles from the object
        'fullFileName' => $refVars->{'fullFileName'},
        'samtools'     => $refVars->{'samtools'},
        'config'        => $refVars->{'config'}
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
    
    (defined($fullFileName)) or croak('filename is not defined');
    
    # check if file is accessable
    (-e $fullFileName) or croak("Cannot access Samfile, or file does not exist '" . $fullFileName . "'");
    
    # if file name has 0 bytes, delete it
    if (-s $fullFileName == 0) {
        croak("Sam file is empty: '" . $fullFileName . "'");   
    }
}

sub checkDefined {
    my $var = shift;
    # print "checking Dir " . $$refDir . "\n";
    defined($var) or croak ("Essential input not defined."); 
}


# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------

# setter

sub setFileName {
    my $self = shift;
    
    my $fullFileName = shift;
    checkSamFile($fullFileName);
}

# getter

sub getFileName {
    my $self = shift;
    
    return ($self->{'fullFileName'});
}



# returns the command for the convertion to bam file

sub toBamfileCommand {
    my $self = shift;
    
    if ($self->headerHasSQLines() == 1) {

        return(Command->new({
            'name'        => 'sam to bamfile',
            'command'     => $self->getSamToolsPath() . " view -bS " . $self->getFileName() . " > " . $self->getBamFileName() . " 2> " . $self->getLogFile(), 
            'inputFiles'  => [ $self->getFileName() ],
            'outputFiles' => [ $self->getBamFileName() ],
            'openthreads' => 1
        }));
    }
    else {
        croak("SQ Headerless Conversion not yet implemented");    
    }
}


# gets the new bam files

sub getBamFile {
    my $self = shift;
    
    # returns a new bam file if succeeded
    return(Bamfile->new({
        'fullFileName' =>  $self->getBamFileName() 
    }))
}
        

# return the sam tools path

sub getSamToolsPath {
    my $self = shift;
    
    return($self->{'samtools'});
}

# returns the BamfileName -

sub getBamFileName {
    my $self = shift;
    
    my $bamFileName = $self->{'fullFileName'};
    $bamFileName =~ s/sam$/bam/;
    return ($bamFileName);
}

sub getShortFileName {
	my $self = shift;
	my $fullname = $self->{'fullFileName'};
	
	if ($fullname =~ /.+\/(.*?)$/) {
		return($1);
	}
	else {
			croak("Getting short name for bam file failed");
	}
}

sub getLogFile {
    my $self = shift;
 
    my $logFileName =  $self->{'config'}->getPath("log") . $self->getShortFileName();
    #$logFileName =~ s/\/.+?$/\/log-/;
    $logFileName .= "-to-bam.txt.log";
    return ($logFileName);   
}


# returns true (1) if sam header has @SQ lines
sub headerHasSQLines {
    my $self = shift;
    
    # TODO!
    return(1);
}


1;

