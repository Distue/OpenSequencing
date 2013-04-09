#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package Command;
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
    
    # TODO CHECK EVERYTHING
    exists($refVars->{'openthreads'})  or $refVars->{'openthreads'} = 0;
    exists($refVars->{'name'})         or $refVars->{'name'}        = "Unspecified Command";
    defined($refVars->{'outputFiles'}) or $refVars->{'outputFiles'} = [];
    exists($refVars->{'inputFiles'})   or $refVars->{'inputFiles'}  = [];

   
    my $self = {  # this variable stores the varibles from the object
        'name'              => $refVars->{'name'},
        'command'           => $refVars->{'command'},
        'dependentCommands' => $refVars->{'dependentCommands'}, # TODO
        'openthreads'       => $refVars->{'openthreads'},
        'outputFiles'       => $refVars->{'outputFiles'},
        'inputFiles'        => $refVars->{'inputFiles'}
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


# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------

sub checkCommand {
    my $command = shift;
    # TODO
    
}


# get all commands as well as dependent commands

sub get {
    my $self = shift;
    
    my $commandline = "";
    
    foreach my $dependency (@{$self->{'dependentCommands'}}) {
        $commandline .= $dependency->get() . "\n";        
    }
    
    $commandline .= $self->{'command'};
    
    return $commandline;
}


# get all commands as well as dependent commands

sub getCommand {
    my $self = shift;
    
    my $commandline = "";
    
    foreach my $dependency (@{$self->{'dependentCommands'}}) {
        $commandline .= $dependency->getCommand() . "\n";        
    }
    
    # only add command, if output files do not exist
    unless($self->outputFilesExist()) {
            $commandline .= $self->{'command'};
    }
    
    return $commandline;
}

sub outputFilesExist {
    my $self  = shift;
    
    # if there is no output file specified, then run the command
    if (scalar(@{$self->{'outputFiles'}}) == 0) {
        return 0;
    }
    
    # else check if the files are specified
    else {
        my $exist = 1;
        foreach my $file (@{$self->{'outputFiles'}}) {
            # if the file does not exist or the file has not non zero file size
            if(!-e $file || !-s $file) {
                $exist = 0;
            }
        }
        
        return($exist);
    }
}


sub set {
    my $self = shift;
    
    my $command = shift;
    
    $self->checkCommand($command);
}

sub clearOutputFiles {
    my $self = shift;
    
    # TODO
}

sub checkInputFiles {
    my $self = shift;
    
    # TODO
}

sub getName {
    my $self = shift;
    
    if (exists($self->{'name'})) {
        return($self->{'name'});
    }
    else {
        return("Unspecified command");
    }
}

sub print {
    my $self = shift;
    # ERROR?
    print $self->{'command'};
}

sub okayForOpenThreads {
    my $self = shift;
    
    return($self->{'openthreads'} == 1);
}



1;
