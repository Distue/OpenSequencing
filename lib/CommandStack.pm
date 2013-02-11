#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package CommandStack;
use strict;
use warnings;
use Carp qw(croak);
use Command;
use Parallel::ForkManager;


# TODO : TEST
use Data::Dumper;                

my $pm = new Parallel::ForkManager(40);
my @stack;

# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    
    # command 
    my $refVars = shift;
    
    my $self = {  # this variable stores the varibles from the object
        'stack' => [],
        'debug' => 0
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

# adds a command to the stack

sub add {
    my $self = shift;
    my $command = shift;
    
    # print "[ADDING COMMAND] " . $command->get()  . "\n"; 

    push(@{$self->{'stack'}}, $command);
}


# runs the commands in the stack

sub run {
    my $self = shift;
                
    # process the commands or if debug (any number of arguments given) print commands
    while (my $command = shift @{$self->{'stack'}}) {
        
        if(defined($command)) {
            # debug
            if ($self->{'debug'} == 1) {          
                print $command->get() . "\n\n";
            } 
            else {
                if($command->okayForOpenThreads()) { 
                    # Forks and returns the pid for the child:
                    my $pid = $pm->start and next; 
                    
                    $self->process($command); 
              
                    $pm->finish; # Terminates the child process
                }
                else { 
                    # have to clear PATH variable if working on shared directory in linux
                    # $ENV{PATH}=""; 
                    $self->process($command); 
                }
            }
        }
    }
    
    # wait for the threaded subprocesses to finish
    $pm->wait_all_children();
}


sub process {
    my $self = shift;
    my $command = shift;
    
    my $commandLine = $command->getCommand();
    
    if ($commandLine =~ /^$/) {
        print " [ SKIPPING  '" . $command->getName() . "']\n " . $command->get() . "\n";
    }
    else {
        print "[ EXECUTING '" . $command->getName() . "']\n" . $command->getCommand() . "\n";
        
        
        if (system($commandLine) != 0) {
            die("Error while proccessing " . "\n " . $command->get() . "\n"); 
        }   
    }

}

# turns on debugging

sub setDebugOn{
    my $self = shift;
    
    $self->setDebug(1);
}

# turns off debuggin

sub setDebugOff{
    my $self = shift;
    
    $self->setDebug(0);
}

# sets debug status

sub setDebug {
    my $self = shift;
    $self->{'debug'} = shift;
}

1;
