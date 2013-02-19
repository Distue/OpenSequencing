#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package OpenSeqConfig;
use strict;
use warnings;
use Carp qw(croak);


# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    
    
    my $self = {  # this variable stores the varibles from the object
        'path' => {},  
        'file' => {
                   'name' => {},
                   'path' => {}  
                  },
        'maxCores' => '1',
        'commands' => {}
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

# checks if the directory is valid and adds '/'

sub checkDir {
    my $refDir = shift;
    # print "checking Dir " . $$refDir . "\n";
    (-d $$refDir) or croak ($$refDir . " cannot be opened"); 
}

sub checkFile {
    my $file = shift;
    # print "checking Dir " . $$refDir . "\n";
    checkDefined($file);
    (-e $file) or croak ($file . " cannot be opened"); 
}

sub checkDefined {
    my $var = shift;
    # print "checking Dir " . $$refDir . "\n";
    defined($var) or croak ("Essential input not defined."); 
}

# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------

sub setCommand {
    my $self = shift;
    my $name = shift;
    my $command = shift;
    $self->{'commands'}{$name} = $command;  
}

sub getCommand {
    my $self = shift;
    my $name = shift; 
    return($self->{'commands'}{$name});   
}

sub setPath {
    my $self = shift;
    my $name = shift;
    my $input = shift;
    $self->{'path'}{$name} = $input;
    
    unless (-d $input) {
        mkdir($input) or croak ("Could not create directory '" . $input . "'");
    }
}
        
sub getPath {
    my $self = shift;
    my $name = shift;
    return($self->{'path'}{$name});
}

sub setFile {
    my $self = shift;
    my $tag = shift;
    my $filename = shift;
    my $path = $_[0];
    
    $self->{'file'}{'name'}{$tag} = $filename;
    $self->setFilePath($tag, $path);
    
    checkDefined($self->getFileFull($tag));
}

sub setFilePath {
    my $self = shift;
    my $tag = shift;
    my $path = $_[0];
    
    if(defined($path)) {
        checkDir(\$path);
        $self->{'file'}{'path'}{$tag} = $path;
    }
}

sub getFileName {
    my $self = shift;
    my $tag = shift;
    
    return($self->{'file'}{'name'}{$tag});
}

sub getFileFull {
    my $self = shift;
    my $tag = shift;
    
    return($self->getFilePath($tag) . $self->getFileName($tag));
}

sub getFilePath {
    my $self = shift;
    my $tag = shift;
    
    if(exists($self->{'file'}{'path'}{$tag})) {
        return($self->{'file'}{'path'}{$tag})
    }
    else {
        return("");
    }
}

sub setMaxNumberCores {
    my $self = shift;
    my $cores = shift;
    $self->{'maxCores'} = $cores;
}

sub getMaxNumberCores {
    my $self = shift;
    return $self->{'maxCores'};
}

1;
