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
        'commands' => {},
        'options'  => {},
        'project'  => "",
        'debug'    => 0
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
    # #print "checking Dir " . $$refDir . "\n";
    (-d $$refDir) or croak ($$refDir . " cannot be opened"); 
}

sub checkFile {
    my $file = shift;
    # #print "checking Dir " . $$refDir . "\n";
    checkDefined($file);
    (-e $file) or croak ($file . " cannot be opened"); 
}

sub checkDefined {
    my $var = shift;
    # #print "checking Dir " . $$refDir . "\n";
    defined($var) or croak ("Essential input not defined."); 
}

sub checkTag {
	my $hash = shift;
	my $tag = shift;
	
	#print "CheckTag: $hash || $tag\n\n";
	defined($hash->{$tag}) or croak ("Requested " . $tag . " does not exist."); 
}

# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------


sub setProject {
	my $self = shift;
	my $name = shift;
	$self->{'project'} = $name;
}

sub getProject {
	my $self = shift;
	
	return($self->{'project'});
}

sub setCommand {
    my $self = shift;
    my $name = shift;
    my $command = shift;
    $self->{'commands'}->{$name} = $command;  
}

sub setOption {
    my $self = shift;
    my $name = shift;
    my $option = shift;
    $self->{'options'}->{$name} = $option;  
}

sub getCommand {
    my $self = shift;
    my $tag = shift; 
    
    checkTag($self->{'commands'}, $tag);
    
    return($self->{'commands'}->{$tag});   
}

sub getOption {
    my $self = shift;
    my $tag = shift; 
    
    checkTag($self->{'options'}, $tag);
    
    return($self->{'options'}->{$tag});   
}


sub setPath {
    my $self = shift;
    my $name = shift;
    my $input = shift;
    
    # if there is no trailing /, add one
    if($input !~ /\/$/) {
    	$input .= "/";	
    }
    
	unless (-d $input) {
        mkdir($input) or croak ("Could not create directory '" . $input . "'");
    }
    
    $self->{'path'}->{$name} = $input;
}
        
sub getPath {
    my $self = shift;
    my $tag = shift;
    
    checkTag($self->{'path'}, $tag);
    
    return($self->{'path'}->{$tag});
}

sub setFile {
    my $self = shift;
    my $tag = shift;
    my $filename = shift;
    my $path = $_[0];
    
    #print "setFile: $tag || $filename || $path\n\n";
    
    $self->setFileName($tag, $filename);
    $self->setFilePath($tag, $path);
    
    checkDefined($self->getFileFull($tag));
}

sub setFileName {
	my $self = shift;
	my $tag = shift;
	my $filename = shift;
	
	$self->{'file'}{'name'}{$tag} = $filename;
}

# path is optional 
sub setFilePath {
    my $self = shift;
    my $tag = shift;
    my $path = $_[0];
    
    #print "setFilePath: $tag\n " . $path . "\n";
    
    if(defined($path)) {
        if($path !~ /\/$/) {
    		$path .= "/";	
    	}
        checkDir(\$path);
        $self->{'file'}->{'path'}->{$tag} = $path;
    }
}

sub getFileName {
    my $self = shift;
    my $tag = shift;
    
    #print "getFileName: $tag\n\n";
    checkTag($self->{'file'}->{'name'}, $tag);
    #print "checkTag done\n";
    return($self->{'file'}->{'name'}->{$tag});
}

sub getFileFull {
    my $self = shift;
    my $tag = shift;
    
    return($self->getFilePath($tag) . $self->getFileName($tag));
}

sub getFilePath {
    my $self = shift;
    my $tag = shift;
    
    #print "getFilePath $tag \n";
    checkTag($self->{'file'}->{'path'}, $tag);
    #print "checkTag Done \n";
    
    if(exists($self->{'file'}->{'path'}->{$tag}) && defined($self->{'file'}->{'path'}->{$tag})) {
        return($self->{'file'}->{'path'}->{$tag})
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


# --------------
# debug options
# --------------

sub isDebugOn {
    my $self = shift;
    return($self->{'debug'} == 1);
}

sub isDebugOff {
    my $self = shift;
    return(!$self->isDebugOn());
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
