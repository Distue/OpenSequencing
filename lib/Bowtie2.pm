#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package Bowtie2;
use strict;
use warnings;
use Carp qw(croak);
use SeqSamples;
use Command;
use CommandStack;
use BamfileHandler;
use OpenSeqConfig;

# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    
    # command 
    my $refVars = shift;
    
    my $config = $refVars->{'config'};     
    checkDefined($refVars->{'stack'});
    checkDefined($refVars->{'seqSamples'});
    checkDefined($config->getPath('log'));

    $refVars->{'logDir'} = $config->getPath('log');

    

    my $self = {  # this variable stores the varibles from the object
        'seqSamples' => $refVars->{'seqSamples'},
        'outputDir'  => $config->getPath("output"),
        'logDir'     => $refVars->{'config'}->getPath("log"),
        'reference'  => $config->getFileName("Bowtie2reference"),
        'command'    => $config->getCommand("bowtie2"),
        'stack'      => $refVars->{'stack'},
        'config'     => $refVars->{'config'},
        'init'       => 0,
        'cores'      => $config->getMaxNumberCores() 
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

sub checkDefined {
    my $var = shift;
    # print "checking Dir " . $$refDir . "\n";
    defined($var) or croak ("Essential input not defined."); 
}

# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------

sub checkInit {
   my $self = shift;
   ($self->{'init'} == 1) or $self->init();
}

sub init {
 	my $self = shift;
 	
#	        # add the actual command
#		$self->{'stack'}->add(Command->new({
#    	    'name'    => "Export Index",
#       	 	'command' => "export BOWTIE2_INDEXES='" . $self->{"config"}->getFilePath("Bowtie2reference") . "'",
#       	 	'openthreads' => 0,
#       		'outputFiles' => [],
#        	'inputFiles' => []   
#	    }));     
	    
	$self->{'init'} = 1;    
}

sub run {
    my $self = shift;
    my $refVars = shift;
    my $options = $refVars->{'options'};
    $self->checkInit();
                
    foreach my $sample ($self->getSamples()) {
    	
    	my $commandline = "export BOWTIE2_INDEXES='" . $self->{"config"}->getFilePath("Bowtie2reference") . "' ";
    	     
     	$commandline .= $self->{'command'}
     	            . " " . $options
     	            . "	--threads " . $self->{'cores'}
	            	. " -x " . $self->{'reference'}  
            		. " ";
      	if($self->{'seqSamples'}->isPairedEnd($sample)) {     
      		$commandline .= " -1 " . $self->{'seqSamples'}->getDirection1FullFileName($sample) . " ";  	
     		$commandline .= " -2 " . $self->{'seqSamples'}->getDirection2FullFileName($sample) . " ";
   		}   
   		else {
   			$commandline .= " -U " . $self->{'seqSamples'}->getDirection1FullFileName($sample) . " ";  	
   		}
   		
   		# define output
   		$commandline .= " -S " . $self->getSamFileName($sample);
   	
   		# log
   		$commandline .= " 2> " . $self->{'logDir'} . $sample . ".bowtie2.log.txt";   
        
        
        
        # create input file array              
        my $inputFiles;
        
        if($self->{'seqSamples'}->isPairedEnd($sample)) { 
        	$inputFiles = [ $self->{'seqSamples'}->getDirection1FullFileName($sample),
        					$self->{'seqSamples'}->getDirection2FullFileName($sample) ]
        }
        else {
        	$inputFiles = [ $self->{'seqSamples'}->getDirection1FullFileName($sample) ]
        }

        				    
        # add the actual command
		$self->{'stack'}->add(Command->new({
    	    'name'    => "Bowtie2",
       	 	'command' => $commandline,
       		'openthreads' => 0,
       		'outputFiles' => [ $self->getSamFileName($sample) ],
        	'inputFiles' => $inputFiles   
	    }));                       
    }                    
                        
                        
                        
    $self->{'stack'}->run();                    
                        
    
    return($self->getSamfilesHandler());
}

sub getSamples {
    my $self = shift;
    return($self->{'seqSamples'}->getSamples());
}


# return the SamfileHandler with all the samples

sub getSamfilesHandler {
    my $self = shift;
    
    my $samfileHandler = SamfileHandler->new({
                    'stack' => $self->{'stack'},
                    'config' => $self->{'config'}
                });
    
    foreach my $sample ($self->getSamples()) {
        $samfileHandler->addFile($sample, $self->getSamFileName($sample));
    }
    
    return ($samfileHandler);
}


sub getSamFileName {
    my $self = shift;
    my $sample = shift;
    
    return($self->{'outputDir'} . $sample . ".bowtie2.sam");
}

1;

 
