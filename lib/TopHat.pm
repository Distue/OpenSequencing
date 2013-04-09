#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package TopHat;
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
        'reference'  => $config->getFileFull("Bowtie2reference"),
        'tophat'     => $config->getCommand("tophat"),
        'stack'      => $refVars->{'stack'},
        'config'     => $refVars->{'config'},
        'init'       => 0
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
 	
	        # add the actual command
		$self->{'stack'}->add(Command->new({
    	    'name'    => "Export Index",
       	 	'command' => "export BOWTIE2_INDEXES='" . $self->{"config"}->getFilePath("Bowtie2reference") . "'",
       	 	'openthreads' => 0,
       		'outputFiles' => [],
        	'inputFiles' => []   
	    }));     
	    
	$self->{'init'} = 1;    
}

sub run {
    my $self = shift;
    my $refVars = shift;
    my $options = $refVars->{'options'};
    $options =~ chomp($options);
	$self->checkInit();
                
    foreach my $sample ($self->getSamples()) {
     	my $commandline = "";
     	
     	$commandline .= $self->{'config'}->getCommand("tophat")   	
				    . " -o " . $self->{'outputDir'} . $sample . "/"
				    . "	-p " . $self->{'config'}->getMaxNumberCores() 
                	. " " . $options
	            	. " " . $self->{'config'}->getFileName("Bowtie2reference");
	            	
        $commandline .= " ";
      	
      	$commandline .= $self->{'seqSamples'}->getDirection1FullFileName($sample) . " ";  	

    	if($self->{'seqSamples'}->isPairedEnd($sample)) {            #
     		$commandline .= $self->{'seqSamples'}->getDirection2FullFileName($sample) . " ";
   		}   
   	
   		$commandline .= " 2> " . $self->{'logDir'} . $sample . ".tophat.log.txt";   
        
        $commandline =~ s/\s+/ /g;
        
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
    	    'name'    => "TopHat",
       	 	'command' => $commandline,
       		'openthreads' => 0,
       		'outputFiles' => [ $self->getBamFileName($sample) ],
        	'inputFiles' => $inputFiles   
	    }));                       
    }                    
         
	$self->{'stack'}->run();   
	                
   	foreach my $sample ($self->getSamples()) {         
		$self->{'stack'}->add(Command->new({
    	    'name'    => "Create symbolic links for TopHat output files",
       	 	'command' => "ln -s " . $self->{'outputDir'} . $sample . "/accepted_hits.bam " . $self->{'outputDir'} . $sample . ".tophat.bam",
       		'openthreads' => 1,
       		'outputFiles' => [ $self->{'outputDir'} .$sample . ".tophat.bam" ],
        	'inputFiles' => [ $self->{'outputDir'} . $sample . "/accepted_hits.bam" ]  
	    }));     
	}               
                    
    $self->{'stack'}->run();                
    
    return($self->getBamfileHandler());
}

sub getSamples {
    my $self = shift;
    return($self->{'seqSamples'}->getSamples());
}


# return the SamfileHandler with all the samples

sub getBamfileHandler {
    my $self = shift;
    
    my $bamfileHandler = BamfileHandler->new({
                    'stack' => $self->{'stack'},
                    'config' => $self->{'config'}
                });
    
    foreach my $sample ($self->getSamples()) {
        $bamfileHandler->addFile($sample, $self->getBamFileName($sample));
    }
    
    return ($bamfileHandler);
}


sub getBamFileName {
    my $self = shift;
    my $sample = shift;
    
    return($self->{'outputDir'} . $sample . ".tophat.bam");
}

1;

 
