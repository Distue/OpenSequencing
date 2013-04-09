#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package BWA;
use strict;
use warnings;
use Carp qw(croak);
use SeqSamples;
use Command;
use CommandStack;
use SamfileHandler;
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
        'BWAreference'  => $config->getFileFull("BWAreference"),
        'bwa'        => $config->getCommand("bwa"),
        'stack'      => $refVars->{'stack'},
        'config'     => $refVars->{'config'}
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
   # ($self->{'init'}) or $self->init();
}


sub createCommandSingleEnds {
    my $self = shift;
    $self->checkInit();
    
    my $inputFile = shift;
    my $tmpFile = shift;
    my $logFile = shift;
   
    my $command = Command->new({
        'name'    => "BWA single direction",
        'command' => $self->{'bwa'} . " aln "
                            # processors
                            . "-t " . $self->{'config'}->getMaxNumberCores() . " "
                            
                            # reference
                            . $self->{'BWAreference'} . " "
                    
                            # in Illumina 1.3+ read format
                            # . "-I "
                    
                            # file
                            . $inputFile . " "
                    
                            # to out file
                            . "> " . $tmpFile . " "
                    
                            # errorlog to file
                            . "2> " . $logFile,
        'openthreads' => 0,
        'outputFiles' => [ $tmpFile ],
        'inputFiles' => [ $inputFile ]   
        
    });
    

    return($command);
}
 

sub createCommandPairedEnd {
    my $self = shift;
    $self->checkInit();
    
    my $direction1Sai = shift;
    my $direction2Sai = shift;    
    my $direction1FileName = shift;
    my $direction2FileName = shift;
    my $samFileName = shift;
    my $logFile = shift;
    
    my $command = Command->new({
            'name'    => "BWA paired ends",
           'command' =>  $self->{'bwa'} . " sampe " . $self->{'BWAreference'} . " "
                            # first alignment file
                            . $direction1Sai . " "
                            . $direction2Sai . " "
                            . $direction1FileName ." "
                            . $direction2FileName . " "
                            . "> "
                            . $samFileName
                            . " 2> "
                            . $logFile,
        'openthreads' => 1,
        'outputFiles' => [ $samFileName ],
        'inputFiles' => [ $direction1Sai , $direction2Sai, $direction1FileName, $direction2FileName ]   
    });
    
    return($command);
}
                          
sub createCommandSingleRead {
    my $self = shift;
    $self->checkInit();
    
    my $direction1Sai = shift; 
    my $direction1FileName = shift;
    my $samFileName = shift;
    my $logFile = shift;
    
    my $command = Command->new({
            'name'    => "BWA single read",
           'command' =>  $self->{'bwa'} . " samse " . $self->{'BWAreference'} . " "
                            # first alignment file
                            . $direction1Sai . " "
                            . $direction1FileName ." "
                            . "> "
                            . $samFileName
                            . " 2> "
                            . $logFile,
        'openthreads' => 1,
        'outputFiles' => [ $samFileName ],
        'inputFiles' => [ $direction1Sai , $direction1FileName ]   
    });
    
    return($command);
}

sub run {
    my $self = shift;
    $self->checkInit();
    
    # for all samples
    foreach my $sample ($self->getSamples()) {
    
        # print("Processing Sample " . $sample . "\n");
       
        # do single read alignment
        if($self->{'seqSamples'}->isSingleRead($sample)) {
            my $tmpDir1File = $self->{'outputDir'} . $sample . "_1.sai";
                        
            $self->{'stack'}->add(
                $self->createCommandSingleEnds(
                    $self->{'seqSamples'}-> getDirection1FullFileName($sample),     # input file
                    $tmpDir1File,                                                   # tmp dir 
                    $self->{'logDir'} . "log-bwa-" . $sample . "_1.txt"             # log file
                )
            );
               
            # combine those alignments through the 'bwa samse' (single-end) command.
            $self->{'stack'}->add(
                $self->createCommandSingleRead (
                    $tmpDir1File,
                    $self->{'seqSamples'}-> getDirection1FullFileName($sample),     # input file
                    $self->getSamFileName($sample),                                 # sam file name
                    $self->{'logDir'} . "log-bwa-" . $sample . ".txt"               # log file
                )
            );
        }
        # do paired end alignment
        elsif($self->{'seqSamples'}->isPairedEnd($sample)) {
            # print "Doing BWA paired end for ". $sample . "\n";
                        
            my $tmpDir1File = $self->{'outputDir'} . $sample . "_1.sai";
            my $tmpDir2File = $self->{'outputDir'} . $sample . "_2.sai";
            
            # create alignments for both ends
            # direction 1
            $self->{'stack'}->add(
                $self->createCommandSingleEnds(
                    $self->{'seqSamples'}-> getDirection1FullFileName($sample),     # input file
                    $tmpDir1File,                                                   # tmp dir 
                    $self->{'logDir'} . "log-bwa-" . $sample . "_1.txt"             # log file
                )
            );
                                  
            # direction 2
            $self->{'stack'}->add(
                $self->createCommandSingleEnds(
                    $self->{'seqSamples'}-> getDirection2FullFileName($sample),     # input file
                    $tmpDir2File,                                                   # tmp dir 
                    $self->{'logDir'} . "log-bwa-" . $sample . "_2.txt"          # log file
                )
            );
            
            
            
            
            # combine those alignments through the 'bwa sampe' (paired-end) command.
            $self->{'stack'}->add(
                $self->createCommandPairedEnd(
                    $tmpDir1File,
                    $tmpDir2File,
                    $self->{'seqSamples'}-> getDirection1FullFileName($sample),     # input file
                    $self->{'seqSamples'}-> getDirection2FullFileName($sample),     # input file
                    $self->getSamFileName($sample),                                 # sam file name
                    $self->{'logDir'} . "log-bwa-" . $sample . ".txt"               # log file
                )
            );
            
        # input files is neither 1 or 2
        } else {
              croak("Sample '" . $sample . "' has strange amount of input files listed: " . $self->{'seqSamples'}->getDirectionCount());
        }
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
    
    return($self->{'outputDir'} . $sample . ".bwa.sam");
}

1;

 
