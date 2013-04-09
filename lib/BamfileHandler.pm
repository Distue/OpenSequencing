#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package BamfileHandler;
use strict;
use warnings;
use Carp qw(croak);
use Bamfile;
use CommandStack;
use SamfileHandler;
# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    my $refVars = shift;
  
    checkDefined($refVars->{'stack'});
    checkDefined($refVars->{'config'});
    
    my $self = {  # this variable stores the varibles from the object
        'files'    => {} ,                  # keeps all the files
        'samtools' => 'samtools',           # todo organize samtools
        'stack'    => $refVars->{'stack'},  # command stack
        'config'   => $refVars->{'config'}
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

sub addFile {
    my $self = shift;
    
    my $sampleName = shift;
    my $fullFileName = shift;
    
    defined($fullFileName) or croak('not a valid file name');
            
    my $samfile = Bamfile->new({
                        'fullFileName' => $fullFileName,
                        'samtools'     => $self->{'samtools'},
                        'config'       => $self->{'config'}
                    });
    
    $self->{'files'}->{$sampleName} = $samfile;
}

# adds all sam files of directory

sub addDirectory {
    my $self = shift;
     
    my $inputDir = shift;
    
    opendir(IMD, $inputDir) or croak("Cannot open directory " . $inputDir);
    
    my @files = readdir(IMD);
    
    foreach my $file (@files) {
        # print $file . "\n";
        if ($file =~ /(.+)(.bwa)?\.bam$/) {
            if ($1 =~ /sorted/) {
                print "file '" . $file .  "' is sorted \n";
            } 
            else { 
                print "# add file " .$1 . "\n";
                $self->addFile($1, $inputDir . $file);
            }   
        }
    }
    
    close(IMD);
}

sub sortAndIndex {
    my $self = shift;
    
    # get all the bamfile objects
    foreach my $sample ($self->getSamples()) {
        my $bamfileName = $self->getBamfile($sample)->get();
        my $sortedname = $self->getSortedBamfileName($sample);
        my $sortedInput = $sortedname;
        $sortedInput =~ s/\.bam$//;  
        my $logfile = $sortedInput;
        my $outputDir = $self->{'config'}->getPath("output");
        my $inputDir = $self->{'config'}->getPath("log");
        $logfile =~ s/$outputDir/$inputDir/;    
 
        my $command1 = Command->new({
                'name'        => "samtools sort",
                'command'     => "samtools sort " . $bamfileName . " " . $sortedInput . " 2> "  . $logfile . ".sort.txt",
                'inputFiles'  => [ $bamfileName ],
                'outputFiles' => [ $sortedname ],
        });        

        $self->addCommand(
            Command->new( {
                'name'             => "samtools index",
                'command'          => "samtools index " . $sortedname . " 2> "  . $logfile . ".index.txt",
                'inputFiles'       => [ $sortedname ],
                'outputFiles'      => [ $sortedname . ".bai" ],
                'openthreads'      => 1,
                'dependentCommands' => [ $command1 ]
            })             
        );
        
        
    }
    
    # return the samhandler
    
    $self->run();
}


sub sortedToSam {
   my $self = shift;
    
    # get all the bamfile objects
    foreach my $sample ($self->getSamples()) {
        my $sortedBamInput = $self->getSortedBamfileName($sample);
        my $sortedSamOutput = $sortedBamInput;
        $sortedSamOutput =~ s/bam$/sam/;  
        my $logfileName = $sortedSamOutput;
        my $outputDir = $self->{'config'}->getPath("output");
        my $logDir = $self->{'config'}->getPath("log");
        $logfileName =~ s/$outputDir/$logDir/;      

        $self->addCommand(
            Command->new({
                'name'        => "conversion sorted bam to sam",
                'command'     => "samtools view -h -o " . $sortedSamOutput ." " . $sortedBamInput . " 2> " . $logfileName . ".log",
                'inputFiles'  => [ $sortedBamInput ],
                'outputFiles' => [ $sortedSamOutput ],
                'openthreads' => 1
        }));        
    }
    
    # do the process
    $self->run();  
    
    # return the handler
    my $samfileHandler = SamfileHandler->new({'config' => $self->{'config'},
    										  'stack'  => $self->{'stack'} });
    foreach my $sample ($self->getSamples()) {
    	$samfileHandler->addFile($sample, $self->getSortedSamfileName($sample));
    }			
    
    return ($samfileHandler);				  								  
}


sub flagstat {
    my $self = shift;
    
    # get all the bamfile objects
    foreach my $sample ($self->getSamples()) {
        my $sortedname = $self->getSortedBamfileName($sample);
        my $flagstatname = $sortedname;
        $flagstatname =~ s/bam$/flagstat.txt/;
        
        $self->addCommand(
            Command->new( {   
                'name'             => "flagstat",
                'command'          => "samtools flagstat " . $sortedname . " > " . $flagstatname,
                'inputFiles'       => [ $sortedname ],
                'outputFiles'      => [ $flagstatname ],
                'openthreads'      => 1,
            })             
        );       
    }
    
    $self->run();
}

sub getBamfile {
    my $self = shift;
    my $sample = shift;
    return($self->{'files'}->{$sample});
}

sub getBamfileName {
    my $self = shift;
    my $sample = shift;
    return($self->getBamfile($sample)->get());
}


sub getSortedBamfileName {
    my $self = shift;
    my $sample = shift;
    
    my $sortedname = $self->getBamfileName($sample);
    $sortedname =~ s/bam$/sorted.bam/;
    return($sortedname);
}

sub getSortedSamfileName {
    my $self = shift;
    my $sample = shift;
    
    my $sortedname = $self->getSortedBamfileName($sample);
    $sortedname =~ s/bam$/sam/;
    return($sortedname);
}


sub getSamples {
    my $self = shift;
    return(sort keys %{$self->{'files'}});
}

sub run {
    my $self = shift;
    $self->{'stack'}->run();
}

sub addCommand {
    my $self = shift;
    my $command = shift;
    $self->{'stack'}->add($command);
}


sub cufflinks {
    my $self = shift;
    
    my $refVars = shift;
    my $config    = $refVars->{'config'};
    
    my $cufflinks = $config->getCommand("cufflinks");
    my $gtf       = $refVars->{'gtf'};
    my $outputdir = $refVars->{'output'};
    my $options   = $refVars->{'options'};
    
    checkDefined($cufflinks);
    checkDefined($gtf);
    checkDefined($outputdir);
     
    unless(-d $outputdir){
        mkdir $outputdir;
    }
    
    # get all the bamfile objects
    foreach my $sample ($self->getSamples()) {
        my $sortedName = $self->getSortedBamfileName($sample);
		my $command = $cufflinks . " " . $options  . " " .  
                                 "-o " . $outputdir . "/" .
                                 $sample . "/ -p " . $config->getMaxNumberCores() .
                                 " -G " . $gtf . " " . $sortedName .
                                 " 2> " . $self->{'config'}->getPath("log") . "/cufflinks-" . $sample . ".txt";

        $self->addCommand(
            Command->new({
                'name'        => 'Cufflinks',
                'command'     => $command,
                'inputFiles'  => [ $sortedName ],
                'outputFiles' => [ $outputdir . "/" . $sample ],
                'openthreads' => 0,
            })         
        );
    }
    
    $self->run(); 
}


1;

