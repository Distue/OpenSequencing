#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package SamfileHandler;
use strict;
use warnings;
use Carp qw(croak);
use Samfile;
use CommandStack;
use BamfileHandler;

# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    my $refVars = shift;
  
    #TODO check stack
    
    my $self = {  # this variable stores the varibles from the object
        'files' => {},             # keeps all the files
        'samtools' => 'samtools',   # todo organize samtools
        'stack' => $refVars->{'stack'} # command stack
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

sub addFile {
    my $self = shift;
    
    my $sampleName = shift;
    my $fullFileName = shift;
    
    defined($fullFileName) && -e $fullFileName or croak('filename not defined or valid ' . $fullFileName);
    defined($sampleName) or croak('sample name not defined or valid');
    
    my $samfile = Samfile->new({
                        'fullFileName' => $fullFileName,
                        'samtools'     => $self->{'samtools'}
                        # configuration?
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
        if ($file =~ /(.+)(.bwa)?\.sam$/) {
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

sub toBamfiles {
    my $self = shift;
    
    
    # get all the samfile objects
    foreach my $sample (sort keys %{$self->{'files'}}) {
        my $samfile = $self->{'files'}->{$sample};
        
        $self->{'stack'}->add(
            $samfile->toBamfileCommand()
        );
    }
    
    $self->{'stack'}->run();
    
    return($self->getBamfileHandler());
}

sub getBamfileHandler {
    my $self = shift;
    
    my $bamfileHandler = BamfileHandler->new({
                    'stack' => $self->{'stack'}
                });
    
    foreach my $sample ($self->getSamples()) {
        #print $sample . "\n";
        $bamfileHandler->addFile(
                             $sample,
                             $self->getBamfileName($sample)
                         );
    }
    
    return ($bamfileHandler);  
}

sub getSamples {
    my $self = shift;
    
    return(sort keys %{$self->{'files'}});
}


# take sam name and replace sam by bam

sub getBamfileName {
    my $self = shift;
    my $sample = shift;

    exists($self->{'files'}->{$sample}) or croak("file for sample '" . $sample . "' does not exist");
    
    my $bamfileName =  $self->{'files'}->{$sample}->getFileName();
    $bamfileName =~ s/sam$/bam/;

    return($bamfileName);
}

sub run {
    my $self = shift;
    $self->{'stack'}->run();
}



sub sort {
    my $self = shift;
    
    # get all the samfile objects
    foreach my $sample ($self->getSamples()) {
        my $samfileName = $self->getSamfile($sample)->getFileName();
        my $sortedName  = $self->getSortedSamfileName($sample);

        $self->addCommand(
            Command->new({
                'command'     => "sort -s -k 1,1 " . $samfileName . " > " . $sortedName . " 2> " . $sortedName . ".txt.log",
                'inputFiles'  => [ $samfileName ],
                'outputFiles' => [ $sortedName ],
                'openthreads' => 0
            })         
        );
    }
    
    $self->run();
}

sub htseq {
    my $self = shift;
    
    my $gtf = shift;
    
    # get all the samfile objects
    foreach my $sample ($self->getSamples()) {
        my $sortedName = $self->getSortedSamfileName($sample);
        my $outputName = $sortedName . ".htseq.count.txt";
        
        $self->addCommand(
            Command->new({
                'command'     => "htseq-count --stranded=no  " . $sortedName . " " . $gtf . " > " . $outputName . " 2> " . $outputName . ".log",
                'inputFiles'  => [ $sortedName ],
                'outputFiles' => [ $outputName ],
                'openthreads' => 1
            })         
        );
    }
    
    $self->run(); 
}

sub getSamfile {
    my $self = shift;
    my $sample = shift;
    return($self->{'files'}->{$sample});
}

sub getSamfileName {
    my $self = shift;
    my $sample = shift;
    return($self->getSamfile($sample)->getFileName());
}


sub getSortedSamfileName {
    my $self = shift;
    my $sample = shift;
    
    my $sortedname = $self->getSamfileName($sample);
    $sortedname =~ s/sam$/sorted.sam/;
    return($sortedname);
}


sub addCommand {
    my $self = shift;
    my $command = shift;
    $self->{'stack'}->add($command);
}



1;

