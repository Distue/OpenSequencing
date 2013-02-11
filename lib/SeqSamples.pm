#!/usr/bin/perl -w
# -----------------------------------------------------------------------------
# Author: Thomas Schwarzl
# Email:  thomas@schwarzl.net
# Info: TODO
# DevInfo: TODO
# -----------------------------------------------------------------------------
# Description: TODO
# -----------------------------------------------------------------------------

package SeqSamples;
use strict;
use warnings;
use Carp qw(croak);

# -----------------------------------------------------------------------------
# Class Methods
# -----------------------------------------------------------------------------

# constructor - is called when object is created

sub new {
    my $class = shift;
    my $refVars = shift;
  
    # TODO check everything
    # stack
 
    checkDir(\$refVars->{'inputDir'});
    checkOrganism($refVars->{'organism'});

    my $self = {  # this variable stores the varibles from the object
        'inputDir' => $refVars->{'inputDir'},
        'organism' => $refVars->{'organism'},
        'samples'  => {},
        'type'     => $refVars->{'type'},
        'stack'    => $refVars->{'stack'}
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

sub checkOrganism {
    my $organism = shift; 
    ($organism =~ /human|mouse/) or croak ("'" . $organism . "' is not recognized");   
}


# -----------------------------------------------------------------------------
# Instance Methods
# -----------------------------------------------------------------------------

sub checkInit {
   my $self = shift;
   ($self->{'init'}) or $self->init();
}

sub checkSample {
    my $self = shift;
    my $sample = shift;
    
    defined($sample) or croak('Function call need sample name as input'); 
}

# initialization
sub init {
    my $self = shift;
    
    my $inputDir = $self->{'inputDir'};
    
    
    opendir(IMD, $inputDir) or croak("Cannot open directory " . $inputDir);
    
    my @theDirs = readdir(IMD);
    
    foreach my $dir (@theDirs) {
       if ($dir =~ "Sample_(.+)") {
        
          #print $dir . " - " .  $1 . "\n";
          my $sample = $1;
          
          # read subdirectory
          my $subDir = $inputDir  . "/" . $dir . "/";
          
          opendir (IMD2, $subDir) or croak("Cannot open directory " . $subDir) ;
               my @theFiles = readdir(IMD2);
               
                foreach my $file (@theFiles) {
                    #print $file . "\n";
                    $self->addFile($sample, $subDir, $file);
                }
                
          closedir(IMD2);
       }
    }
    
   closedir(IMD);
   
   $self->{'init'} = 1;
   $self->run(); 
}

# run the stack

sub run { 
    my $self = shift;
    $self->checkInit();
    $self->{'stack'}->run();
}

# adds file


##################################################################
# at the moment only for fastq! TODO for all different file types, generic etc.
##################################################################

sub addFile {
    my $self = shift;
    my $sample = shift;
    my $subDir = shift; 
    my $fileName = shift;

    # has ending fastq or fastq.gz
    if($fileName =~ /(\.fastq)(\.gz)?$/) {
        
        # if there are zipped files, unzip them if not already unzipped.
        my $unzippedFileName = $fileName; 
        if($fileName =~ /(.+)\.gz/) {
            $unzippedFileName = $1;
            unless(-e $unzippedFileName) {
		 $self->{'stack'}->add(Command->new( {
                        'command' => "gunzip " . $subDir . $fileName,
                        'openthreads'=> 1,
                        'inputFiles' =>  [ $subDir . $fileName ],
                        'name'      => "Gunzipping compressed file"
                }));
            }
        }
       
        # get read direction 1
        if ($unzippedFileName =~ /_R1_/) {
            $self->{'samples'}->{$sample}->{1} = $subDir . $unzippedFileName;
        }
     
        # get read direction 2
        elsif ($unzippedFileName =~ /_R2_/) { 
            $self->{'samples'}->{$sample}->{2} = $subDir . $unzippedFileName;
        }
    }
}         

# returns a list with all the samples

sub getSamples {
    my $self = shift;
    $self->checkInit();

    return(sort keys %{$self->{'samples'}});
}

# prints sample list 

sub printSamples {
    my $self = shift;
    $self->checkInit();
     
    print "[Print Samples]\n";
    foreach my $sample ($self->getSamples()) {
        print "- " . $sample . "\n";
    }
}

# print the sample count

sub printSampleCount {
    my $self = shift;
    $self->checkInit();

    print $self->getSampleCount() . "\n";
}

sub printSampleDetails {
    my $self = shift;
    $self->checkInit();

    # TODO
}

# returns the number of samples

sub getSampleCount {
    my $self = shift;
    $self->checkInit();

    my @samples = $self->getSamples();
    return(scalar(@samples));
}

# returns the number of direction files

sub getDirectionCount {
    my $self = shift;
    $self->checkInit();
 
    my $sample = shift;
    $self->checkSample($sample);
    
    return(scalar(keys %{$self->{'samples'}->{$sample}}));
}

# returns 1 if sample is single read, else 0

sub isSingleRead {
    my $self = shift;
    $self->checkInit();
 
    my $sample = shift;
    $self->checkSample($sample);
    
    # tests if there is only one read direction
    if ($self->getDirectionCount($sample) == 1) {
        return(1)
    }
    else {
        return(0)
    }
}

# returns 1 if sample is paired end, else 0

sub isPairedEnd {
    my $self = shift;
    $self->checkInit();

    my $sample = shift;
    $self->checkSample($sample);
    
    # tests if there are two read direction files
    if ($self->getDirectionCount($sample) == 2) {
        return(1)
    }
    else {
        return(0)
    }
}

# returns file name of direction 1
# TODO delete
sub getDirection1FileName {
    my $self = shift;
    $self->checkInit();

    my $sample = shift;
    $self->checkSample($sample);
    
    return($self->{'samples'}->{$sample}->{1})
}

# returns full file name of direction 1

sub getDirection1FullFileName {
    my $self = shift;
    $self->checkInit();

    my $sample = shift;
    $self->checkSample($sample);
    
    return($self->{'samples'}->{$sample}->{1});
}


# returns file name of direction 2
# TODO Delete
sub getDirection2FileName {
    my $self = shift;
    $self->checkInit();
 
    my $sample = shift;
    $self->checkSample($sample);
    
    return($self->{'samples'}->{$sample}->{2})
}


# returns full file name of direction 1

sub getDirection2FullFileName {
    my $self = shift;
    $self->checkInit();

    my $sample = shift;
    $self->checkSample($sample);
    
    return($self->{'samples'}->{$sample}->{2});
}


# get organism

sub getOrganism {
    my $self = shift;
    $self->checkInit();
 
    return ($self->{'organism'});
}

1;

 
