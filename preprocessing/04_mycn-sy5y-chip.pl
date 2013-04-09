#!/bin/perl
# ---------------------------------------------------------
# @REV: 2
# Thomas Schwarzl <thomas@schwarzl.net>
# ---------------------------------------------------------
# Script for doing BWA alignments on paired-end Illumina
# reads
# ---------------------------------------------------------

use warnings;
use strict;


use lib 'lib';
use SeqSamples;
use CommandStack;

use SamfileHandler;
use BamfileHandler;
use config;

use BWA;
use Bowtie2;

# ---------------------------------------------------------
# Config
# ---------------------------------------------------------

my $config = config->getConfig();
$config->setProject("mycn-sy5y-chip");
$config->setPath("base",            "/share/massstorage/projects/" . $config->getProject() . "/");
$config->setPath("input",           $config->getPath("base")   . "raw/9nov/");
$config->setPath("output",          $config->getPath("base")   . "alignment/9nov-b2/");
$config->setPath("log",             $config->getPath("output") . "log");
$config->setOption("sampleInput", "files");
    
$config->setMaxNumberCores(5);
#$config->setDebugOn();

# ---------------------------------------------------------
# Initiate stack
# ---------------------------------------------------------

print "Start Command Stack \n";

my $commandstack = CommandStack->new({'config' => $config});

# ---------------------------------------------------------
# Read samples
# ---------------------------------------------------------

print "Seq Samples \n";
my $seqSamples = SeqSamples->new( { 'config'   => $config,
                                    'stack'    => $commandstack } );
$seqSamples->init();
$seqSamples->printSampleCount();
$seqSamples->printSamples(); 


# ---------------------------------------------------------
# Bowtie2
# ---------------------------------------------------------

print "Start Bowtie2 \n"; 

my $bowtie2 = Bowtie2->new( {'seqSamples' => $seqSamples,
                            'config'  => $config,
                            'stack'   => $commandstack
                            });

my $samfileHandler_bowtie2 = $bowtie2->run({'options' =>  " "});

my $bamfileHandler_bowtie2 = $samfileHandler_bowtie2->toBamfiles();
$bamfileHandler_bowtie2->sortAndIndex();
$bamfileHandler_bowtie2->flagstat();

# ---------------------------------------------------------
# BWA
# ---------------------------------------------------------

die();
$config->setPath("input",           $config->getPath("base")   . "alignment/9nov-bwa/");

print "Start BWA \n"; 

my $bwa = BWA->new(  {  # this variable stores the varibles from the object
               'seqSamples' => $seqSamples,
               'config'     => $config,
               'stack'      => $commandstack
           } );

my $samfileHandler_bwa = $bwa->run();

my $bamfileHandler_bwa  = $samfileHandler_bwa->toBamfiles();
$bamfileHandler_bwa->sortAndIndex();
$bamfileHandler_bwa->flagstat();

#$bamfileHandler->sortedToSam();




