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
use BWA;
use SamfileHandler;
use BamfileHandler;
use config;

# ---------------------------------------------------------
# Config
# ---------------------------------------------------------

my $config = config->getConfig();
$config->setPath("base",            "/share/massstorage/projects/mcf10a-diff/");
$config->setPath("input",           $config->getPath("base")   . "raw/");
$config->setPath("output",          $config->getPath("base")   . "alignment/run1/");
$config->setPath("cufflinksOutput", $config->getPath("output") . "cufflinks202");
$config->setPath("log",             $config->getPath("output") . "log");
$config->setMaxNumberCores(5);


# ---------------------------------------------------------
# Initiate stack
# ---------------------------------------------------------

print "Start Command Stack \n";

my $commandstack = CommandStack->new({'config' => $config});

# $commandstack->setDebugOn();

# ---------------------------------------------------------
# Read samples
# ---------------------------------------------------------

print "Seq Samples \n";
my $seqSamples = SeqSamples->new( {
                       'organism' => "human",
                       'inputDir' => $config->getPath("input"),
                       'format'   => "fastq",
                       'stack'    => $commandstack } );
$seqSamples->init();
$seqSamples->printSampleCount();
$seqSamples->printSamples(); 


# ---------------------------------------------------------
# Align
# ---------------------------------------------------------

print "Start BWA \n"; 

my $bwa = BWA->new(  {  # this variable stores the varibles from the object
               'seqSamples' => $seqSamples,
               'config'     => $config,
               'stack'      => $commandstack
           } );

my $samfileHandler = $bwa->run();

# ---------------------------------------------------------
# Preparation
# ---------------------------------------------------------

my $bamfileHandler = $samfileHandler->toBamfiles();

$bamfileHandler->sortAndIndex();
$bamfileHandler->flagstat();
$bamfileHandler->sortedToSam();

# ---------------------------------------------------------
# Cufflinks
# ---------------------------------------------------------

$bamfileHandler->cufflinks({ 'config' => $config,
                             'gtf'    => $config->getFileFull("gtfEnsembl"),
                             'output' => $config->getPath("cufflinksOutput")});

# ---------------------------------------------------------
# HTSeq
# ---------------------------------------------------------

$samfileHandler->htseqCount({'gtf'    => $config->getFileFull("gtfEnsembl"),
                             'output' => $config->getPath("output")  });


