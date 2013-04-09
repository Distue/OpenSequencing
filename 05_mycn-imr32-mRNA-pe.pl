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
use TopHat;
use Bowtie2;

# ---------------------------------------------------------
# Config
# ---------------------------------------------------------

my $config = config->getConfig();

$config->setProject("mycn-imr32-mRNA-pe");
$config->setPath("base",            "/share/massstorage/projects/" . $config->getProject() . "/");
$config->setPath("input",           $config->getPath("base")   . "raw/");
$config->setPath("output",          $config->getPath("base")   . "alignment/run2/");
$config->setPath("cufflinksOutput", $config->getPath("output") . "cufflinks202quant");
$config->setPath("log",             $config->getPath("output") . "log");
$config->setMaxNumberCores(30);

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
# Align
# ---------------------------------------------------------

my $tophat = TopHat->new( {'seqSamples' => $seqSamples,
                            'config'  => $config,
                            'stack'   => $commandstack
                            });

my $bamfileHandler = $tophat->run({'options' =>
								" -G /share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
							 # "‐-GTF " . $config->getFileFull("gtf")
							 # 434 - 2 * 42 - 120 = 230
				             . " --mate-inner-dist 230 "        # This is the expected (mean) inner distance between mate pairs. For, example, for paired end runs with fragments selected at 300bp, where each end is 50bp, you should set -r to be 200. The default is 50bp.
					         . " --mate-std-dev 150"  # The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.
				            });





#print "Start BWA \n"; 

#my $bwa = BWA->new(  {  # this variable stores the varibles from the object
#               'seqSamples' => $seqSamples,
#               'config'     => $config,
#               'stack'      => $commandstack
#           } );

#my $samfileHandler = $bwa->run();

# ---------------------------------------------------------
# Preparation
# ---------------------------------------------------------


$bamfileHandler->sortAndIndex();
$bamfileHandler->flagstat();
my $samfileHandler = $bamfileHandler->sortedToSam();

# ---------------------------------------------------------
# Cufflinks
# ---------------------------------------------------------

$bamfileHandler->cufflinks({ 'config' => $config,
							'options' => "-N",
							'gtf'    => "/share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                            # 'gtf'    => "/share/massstorage/seq-genomes/Homo_sapiens.Ensembl.GRCh37.converted.gtf",
                             #'gtf'   => $config->getFileFull("gtfAnnotation"),
                             'output' => $config->getPath("cufflinksOutput")});

# ---------------------------------------------------------
# HTSeq
# ---------------------------------------------------------

$samfileHandler->htseqCount({#'gtf'    => $config->getFileFull("gtf"),
			    'gtf'    => "/share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf",
                             'output' => $config->getPath("output")  });


