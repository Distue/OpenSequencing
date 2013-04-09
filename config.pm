#!/bin/perl
# ---------------------------------------------------------
# Thomas Schwarzl <thomas@schwarzl.net>
# ---------------------------------------------------------

use warnings;
use strict;

package config;

use lib 'lib';
use OpenSeqConfig;


sub getConfig {
    our $config = OpenSeqConfig->new();
    
    $config->setCommand("bwa",     "bwa");
    $config->setCommand("bowtie2", "bowtie2");
    $config->setCommand("tophat",  "tophat");

    $config->setFile("Bowtie2reference",  "genome",        "/share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/");
    $config->setFile("BWAreference",      "genome",        "/share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/");

    $config->setPath("cufflinks",    "/share/massstorage/seq-tools/cufflinks/current/");
    $config->setCommand("cufflinks", $config->getPath("cufflinks") . "cufflinks");
    $config->setCommand("cuffdiff",  $config->getPath("cufflinks") . "cuffdiff");
	$config->setCommand("gunzip",     "gunzip");
    
    $config->setPath("cufflinksOutput", "/share/massstorage/seq-tools/cufflinks/current/");
    $config->setMaxNumberCores(42);
    
    $config->setFile("gtf",  "genes.gtf", "/share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/");
    $config->setFile("gtfAnnotation", "Homo_sapiens.Ensembl.GRCh37.converted.gtf", "/share/massstorage/seq-genomes/");
    
    # parameters
    $config->setOption("sampleInput", "subdirs");
    $config->setOption("direction1_filename_marker", "_R1_");
    $config->setOption("direction2_filename_marker", "_R2_");
    
    return($config);
}


1;





    # just debug:
    # $config->setFile("BWAreference", "hg19bwaidx", "/share/massstorage/seq-tools/bwa/reference/"); 
    # $config->setFile("gtfNCBI",    "Homo_sapiens.NCBI.build37.2.converted.gtf",    $config->getPath("gtf"));
    # $config->setFile("gtfUCSC",    "Homo_sapiens.UCSC.hg19.gtf",                   $config->getPath("gtf"));
   # $config->setFile("gtf",  "Homo_sapiens.UCSC.hg19.gtf", "/share/massstorage/seq-genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/");   
    


#my $bwacommand = "bwa";
#
#my $inputDir = "/share/massstorage/projects/mycn-lines-mRNA-pe/raw/";
#my $outputDir = "/share/massstorage/projects/mycn-lines-mRNA-pe/alignment/run1/";
#
#my $referenceDir  = "/share/massstorage/seq-tools/bwa/reference/";
#my $referenceName = "hg19bwaidx";
#
#
##my $gtfFile = "/share/massstorage/seq-genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf";
#my $gtfFile = "/share/massstorage/seq-genomes/Homo_sapiens.Ensembl.GRCh37.65.converted.gtf";
#my $cufflinksDir = "/share/massstorage/seq-tools/cufflinks/cufflinks-2.0.2.Linux_x86_64/";
#my $cufflinks = $cufflinksDir . "cufflinks";
#my $cuffdiff  = $cufflinksDir . "cuffdiff";
#my $cufflinksOutputDir = $outputDir . "cufflinks202";
#
#
#
#mkdir($outputDir);