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
    
    $config->setCommand("bwa", "bwa");
    $config->setFile("reference", "hg19bwaidx", "/share/massstorage/seq-tools/bwa/reference/");
    
    $config->setPath("gtf",        "/share/massstorage/seq-genomes/");
    $config->setFile("gtfNCBI",    "Homo_sapiens.NCBI.build37.2.converted.gtf",    $config->getPath("gtf"));
    $config->setFile("gtfUCSC",    "Homo_sapiens.UCSC.hg19.gtf",                   $config->getPath("gtf"));
    $config->setFile("gtfEnsembl", "Homo_sapiens.Ensembl.GRCh37.65.converted.gtf", $config->getPath("gtf"));
    
    $config->setPath("cufflinks", "/share/massstorage/seq-tools/cufflinks/current/");
    $config->setCommand("cufflinks", $config->getPath("cufflinks") . "cufflinks");
    $config->setCommand("cuffdiff",  $config->getPath("cufflinks") . "cuffdiff");
    
    $config->setPath("cufflinksOutput", "/share/massstorage/seq-tools/cufflinks/current/");
    $config->setFile("gtfUCSC",         "Homo_sapiens.UCSC.hg19.gtf");
    $config->setFile("gtfEnsembl",      "Homo_sapiens.Ensembl.GRCh37.65.converted.gtf");
    $config->setMaxNumberCores(42);
    
    
    return($config);
}


1;


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