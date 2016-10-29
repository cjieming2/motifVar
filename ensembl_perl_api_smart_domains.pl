#!/usr/bin/perl

#####!/usr/local/ensembl/bin/perl


use strict;
use warnings;
use lib qw(ensembl/modules);

use Bio::EnsEMBL::DBSQL::DBAdaptor;;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
       -host =>'ensembldb.ensembl.org',
       -user => 'anonymous'
);

#Containers for the strings to be synthesized and then written to a file each. Each of the files will be a table.
my $strAllChromosomes="ChromosomeID\n";
my $strAllPrDomains="Pfam \t \t Start-End \t \t translationID\n";
###

my @chromosomes = ("Y");
my @db_adaptors = @ { $registry->get_all_DBAdaptors()};
my $slice_adaptor = $registry->get_adaptor('Human', 'Core', 'Slice');

#Get genes, transcripts and domains for each chromosome

foreach my $chr(@chromosomes){
       $strAllChromosomes .= "$chr\n";
       print "Chromosome: $chr...\n";

       #Define a slice containing a chromosome
       my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr);

       #Get all genes in chromosome
       my $AllGenes = $slice->get_all_Genes;

       foreach my $gene(@$AllGenes){

                       #Get transcripts of each gene...and then their translation and domain features..
                       foreach my $transcript(@{$gene->get_all_Transcripts()}){
                               my $translation = $transcript->translation();

                               if (not defined $translation ){
                               next ;}

                               my $domainFeatures = $translation->get_all_DomainFeatures();


                               foreach my $pf (@$domainFeatures) {
                               my $analysis = $pf->analysis();

                               if (defined($analysis) && lc($analysis->logic_name() ) eq "smart" ) {
                               $strAllPrDomains .=  $pf->hseqname() ."\t". $pf->start() .'-'. $pf->end() ." \t" . "(" .$translation->stable_id().")". "\n";
                               }
                               }

                       }
               }
       }

 print "$strAllPrDomains\n";
