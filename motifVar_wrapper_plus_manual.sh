#!/bin/bash

## This script tries to integrate the automated pipeline motifVar.sh and the manual portions of it
## the aim is to include all the manual commands that you have already figured out
## and make the run with only a single domain argument so that this can be looped
## Hence, many of the 'manual' parameters are hardcoded into this wrapper script.
## The logic and comments should be able to guide you through how to adapt this for your own use.
## you can loop it like this:
## for i in TPR 
##   do bsub-make-plus.sh motifvar-wrapper-"$i" "motifVar_wrapper_plus_manual.sh $i"
##   bsub -q gerstein < bsub-script-rdy-motifvar-wrapper-"$i".sh 
## done


if [[ "$#" -ne 1 && "$1" -ne -1 ]] ; then
	echo "==============================="
	echo "== USAGE ======================"
	echo "==============================="
	echo "motifVar_wrapper_plus_manual.sh <smart domain name>" 
	echo "e.g. motifVar_wrapper_plus_manual.sh TPR"

	exit 1
fi

## make directory of domain 
## enter domain
#mkdir $1
cd $1

#######################################################
####### protPos2gPos; 5 args + m
## e.g. motifvar.sh protPos2gPos /ens/path/ensembl2coding_ens.noErr.txt /ens/path/allchr.ens73.noErr.tsv /smart/path/smart_domains_1158_all_131025.txt 73
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module and modules 1 and 2
####### they use the same files, so it seems they can't be run concurrently
####### solution is to run 0,1,2 first in series on a separate script 
####### then run a second script to run the rest in parallel
 
#motifVar.sh protPos2gPos /gpfs/scratch/fas/gerstein/jc2296/ensembl/ensembl73/ensembl2coding_ensembl73.proteinIDs.genomicPos.chrs.strand.noErr.txt /gpfs/scratch/fas/gerstein/jc2296/ensembl/ensembl73/allchromosomes.ens73.alldomainfeatures.smart.mod.noErr.tsv /gpfs/scratch/fas/gerstein/jc2296/smart/131025/smart_domains_1158_all_131025.txt 73
#
### manual protPos2gPos
#cd 0-motifVar_protPos2gPos-m
#
#REMOVE_CDS_FILE="/gpfs/scratch/fas/gerstein/jc2296/gencode/gencode.17.cds_start_end_NF.ensembl2coding_ensembl73.proteinIDs.txt"
#PROTPOS2GPOS_OFILE="motifVar_protPos2gPos.ens73.smartdomains.txt"
#PROTPOS2GPOS_OFILE_OLD="motifVar_protPos2gPos.ens73.smartdomains.unprocessed.txt"
#
## save original file
#mv ${PROTPOS2GPOS_OFILE} ${PROTPOS2GPOS_OFILE_OLD}
#
## remove incomplete/truncated CDS sequences
#fsieve2 -s <(cut -f3 ${REMOVE_CDS_FILE}) -m <(awk '{OFS="\t"}{FS="\t"}{print $8,$0}' ${PROTPOS2GPOS_OFILE_OLD} | sed 's/EnsemblProtID/EnsemblProtID1/1') | cut -f2- > ${PROTPOS2GPOS_OFILE}
#
#
## get out of folder
#cd ..


#######################################################
####### 1 fasta2prot; 9 args
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module

#SMART_FASTA_PATH="/scratch/fas/gerstein/jc2296/smart/131025/fasta/HUMAN_$1.fa"
#motifVar.sh 1 $1 ${SMART_FASTA_PATH} 73 - - - - -
#
########################################################
######## 2 domain info (domain2info); 9 args
######## the columns to grab are hardcoded into the motifVar.sh 
######## for this module
#
#motifVar.sh 2 $1 - 73 - - - - -

#######################################################
####### 3  add domain sequence (info2seq); 9 args
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module

motifVar.sh 3 $1 - 73 - - - - -

#######################################################
####### 4 domain2codon; 9 args

ENSEMBL2CODING_PATH="/gpfs/scratch/fas/gerstein/jc2296/ensembl/ensembl73/ensembl2coding_ensembl73.proteinIDs.genomicPos.chrs.strand.noErr.txt"
motifVar.sh 4 $1 - 73 ${ENSEMBL2CODING_PATH} - - - -

## manual 4 domain2codon
SNV_BED_PATH="/gpfs/scratch/fas/gerstein/jc2296/exomes/exac/ExAC.r0.3.sites.vep.snps.pass.mod.bed"
cd 4-domain2codon-$1-m

ln -s ${SNV_BED_PATH}

cd ..


#######################################################
####### 5 codonIntersectSnv; 9 args

motifVar.sh 5 $1 - 73 - ExAC.r0.3.sites.vep.snps.pass.mod.bed ExAC.r0.3 - -


#######################################################
####### 6 vep; 9 args
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module

motifVar.sh 6 $1 - 73 - ExAC.r0.3.sites.vep.snps.pass.mod.bed ExAC.r0.3 /home/fas/gerstein/jc2296/software/vep_v73/ensembl-tools-release-73/scripts/variant_effect_predictor/variant_effect_predictor.pl -


#######################################################
####### 7 SNVmerge; 9 args
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module

motifVar.sh 7 $1 - 73 - - ExAC.r0.3 - -


#######################################################
####### 8 SNVprofiles; 9 args
####### by default, this module grabs module 7's auto file
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module
####### col63 for condel/sift/polyphen

motifVar.sh 8 $1 - 73 - - ExAC.r0.3 - -

#######################################################
####### 9 delta DAF; 9 args
####### by default, this module grabs module 7's auto file
####### the columns to grab are hardcoded into the motifVar.sh 
####### for this module
####### col63 for condel/sift/polyphen
AA_PATH="/scratch/fas/gerstein/jc2296/1KG/1KG_phase1_hg19/3-daf/all/nonmono.ALL.phase1_release_v3.20101123.allcov.35233368snps.txt.anc.allele.bed"
motifVar.sh 9 $1 - 73 - - ExAC.r0.3 - ${AA_PATH}

#######################################################
####### packaging for website

mkdir $1-motifvar
cd $1-motifvar

## aa freq tables and seq logo

for i in ../3-domaininfo2seq-$1/*.pdf ../3-domaininfo2seq-$1/*.aamat ../3-domaininfo2seq-$1/*.global ../3-domaininfo2seq-$1/*.re ../3-domaininfo2seq-$1/*.rfreq
	do
		ln -s "$i"
	done


## sift, r/c and ns/s
for i in ../8-SNVprofiles-$1/*.sift.bed ../8-SNVprofiles-$1/*.sift.enrich
	do
		ln -s "$i"
	done 
	
	
## delta daf
for i in ../9-deltaDAF-$1/*.ddaf.txt
	do
		ln -s "$i"
	done 


#######################################################
####### exit domain
cd ..