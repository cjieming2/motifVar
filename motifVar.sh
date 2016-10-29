#!/bin/bash

#########################################################
## motifVar.sh shell script
## A semi-automated pipeline
## Please refer to README for more detailed description
#########################################################
 
if [[ "$#" -ne 9 && "$1" -ne -1 && ( "$1" -ne "protPos2gPos" && "$#" -ne 5 ) ]] ; then
	echo "==============================="
	echo "== USAGE ======================"
	echo "==============================="
	echo "motifvar.sh <flag> <smart domain name> <Ensembl version> <Ensembl fasta file> <ensembl2coding file> <SNV BED file> <SNV catalog name> <VEP path> <anc allele BED file>" 
	echo "e.g. motifvar.sh 12 TPR /path/HUMAN_TPR.fa 73 /path/to/ensembl2coding.file ExAC.pass.bed ExAC.r0.3 /path/to/variant_effect_predictor.pl /path/to/anc.allele.bed > motifvar-160424.log"
  echo "--note for paths do not end with '/' pls"
  echo "this script needs 9 arguments"
  echo "--if a file is module-specific, you can use '-' to replace"
  echo "arg1:flag (refer to flag code) "
  echo "arg2:domain name from SMART database, e.g. TPR"
  echo "arg3:module 1; full path of fasta file from SMART database, with header >smart|TPR-ensembl|ENSP00000245105|ENSP00000245105/560-593 no description [Homo sapiens]"
  echo "arg4:Ensembl version e.g. 73"
  echo "arg5:module 4; full path of ensembl2coding file (adapted from Ensembl BioMart); required by module 'protPos2gPos' as well"
  echo "      --file header: EnsemblGeneID   EnsemblTranscriptID"
	echo "        EnsemblProteinID    EnsemblExonID  chr      strand"  echo "genomicCodingStart (1-based)      genomicCodingEnd"
	echo "arg6:modules 5, 6 and 7; SNV BED filename (should be in folder 4)"
	echo "arg7:modules 5 and 6; SNV catalog name, e.g. ExAC.r0.3"
	echo "arg8:module 6; path of VEP 73 tool"
	echo "arg9:module 9; path of ancestral allele BED file with only 4 cols and the 4th col is the AA (not derived allele)"
  echo -e "\n\n"
  echo "==============================="
  echo "== FLAG CODE =================="
  echo "==============================="
  echo "\"protPos2gPos\": requires 5 args; run only the module 'protPos2gPos' (SmartAApos2genomePos), which converts protein positions in SMART domains to genomic positions using Ensembl IDs of proteins with SMART domains"
  echo ""
  echo "--- ALL FLAG CODES below require 9 args ---"
  echo "--- optional args can be replaced by '-' ---"
  echo "1:fasta2prot; requires args 1, 2, 3, 4; convert fasta files from SMART to tab-delimited files, using the fasta headers"
  echo "2:domain2info; requires args 1, 2, 4; requires \"protPos2gPos\", obtain domain info from protPos2gPos info master file"
  echo "3:info2seq; requires args 1, 2, 4; requires modules 1 and 2, currently grabs fasta sequence remotely from Ensembl 73; requires installation of WebLogo and module Bio::EnsEMBL::Registry in Perl installation, with .bash_profile modified to run WebLogo (download from GitHub), grabs protein sequences and make a sequence logo using WebLogo"
  echo "4:domain2codon; requires args 1, 2, 4, 5; requires modules 1 and 3, outputs a BED file; splice up the DNA sequences into codons using ensembl2coding file"
  echo "5:codonIntersectSnv; requires args 1, 2, 4, 6, 7; requires modules 1 and 4, output a BED file; intersects codon BED file with SNV BED file, DOES NOT remove chrX, chrY, chrM or singletons (allele count 1)"
  echo "6:VEP and retains only coding, canonical transcripts, output a BED-like TXT file (sed 1d to convert to BED); requires args 1, 2, 4, 6, 7, 8; requires module 5 and VEP 73, Ensembl v73 cache and API to be installed. Do NOT do this on the head node. It assumes also that the SNV BED file has the following columns: \n col1:chr, col3:end position 1 based, col5:ref allele, col6:alt allele \n NOTE: if SNV file is not in the format above, and not VEP 73, this module needs to be manually edited in this script. "
  echo "7:SNVmerge, merges the SNV BED+codon file using the VEP output in module 6; requires args 1, 2, 4, 7; requires modules 6; produces 2 output files, one with and without sex chromosomes"
  echo "8:SNVprofiles, provides files for SIFT, R/C, NS/S calculations; requires args 1, 2, 4, 7; requires modules 7; uses autosomal file from module 7 to produce SNV profiles"
  echo "9:deltaDAF, provides results for deltaDAF; requires args 1, 2, 4, 7, 9; requires modules 7; uses autosomal file from module 7 to produce SNV profiles"
  echo "a combination of modules:e.g. 12, will run modules 1 and 2; 123, runs modules 1,2 and 3"
  echo "-1: if $1 -eq -1, clean up everything; creates a folder 'trash'"
  echo "e.g. motifvar.sh -1"
  echo -e "\n\n"
	echo "==============================="
	echo "== 'protPos2gPos' module ======"
	echo "==============================="
	echo "arg1: protPos2gPos"
	echo "arg2: ensembl2coding file path (adapted from Ensembl BioMart)"
	echo "      --file header: EnsemblGeneID   EnsemblTranscriptID"
	echo "        EnsemblProteinID    EnsemblExonID  chr      strand"  echo "genomicCodingStart (1-based)      genomicCodingEnd"
	echo "arg3: SMART domain file path with EnsemblProtID (converted from"
	echo "      motifVar_smartAApos2tsv after obtaining SMART domain info from" 
	echo "Ensembl using perl_api_smart_domains_suganthi.pl)"
	echo "      --file header: chr     smart   protaastart     protaaend"
	echo "        EnsemblProtID"
	echo "arg4: SMART domain file path directly from SMART database"
	echo "      --file header: DOMAIN  ACC     DEFINITION      DESCRIPTION"
	echo "arg5: ensembl version"
	echo "e.g. motifvar.sh protPos2gPos /ens/path/ensembl2coding_ens.noErr.txt /ens/path/allchr.ens73.noErr.tsv /smart/path/smart_domains_1158_all_131025.txt 73"

	exit 1
fi


#######################################################
## this causes the script to exit if there's any error in any part of the pipeline 
## The shell does NOT exit if the command that fails is part of the command list 
## immediately following a while or until keyword, part of the test in an if statement, 
## part of a && or || list, or if the command's return value is being inverted via !
#set -e

#############################################################
## -1 clean up everything
#############################################################
if [[ $1 -eq -1 ]] ; then

	mkdir trash
	mv motifVar_protPos2gPos [12]-* trash
	mv *.err *.log trash
	
	exit 0

fi

#############################################################
## protPos2gPos - converts protein positions in SMART domains to genomic positions using EnsemblProtIDs and SMART domains
#############################################################
if [[ $1 == "protPos2gPos" ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 0-motifVar_protPos2gPos-m
	cd 0-motifVar_protPos2gPos-m
	
	ENSEMBLFILE=$2
	DOMAINENSFILE=$3
	DOMAINFILE=$4
	ENS_VER=$5
	
	LFILE="0-motifVar_protPos2gPos.log"
	
	ln -s ${ENSEMBLFILE}
	ln -s ${DOMAINENSFILE}
	ln -s ${DOMAINFILE}

	## start log print
	echo "#######################" > ${LFILE}
	echo "## protPos2gPos #######" >> ${LFILE}
	echo "#######################" >> ${LFILE}
	 
	echo "Parameters:" >> ${LFILE}
	echo "ENSEMBL FILE         =${ENSEMBLFILE}" >> ${LFILE}
	echo "DOMAIN ENSEMBL FILE  =${DOMAINENSFILE}" >> ${LFILE}
	echo "DOMAIN FILE          =${DOMAINFILE}" >> ${LFILE}
	echo "ENSEMBL DB VERSION   =${ENS_VER}" >> ${LFILE}
	
	echo -e "\n\n" >> ${LFILE}
	
	date >> ${LFILE}
	
	## convert protein positions of SMART domains to genomic positions using Ensembl
	motifVar_smartAApos2genomePos -e ${ENSEMBLFILE} ${DOMAINENSFILE} | awk 'NR == 1; NR > 1 {if($2>$3){print $1"\t"$3"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}else{print $0}}' > temp.txt
	
	echo "Converting protein positions in SMART domains to genomic positions using EnsemblProtIDs and SMART domains... Done." >> ${LFILE}
	
	## integrate SMART domain information
	fselect a.+,b.DOMAIN,b.DEFINITION from temp.txt, $DOMAINFILE where a.smart=b.ACC | sed 's/ /_/g' | sed -e 's/\t\t/\tNA\t/g' -e 's/\t\n/\tNA\n/g' > motifVar_protPos2gPos.ens${ENS_VER}.smartdomains.txt
	
	## housekeeping	
	echo -e "motifVar_protPos2gPos.ens${ENS_VER}.smartdomains.txt created."
	echo ""
	echo ""
	echo "Program ran successfully" >> ${LFILE}
	echo "[MANUAL WORK] Pls note: post-processing on the smartDomain2gPos output file might be needed - e.g. ID errors, missing data, linesWprobs, CDS not fully annotated etc. This script doesn't check for these. Check Word documentation for details." >> ${LFILE}
	date >> ${LFILE}
	
	rm temp.txt
	
	cd ..
	
	exit 0
	
fi
 

#############################

FLAG=$1
DOMAIN=$2  
SMARTFASTA_PATH=$3
ENS_VER=$4
#ENSEMFASTA_PATH=$5
ENSEMBFILE_PATH=$5
SNV_BED=$6
SNV_NAME=$7
VEP_PATH=$8
AA_PATH=$9


#############################
## print parameters to log
#############################
echo "FLAG               =${FLAG}"
echo "DOMAIN             =${DOMAIN}"
echo "SMARTFASTA_PATH    =${SMARTFASTA_PATH}"  
echo "ENS_VER            =${ENS_VER}"
#echo "ENSEMFASTA_PATH    =${ENSEMFASTA_PATH}"
echo "ENSEMBFILE_PATH    =${ENSEMBFILE_PATH}"
echo "SNV_BED            =${SNV_BED}"
echo "SNV_CATALOG_NAME   =${SNV_NAME}"
echo "VEP_PATH           =${VEP_PATH}"
echo "AA_BED_PATH        =${AA_PATH}"
echo -e "\n"
#############################################################
## 1 fasta2prot
#############################################################
if [[ ${FLAG} =~ 1 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 1-fasta2prot-${DOMAIN}
	cd 1-fasta2prot-${DOMAIN}
	
	IFILE="${DOMAIN}.fasta"
	OFILE_prot="${DOMAIN}.prot"
	OFILE_indi="${DOMAIN}.prot.indiv"
	LFILE="1-fasta2prot-${DOMAIN}.log"
	SFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ${SMARTFASTA_PATH} ${IFILE}
	
	## start log print
	echo "######################" > ${LFILE}
	echo "## 1-fasta2prot ######" >> ${LFILE}
	echo "######################" >> ${LFILE}
	echo "INPUT FILE = ${IFILE}" >> ${LFILE}
	echo "OUTPUT FILE-prot = ${OFILE_prot}" >> ${LFILE}
	echo "OUTPUT FILE-indi = ${OFILE_indi}" >> ${LFILE}
	echo "motif length file= ${SFILE}" >> ${LFILE}
	
	date >> ${LFILE}
	
	## converts smart fasta to tab-delimited prot grouped by proteins
	motifVar_fastaSmart2tsv4lyn ${IFILE} -o ${OFILE_prot}

	## converts smart fasta to tab-delimited prot split by individual motif
	motifVar_fastaSmart2tsv4lyn ${IFILE} -o ${OFILE_indi} -i 1
	
	echo "Converted SMART fasta to prot and prot.indiv... Done." >> ${LFILE}
	
	## obtain most common motif length
	distinct -kc$(head -2 ${OFILE_indi} | ftranspose | grep -n '' | sed 's/:/\t/g' | grep length | cut -f1) ${OFILE_indi} | sort -nrk2 | head -1 | cut -f1 > ${SFILE}
	
	echo -e "\n\n" >> ${LFILE}
	distinct -kc$(head -2 ${OFILE_indi} | ftranspose | grep -n '' | sed 's/:/\t/g' | grep length | cut -f1) ${OFILE_indi} | sort -nrk2 >> ${LFILE}
	
	echo "Obtained most common motif length... Done." >> ${LFILE}
	
	## housekeeping
	date >> ${LFILE}
	
	cd ..
fi

##############################################################
### 2 domain info (domain2info)
##############################################################
if [[ ${FLAG} =~ 2 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 2-domain2info-${DOMAIN}
	cd 2-domain2info-${DOMAIN}
	
	IFILE="motifVar_protPos2gPos.ens${ENS_VER}.smartdomains.txt"
	OFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.txt"
	LFILE="2-domain2info-${DOMAIN}.log"
	SFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.ensProteinIDList"
	ln -s ../0-motifVar_protPos2gPos-m/${IFILE}
	
	## start log print
	echo "#########################" > ${LFILE}
  echo "## 2-domain2info ########" >> ${LFILE}
  echo "#########################" >> ${LFILE}
  echo "INPUT FILE  = ${IFILE}" >> ${LFILE}
  echo "OUTPUT FILE = ${OFILE}" >> ${LFILE}
  echo "OUTPUT FILE2= ${SFILE}" >> ${LFILE}
  
  date >> ${LFILE}
	
	## get domain info file from master file created in protPos2gPos
	d=$(echo ${DOMAIN})
	awk '{OFS="\t"}{FS="\t"}{if($12 == "DOMAIN" || $12 == domain){print $0}}' domain=$d ${IFILE} > ${OFILE}
	
	echo "cut -f5 | uniq -c" >> ${LFILE}
	cut -f5 ${OFILE} | uniq -c >> ${LFILE}
	
	## produce unique Ensembl protein IDs with this domain
	cut -f8 ${OFILE} | sed 1d | sort | uniq > ${SFILE}
	
	## housekeeping
	date >> ${LFILE}

	cd ..
fi



##############################################################
### 3  add domain sequence (info2seq)
##############################################################
if [[ ${FLAG} =~ 3 ]] ; then
	
	## set up
	mkdir 3-domaininfo2seq-${DOMAIN}
	cd 3-domaininfo2seq-${DOMAIN}
	
	IFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.txt"
	IFILE2="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.ensProteinIDList"
	OFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.txt"
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../1-fasta2prot-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
  OFILE2="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.${CLEN}aa.txt"
  OFILE2_NAME="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.${CLEN}aa"
  OFILE3="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.ensProteinIDList.fasta"
	LFILE="3-domaininfo2seq-${DOMAIN}.log"
	SFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.${CLEN}aa.1${DOMAIN}"
	ln -s ../2-domain2info-${DOMAIN}/${IFILE}
	ln -s ../2-domain2info-${DOMAIN}/${IFILE2}
	
	## start log printing	
	echo "##########################" > ${LFILE}
  echo "## 3-info2seq ############" >> ${LFILE}
  echo "##########################" >> ${LFILE}
  echo "INPUT FILE  = ${IFILE}" >> ${LFILE}
  echo "INPUT FILE2 = ${IFILE2}" >> ${LFILE}
  echo "OUTPUT FILE = ${OFILE}" >> ${LFILE}
  echo "OUTPUT FILE2= ${OFILE2}" >> ${LFILE}
  echo "ENSEM FASTA = ${OFILE3}" >> ${LFILE}
  
  date >> ${LFILE}
  
  ## grab fasta file from Ensembl 73
  echo "Grab fasta file remotely from Ensembl 73..."
  ensembl_fetch_protein_fasta_sequences_motifvar ${IFILE2} ${OFILE3} >> ${LFILE}

	## integrate sequence from fasta file from BioMart
	motifVar_fastaGrab -f ${OFILE3} ${IFILE} > ${OFILE}
  
  ## grab info and sequences of the most common motif length
  awk '{OFS="\t"}{FS="\t"}{if($11 == commonlength || $11 == "motifSize"){print $0}}' commonlength=${CLEN} ${OFILE} > ${OFILE2}

	## create sequence logos of most common motif length using WebLogo
	cut -f12 ${OFILE2} | sed 1d | sort | uniq > ${SFILE}
	weblogo -f ${SFILE} -o ${SFILE}.pdf -F pdf -A protein -U bits --composition "{'L':9.975,'A':7.013,'S':8.326,'V':5.961,'G':6.577,'K':5.723,'T':5.346,'I':4.332','E':7.096,'P':6.316,'R':5.650,'D':4.728,'F':3.658,'Q':4.758,'N':3.586,'Y':2.653,'C':2.307,'H':2.639,'M':2.131,'W':1.216}"  -n ${CLEN} -c chemistry --stack-width 25 --errorbars no
	
	## get frequency tables
	msa2consensus4lyn -c 12 -i 1 -n ${OFILE2_NAME} ${OFILE2}

  ## housekeeping
	date >> ${LFILE}
	
	cd ..
fi


##############################################################
### 4 domain2codon
##############################################################
if [[ ${FLAG} =~ 4 ]] ; then
	## set up
	mkdir 4-domain2codon-${DOMAIN}-m
	cd 4-domain2codon-${DOMAIN}-m 
	
	LFILE="4-domain2codon-${DOMAIN}.log"
	EFILE=${ENSEMBFILE_PATH}
	ln -s ${EFILE}
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../1-fasta2prot-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
  IFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.${CLEN}aa.txt"
  ln -s ../3-domaininfo2seq-${DOMAIN}/${IFILE}
  OFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.${CLEN}aa.codon.bed"
	
	## start log printing
	echo "############################" > ${LFILE}
  echo "## 4-domain2codon ##########" >> ${LFILE}
  echo "############################" >> ${LFILE}
  echo "INPUT FILE   = ${IFILE}" >> ${LFILE}
  echo "ENSEM FILE   = ${EFILE}" >> ${LFILE}
  echo "OUTPUT FILE  = ${OFILE}" >> ${LFILE}
   
  date >> ${LFILE}
  
	## converting to codons
	motifVar_Domain2ResidueBed -e ${EFILE} ${IFILE} > ${OFILE}

	date >> ${LFILE}
	echo "[MANUAL WORK] This module prepares the domains, by converting them into codons. The next file needs to be prepared manually: the SNV BED file from a catalog, e.g. ExAC. Since the last info column of this BED file can vary depending on the catalog source and since there isn't header in a BED file, the SNV BED file needs to be prepared separately and manually. Put this file in this 4-~ folder. e.g. ExAC.r0.3.sites.vep.snps.pass.mod.bed - you do not need to remove singletons" >> ${LFILE}
	
	cd ..
fi


##############################################################
### 5 codonIntersectSnv
##############################################################
if [[ ${FLAG} =~ 5 ]] ; then

	## set up
	mkdir 5-codonIntersectSnv-${DOMAIN}
	cd 5-codonIntersectSnv-${DOMAIN} 
	
	LFILE="5-codonIntersectSnv-${DOMAIN}.log"
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../1-fasta2prot-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
  IFILE="motifVar_protPos2gPos.ens${ENS_VER}.${DOMAIN}.seq.${CLEN}aa.codon.bed"
  ln -s ../4-domain2codon-${DOMAIN}-m/${IFILE}
  ln -s ../4-domain2codon-${DOMAIN}-m/${SNV_BED}
  OFILE="motifVar_${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.bed"
	
	## start log printing
	echo "############################" > ${LFILE}
  echo "## 5-codonIntersectSnv ##########" >> ${LFILE}
  echo "############################" >> ${LFILE}
  echo "INPUT FILE     = ${IFILE}" >> ${LFILE}
  echo "SNV BED FILE   = ${SNV_BED}" >> ${LFILE}
  echo "OUTPUT FILE    = ${OFILE}" >> ${LFILE}
  
  date >> ${LFILE}
  
	## intersect codon with SNVs
	intersectBed -a ${SNV_BED} -b ${IFILE} -wa -wb > ${OFILE}
	
	
	#col=${AC_COL} ## ## remove chrX chrY chrM ; removes singletons
	#intersectBed -a ${SNV_BED} -b ${IFILE} -wa -wb | awk '{OFS="\t"}{FS="\t"}{if($1 != "chrX" && $1 != "chrY" && $1 != "chrM"){print $0}}' | awk '{OFS="\t"}{FS="\t"}{if($col > 1){print $0}}' col=$col > ${OFILE}
		
	## count how many SNVs
	wc -l ${OFILE} >> ${LFILE}
	
	## count how many unique SNVs
	echo "Number of unique SNVs = " $(cut -f1-3 ${OFILE} | sortByChr.sh - | uniq | wc -l) >> ${LFILE}

	echo "NOTE: ${OFILE} can be redundant because an SNV can affect codons from different proteins from the same gene." >> ${LFILE}
	date >> ${LFILE}
	
fi


##############################################################
### 6 vep
##############################################################
## note that this module uses the SNV file to run VEP
##############################################################

if [[ ${FLAG} =~ 6 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 6-vep-${DOMAIN}
	cd 6-vep-${DOMAIN}
	
	## set up
	LFILE="6-vep-${DOMAIN}.log"
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../5-codonIntersectSnv-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
	IFILE="motifVar_${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.bed"
	VIPUT="motifVar_${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepinput"
	ln -s ../5-codonIntersectSnv-${DOMAIN}/${IFILE}
	ln -s ../5-codonIntersectSnv-${DOMAIN}/${VIPUT}
	ln -s ../5-codonIntersectSnv-${DOMAIN}/${SNV_BED}
	VOPUT="motifVar_${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput"

	## start log printing
	echo "############################" > ${LFILE}
  echo "## 6-vep-${DOMAIN}" >> ${LFILE}
  echo "############################" >> ${LFILE}
  echo "INPUT FILE   = ${IFILE}" >> ${LFILE}
  echo "VEP INPUT    = ${VIPUT}" >> ${LFILE}
  echo "VEP OUTPUT   = ${VOPUT}" >> ${LFILE}
	
  date >> ${LFILE}
 	
 	echo "input file = ${IFILE}" >> ${LFILE}
 	
 	## make vepinput
 	awk '{OFS="\t"}{FS="\t"}{print $1,$3,$3,$5"/"$6,"+"}' ${IFILE} | sortByChr.sh - | sed 's/^chr//g' | uniq > ${VIPUT}
 	
 	## run vep by chromosome
 	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
 		do 
 			awk '{OFS="\t"}{FS="\t"}{if($1 == "'$i'"){print $0}}' ${VIPUT} > vepinput."$i".txt 
 			perl ${VEP_PATH} -i vepinput."$i".txt -o vepoutput."$i".txt --force -sift=b -polyphen=b --ccds --numbers --domains -regulatory --canonical --protein --symbol --cache >& vepinput."$i".log &
 		done
 	
 	wait
 	
 	## process and combine vepoutput
 	grep "#" vepoutput.1.txt > header.vepoutput
 	
 	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do grep -v "#" vepoutput."$i".txt ; done | cat header.vepoutput - > ${VOPUT}
 	
 	## produce unique SNV file
 	## grepping for coding regions has redundancy because of multiple transcripts
 	## grepping for canonical to get 'unique' canonical transcript
 	# not all unique though even with canonical, due to overlapping gene annotation
 	# need to check by comparing with vepinput
 	grep -v "##" ${VOPUT} | awk '{OFS="\t"}{FS="\t"}{if($6 == "Feature_type" || $6 == "Transcript"){print $0}}' | awk '{OFS="\t"}{FS="\t"}{if($7 == "Consequence" || $7 == "missense_variant" || $7 == "missense_variant,splice_region_variant" || $7 == "stop_gained" || $7 == "stop_lost" || $7 == "synonymous_variant" || $7 == "stop_gained,splice_region_variant" || $7 == "initiator_codon_variant" || $7 == "coding_sequence_variant" || $7 == "splice_acceptor_variant" || $7 == "splice_donor_variant" || $7 == "splice_region_variant,synonymous_variant"){print $0}}' | egrep "Feature|CANONICAL=YES" | egrep "SYMBOL|Feature" > ${VOPUT}.coding.canonical
 	
 	echo "distinct -kc1 ${VOPUT}.coding.canonical | sort -nrk2 | distinct -kc2 " >> ${LFILE}
 	distinct -kc1 ${VOPUT}.coding.canonical | sort -nrk2 | distinct -kc2 >> ${LFILE}
 	
 	echo "If there exists values of 2, it means that not all transcripts are unique, even with canonical." >> ${LFILE}
 	echo "Need to check by comparing with vepinput." >> ${LFILE}
 	echo "Note also that VEP gives an \"undefined value\" error when the vepinput file is empty." >> ${LFILE}
 	echo "The script will proceed with merging duplicate lines using mergeDuplicates.py and column 1, separating different elements in duplicated rows by ';'. Use coding.canonical file to reformat file if merging is undesirable." >> ${LFILE}
 	#echo "Please edit the file header.${SNV_BED} manually, and find the EnsemblTranscriptID for use in module 7." >> ${LFILE}
 	#echo "create the transID, by combining the header with the edited version of the SNV BED to show the corresponding fields with EnsemblTranscriptID. E.g. sed -e 's/|/\t/1' -e 's/|/\t/1' -e 's/|/\t/1' ExAC. r0.3.sites.vep.snps.pass.mod.bed | cat header.ExAC.r0.3.sites.vep.snps.pass.mod.bed - > ExAC.r0.3.sites.vep.snps.pass.mod.bed.transID"
 	
 	## header of SNV BED
 	head -1 ${SNV_BED} > header.${SNV_BED}
 	
	## convert to BED-like TXT file (with header)
	mergeDuplicates.py ${VOPUT}.coding.canonical 1 | sed -e 's/#Uploaded_variation/chr_position_allele/1' -e 's/_/\t/1' -e 's/_/\t/1' | awk '{OFS="\t"}{FS="\t"}{print "chr"$1, $2-1,$2,$0}' | cut -f1-3,6- > ${VOPUT}.coding.canonical.mergeDup.txt
	
	echo "awk '{OFS=\"\t\"}{FS=\"\t\"}{print \$1\"_\"\$2\"_\"\$3}' ${VOPUT}.coding.canonical.mergeDup.txt | distinct -kc1 | sort -nrk2 | distinct -kc2 " >> ${LFILE}
 	awk '{OFS="\t"}{FS="\t"}{print $1"_"$2"_"$3}' ${VOPUT}.coding.canonical.mergeDup.txt | distinct -kc1 | sort -nrk2 | distinct -kc2 >> ${LFILE}
 	
 	## clean up
 	mkdir trash
 	mv vepinput.*.txt vepinput.*.log vepoutput.*.txt vepoutput.*.html trash
	
 	
	date >> ${LFILE}
	
fi


##############################################################
### 7 SNVmerge
##############################################################
## want infor from (1) SNV BED with (2) domain codon info and (3) VEP output
## merge #1 input SNV+codon file with #2 VEP output
##############################################################

if [[ ${FLAG} =~ 7 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 7-SNVmerge-${DOMAIN}
	cd 7-SNVmerge-${DOMAIN}
	
	## set up
	LFILE="7-SNVmerge-${DOMAIN}.log"
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../6-vep-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
	IFILE="motifVar_${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.bed" #1 
	ln -s ../6-vep-${DOMAIN}/${IFILE}
	IFILE2="motifVar_${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.txt" #2
	ln -s ../6-vep-${DOMAIN}/${IFILE2}
	OFILE="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.bed"
	OFILE_AUTO="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.auto.bed"

	## start log printing
	echo "############################" > ${LFILE}
  echo "## 7-SNVmerge-${DOMAIN}" >> ${LFILE}
  echo "############################" >> ${LFILE}
	
  date >> ${LFILE}
 	
 	echo "SNV+CODON = ${IFILE}" >> ${LFILE}
 	echo "VEPOUTPUT = ${IFILE2}" >> ${LFILE}
 	echo "merged OUTPUT = ${OFILE}" >> ${LFILE}
 	echo "merged OUTPUT = ${OFILE_AUTO}" >> ${LFILE}
 	
 	## merge
	sed 1d ${IFILE2} | intersectBed -a ${IFILE} -b - -wa -wb | sed 's/#/\t/g' | awk '{OFS="\t"}{FS="\t"}{if($38==$54){print $0}}' | sortByChr.sh - | tee ${OFILE} | awk '{OFS="\t"}{FS="\t"}{if($1 != "chrX" && $1 != "chrY" && $1 != "chrM"){print $0}}' > ${OFILE_AUTO}
 	
 	echo "distinct -kc1 ${OFILE}" >> ${LFILE}
 	distinct -kc1 ${OFILE} >> ${LFILE}
 	echo "" >> ${LFILE}
 	echo "distinct -kc1 ${OFILE_AUTO}" >> ${LFILE}
 	distinct -kc1 ${OFILE_AUTO} >> ${LFILE}
	
 	
	date >> ${LFILE}
	
fi


##############################################################
### 8 SNVprofiles
##############################################################
## compute SIFT, R/C, NS/S
##############################################################

if [[ ${FLAG} =~ 8 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 8-SNVprofiles-${DOMAIN}
	cd 8-SNVprofiles-${DOMAIN}
	
	## set up
	LFILE="8-SNVprofiles-${DOMAIN}.log"
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../7-SNVmerge-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
	IFILE="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.auto.bed"
	ln -s ../7-SNVmerge-${DOMAIN}/${IFILE}
	OFILE="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.auto.sift.bed"
	OFILE_ENRICH="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.auto.sift.enrich"
	
	## start log printing
	echo "############################" > ${LFILE}
  echo "## 8-SNVprofiles-${DOMAIN}" >> ${LFILE}
  echo "############################" >> ${LFILE}
  echo "INPUT FILE    = ${IFILE}" >> ${LFILE}
  echo "WITH SIFT INFO= ${OFILE}" >> ${LFILE}
  echo "WITH R,C,NS,S = ${OFILE_ENRICH}" >> ${LFILE}
	
  date >> ${LFILE}
 	
 	## condel, sift and polyphen
 	motifVar_SiftPolyCondel -c 63 ${IFILE} > ${OFILE}
 	
 	## R/C and NS/S
 	motifVar_CalcRareNS -f 0.005 -c AF=20,NS=56,AC=8,ID=39,resNum=46 ${OFILE} > ${OFILE_ENRICH}
 	
	date >> ${LFILE}
	
fi

##############################################################
### 9 deltaDAF
##############################################################
## compute delta DAF between populations
##############################################################

if [[ ${FLAG} =~ 9 ]] ; then
	## setting up
	## make directories, go in directory, set up links
	mkdir 9-deltaDAF-${DOMAIN}
	cd 9-deltaDAF-${DOMAIN}
	
	## set up
	LFILE="9-deltaDAF-${DOMAIN}.log"
	NFILE="${DOMAIN}_mostCommonMotifLen.txt"
	ln -s ../7-SNVmerge-${DOMAIN}/${NFILE}
	CLEN=$(cat ${NFILE})
	IFILE="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.auto.bed"
	ln -s ../7-SNVmerge-${DOMAIN}/${IFILE}
	PREFIX="final.merged.${SNV_NAME}.ens${ENS_VER}.${DOMAIN}.${CLEN}aa.codon.vepoutput.coding.canonical.mergeDup.sorted.auto"
	
	## post-intersectBed
	OFILE_IB="${PREFIX}.intersected.bed"
	OFILE_DAF="${PREFIX}.daf.txt" # this will have a header
	OFILE_DDAF="${PREFIX}.ddaf.txt" # this will have a header
	
	## start log printing
	echo "############################" > ${LFILE}
  echo "## 9-deltaDAF-${DOMAIN}" >> ${LFILE}
  echo "############################" >> ${LFILE}
  echo "INPUT FILE         = ${IFILE}" >> ${LFILE}
  echo "PREFIX             = ${PREFIX}" >> ${LFILE}
  echo "FILE INTERSECTED   = ${OFILE_IB}" >> ${LFILE}
  echo "FILE WITH DAF      = ${OFILE_DAF}" >> ${LFILE}
	
  date >> ${LFILE}
 	
 	## delta DAF 
 	## for exac
 	## if $67 AA == $5 REF allele, then AC/AN since AC=alternate counts
	## order of pops: AFR AMR EAS FIN NFE OTH SAS
	## AFR=african-american/african, $9 and $22 (AC and AN)
	## AMR=american, $10 and $23
	## EAS=east asian, $12 and $25
	## FIN=finnish $13 and $26 + NFE=non-finnish european $17 and $27 = EUR=european
	## OTH=other $18 and $28
	## SAS=south asian $19 and $29
 
	## print intermediate files for possible debugging, tweaking of columns (since it's hardcoded here) or workarounds
	## there are some cases of 0 in the AN of any of the populations
	## this line produces 2 files: 1) BED file directly after intersection, to see which cols to take
	## 2) BED-like txt file with header for AC and AN counts only
	
	intersectBed -a ${IFILE} -b ${AA_PATH} -wa -wb > ${OFILE_IB} 
	
	## calculates DDAF
	PSEUDOCOUNT=1
	p=$(echo ${PSEUDOCOUNT})
	awk '{OFS="\t"}{FS="\t"}{if($67 == $5){print $1,$2,$3,$5,$6,$67,$46,$8,$56,$9,$22,$9/($22+p),$10,$23,$10/($23+p),$12,$25,$12/($25+p),$13+$17,$26+$27,($13+$17)/($26+$27+p),$19,$29,$19/($29+p)}else{print $1,$2,$3,$5,$6,$67,$46,$8,$56,$9,$22,1-($9/($22+p)),$10,$23,1-($10/($23+p)),$12,$25,1-($12/($25+p)),$13+$17,$26+$27,1-(($13+$17)/($26+$27+p)),$19,$29,1-($19/($29+p))}}' p=$p ${OFILE_IB} | awk '{OFS="\t"}{FS="\t"}function abs(x){return ((x < 0.0) ? -x : x)}{print $0,abs($12-$15),abs($12-$18),abs($12-$21),abs($12-$24),abs($15-$18),abs($15-$21),abs($15-$24),abs($18-$21),abs($18-$24),abs($21-$24)}' | cat <(echo -e "chr\tstart\tend\tref\talt\tAA\tresNum\tAC\tNS\tAFR_AC\tAFR_AN\tAFR_DAF\tAMR_AC\tAMR_AN\tAMR_DAF\tEAS_AC\tEAS_AN\tEAS_DAF\tEUR_AC\tEUR_AN\tEUR_DAF\tSAS_AC\tSAS_AN\tSAS_DAF\tdDAF-AFR_AMR\tdDAF-AFR_EAS\tdDAF-AFR_EUR\tdDAF-AFR_SAS\tdDAF-AMR_EAS\tdDAF-AMR_EUR\tdDAF-AMR_SAS\tdDAF-EAS_EUR\tdDAF-EAS_SAS\tdDAF-EUR_SAS") - > ${OFILE_DDAF}
	
	date >> ${LFILE}
	
fi