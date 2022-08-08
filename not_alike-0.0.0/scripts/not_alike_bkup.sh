#!/usr/bin/env bash

#
#	not_alike pipeline.
#	Finds dissimilar sequences of a query genome
#	Comparing it to a database of determined size.
#	This database could have hugh size.
#

#
#	Author: Javier Montalvo-Arredondo.
#	Versin: 3.0.0.
#	Contact: Calzada Antonio Narro 1923
#	Buenavista Saltillo Coahuila México.
#	Basic Sciences division.
#


while getopts "i:s:S:d:t:h" opt;
do

	case $opt in

		h) echo "Usage: ...";;
		i) sequences=$OPTARG;;
		s) size=$OPTARG;;
		S) step=$OPTARG;;
		d) dbase=$OPTARG;;
		t) task=$OPTARG;;
		*) echo "Usage ...";;

	esac

done
shift $((OPTIND-1))

PID=$RANDOM$RANDOM

#	NOT_ALIKE BEGINS.

echo -e "\n[NOT_ALIKE][BEGINS]\n"

echo "Sequences file name $sequences"
echo "Size value = $size"
echo "Step size value = $step"
split_seqs_file=$(basename ${sequences%.*}_split.$size.$step.fasta)
echo "Output file name = $split_seqs_file"
echo "DataBase name is = $dbase"
echo "Task = $task"
echo "Process PID = $PID"
dbase_dirname=$(dirname $dbase)

#	MANAGER :: SPLIT_GENOME BEGINS.

echo -e "\n[MANAGER::SPLIT_GENOME][BEGINS]\n"

if [ ! -d split_out ];
then
	mkdir split_out
else
	echo "split_out directory exists!"
fi

if [ ! -f split_out/${split_seqs_file}  ];
then
	manager split_genome $sequences $size $step split_out/${split_seqs_file}
else
	echo "split_out/${split_seqs_file} exists!"
fi



echo -e "\n[MANAGER::SPLIT_GENOME][ENDS]\n"


if [ ! -d blast_out ];
then
	mkdir blast_out
else
	echo "blast_out directory exists!"
fi

#	MEGABLAST BEGINS.

echo "Copy split_out/${split_seqs_file} to split_out/input_seqs.${PID}.txt"
cp split_out/${split_seqs_file} split_out/input_seqs.${PID}.txt

echo -e "\n[BLAST-SEARCHING][BEGINS]\n"

cat $dbase | \
	while read -r line;
	do
	echo "BLAST split_out/input_seqs.txt against $dbase_dirname/${line%.*}.db"
	blastn -query split_out/input_seqs.${PID}.txt \
		-db $dbase_dirname/${line%.*}.db \
		-task $task \
		-out blast_out/$(basename ${sequences%.*}).blast \
		-evalue 1 \
		-outfmt "6 qseqid" \
		-max_target_seqs 1

##		FOR EACH BLAST HITS FILE.
##		MANAGER :: HIDE_MATCHED_SEQS BEGINS.
	
	NUM_OF_HITS=$(wc -l blast_out/$(basename ${sequences%.*}).blast | cut -d" " -f1)
	echo -e "Number of hits = $NUM_OF_HITS\n"
		
	if [ $NUM_OF_HITS -gt 2 ];
	then
		cat blast_out/$(basename ${sequences%.*}).blast | \
			while read line; 
			do 
				echo ">${line}" >> blast_out/tmp
			done

		uniq blast_out/tmp > blast_out/tmp2
		rm blast_out/tmp

	#### I introduced the line tagged below because the *.blast file is corrupted in somehow.
	#### I do not know if it came from cat, from while loop or from uniq linux command.
	#### But such corruption also corrupts the memory allocation in fopen function in
	#### write_seqs function in libio.c.
	#### Adding a last return character to the *.blast file reverts the corruption.
	#### Do not ask me why, please.
	
		echo "" >> blast_out/tmp2 ### This is the line I introduced to the code.
	
		mv blast_out/tmp2 blast_out/$(basename ${sequences%.*}).blast
		tail blast_out/$(basename ${sequences%.*}).blast
	
		echo -e "\n[MANAGER::HIDE_MATCHED_SEQS][BEGINS]\n"
		manager hide_matched_seqs split_out/input_seqs.${PID}.txt blast_out/$(basename ${sequences%.*}).blast DNA split_out/tmp.${PID}.txt
		if [ ! -s split_out/tmp.${PID}.txt ];
		then
			echo -e "ERROR!\n"
			echo -e "split_out/tmp.${PID}.txt does not exist!\n"
			break
		fi
		echo -e "\n[MANAGER::HIDE_MATCHED_SEQS][ENDS]\n"
		echo -e "Updating split_out/input_seqs.${PID}.txt\n"
		mv split_out/tmp.${PID}.txt split_out/input_seqs.${PID}.txt
	else
		echo "\nNo BLAST hits, skipping updating split_out/input_seqs.txt\n"
	fi

	done

echo -e "\n[BLAST-SEARCHING][ENDS]\n"

echo -e "\n[HISAT2-BUILD][BEGINS]\n"

if [ -f ht2-idx/$(basename ${sequences%.*}).1.ht2 ];
then

	echo -e "I found $(basename ${sequences%.*}).1.ht2 in ht2-idx folder\n"
	echo -e "Skipping ht2-idx building\n"

else
	mkdir -p ht2-idx
	echo -e "building ht2-idx of $sequences\n"
	hisat2-build $sequences ht2-idx/$(basename ${sequences%.*})
fi

echo -e "\n[HISAT2-BUILD][ENDS]\n"

echo -e "\n[HISAT2-ALIGN][BEGINS]\n"


if [ -d sam_out ];
then
	echo -e "sam_out directory exists!\n"
else
	mkdir sam_out
fi

hisat2 -x ht2-idx/$(basename ${sequences%.*}) -f -U split_out/input_seqs.${PID}.txt -S sam_out/nal_frags.${PID}.sam


echo -e "\n[HISAT2-ALIGN][ENDS]\n"

echo -e "\n[SAMTOOLS][BEGINS]\n"

samtools view -b -h -f 0 -o sam_out/nal_frags.${PID}.bam sam_out/nal_frags.${PID}.sam
samtools sort -o sam_out/nal_frags.${PID}.sort.bam -O BAM --reference $sequences sam_out/nal_frags.${PID}.bam

echo -e "\n[SAMTOOLS][ENDs]\n"

echo -e "\n[STRINGTIE][BEGINS]\n"

if [ -d gtfs ];
then
	echo -e "gtfs folder exists!"
else
	mkdir -p gtfs
fi

stringtie sam_out/nal_frags.${PID}.sort.bam -L -s 1.5 -g 50 -o gtfs/nal_frags.${PID}.gtf 

echo -e "\n[STRINGTIE][ENDS]\n"

echo -e "\n[NOT_ALIKE][ENDS]\n"

echo -e "My best wishes!\n"
echo -e "Sincerely Javier Montalvo."

