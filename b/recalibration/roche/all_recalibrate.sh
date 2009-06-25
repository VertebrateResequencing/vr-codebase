DIR=`pwd`
SSAHA=$G1K/bin/ssaha2
MALE_REF_PREFIX=$G1K/ref/human_b36_male
FEMALE_REF_PREFIX=$G1K/ref/human_b36_female
MALE_REF_FASTA=$G1K/ref/human_b36_male.fa
FEMALE_REF_FASTA=$G1K/ref/human_b36_female.fa
DB_SNP=$G1K/misc/dbSNP/snp126.snp

if [ $# -ne 3 ]
then
    echo "Syntax: $0 root_dir projects_list_file sample_size"
    exit
fi

ROOT=$1
PROJECTS_LIST_FILE=$2
SAMPLE_SIZE=$3

if [ ! -f $PROJECTS_LIST_FILE ]
then
	echo "Cant find projects list file"
	exit
fi

cd $ROOT
DIR=`pwd`
echo "Starting in $DIR"

#populations
while read i
do
	if [ ! -d $i ]
	then
		continue
	fi
	
        echo "In $i"
        cd $i

	#individuals
        for j in *
        do
                echo "In $j"
                cd $j
		
		if [ ! -d 454 ]
		then
			cd ..
			continue
		fi
		
		#only for 454 reads
		cd 454
		
		#libraries
                for k in *
                do
                        echo "In $k"
                        cd $k
			
			#lanes
			for l in *
			do
				cd $l
				echo "In $l"
				
				if [ ! -d sample ]
				then
					mkdir sample
				fi
				
				#check if the sampling and mapping already done
				if [ -f sample/sample.table ] && [ `ls -l sample/sample.table |awk '{print $5}'` -gt 10 ]
				then
					echo "Already finished lane"
					cd ..
					continue
				fi
				
				#delete any old files
				for m in sample/*
				do
					if [ -f $m ]
					then
						rm $m
					fi
				done
				
				gender=`perl -w -e "use G1KUtilities;G1KUtilities::path2Gender('/lustre/sf4/1kgenomes/ref/genders.txt');"`
				
				pwd
				
				#cat the fastq files together and sample from it
				cd sample
				laneID=$RANDOM
				bsub -J cat.$laneID -q small -o recalibrate.o -e recalibrate.e "zcat ../*.fastq.gz > combined.fastq"
				
				bsub -J sample.$laneID -w "done(cat.$laneID)" -q small -o recalibrate.o -e recalibrate.e perl -w -e "use Recalibration;Recalibration::sampleSubsetBasesFastq( 'combined.fastq', '.', 'sample.fastq', '$SAMPLE_SIZE' );"
				
				#align the reads using ssaha2
				if [ $gender == "male" ] || [ $gender == "unknown" ]
				then
					bsub -J ssaha.$laneID -w "done(sample.$laneID)" -q normal -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000000 -o recalibrate.o -e recalibrate.e "rm combined.fastq;$SSAHA -454 -output cigar -save $MALE_REF_PREFIX sample.fastq | gzip -c > sample.cigar.gz"
				elif [ $gender == "female" ]
				then
					bsub -J ssaha.$laneID -w "done(sample.$laneID)" -q normal -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000000 -o recalibrate.o -e recalibrate.e "rm combined.fastq;$SSAHA -454 -output cigar -save $FEMALE_REF_PREFIX sample.fastq | gzip -c > sample.cigar.gz"
				fi
				
				bsub -J unique.$laneID -w "done(ssaha.$laneID)" -q small -o recalibrate.o -e recalibrate.e perl -w -e "use AssemblyTools;AssemblyTools::printSingleHitReadsCigar( 'sample.cigar.gz', 'sample.unique.cigar' );"
				
				if [ $gender == "male" ] || [ $gender == "unknown" ]
				then
					bsub -J table.$laneID -w "done(unique.$laneID)" -q normal -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000000 -o recalibrate.o -e recalibrate.e perl -w -e "use Recalibration;Recalibration::cigar2MatchMismatchTable( 'sample.fastq', '$MALE_REF_FASTA', 'sample.unique.cigar', '$DB_SNP', 'sample.table' );"
				elif [ $gender == "female" ]
				then
					bsub -J table.$laneID -w "done(unique.$laneID)" -q normal -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000000 -o recalibrate.o -e recalibrate.e perl -w -e "use Recalibration;Recalibration::cigar2MatchMismatchTable( 'sample.fastq', '$FEMALE_REF_FASTA', 'sample.unique.cigar', '$DB_SNP', 'sample.table' );"
				fi
				
				cd ..
				
				#apply the recalibration table
				count=1
				for original in `ls *.fastq.gz | grep -v recal`
				do
					newName=`echo $original | awk -F. '{print $1".recal.fastq.gz"}'`
					bsub -J apply.$laneID.$count -w "done(table.$laneID)" -q small -o sample/recalibrate.o -e sample/recalibrate.e perl -w -e "use Recalibration;Recalibration::applyMatchMismatchTable( '$original', 'sample/sample.table', '$newName', 'sample/sample.recal' );"
					((count++))
				done
							
				cd ..
                	done
			cd ..
		done
                cd ..
		cd .. #go up technology
        done
        cd ..
done < $PROJECTS_LIST_FILE
