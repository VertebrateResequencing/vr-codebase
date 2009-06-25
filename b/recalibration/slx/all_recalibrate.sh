MAQ=$G1K/bin/maq
MALE_REF=/lustre/sf4/1kgenomes/ref/human_b36_male.bfa
FEMALE_REF=/lustre/sf4/1kgenomes/ref/human_b36_female.bfa
DB_SNP=/lustre/sf4/1kgenomes/misc/dbSNP/snp126.snp
OUT_FILE=recalibrate.o
ERROR_FILE=recalibrate.e
MIN_MAPPING_QUALITY=0

if [ $# -ne 4 ]
then
    echo "Syntax: $0 root_dir projects_list_file sample_size LSF_queue"
    exit
fi

ROOT=$1
PROJECTS_LIST_FILE=$2
SAMPLE_SIZE=$3
LSF_QUEUE=$4

if [ ! -f $PROJECTS_LIST_FILE ]
then
	echo "Cant find projects list file"
	exit
fi 

if [ $LSF_QUEUE != "normal" ] && [ $LSF_QUEUE != "long" ] && [ $LSF_QUEUE != "basement" ] && [ $LSF_QUEUE != "yesterday" ]
then
	echo "LSF queue must be: long, basement, yesterday or normal"
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
        for individual in *
        do
		if [ ! -d $individual ]
		then
			continue
		fi
                echo "In $individual"
                cd $individual
		
		if [ ! -d SLX ]
		then
			cd ..
			continue
		fi
		
		#only for solexa reads
		cd SLX
		
		#libraries
                for library in *
                do
			if [ ! -d $library ]
			then
				continue
			fi
                        echo "In $library"
                        cd $library
			
			libID=$RANDOM
			
			#lanes
			read1Name=""
			read2Name=""
			for lane in *
			do
				if [ ! -d $lane ]
				then
					continue
				fi
				cd $lane
				echo "In $lane"
				
				if [ ! -d sample ]
				then
					mkdir sample
				fi
				
				laneID=$RANDOM
				
				#check if the sampling and mapping already done
				if [ -f sample/sample.map ]
				then
					size=`ls -la sample/sample.map | awk '{print $5}'`
					
					if [ $size -gt 1000 ]
					then
						echo "Already finished lane"
						cd ..
						continue
					fi
				fi
				
				#delete any old files
				for m in sample/*
				do
					if [ -f $m ]
					then
						rm $m
					fi
				done
				
				lane_read0="na"
				lane_read1="na"
				lane_read2="na"
				expectedInsert="na"
				if [ -f meta.info ] || [ -h meta.info ]
				then
					lane_read0=`grep read0 meta.info | head -1 | awk -F: '{print $2}'`
					lane_read1=`grep read1 meta.info | head -1 | awk -F: '{print $2}'`
					lane_read2=`grep read2 meta.info | head -1 | awk -F: '{print $2}'`
				fi
				
				gender=`perl -w -e "use G1KUtilities;G1KUtilities::path2Gender('/lustre/sf4/1kgenomes/ref/genders.txt');"`
				
				pwd
				
				#sample 30million bases from the fastq files (paired OR unpaired)
				paired=0
				numFiles=`ls *.fastq.gz | grep -v "recal" | wc -l`
				echo "Found $numFiles uncalibrated fastq files"
				
				if [ $lane_read0 ] && [ ! $lane_read1 ] && [ ! $lane_read2 ]
				then
					bsub -J sample.$libID.$laneID -q normal -o sample/$OUT_FILE -e sample/$ERROR_FILE perl -w -e "use Recalibration;Recalibration::sampleSubsetBasesLaneSolexa( '.', '$SAMPLE_SIZE', '$lane_read0', 'sample/sample_0.fastq' );"
						
						cd sample
						if [ $gender == "male" ]
						then
							bsub -J map.$libID.$laneID -w "done(sample.$libID.$laneID)" -q normal -o $OUT_FILE -e $ERROR_FILE "$MAQ fastq2bfq sample_0.fastq sample_0.bfq; $MAQ map sample.map $MALE_REF sample_0.bfq"
						elif [ $gender == "female" ]
						then
							bsub -J map.$libID.$laneID -w "done(sample.$libID.$laneID)" -q normal -o $OUT_FILE -e $ERROR_FILE "$MAQ fastq2bfq sample_0.fastq sample_0.bfq; $MAQ map sample.map $FEMALE_REF sample_0.bfq"
						fi
				elif [ $lane_read1 ] && [ $lane_read2 ]
				then
					bsub -J sample.$libID.$laneID -q normal -o sample/$OUT_FILE -e sample/$ERROR_FILE perl -w -e "use Recalibration;Recalibration::sampleSubsetBasesLaneSolexa( '.', '$SAMPLE_SIZE', '$lane_read1', 'sample/sample_1.fastq', '$lane_read2', 'sample/sample_2.fastq' );"
							
					cd sample
					if [ $gender == "male" ]
					then
						bsub -J map.$libID.$laneID -w "done(sample.$libID.$laneID)" -q normal -o $OUT_FILE -e $ERROR_FILE "$MAQ fastq2bfq sample_1.fastq sample_1.bfq; $MAQ fastq2bfq sample_2.fastq sample_2.bfq; $MAQ map -a 1000 sample.map $MALE_REF sample_1.bfq sample_2.bfq"
					elif [ $gender == "female" ]
					then
						bsub -J map.$libID.$laneID -w "done(sample.$libID.$laneID)" -q normal -o $OUT_FILE -e $ERROR_FILE "$MAQ fastq2bfq sample_1.fastq sample_1.bfq; $MAQ fastq2bfq sample_2.fastq sample_2.bfq; $MAQ map -a 1000 sample.map $FEMALE_REF sample_1.bfq sample_2.bfq"
					fi
					
					paired=1
				fi
				
				#see if at least 100,000 reads were mapped - create a mapview file
				echo -e 'numMapped=`'$MAQ' mapstat sample.map | grep "Total number of reads" | awk \047{print $6}\047`' > commands.sh
				echo 'if [ $numMapped -gt 1000 ]
				then' >> commands.sh
				
				if [ $MIN_MAPPING_QUALITY -gt 0 ]
				then
					echo "mv sample.map sample.all.map" >> commands.sh
					echo "$MAQ submap -q $MIN_MAPPING_QUALITY sample.map sample.all.map" >> commands.sh
				fi
				
				if [ $gender == "male" ]
				then
					#run the map2qmap command with the dbSNP pos masked out
					if [ $paired -eq 1 ]
					then
						echo "$G1K/bin/map2qMapFile.pl -p $DB_SNP sample.map $MALE_REF > qualmapLH.txt" >> commands.sh
					else
						echo "$G1K/bin/map2qMapFile.pl -s -p $DB_SNP sample.map $MALE_REF > qualmapLH.txt" >> commands.sh
					fi
					
					echo "$G1K/bin/makeQualitiesBayesian.pl qualmapLH.txt > qualmapBayesian.txt" >> commands.sh
					echo "fi" >> commands.sh
					bsub -J recal.$libID.$laneID -w "done(map.$libID.$laneID)" -q normal -o $OUT_FILE -e $ERROR_FILE sh commands.sh
					
				elif [ $gender == "female" ]
				then
					#run the map2qmap command with the dbSNP pos masked out
					if [ $paired -eq 1 ]
					then
						echo "$G1K/bin/map2qMapFile.pl -p $DB_SNP sample.map $FEMALE_REF > qualmapLH.txt" >> commands.sh
					else
						echo "$G1K/bin/map2qMapFile.pl -s -p $DB_SNP sample.map $FEMALE_REF > qualmapLH.txt" >> commands.sh
					fi
					
					echo "$G1K/bin/makeQualitiesBayesian.pl qualmapLH.txt > qualmapBayesian.txt" >> commands.sh
					echo "fi" >> commands.sh
					bsub -J recal.$libID.$laneID -w "done(map.$libID.$laneID)" -q normal -o $OUT_FILE -e $ERROR_FILE sh commands.sh
				fi
				
				#create the recalibrated files
				cd ..
				dir=`pwd`
				if [ $paired -eq 1 ]
				then
					if [ $lane_read0 ] && [ -f $lane_read0 ]
					then
						recalread0=`echo $lane_read0 | sed 's/\(.*\)_1.fastq.gz/\1_1.recal.fastq.gz/'`
						bsub -J apply.$libID.$laneID -w "done(recal.$libID.$laneID)" -q normal -o sample/$OUT_FILE -e sample/$ERROR_FILE perl -w -e "use Recalibration;Recalibration::applyPosQualMap( 'sample/qualmapBayesian.txt', '$dir', '$lane_read0', '$recalread0' );"
					fi
					recalread1=`echo $lane_read1 | sed 's/\(.*\)_1.fastq.gz/\1_1.recal.fastq.gz/'`
					recalread2=`echo $lane_read2 | sed 's/\(.*\)_2.fastq.gz/\1_2.recal.fastq.gz/'`
					bsub -J apply.$libID.$laneID -w "done(recal.$libID.$laneID)" -q normal -o sample/$OUT_FILE -e sample/$ERROR_FILE perl -w -e "use Recalibration;Recalibration::applyPosQualMap( 'sample/qualmapBayesian.txt', '$dir', '$lane_read1', '$recalread1', '$lane_read2', '$recalread2' );"
				else
					recalread0=`echo $lane_read0 | sed 's/\(.*\).fastq.gz/\1.recal.fastq.gz/'`
					bsub -J apply.$libID.$laneID -w "done(recal.$libID.$laneID)" -q normal -o sample/$OUT_FILE -e sample/$ERROR_FILE perl -w -e "use Recalibration;Recalibration::applyPosQualMap( 'sample/qualmapBayesian.txt', '$dir', '$lane_read0', '$recalread0' );"
				fi				
				cd ..
                	done
			cd ..
		done
                cd ..
		cd .. #go up technology
        done
        cd ..
done < $PROJECTS_LIST_FILE
