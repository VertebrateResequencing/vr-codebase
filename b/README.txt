                    DMC's 1000-genomes perl scripts
                    ===============================
                     David Carter, dmc@sanger.ac.uk
                              July 2008

If anything is missing, look in ~dmc/perl (or ~dmc/bin).

Many of these scripts expect the environment variable G1K to be set
appropriately. See ~dmc/.cshrc for how to do this.

THE MAIN PIPELINE SCRIPT

g1k.pl
  Run without arguments to remind yourself of the options.

COLLECTING THE DATA FROM NPG

updateDataDirectory.pl
  Query NPG to construct the $G1K/DATA hierarchy. The shell script
  updateDataDirectory10 calls this 10 times in parallel for speed.

collectSlxFastq.pl
  Collect a list of all likely-looking fastq files in repository directories.
  Saves making lots of slow calls to repository in updateDataDirectory.pl

RECALIBRATION

pairUpFastqFiles.pl
  Merge separate files for the two ends in e.g. Broad fastq files into single
  paired-end files (Sanger does this anyway and g1k.pl requires it). You
  need to do this whenever new Broad (and other?) fastq data is shipped in. 

makeQualitiesBayesian.pl
  Filter for quality files to decide qualities using uniform Bayesian prior

applyQualityMap.pl
  Called from g1k.pl cal to apply qualmapBayesian.txt to an uncalibrated fastq file 
  to create recalibrated.fastq.gz.

linkToRelatedQualmaps.pl
  When quality maps don't get made because of insufficient good data in a lane,
  this script tries to rescue the situation by hard-linking to a quality map
  for another lane in the same run. Lanes in same run tend to have very
  similar quality behaviour. (You will need to "sh" the output).

combineQualmaps.pl
  Combine lots of quality maps into one to get an impression of how a set of
  lanes (e.g. all those from a given centre) are behaving. Diagnostic only --
  not part of the main pipeline.

copyRecalibratedToRepository.pl
  Output commands to copy recalibrated.fastq.gz to /nfs/repository/.../...recal.fastq.gz.
  Then someone (e.g. Tom) actually does the copying. There has to be this
  manual intervention because directories can be full.

qualityAndMappingTable.pl
  Generate a spreadsheet for Richard to pass around the DATARELEASE group.

extractEntropyAndQualities.pl
fastqEntropy.pl
  More quality-checking code (looks at quality distributions and also
  calculates "discounted entropy" as defined in my mail to DATARELEASE).  






