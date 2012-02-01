package Pathogens::RNASeq::Exceptions;

use Exception::Class (
    Pathogens::RNASeq::Exceptions::FailedToOpenAlignmentSlice => { description => 'Couldnt get reads from alignment slice. Error with Samtools or BAM' },
    Pathogens::RNASeq::Exceptions::FailedToOpenExpressionResultsSpreadsheetForWriting => { description => 'Couldnt write out the results for expression' },
		Pathogens::RNASeq::Exceptions::InvalidInputFiles => { description => 'Invalid inputs, sequence names or lengths are incorrect' },
		Pathogens::RNASeq::Exceptions::FailedToCreateNewBAM => { description => 'Couldnt create a new bam file' },
		Pathogens::RNASeq::Exceptions::FailedToCreateMpileup => { description => 'Couldnt create an mpileup' }
);

1;
