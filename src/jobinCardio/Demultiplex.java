package jobinCardio;

import java.io.File;

public class Demultiplex
{
	public static File forwardSeqs = new File(
			"/projects/afodor_research/JobinCardio/*/Sample_1-26241228/Data/Intensities/BaseCalls/Sample-Name-1_S1_L001_R1_001.fastq.gz");
	
	public static File ReverseSeqs = new File(
			"/projects/afodor_research/JobinCardio/*/Sample_1-26241228/Data/Intensities/BaseCalls/Sample-Name-1_S1_L001_R1_001.fastq.gz");
	
	public static File keyFile = new File(
			"/projects/afodor_research/JobinCardio/keys/barcode Run2 5-30-2015.txt");
	
}
