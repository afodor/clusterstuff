package jobinCardio;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.zip.GZIPInputStream;

import parsers.FastQ;
import utils.Translate;

public class Demultiplex
{
	public static File forwardSeqs = new File(
			"/projects/afodor_research/JobinCardio/*/Sample_1-26241228/Data/Intensities/BaseCalls/Sample-Name-1_S1_L001_R1_001.fastq.gz");
	
	public static File reverseSeqs = new File(
			"/projects/afodor_research/JobinCardio/*/Sample_1-26241228/Data/Intensities/BaseCalls/Sample-Name-1_S1_L001_R1_002.fastq.gz");
	
	public static File keyFile = new File(
			"/projects/afodor_research/JobinCardio/keys/barcode Run2 5-30-2015.txt");
	
	public static void main(String[] args) throws Exception
	{
		List<PrimerKeyLine> primerList = PrimerKeyLine.getList(keyFile.getAbsolutePath());
		
		BufferedReader forwardReader =
				new BufferedReader(new InputStreamReader( 
						new GZIPInputStream( new FileInputStream( forwardSeqs))));
		
		BufferedReader backwardsReader =
				new BufferedReader(new InputStreamReader( 
						new GZIPInputStream( new FileInputStream( reverseSeqs))));		
				
		int numForwardMatch =0;
		int numExamined =0;
		for( FastQ forward = FastQ.readOneOrNull(forwardReader); forward != null;
								forward = FastQ.readOneOrNull(forwardReader))
		{
			FastQ back = FastQ.readOneOrNull(backwardsReader);
			
			if(! forward.getFirstTokenOfHeader().equals(back.getFirstTokenOfHeader()))
				throw new Exception("No " + forward.getFirstTokenOfHeader() + " " + 
							back.getFirstTokenOfHeader());
			
			String revSeq = Translate.safeReverseTranscribe(back.getSequence());
			
			boolean gotOne = false;
			
			for( PrimerKeyLine pkl : primerList)
			{
				if( pkl.matchesForward(forward.getSequence()) &&
							pkl.matchesReverse(revSeq))
					gotOne = true;
			}
			
			numExamined++;
			
			if( gotOne)
				numForwardMatch++;
			
			System.out.println(numForwardMatch + " " + numExamined);
			
			if( numExamined == 10000)
				System.exit(1);
		}
	}
	
}
