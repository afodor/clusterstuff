package jobinCardioRun1;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.List;
import java.util.zip.GZIPInputStream;
import parsers.FastQ;

public class Demultiplex
{
	public static File forwardSeqs = new File(
			"/projects/afodor_research/JobinCardio/*/Sample_1-26241228/Data/Intensities/BaseCalls/Sample-Name-1_S1_L001_R1_001.fastq.gz");
	
	public static File reverseSeqs = new File(
			"/projects/afodor_research/JobinCardio/*/Sample_1-26241228/Data/Intensities/BaseCalls/Sample-Name-1_S1_L001_R2_001.fastq.gz");
	
	public static File keyFile = new File(
			"/projects/afodor_research/JobinCardio/keys/barcode Run2 5-30-2015.txt");
	
	public static File fastaOutDir = new File("/projects/afodor_research/JobinCardio/fastaOut");
	
	public static void main(String[] args) throws Exception
	{
		HashMap<Integer, Holder> outMap = new HashMap<Integer, Holder>();
		
		List<PrimerKeyLine> primerList = PrimerKeyLine.getList(keyFile.getAbsolutePath());
		
		BufferedReader forwardReader =
				new BufferedReader(new InputStreamReader( 
						new GZIPInputStream( new FileInputStream( forwardSeqs))));
		
		BufferedReader backwardsReader =
				new BufferedReader(new InputStreamReader( 
						new GZIPInputStream( new FileInputStream( reverseSeqs))));		
				
		long numMatched=0;
		long numExamined =0;
		for( FastQ forward = FastQ.readOneOrNull(forwardReader); forward != null;
								forward = FastQ.readOneOrNull(forwardReader))
		{
			FastQ back = FastQ.readOneOrNull(backwardsReader);
			
			if(! forward.getFirstTokenOfHeader().equals(back.getFirstTokenOfHeader()))
				throw new Exception("No " + forward.getFirstTokenOfHeader() + " " + 
							back.getFirstTokenOfHeader());
			
			//String revSeq = Translate.safeReverseTranscribe(back.getSequence());
			
			PrimerKeyLine matchedPkl = null;
			
			for( PrimerKeyLine pkl : primerList)
			{
				if( pkl.matchesForward(forward.getSequence()) &&
							pkl.matchesReverse(back.getSequence()))
				{
					if (matchedPkl != null)
						throw new Exception("Double match! " + forward.getFirstTokenOfHeader());
					
					matchedPkl = pkl;
				}
			}
			
			numExamined++;
			
			if( matchedPkl != null)
			{
				numMatched++;
				
				Holder h = outMap.get(matchedPkl.getSampleIndex());
				
				if( h==  null)
				{
					h = new Holder(matchedPkl.getSampleIndex());
					h.initWriters();
					outMap.put(matchedPkl.getSampleIndex(), h);
				}
				
				h.forwardWriter.write(">index_" + matchedPkl.getSampleIndex() + "_R1_" + numExamined + "\n");
				h.forwardWriter.write(forward.getSequence() + "\n");
				h.backwardsWriter.write(">index_" + matchedPkl.getSampleIndex() + "_R2_" + numExamined + "\n");
				h.backwardsWriter.write(back.getSequence() + "\n");
			}
				
			if( numExamined % 10000 == 0 )
				System.out.println(numExamined + " " + numMatched);
		}
		
		for( Holder h : outMap.values())
			h.closeWriters();
		
		System.out.println("finished " + numExamined + " " + numMatched);
	}
	
	private static class Holder
	{
		BufferedWriter forwardWriter ;
		BufferedWriter backwardsWriter ;
		
		final int experimentNum;
		
		Holder(int experimentNum)
		{
			this.experimentNum = experimentNum;
		}
		
		void initWriters() throws Exception
		{
			this.forwardWriter = new BufferedWriter(new FileWriter(new File(
					fastaOutDir.getAbsolutePath() + File.separator + "sample_" + experimentNum + "_R1.txt")));
			
			this.backwardsWriter= new BufferedWriter(new FileWriter(new File(
					fastaOutDir.getAbsolutePath() + File.separator + "sample_" + experimentNum + "_R2.txt")));
	
		}
		
		void closeWriters() throws Exception
		{
			this.forwardWriter.flush();
			this.forwardWriter.close();
			this.backwardsWriter.flush();
			this.backwardsWriter.close();
		}
	}
}
