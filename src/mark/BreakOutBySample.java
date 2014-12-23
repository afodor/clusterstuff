package mark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.StringTokenizer;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class BreakOutBySample
{
	private static final File IN_SEQUENCE_DIR = 
			new File("/projects/afodor_research/mark/sequences/mrg_take2_l25_r25_CR");
	
	private static final File OUT_SEQUENCE_DIR = 
			new File("/projects/afodor_research/mark/fastaBySample");
	
	
	public static void main(String[] args) throws Exception
	{
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime( 
			IN_SEQUENCE_DIR.getAbsolutePath() + File.separator + "mrg_j200p5_q20_take2_l25_r25.fna"	);
		
		HashMap<String, BufferedWriter> fileMap = new HashMap<String, BufferedWriter>();
		
		int numDone =0;
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			String sampleID = new StringTokenizer(fs.getFirstTokenOfHeader(), "_").nextToken();
			
			BufferedWriter writer = fileMap.get(sampleID);
			
			if( writer == null)
			{
				writer = new BufferedWriter(new FileWriter(OUT_SEQUENCE_DIR.getAbsolutePath() + File.separator + 
						sampleID));
				fileMap.put(sampleID, writer);
				
			}
			
			writer.write(">" + fs.getFirstTokenOfHeader() + "\n");
			writer.write(fs.getSequence() + "\n");
			
			numDone++;
			
			if( numDone % 10000 == 0 )
				System.out.println(numDone + " " + fileMap.size());
		}
		
		for(BufferedWriter w: fileMap.values())
		{
			w.flush(); w.close();
		}
	}
}
