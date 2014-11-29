package urbanVsRural;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class WriteSequenceToSampleFile
{
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer =new BufferedWriter(new FileWriter("/projects/afodor/ChinaSequences/forwardReadToSample"));
		File topDir = new File("/projects/afodor/ChinaSequences/rdpResults");
		
		for(String s : topDir.list())
		{
			if( s.endsWith(".fasta.gz") && CountSeqs.getReadNumber(s) == 1)
			{
				System.out.println(s);
				String sampleId = s.substring(s.indexOf("gz") + 1);
				sampleId = sampleId.replaceAll(".fasta.gz", "").replaceAll("filepart_", "");
				
				FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(
						new File(topDir.getAbsolutePath() + File.separator + s), true);
				
				for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
				{
					writer.write(fs.getFirstTokenOfHeader() + "\t" + sampleId + "\n");
				}
				
				writer.flush();  
			}
		}
		
		writer.flush(); writer.close();
	}
}
