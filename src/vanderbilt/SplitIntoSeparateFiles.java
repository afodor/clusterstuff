package vanderbilt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class SplitIntoSeparateFiles
{
	static final File outDir = new File( "/projects/afodor/vanderbilt/split16S");
	
	public static void main(String[] args) throws Exception
	{
		HashMap<String, BufferedWriter> fileMap = new HashMap<String, BufferedWriter>();
		split("/projects/afodor/vanderbilt/VanderbiltSequences_Dec52014/MSHR1/seqs_all1.fna", fileMap, "all1");
		split("/projects/afodor/vanderbilt/VanderbiltSequences_Dec52014/MSHR2/seqs_all2.fna", fileMap, "all2");
		
		for( BufferedWriter writer : fileMap.values())
		{
			writer.flush();  writer.close();
		}
	}
	
	private static void split(String fileToSplit, HashMap<String, BufferedWriter> map, String suffix) throws Exception
	{
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(fileToSplit);
		
		int numDone =0;
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence() )
		{
			String sampleID = fs.getFirstTokenOfHeader().split("_")[0] + "_" + suffix;
			
			BufferedWriter writer= map.get(sampleID);
			
			if( writer== null)
			{
				writer= new BufferedWriter(new FileWriter(outDir + File.separator+ sampleID + ".fasta" ));
				map.put(sampleID, writer);
			}
			
			writer.write(fs.getHeader() + "\n");
			writer.write(fs.getSequence() + "\n");
			
			numDone++;
			
			if( numDone % 100000 == 0)
			{
				System.out.println(numDone + " " + map.size());
			}
		}
	}
}
