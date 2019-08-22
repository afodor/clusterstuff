package tb_Dada2Map;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class BreakUpFastaFile
{
	public static String IN_FILE_PATH = "/nobackup/afodor_research/datasets/maigaTB_June2019/dada2Fasta/dada2_fastaSeqs.txt";
	
	public static String OUT_FILE_DIR = "/nobackup/afodor_research/datasets/maigaTB_June2019/dada2Fasta/queryFiles";
	
	public static void main(String[] args) throws Exception
	{
		HashMap<Integer, BufferedWriter> writerMap = new HashMap<Integer, BufferedWriter>();
		
		for( int x=0; x < 100; x++)
		{
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				OUT_FILE_DIR + File.separator + "out_" + x + ".fasta"	)));
			
			writerMap.put(x, writer);
		}
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(IN_FILE_PATH);
		
		int index =0;
		
		for(FastaSequence fs= fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			String s= fs.getHeader();
			int count = Integer.parseInt(s.substring(s.indexOf("=") +1, s.length()));
			
			if(count >= 100)
			{
				BufferedWriter writer = writerMap.get(index);
				
				writer.write(">" + fs.getHeader() + "\n");
				writer.write(fs.getSequence() + "\n");
				
				index++;
				
				if( index == 100)
					index =0;
			}
		}
		
		for( BufferedWriter writer : writerMap.values())
		{
			writer.flush(); writer.close();
		}
	}
}
