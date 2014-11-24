package urbanVsRural;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

public class CountSeqs
{
	
	public static int getReadNumber(String s) throws Exception
	{
		int readNum = Integer.parseInt( s.charAt(s.indexOf("fq")-2) + "");
		
		if (readNum < 1 || readNum > 2)
			throw new Exception("No");
		
		return readNum;
		
	}
	
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File( 
				"/projects/afodor/ChinaSequences/counts.txt")));
		
		File topDir = new File("/projects/afodor/ChinaSequences/rdpResults");
		
		String[] list = topDir.list();
		
		for(String s : list)
		{
			if( s.endsWith(".fasta"))
			{
				List<FastaSequence> fastaList= 
						FastaSequence.readFastaFile(topDir.getAbsolutePath() + File.separator + 
								s);
				
				writer.write(s + "\t" + getReadNumber(s) + "\t" + fastaList.size() + "\n");
			}
		}
		
		writer.flush();  writer.close();
	}
}
