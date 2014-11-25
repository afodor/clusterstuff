package urbanVsRural;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import parsers.FastaSequence;

public class RarifyDown
{
	private static final int RARIFICATION_NUMBER = 58977;
	private static final Random RANDOM = new Random(323021);
	
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
				"/projects/afodor/ChinaSequences/rariedForward.txt")));
		
		File topDir = new File("/projects/afodor/ChinaSequences/rdpResults");
		
		String[] list = topDir.list();
		
		for(String s : list)
		{
			if( s.endsWith(".fasta") && getReadNumber(s) == 1)
			{
				List<FastaSequence> fastaList= 
						FastaSequence.readFastaFile(topDir.getAbsolutePath() + File.separator + 
								s);
				
				
				Collections.shuffle(fastaList,RANDOM);
				
				for(int x=0; x < RARIFICATION_NUMBER; x++)
				{
					FastaSequence fs = fastaList.get(x);
					writer.write(fs.getHeader() + "\n");
					writer.write(fs.getSequence() + "\n");
				}
				writer.flush();
				System.out.println(s);
			}
			
			writer.flush();  
		}
		
		writer.flush();  writer.close();
	}
}
