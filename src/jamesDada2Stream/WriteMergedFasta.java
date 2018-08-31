package jamesDada2Stream;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.StringTokenizer;


public class WriteMergedFasta
{
	private static String TOP_DIR = "/scratch/afodor_research/datasets/uegp/raw_Stream_16S/forward_filtered/APR";

	public static void main(String[] args) throws Exception
	{
		HashSet<String> set = getSet();
		
		System.out.println("Got " + set.size() + " unique sequences");
		
		
		BufferedWriter writer = new BufferedWriter(new  FileWriter(
				new File("/scratch/afodor_research/datasets/uegp/dada2ToBlast/collectedSeqVariants.txt")));
		
		int index =1;
		
		for(String s: set)
		{
			writer.write(">" + index + "\n");
			index++;
			writer.write(s + "\n");
		}
		
		writer.flush(); writer.close();
	}
	
	public static HashSet<String> getSet() throws Exception
	{
		HashSet<String> set =new LinkedHashSet<String>();
		
		File topDir = new File(TOP_DIR);
		
		for(String nextDir :topDir.list() )
		{
			File bottomDir = new File(topDir.getAbsolutePath() + File.separator + nextDir);
			
			if( bottomDir.isDirectory())
			{
				String[] list2 = bottomDir.list();
				
				for(String fileName : list2)
				{
					if( fileName.endsWith("_F.txt"))
					{
						File aFile = new File(bottomDir.getAbsolutePath() +File.separator + fileName);
						addToSet(set, aFile);
					}
				}
			}
		}
		
		return set;
	}
	
	private static void addToSet(HashSet<String> set, File file) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String firstLine = reader.readLine();
		StringTokenizer sToken = new StringTokenizer(firstLine);
		
		while(sToken.hasMoreTokens())
			set.add(sToken.nextToken().replaceAll("\"", ""));
		
		reader.close();
	}
}
