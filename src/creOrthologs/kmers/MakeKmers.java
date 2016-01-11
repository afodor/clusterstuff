package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class MakeKmers
{
	public static final int KMER_LENGTH = 12;
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 2)
		{
			System.out.println("Usage inFile outFile");
			System.exit(1);
		}
		
		HashMap<String, Integer> map = breakIntoKmers(new File(args[0]));
		writeResults(new File(args[1]), map);
		
	}
	
	private static boolean isACGT(String s)
	{
		
		for( int x=0; x < s.length(); x++)
		{
			char c = s.charAt(x);
			
			if( c != 'A' && c != 'C' && c != 'G' && c != 'T')
				return false;
		}
		
		return true;
	}
	
	private static void writeResults( File outFile, HashMap<String, Integer> map) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		
		for(String s : map.keySet())
			writer.write(s + "\t" + map.get(s) + "\n");
		
		writer.flush();  writer.close();
	}
	
	private static HashMap<String, Integer> breakIntoKmers(File inFile) throws Exception
	{
		HashMap<String, Integer> map = new HashMap<String,Integer>();
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(inFile);
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null;
						fs = fsoat.getNextSequence())
		{
			String seq =fs.getSequence();
			
			for( int x=0; x < seq.length()- KMER_LENGTH; x++)
			{
				String sub = seq.substring(x, x + KMER_LENGTH);
				
				if( isACGT(sub))
				{
					Integer count = map.get(sub);
					
					if( count == null)
						count =0;
					
					count++;
					
					map.put(sub,count);
				}
			}
		}
	
		fsoat.close();
		return map;
		
	}
}
