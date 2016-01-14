package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

import creOrthologs.RunBlastAll;
import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;
import utils.Translate;

public class MakeKmers
{
	public static final int KMER_LENGTH = 12;
	public static final File KMER_DIR = new File("/nobackup/afodor_research/af_broad/kmers");
	
	public static void main(String[] args) throws Exception
	{
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/nobackup/afodor_research/af_broad" + File.separator + d);
			
			String[] list = genomeDir.list();

			for( String s : list)
			{	
				if( s.endsWith("fasta"))
				{
					System.out.println(s);
					File inSeqs= new File( genomeDir.getAbsolutePath() + File.separator + s);
					HashMap<String, Integer> map = breakIntoKmers(inSeqs);
				 	
					File outFile = new File( KMER_DIR + File.separator + 
								s.replace(".scaffolds.fasta", "") + "_kmers.txt");
					
					writeResults(outFile, map);
				}
			}
		}
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
					{
						String reverse = Translate.reverseTranscribe(sub);
						
						count = map.get(reverse);
						
						if( count != null)
							sub = reverse;
					}
					
					if( count == null)
						count =0;
					
					count++;
					
					map.put(sub,count);
				}
			}
		}
	
		fsoat.close();
		
		for(String s : map.keySet())
			if( map.containsKey(Translate.reverseTranscribe(s)))
				throw new Exception("Logic error " + s + " " + Translate.reverseTranscribe(s));
		
		return map;
		
	}
}
