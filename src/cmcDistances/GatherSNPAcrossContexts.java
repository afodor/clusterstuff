package cmcDistances;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

import coPhylog.ContextCount;
import utils.Translate;

public class GatherSNPAcrossContexts
{
	public static final int CUTOFF = 10;
	
	private static final File RESULTS_SPREADSHEET 
		= new File("/nobackup/afodor_research/fromNury52016/bothSnpsAt" + CUTOFF + ".txt");
			
	private static void writeResults( HashSet<String> set ) throws Exception
	{
		List<String> list = new ArrayList<String>(set);
		Collections.sort(list);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(RESULTS_SPREADSHEET));
		
		writer.write("genome");
		
		for(String s : list)
			writer.write("\t" + s);
		
		writer.write("\n");
		
		for(String s : WriteContextScripts.OUTPUT_DIRECTORY.list())
			if( s.endsWith("context"))
			{
				writer.write(s.replaceAll(".context", ""));
				
				File inFile = new File(WriteContextScripts.OUTPUT_DIRECTORY + 
							File.separator + s);
				System.out.println("reading " + inFile.getAbsolutePath());
				HashMap<Long, ContextCount> map = 
						WriteSNPFile.readFileRequireMin(inFile, 0);
				
				for(String s2 : list)
				{
					Long aLong = Encode.makeLong(s2);
					ContextCount cc = map.get(aLong);
					
					if( cc == null)
					{
						cc = map.get(Encode.makeLong(Translate.reverseTranscribe(s2)));
						
						if( cc != null)
							cc.reverseTranscribe();
					}
						
					if( cc== null)
						writer.write("\t-");
					else
						writer.write("\t" + cc.toString());
					
				}
				
				writer.write("\n");
				writer.flush();
			}
		
		writer.flush();  writer.close();
	}
	
	
	public static void main(String[] args) throws Exception
	{
		HashSet<String> set = gatherSNPs();
		
		for(String s : set)
			System.out.println(s);
		
		writeResults(set);
	}
	
	
	private static ContextCount parseToken( String s) throws Exception
	{
		s = s.replace("[", "").replace("]", "");
		
		StringTokenizer sToken = new StringTokenizer(s, ",");
		
		ContextCount cc = new ContextCount( Integer.parseInt(sToken.nextToken()),
				Integer.parseInt(sToken.nextToken()),Integer.parseInt(sToken.nextToken()),
				Integer.parseInt(sToken.nextToken()));
		
		if( sToken.hasMoreTokens())
			throw new Exception("Parsing error " + s);
		
		return cc;
	}
	
	private static HashSet<String> gatherSNPs() throws Exception
	{
		HashSet<String> set = new HashSet<String>();
		
		String[] files = WriteDistanceScripts.SNP_DIRECTORY.list();
		
		for(String s : files)
			if( s.endsWith(".txt") && s.indexOf("_023.txt") == -1)
			{
				BufferedReader reader = new BufferedReader(new FileReader(
					new File(WriteDistanceScripts.SNP_DIRECTORY.getAbsolutePath() + 
							File.separator + s)));
				
				reader.readLine();  reader.readLine();  reader.readLine();
				
				for(String s2 = reader.readLine(); s2 != null; s2=reader.readLine())
				{
					String[] splits =s2.split("\t");
					
					if( splits.length != 4)
						throw new Exception("Unexpected input " + s);
					
					ContextCount cc1 = parseToken(splits[1]);
					ContextCount cc2 = parseToken(splits[2]);
					
					if( cc1.getMax() >= CUTOFF && 
							cc2.getMax() >= CUTOFF)
					{
						String seq = Encode.getKmer( Long.parseLong(splits[0]) , 
								WriteSNPFile.KMER_SIZE);
						
						String flip = Translate.reverseTranscribe(seq);
						
						if( ! set.contains(seq) && ! set.contains(flip) )
						{
							set.add(seq);
						}
					}
				}
				
				reader.close();
			}
		
		return set;
	}
}
