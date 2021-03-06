package swarm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.StringTokenizer;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;
import parsers.PivotOTUs;

public class PivotSwarmResults
{
	private static final int MIN_SEQUENCES = 0;
	
	public static void main(String[] args) throws Exception
	{
		HashMap<String, Integer> sequenceIDtoOTUMap = getSequenceIdToOTUMap(
				"/nobackup/afodor_research/topeOneAtATime/swarmOutWithSingles.txt",MIN_SEQUENCES,
				"/nobackup/afodor_research/topeOneAtATime/OTUPlusSeqIDsWithSingles.txt");
		
		System.out.println("Got " + sequenceIDtoOTUMap.size() + " OTUs");
		
		HashMap<String, Integer> rawSequenceToIDMap = 
				getRawSequenceToIDMap(sequenceIDtoOTUMap , 
						"/nobackup/afodor_research/topeOneAtATime/mergedForSwarmIncludingSingletons.txt");
		
		System.out.println("Got " + rawSequenceToIDMap.size() + " sequences ");
		
		// ok to garbage collect
		sequenceIDtoOTUMap = null;
		
		HashMap<String, HashMap<String,Integer>> countMap = 
				new HashMap<String, HashMap<String,Integer>>();
		
		addToCountMap( new File( "/nobackup/afodor_research/topeOneAtATime/file3/mergedOut"),
					countMap, rawSequenceToIDMap);
		
		addToCountMap( new File( "/nobackup/afodor_research/topeOneAtATime/file4/mergedOut"),
					countMap, rawSequenceToIDMap);
		
		PivotOTUs.writeResults(countMap, 
				"/nobackup/afodor_research/topeOneAtATime/swarmOTUsAsColumnsWithSingles.txt");
	}
	
	public static void addToCountMap(File directory, 
			HashMap<String, HashMap<String,Integer>> countMap,
			HashMap<String, Integer> rawSequenceToIDMap 
			) throws Exception
	{
		String[] files = directory.list();
		
		for(String s : files)
		{
			if( countMap.containsKey(s))
				throw new Exception("Parsing error");
			
			FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(directory.getAbsolutePath() + 
					File.separator + s );
			
			HashMap<String, Integer> innerMap = new HashMap<String,Integer>();
			
			for(FastaSequence fs = fsoat.getNextSequence(); fs != null; 
							fs = fsoat.getNextSequence())
			{
				Integer otuID = rawSequenceToIDMap.get(fs.getSequence());
				
				if( otuID != null)
				{
					String otuString = "OTU_" + otuID;
					
					Integer count = innerMap.get(otuString);
					
					if( count == null)
						count = 0;
					
					count++;
					
					innerMap.put(otuString, count);
				}
			}
			
			if( innerMap.size() > 0 )
				countMap.put(s, innerMap);
		}
	}
	
	public static HashMap<String, Integer> getRawSequenceToIDMap(
			HashMap<String, Integer> sequenceIdToOTUMap,
			String mergedFile) throws Exception
	{
		HashMap<String, Integer> map = new HashMap<String,Integer>();
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(mergedFile);
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			String key = fs.getFirstTokenOfHeader();
			
			Integer val = sequenceIdToOTUMap.get(key);
			
			if( val != null)
				map.put(fs.getSequence(), val);
		}
		
		return map;
	}
	
	private static int getCount(String line) throws Exception
	{
		int count =0;
		
		StringTokenizer sToken = new StringTokenizer(line);
		
		while(sToken.hasMoreTokens())
		{
			String[] splits = sToken.nextToken().split("_");
			
			if( splits.length != 2)
				throw new Exception("Parsing error " + line);
			
			count += Integer.parseInt(splits[1]);
		}
		
		return count;
	}
	
	public static HashMap<String, Integer> getSequenceIdToOTUMap(
			String swarmOutFile, int minNumSequences, String namedOtuFile) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(namedOtuFile));
		
		int index =1;
		HashMap<String, Integer> map =new HashMap<String,Integer>();
		
		BufferedReader reader = new BufferedReader(new FileReader(swarmOutFile));
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			if( getCount(s) >= minNumSequences )
			{
				writer.write(index + "");
				
				StringTokenizer sToken =new StringTokenizer(s);
				
				while( sToken.hasMoreTokens())
				{
					String seqId = sToken.nextToken();
					map.put(seqId, index);
					writer.write("\t" + seqId);
				}
					
				writer.write("\n");
				
				index++;
			}
		}
		
		reader.close();
		
		writer.flush();  writer.close();
		
		return map;
	}
	
}
