package abundantOTU;

import java.io.File;
import java.util.HashMap;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;
import parsers.PivotOTUs;

public class PivotAbundantOTUResults
{
	public static void main(String[] args) throws Exception
	{
		HashMap<String, String> sequenceIDtoOTUMap = 
				AbundOTUClustParser.getSequenceToOtuMap(
						"/nobackup/afodor_research/topeOneAtATime/abundantOTU/clust.clust");
		
		System.out.println("Got " + sequenceIDtoOTUMap.size() + " OTUs");
		
		HashMap<String, String> rawSequenceToIDMap = 
				getRawSequenceToIDMap(sequenceIDtoOTUMap , 
						"/nobackup/afodor_research/topeOneAtATime/mergedForSwarm.txt");
		
		System.out.println("Got " + rawSequenceToIDMap.size() + " sequences ");
		
		// ok to garbage collect
		sequenceIDtoOTUMap = null;
		
		HashMap<String, HashMap<String,Integer>> countMap = 
				new HashMap<String, HashMap<String,Integer>>();
		
		addToCountMap( new File( "/nobackup/afodor_research/topeOneAtATime/file3/mergedOut"),
					countMap, rawSequenceToIDMap);
		
		addToCountMap( new File( "/nobackup/afodor_research/topeOneAtATime/file4/mergedOut"),
					countMap, rawSequenceToIDMap);
		
		PivotOTUs.writeResults(countMap, "/nobackup/afodor_research/topeOneAtATime/abundantOTU/abundantOTUsAsColumns.txt");
	}
	
	public static void addToCountMap(File directory, 
			HashMap<String, HashMap<String,Integer>> countMap,
			HashMap<String, String> rawSequenceToIDMap 
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
				String otuID = rawSequenceToIDMap.get(fs.getSequence());
				
				if( otuID != null)
				{
					Integer count = innerMap.get(otuID);
					
					if( count == null)
						count = 0;
					
					count++;
					
					innerMap.put(otuID, count);
				}
			}
			
			if( innerMap.size() > 0 )
				countMap.put(s, innerMap);
		}
	}
	
	public static HashMap<String, String> getRawSequenceToIDMap(
			HashMap<String, String> sequenceIdToOTUMap,
			String mergedFile) throws Exception
	{
		HashMap<String, String> map = new HashMap<String,String>();
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(mergedFile);
		
		for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			String key = fs.getFirstTokenOfHeader();
			
			String val = sequenceIdToOTUMap.get(key);
			
			if( val != null)
				map.put(fs.getSequence(), val);
		}
		
		return map;
	}
}
