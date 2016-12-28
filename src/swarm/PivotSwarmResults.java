package swarm;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.StringTokenizer;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class PivotSwarmResults
{
	public static void main(String[] args) throws Exception
	{
		HashMap<String, Integer> sequenceIDtoOTUMap = getSequenceIdToOTUMap(
				"/nobackup/afodor_research/topeOneAtATime/swarmOut.txt",10);
		
		System.out.println("Got " + sequenceIDtoOTUMap.size() + " OTUs");
		
		HashMap<String, Integer> rawSequenceToIDMap = 
				getRawSequenceToIDMap(sequenceIDtoOTUMap , 
						"/nobackup/afodor_research/topeOneAtATime/mergedForSwarm.txt");
		
		System.out.println("Got " + rawSequenceToIDMap.size() + " sequences ");
		
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
	
	public static HashMap<String, Integer> getSequenceIdToOTUMap(
			String swarmOutFile, int minNumSequences) throws Exception
	{
		
		int index =1;
		HashMap<String, Integer> map =new HashMap<String,Integer>();
		
		BufferedReader reader = new BufferedReader(new FileReader(swarmOutFile));
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			StringTokenizer sToken =new StringTokenizer(s);
			
			if( sToken.countTokens() >= minNumSequences )
			{
				while( sToken.hasMoreTokens())
					map.put(sToken.nextToken(), index);
				
				index++;
			}
		}
		
		reader.close();
		
		return map;
	}
	
}
