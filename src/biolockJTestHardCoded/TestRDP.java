package biolockJTestHardCoded;

import java.io.File;
import java.util.HashMap;
import java.util.StringTokenizer;

import parsers.NewRDPNode;
import parsers.NewRDPParserFileLine;

public class TestRDP
{
	public static final String IN_SEQUENCE_DIRECTORY = "/scratch/afodor_research/blj_pipeline/ieclabau/rhizosphere_5_2019Jun02/01_PearMergeReads/output";
	
	public static final String IN_RDP_DIRECTORY = "/scratch/afodor_research/blj_pipeline/ieclabau/rhizosphere_5_2019Jun02/02_RdpClassifier/output";
	
	public static final String BULK_RDP_DIRECTORY = "/scratch/afodor_research/blj_pipeline/ieclabau/rhizosphere_5_2019Jun02/03_RdpParser/output";
	
	public static final int THRESHOLD = 80;
	
	public static void main(String[] args) throws Exception
	{
		getCountsForDirectory(new File(IN_RDP_DIRECTORY), "family");
	}
	
	private static HashMap<String, HashMap<String,Long>> getCountsForDirectory(File inDir, String level) throws Exception
	{
		HashMap<String, HashMap<String,Long>>  map = new HashMap<>();
		
		for(String fileName : inDir.list())
		{
			File file = new File( inDir.getAbsoluteFile() + File.separator + fileName);
			
			HashMap<String, Long> innerMap = getCountAtLevel(file, level);
			
			String aName = file.getName();
			
			aName = new StringTokenizer(aName, "_").nextToken();
			
			System.out.println("Got " + aName + " "+  innerMap.size() + " "+  level);
			
			
			map.put(aName, innerMap);
		}
		
		return map;
	}
	
	private static HashMap<String, Long> getCountAtLevel(File file, String level) throws Exception
	{
		HashMap<String, Long>  map = new HashMap<>();
		
		HashMap<String, NewRDPParserFileLine> rdpMap = 
		NewRDPParserFileLine.getAsMapFromSingleThread(file);
		
		for(NewRDPParserFileLine rdpLine : rdpMap.values())
		{
			NewRDPNode node = rdpLine.getTaxaMap().get(level);
			
			if( node != null &&  node.getScore() >= THRESHOLD )
			{
				Long count =map.get(node.getTaxaName());
				
				if( count == null)
					count =0l;
				
				count = count +1;
				
				map.put(node.getTaxaName(), count);
			}
		}
		
		return map;
		
	}
	
}
