package biolockJTestHardCoded;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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
		HashMap<String, HashMap<String,Long>> map = 
				getCountsForDirectory(new File(IN_RDP_DIRECTORY), "family");
		
		checkForDirectory(new File(BULK_RDP_DIRECTORY), map);
	}
	
	private static void checkForDirectory(File inDir,HashMap<String, HashMap<String,Long>> map ) throws Exception
	{
		String[] files = inDir.list();
		
		for(String s : files )
		{
			if( s.endsWith(".tsv"))
			{
				File inFile =new File(inDir.getAbsoluteFile() + File.separator + s);
				
				String aName = inFile.getName();
				
				aName =aName.substring(aName.lastIndexOf("_")+1, aName.length()-4);
				
				if( ! map.containsKey(aName) )
					throw new Exception("Could not find " + aName + " in " + map.keySet());
				
				HashMap<String,Long> innerMap = map.get(aName);
				assertEqual(inFile, innerMap);
			}
		}
	}
	
	private static void assertEqual( File inFile, HashMap<String,Long> innerMap ) throws Exception
	{
		System.out.println("check " + inFile.getAbsolutePath());
		BufferedReader reader = new BufferedReader(new FileReader(inFile));
		
		for(String s= reader.readLine(); s!= null; s= reader.readLine() )
		{
			StringTokenizer sToken = new StringTokenizer(s, "\t");
			
			if( sToken.countTokens() != 2)
				throw new Exception("Expecting two tokens");
		}
		
		System.out.println("Ok");
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
				if( node.getTaxaName().toLowerCase().indexOf("unclassified") == -1)  
				{
					Long count =map.get(node.getTaxaName());
					
					if( count == null)
						count =0l;
					
					count = count +1;
					
					map.put(node.getTaxaName(), count);
				}
				
			}
		}
		
		return map;
		
	}
	
}
