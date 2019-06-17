package farnazKrakenExport;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.StringTokenizer;


public class ExportKrakenTable
{
	public static final String[] LEVELS = { "phylum", "class", "order", "family", "genus", "species" };
	
	
	private static final String[] IN_DIRS = 
		{
				"/scratch/afodor_research/blj_pipeline/ffouladi/emilyNovSeq_2019May31/02_Kraken2Classifier/output",
				"/scratch/afodor_research/blj_pipeline/ffouladi/bariatricSurgery_March19_2019Mar09/02_Kraken2Classifier/output",
				"/scratch/afodor_research/blj_pipeline/ffouladi/nrap2_gut_wgs_Kraken2_2019Jan29/3_Kraken2Classifier/output"	
		};
	
	private static final String[] PROJECT_NAMES = 
		{
				"emily", "surgery", "tanya"
		};
	
	private static final String OUTPUT_DIR = "/users/afodor/farnazOutput";
	
	public static void main(String[] args) throws Exception
	{
		for( String level : LEVELS)
		{
			for( int x=0; x < IN_DIRS.length; x++)
			{
				System.out.println(level +  " " + PROJECT_NAMES[x]);
				File krakenOutoutDir = new File( IN_DIRS[x]);
				
				File outputDirectory = new File( OUTPUT_DIR + File.separator + PROJECT_NAMES[x] );
				
				outputDirectory.mkdirs();
				
				HashMap<String, HashMap<String, Long>> map = getMapForLevel(krakenOutoutDir, level);
				 
				File pivotedFile = new File(outputDirectory + File.separator + PROJECT_NAMES[x] + "_kraken_" + level + ".txt");
				 
				 writeResults(map, pivotedFile.getAbsolutePath());
			}
		}
	}
	
	private static HashMap<String, HashMap<String, Long>> getMapForLevel( File aDir, String level) throws Exception
	{
		HashMap<String, HashMap<String, Long>> map = new LinkedHashMap<>();
		
		String[] files = aDir.list();
		
		for(String fileName: files)
		{
			if(fileName.endsWith("_reported.tsv"))
			{
				String sampleName = fileName.replace("_reported.tsv", "");
				
				if( map.containsKey(sampleName))
					throw new Exception("No");
				
				File file = new File(aDir.getAbsolutePath() + File.separator + fileName);
				
				HashMap<String, Long> innerMap =getExpectedAtLevel(file, "" + level.charAt(0));
				
				map.put(sampleName, innerMap);
			}
		}
		
		return map;
	}
	
	public static HashMap<String, Long> getExpectedAtLevel(File file, String level) throws Exception
	{
		HashMap<String, Long> map = new LinkedHashMap<>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		for(String s = reader.readLine();  s != null; s = reader.readLine())
		{
			String[] splits =s.split("\t");
			
			if( splits.length != 2)
			{
				throw new Exception("Expecting 2");
			}
			
			String taxa = splits[0];
			taxa = taxa.substring(taxa.lastIndexOf("|") +1, taxa.length());
			
			if( taxa.startsWith(level + "__"))
			{
				taxa = taxa.replace(level + "__", "");
				
				if( map.containsKey(taxa))
					throw new Exception("Duplicate");
				
				map.put(taxa, Long.parseLong(splits[1]));
			}
		}
		
		return map;
	}
	
	public static void writeResults(HashMap<String, HashMap<String, Long>>  map, String filepath ) 
			throws Exception
{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(filepath)));
		
		writer.write("sample");
		List<String> otuList = getOTUSAtThreshold(map, 0);
		
		for( String s : otuList)
		writer.write("\t" +  s);
		
		writer.write("\n");
		
		for( String s : map.keySet())
		{
		//String expandedString = PivotRDPs.getExpandedString( s);
		//writer.write( expandedString );
		writer.write(s);
		
		HashMap<String, Long> innerMap = map.get(s);
		
		for( String otu : otuList)
		{
		Long aVal = innerMap.get(otu);
			
		if(aVal== null)
			aVal = 0l;
			
		writer.write("\t" + aVal );
		}
		
		writer.write("\n");
		}
		
		writer.flush();  writer.close();
		}

	

	private static List<String> getOTUSAtThreshold(
			HashMap<String, HashMap<String, Long>>  map,
									int threshold) throws Exception
	{
		
		HashMap<String, Long> countMap = new HashMap<String, Long>();
		
		for( String s: map.keySet() )
		{
			HashMap<String, Long> innerMap = map.get(s);
				
			for(String possibleOtu : innerMap.keySet())
			{
				Long oldCount = countMap.get(possibleOtu);
					
				if(oldCount == null)
						oldCount = 0l;
					
				oldCount += innerMap.get(possibleOtu);
					
				countMap.put(possibleOtu, oldCount);
			}
		}
			
		List<String> otuList= new ArrayList<String>();
		
		for( String s : countMap.keySet() )
			if( countMap.get(s) >= threshold )
				otuList.add(s);
		
		return otuList;
	
	}
	

	
}
