package creOrthologs.automatedDistanceMatrix;

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

public class PivotDistanceMatrices
{
	private static File distanceMatrix = new File("/nobackup/afodor_research/af_broad/distanceMatrices");
	
	public static void main(String[] args) throws Exception
	{
		HashMap<String, HashMap<String,Integer>>  map = parseAll();
		writeResults(map);
	}
	
	private static void writeResults(HashMap<String, HashMap<String,Integer>>  map) throws Exception
	{
		System.out.println("Writing");
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/pivotedTo11.txt")));
		
		HashSet<String> innerKeys = new HashSet<String>();
		
		for(HashMap<String, Integer> innerMap : map.values())
		{
			innerKeys.addAll(innerMap.keySet());
		}
		
		List<String> keyList = new ArrayList<String>(innerKeys);
		Collections.sort(keyList);
		
		writer.write("region");
		
		for(String s : keyList)
			writer.write("\t" + s);
		
		writer.write("\n");
		
		for(String s : map.keySet())
		{
			writer.write(s);
			
			HashMap<String, Integer> innerMap = map.get(s);
			
			for(String s2 : keyList)
			{
				Integer val = innerMap.get(s2);
				
				if( val == null)
					writer.write("\tNA");
				else
					writer.write("\t" + val);
			}
			
			writer.write("\n");
		}
		
		writer.flush();  writer.close();
	}
	
	// outer key is region (e.g. 7000000220927531_3367_4032 ) 
	// inner key is genome1@genome2 ; inner value is count
	private static HashMap<String, HashMap<String,Integer>> parseAll() throws Exception
	{
		HashMap<String, HashMap<String,Integer>> map = 
				new HashMap<String, HashMap<String,Integer>>();
		
		String[] files = distanceMatrix.list();
		
		int numDone =1;
		for(String s : files)
		{
			numDone++;
			System.out.println("reading " +  numDone + " "+ files.length + " "+  s);
			if( s.startsWith("klebsiella_pneumoniae_chs_11.0.scaffolds.fasta_"))
			{
				String key = s.replace("klebsiella_pneumoniae_chs_11.0.scaffolds.fasta_", "");
				
				if( map.containsKey(key))
					throw new Exception("No");
				
				HashMap<String, Integer> innerMap = new HashMap<String,Integer>();
				
				map.put(key,innerMap);
				
				BufferedReader reader = new BufferedReader(new FileReader(new File(
					distanceMatrix.getAbsolutePath() + File.separator + s	)));
				
				reader.readLine();
				
				for(String s2 = reader.readLine(); s2 != null; s2 =reader.readLine())
				{
					String[] splits =s2.split("\t");
					
					if( s2.length() != 3)
						throw new Exception("No");
					
					List<String> list = new ArrayList<String>();
					list.add(splits[0]);
					list.add(splits[1]);
					Collections.sort(list);
					
					String innerKey = list.get(0) + "@" + list.get(1);
					
					if( innerMap.containsKey(innerKey))
						throw new Exception("No");
					
					innerMap.put(innerKey, Integer.parseInt(splits[2]));
				}
				
				reader.close();
			}
		}
		
		return map;
	}
}
