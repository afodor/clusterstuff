package creOrthologs.kmers;

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

public class GatherDistanceMatrix
{
	public static final File GATHERED_DIR = new File( "/nobackup/afodor_research/af_broad/gatheredKmerMatrices");
	
	public static void main(String[] args) throws Exception
	{
		writeResults(getDistances());
	}
	
	private static void writeResults(HashMap<String, Double> map) throws Exception
	{
		HashSet<String> set = new HashSet<String>();
		for(String s : map.keySet() )
		{
			String[] splits = s.split("@");
			
			if( splits.length != 2)
				throw new Exception("No");
			
			set.add(splits[0]);  set.add(splits[1]);
		}
		
		List<String> list = new ArrayList<String>(set);
		
		System.out.println("Finish with " + list.size());
		
		Collections.sort(list);
		
		
		int index = 0;
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(GATHERED_DIR.getAbsolutePath() + 
				File.separator + "allDist.txt"));
		
		BufferedWriter keyWriter = new BufferedWriter(new FileWriter(GATHERED_DIR.getAbsolutePath() + 
				File.separator + "allKey.txt"));
		
		String name = "BUG_" + index;
		
		while( name.length() < 10)
			name = name + "_";
		
		writer.write(list.size() + "\n");
		
		for( int x=0; x < list.size(); x++)
		{
			writer.write(name);
			keyWriter.write(name + " " + list.get(x) + "\n");
			
			for( int y=0; y < list.size(); y++)
			{
				String key = getKey(list.get(x), list.get(y));
				
				Double val = map.get(key);
				
				// some of the jobs failed
				// we picked up contrasts from duplicate runs (a is compared to b and b to a )
				// but this will miss some identities (where a is compared to a)
				// so this hack fixes those since the identities have distances of 0
				if( val == null)
				{
					if( list.get(x).equals(list.get(y)))
							val=0.0;
					else
						throw new Exception("No " + list.get(x) + " " + list.get(y));
				}
				
				writer.write(" " + val);
				//System.out.println(list.get(x) + " "+ list.get(y) + " " + val);
			}
			
			writer.write("\n");
		}
		
		writer.flush();  writer.close();
		keyWriter.flush();  keyWriter.close();
	}
	
	private static HashMap<String, Double> getDistances() throws Exception
	{
		String[] files = WriteDistanceMatrixOneVsAll.KMER_DISTANCE_MATRIX_DIR.list();
		
		HashMap<String, Double> map = new HashMap<String, Double>();
		
		for(String s : files)
		{
			BufferedReader reader = new BufferedReader(new FileReader(
					WriteDistanceMatrixOneVsAll.KMER_DISTANCE_MATRIX_DIR.getAbsolutePath() + 
					File.separator + s));
			
			for(String s2= reader.readLine(); s2 != null; s2= reader.readLine())
			{
				String[] splits = s2.split("\t");
				
				if( splits.length != 3)
					throw new Exception("No");
				
				String key = getKey(splits[0], splits[1]);
				
				Double oldVal = map.get(key);
				
				if( oldVal == null)
				{
					map.put(key, Double.parseDouble(splits[2]));
				}
				else
				{
					Double newVal = Double.parseDouble(splits[2]);
					
					if( Math.abs(newVal - oldVal) > 0.001)
						throw new Exception("No " + newVal + " " + oldVal + " " + s + " " + key);
				}
					
			}
			
			reader.close();
		
		}
		
		return map;
	}
	
	private static String getKey(String s1, String s2)
	{
		List<String> list = new ArrayList<String>();
		list.add(s1);
		list.add(s2);
		Collections.sort(list);
		
		return list.get(0) + "@" + list.get(1);
	}
}
