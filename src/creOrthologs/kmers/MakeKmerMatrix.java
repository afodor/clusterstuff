package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class MakeKmerMatrix
{
	private static List<String> getAllNames() throws Exception
	{
		List<String> list = new ArrayList<String>();
		String[] names = MakeKmers.KMER_DIR.list();
		
		for(String s : names)
		{
			if( s.endsWith("_kmers.txt"))
				list.add(s);
		}
		
		return list;
	}
	
	/*
	 * This currently requires too much RAM...
	 */
	public static void main(String[] args) throws Exception
	{
		HashMap<String, Integer[]> countMap = getCountMap();
		writeResults(countMap);
	}
	
	private static void writeResults( HashMap<String, Integer[]> map) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(
				"/nobackup/afodor_research/af_broad/allKmers.txt"));
		
		List<String> list = getAllNames();
		
		writer.write("kmer");
		
		for(String s  : list )
			writer.write("\t" + s);
		
		writer.write("\n");
		
		for(String s : map.keySet())
		{
			writer.write(s);
			
			for( Integer i : map.get(s))
				writer.write("\t" + i);
			
			writer.write("\n");
		}
		
		writer.flush();  writer.close();
				
	}
	
	// key is the kmer
	private static HashMap<String, Integer[]> getCountMap() throws Exception
	{
		HashMap<String, Integer[]> map = new HashMap<String,Integer[]>();
		int listLength= getAllNames().size();
		int index=0;
		
		String[] names = MakeKmers.KMER_DIR.list();
		
		for(String s : names)
		{
			if( s.endsWith("_kmers.txt"))
			{
				BufferedReader reader = new BufferedReader(new FileReader(new File(
					MakeKmers.KMER_DIR.getAbsolutePath() + File.separator + s	)));
				
				for(String s2 = reader.readLine(); s2 != null; s2 = reader.readLine())
				{
					String[] splits = s2.split("\t");
					if( splits.length != 2)
						throw new Exception("No");
					
					Integer[] counts = map.get(splits[0]);
					
					if( counts == null)
					{
						counts =  new Integer[listLength];
						
						for(int x=0; x < listLength; x++)
							counts[x] = 0;
						
						map.put(splits[0], counts);
					}
					
					counts[index] = counts[index] + 1;
				}
				
				reader.close();
				
				index++;
			}
				
		}
		
		return map;
		
	}
}
