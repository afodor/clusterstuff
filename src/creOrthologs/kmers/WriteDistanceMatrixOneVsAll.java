package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

public class WriteDistanceMatrixOneVsAll
{
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("Usage genome name");
			
			HashMap<String, Integer> mainCounts = getCounts(args[0]);
		}
	}
	
	/*
	private double getSumSquare(HashMap<String, Integer> counts)
	{
		
	}*/
	
	private static HashMap<String, Integer> getCounts(String genomeName) throws Exception
	{
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		
		File inFile = new File(MakeKmers.KMER_DIR.getAbsolutePath() + File.separator + 
						genomeName);
		
		if( ! inFile.exists())
			throw new Exception("Could not find " + inFile.getAbsolutePath() );
		
		BufferedReader reader = new BufferedReader(new FileReader(inFile));
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			if( splits.length != 2)
				throw new Exception("no");
			
			if( map.containsKey(splits[0]))
				throw new Exception("No");
			
			map.put(splits[0], Integer.parseInt(splits[1]));
		}
		
		reader.close();
				
		return map;
	}
}
