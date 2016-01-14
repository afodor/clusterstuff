package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

import utils.Translate;

public class WriteDistanceMatrixOneVsAll
{
	public static final File KMER_DISTANCE_MATRIX_DIR = 	
			new File("/nobackup/afodor_research/af_broad/kmerDistanceMatrices");
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("Usage genome name");
			System.exit(1);
		}
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				KMER_DISTANCE_MATRIX_DIR + File.separator 
					+ args[0].replaceAll("_kmers.txt", "") +"_toAll.txt"	)));
			
		HashMap<String, Integer> aCounts = getCounts(args[0]);
		long aSumSquared = getSumSquare(aCounts);
			
		for(String s : MakeKmers.KMER_DIR.list())
		{
			if( s.endsWith("_kmers.txt"))
			{
				HashMap<String, Integer> bMap = getCounts(s);
				double distance = getDistance(aCounts, bMap, aSumSquared);
				writer.write(args[0].replaceAll("_kmers.txt", "") + "\t");
				writer.write(s.replaceAll("_kmers.txt", "") + "\t");
				writer.write(distance + "\n");			
				writer.flush();
			}
		}
			
		writer.flush();  writer.close();
	}
	
	private static long getSumSquare(HashMap<String, Integer> counts)
	{
		long sum =0;
		
		for( Integer i : counts.values())
			sum = sum + i * i;
		
		return sum;
	}
	
	private static double getDistance(HashMap<String, Integer> aMap, HashMap<String, Integer> bMap,
						long sumASquared) throws Exception
	{
		long sumBSquared = getSumSquare(bMap);
		
		long topSum = 0;
		
		for( String s : aMap.keySet() )
		{
			if( bMap.containsKey(s))
			{
				topSum += aMap.get(s) * bMap.get(s);
			}
			else
			{
				String reverse = Translate.reverseTranscribe(s);
				
				if( bMap.containsKey(reverse))
				{
					topSum += aMap.get(s) * bMap.get(reverse);
				}
			}
				
		}
			
		return 1- topSum / Math.sqrt(sumASquared * sumBSquared);
	}
	
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
