package kw_creOrthologs;

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

import creOrthologs.kmers.MakeKmers;
import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;
import utils.Translate;

public class ConstrainKMersToOrtholog
{
	private static final File GATHERED_DIR =
				new File("/nobackup/afodor_research/af_broad/orthologs/gatheredKmerDistanceMatrices");
	
	private static double getSumSquared( HashMap<String, HashMap<String,Integer>> bigMap, String s )
	{
		long count =0;
		
		HashMap<String, Integer> innerMap = bigMap.get(s);
		
		for( Integer i : innerMap.values())
			count = count + i * i;
		
		return count;
	}
	
	private static long getSum( HashMap<String, HashMap<String,Integer>> bigMap, String s )
	{
		long count =0;
		
		HashMap<String, Integer> innerMap = bigMap.get(s);
		
		for( Integer i : innerMap.values())
			count = count + i;
		
		return count;
	}
	
	private static double getDistance(HashMap<String, HashMap<String,Integer>> bigMap,
					String s1, String s2) throws Exception
	{
		double sumASquared = getSumSquared(bigMap, s1);
		double sumBSquared = getSumSquared(bigMap, s2);
		
		if( sumASquared == 0 || sumBSquared == 0 )
			return 1;
		
		long topSum = 0;
		
		HashMap<String, Integer> aMap = bigMap.get(s1);
		HashMap<String, Integer> bMap = bigMap.get(s2);
		
		for( String s : aMap.keySet() )
		{
			if( bMap.containsKey(s))
			{
				topSum += aMap.get(s) * bMap.get(s);
			}	
		}

		return 1- topSum / Math.sqrt(sumASquared * sumBSquared);
	}
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("usage orthologFilePath");
			System.exit(1);
		}
		
		HashMap<String, HashMap<String,Integer>> bigMap = getBigMap( args[0]);
		
		String outFileBase = args[0];
		
		writeDistanceMatrix(bigMap, outFileBase);
	}
	
	private static void writeDistanceMatrix( HashMap<String, HashMap<String,Integer>> bigMap,
			String outFileBase)
		throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			GATHERED_DIR.getAbsolutePath() + File.separator + outFileBase + "_dist.txt"	)));
		
		BufferedWriter keyWriter = new BufferedWriter(new FileWriter(new File(
				GATHERED_DIR.getAbsolutePath() + File.separator + outFileBase + "_Key.txt"	)));
		keyWriter.write("shortName\tlongName\tnumberOfKmers\n");
		
		List<String> list = new ArrayList<String>(bigMap.keySet());
		Collections.sort(list);
	
		writer.write(list.size() + "\n");
		
		for( int x=0; x < list.size(); x++)
		{
			String s1 = list.get(x);
			String name = "BUG_" + x;
			
			while( name.length() < 10)
				name = name + "_";
			
			writer.write(name);
			keyWriter.write(name + " " + list.get(x) + "\t" + getSum(bigMap, list.get(x)) + "\n");
			
			for( int y=0; y < list.size(); y++)
			{
				writer.write(" " + getDistance(bigMap, s1, list.get(y)));
				//System.out.println(list.get(x) + " "+ list.get(y) + " " + val);
			}
			
			writer.write("\n");
		}
		
		writer.flush();  writer.close();
		keyWriter.flush();  keyWriter.close();
	}
	
	// outer key is the genome may; inner key is the k-mer
	private static HashMap<String, HashMap<String,Integer>> getBigMap(String filepath) throws Exception
	{
		HashMap<String, HashMap<String,Integer>> map = new HashMap<String, HashMap<String,Integer>>();
		HashSet<String> constrainingSet = getConstrainingSet(new File(filepath));
		System.out.println("Got constraining set with " + constrainingSet.size());
		
		for(String s : MakeKmers.KMER_DIR.list() )
			if( s.endsWith("_kmers.txt"))
			{
				String shortName = s.replace("_kmers.txt", "");
				HashMap<String, Integer> innerMap = new HashMap<String,Integer>();
				map.put(shortName,innerMap);
				
				BufferedReader reader = new BufferedReader(new FileReader(
						new File(MakeKmers.KMER_DIR.getAbsolutePath() + File.separator + s)));
				
				for(String s2= reader.readLine(); s2 != null; s2= reader.readLine())
				{
					String[] splits = s2.split("\t");
					
					if( splits.length !=2 )
						throw new Exception("No");
					
					String key = null;
					
					if( constrainingSet.contains(splits[0]) )
					{
						key = splits[0];
					}
					else
					{
						String reverse = Translate.reverseTranscribe(splits[0]);
						
						if( constrainingSet.contains(reverse))
						{
							key = reverse;
						}
					}
					
					if( key != null)
					{
						if( innerMap.containsKey(key))
							throw new Exception("No");
						
						innerMap.put(key, Integer.parseInt(splits[1]));
					}
				}
				
				reader.close();
			}
		
		return map;
	}
	
	private static HashSet<String> getConstrainingSet(File ortholgFile) throws Exception
	{	
		HashMap<String, Integer> map = new HashMap<String,Integer>();
		
		
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(ortholgFile);
		
		for( FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
		{
			String seq = fs.getSequence();
			
			for( int x=0; x < seq.length()- MakeKmers.KMER_LENGTH; x++)
			{
				String sub = seq.substring(x, x +  MakeKmers.KMER_LENGTH);
				
				if( MakeKmers.isACGT(sub))
				{
					MakeKmers.addToMap(map, seq);
				}
			}
		}
		
		MakeKmers.checkForNoReverse(map);
		
		return new HashSet<String>(map.keySet());
	}
}
