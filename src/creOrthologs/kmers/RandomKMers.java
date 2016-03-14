package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

import utils.Translate;

public class RandomKMers
{
	public static void main(String[] args) throws Exception
	{
		if( args.length != 4)
		{
			System.out.println("Usage numKMers inKmerPath distFileOut keyFileOut");
			System.exit(1);
		}
		
		List<String> kmerList = getKmerList(args[1]);
		
		System.out.println("Got " + kmerList.size());
		
		Collections.shuffle(kmerList);
		
		System.out.println("Shuffled");
		
		HashSet<String> set = new HashSet<String>();
		
		int index=0;
		
		int numKmers = Integer.parseInt(args[0]);
		
		// will throw if there an insufficient number of kmers
		while( set.size() < numKmers )
		{
			set.add(kmerList.get(index));
			index++;
		}
		
		HashMap<String, HashMap<String,Integer>> bigMap = getBigMap(set);
		ConstrainKMersToRegion.writeDistanceMatrix(bigMap, new File( args[2]), new File(args[3]));
	}
	
	// outer key is the genome name; inner key is the k-mer
	private static HashMap<String, HashMap<String,Integer>> getBigMap(HashSet<String> constrainingSet) throws Exception
	{
		HashMap<String, HashMap<String,Integer>> map = new HashMap<String, HashMap<String,Integer>>();
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
	
	private static List<String> getKmerList(String inKmerFilePath) throws Exception
	{
		List<String> list = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(
				inKmerFilePath)));

		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			
			String seq = sToken.nextToken();
			int aNum = Integer.parseInt(sToken.nextToken());
			
			if( sToken.hasMoreTokens())
				throw new Exception("No");
			
			for(int x=0; x < aNum; x++)
				list.add(seq);
		}
		
		reader.close();
		
		return list;
	}
}
