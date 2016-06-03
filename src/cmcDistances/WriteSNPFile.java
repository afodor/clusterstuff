package cmcDistances;

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

import coPhylog.ContextCount;

public class WriteSNPFile
{
	private static int MIN_NUMBER_READS = 5;
	
	public static HashMap<Long, ContextCount> readFileRequireMin(File file, int minRequiredReads) 
			throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();

		for(String s = reader.readLine(); s != null; s= reader.readLine())
		{
			String[] splits = s.split("\t");
			
			if( splits.length != 5)
				throw new Exception("Unexpectd input");
			
			Long longKey = Encode.makeLong(splits[0]);
			
			
			if( map.containsKey(longKey))
				throw new Exception("Duplicate");
			

			ContextCount cc = new ContextCount(Byte.parseByte(splits[1]), 
					Byte.parseByte(splits[2]), Byte.parseByte(splits[3]), Byte.parseByte(splits[4]));
			
			if( cc.getMax() >= minRequiredReads)
				map.put(longKey,cc);
		}
		
		reader.close(); 
		return map;	
	}

	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 3)
		{
			System.out.println("Usage context1 context2 outFile");
			System.exit(1);
		}
		
		File outFile = new File(args[2]);
		
		if( outFile.exists())
			throw new Exception(outFile.getAbsolutePath() + " already exists ");
		
		HashMap<Long, ContextCount> map1 = 
				readFileRequireMin(new File(args[0]), MIN_NUMBER_READS);
		
		HashMap<Long, ContextCount> map2 = 
				readFileRequireMin(new File(args[1]), MIN_NUMBER_READS);
		
		System.out.println("Comparing " + map1.size() + " " + map2.size());
		
		List<Holder> snpList = new ArrayList<WriteSNPFile.Holder>();
		
		double numMatch =0;
		
		for( Long l : map1.keySet() )
		{
			if( map2.containsKey(l))
			{
				ContextCount cc1 = map1.get(l);
				ContextCount cc2 = map2.get(l);
				
				if( cc1.isDifferentInHighest(cc2))
				{
					HashSet<Character>  high1 = cc1.getHighest();
					
					boolean allZeros = true;
					
					for(Character c : high1)
					{
						if( c == 'A' && cc2.getNumA() > 0)
							allZeros = false;
						
						if( c == 'C' && cc2.getNumC() > 0)
							allZeros = false;

						if( c == 'G' && cc2.getNumG() > 0)
							allZeros = false;
						
						if( c == 'T' && cc2.getNumT() > 0)
							allZeros = false;
						
					}
					
					HashSet<Character> high2 = cc2.getHighest();
					
					for(Character c : high2)
					{
						if( c == 'A' && cc1.getNumA() > 0)
							allZeros = false;
						
						if( c == 'C' && cc1.getNumC() > 0)
							allZeros = false;

						if( c == 'G' && cc1.getNumG() > 0)
							allZeros = false;
						
						if( c == 'T' && cc1.getNumT() > 0)
							allZeros = false;
						
					}
					
					if ( allZeros)
					{
						Holder h = new Holder();
						h.id = l;
						h.cc1 = cc1;
						h.cc2 = cc2;
						h.distance = cc1.getRawDistance(cc2);
						snpList.add(h);
		
					}
				}
				else
				{
					numMatch++;
				}
			}
			
		}
			
		Collections.sort(snpList);
		System.out.println("Found " + snpList.size() + " out of " + map1.size() + " " + map2.size() + " " + numMatch + " " + 
		"fractionMatch = " + (numMatch/map1.size()));
		writeResults(outFile, snpList);
	}
	
	private static void writeResults(File outFile, List<Holder> snpList) throws Exception
	{
		if( outFile.exists())
			throw new Exception(outFile.getAbsolutePath() + " already exists ");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		writer.write("sequence\tcontext1\tcontext2\tdistance\n");
		
		for( Holder h : snpList)
		{
			writer.write( h.id + "\t");
			writer.write(h.cc1 + "\t");
			writer.write(h.cc2 + "\t");
			writer.write( h.distance + "\n");
		}
		
		writer.flush();  writer.close();
	}
	
	private static class Holder implements Comparable<Holder>
	{
		long id;
		ContextCount cc1;
		ContextCount cc2;
		double distance;
		
		@Override
		public int compareTo(Holder o)
		{
			return Double.compare(o.distance, this.distance);
		}
	}
}
