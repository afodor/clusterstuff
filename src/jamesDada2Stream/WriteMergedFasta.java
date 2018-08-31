package jamesDada2Stream;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.StringTokenizer;

public class WriteMergedFasta
{
	private static String TOP_DIR = "/scratch/afodor_research/datasets/uegp/raw_Stream_16S/forward_filtered";

	private static class Holder implements Comparable<Holder>
	{
		String seq;
		Long count;
		
		@Override
		public int compareTo(Holder other)
		{
			return other.count.compareTo(this.count);
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		HashMap<String, Long>  map = getMap();
		
		System.out.println("Got " + map.size() + " unique sequences");
		
		List<Holder> list = new ArrayList<Holder>();
		
		for(String s : map.keySet())
		{
			Holder h = new Holder();
			h.count = map.get(s);
			h.seq = s;
			list.add(h);
		}
		
		Collections.sort(list);
		
		BufferedWriter writer = new BufferedWriter(new  FileWriter(
				new File("/scratch/afodor_research/datasets/uegp/dada2ToBlast/collectedSeqVariants.txt")));
		
		BufferedWriter writer2 = new BufferedWriter(new  FileWriter(
				new File("/scratch/afodor_research/datasets/uegp/dada2ToBlast/rankAbundance.txt")));
		
		writer2.write("rank\tcount\tsequence\n");
		
		int index = 1;
		for(Holder h : list)
		{
			writer.write(">Seq" + index + " counts=" + h.count + "\n");
			writer.write(h.seq + "\n");
			
			writer2.write(index + "\t" + h.count + "\t" + h.seq + "\n");
		}
		
		writer.flush(); writer.close();
		writer2.flush();  writer2.close();	
	}
	
	public static HashMap<String, Long> getMap() throws Exception
	{
		HashMap<String,Long> map=new LinkedHashMap<String,Long>();
		
		File topDir = new File(TOP_DIR);
		
		for(String nextDir :topDir.list() )
		{
			File bottomDir = new File(topDir.getAbsolutePath() + File.separator + nextDir);
			
			if( bottomDir.isDirectory())
			{
				String[] list2 = bottomDir.list();
				
				for(String fileName : list2)
				{
					if( fileName.endsWith("_F.txt"))
					{
						File aFile = new File(bottomDir.getAbsolutePath() +File.separator + fileName);
						System.out.println(aFile.getAbsolutePath());
						addToMap(map, aFile);
					}
				}
			}
		}
		
		return map;
	}
	
	private static void addToMap(HashMap<String,Long> map, File file) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		String firstLine = reader.readLine();
		String secondLine = reader.readLine();
		StringTokenizer sToken = new StringTokenizer(firstLine);
		
		StringTokenizer sToken2 = new StringTokenizer(secondLine);
		
		sToken2.nextToken();
		
		while(sToken.hasMoreTokens())
		{
			String seq = sToken.nextToken().replaceAll("\"", "");
			long count = Long.parseLong(sToken2.nextToken().replaceAll("\"", ""));
			
			Long oldCount = map.get(seq);
			
			if( oldCount == null)
				oldCount = 0l;
			
			oldCount += count;
			
			map.put(seq, oldCount);
		}
		
		if( sToken2.hasMoreTokens())
			throw new Exception("Parsing error");
		
		reader.close();
	}
}
