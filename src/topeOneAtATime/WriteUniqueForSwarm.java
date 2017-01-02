package topeOneAtATime;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class WriteUniqueForSwarm
{
	public static void main(String[] args) throws Exception
	{
		HashMap<String, Integer> map = new HashMap<String,Integer>();
		addDirectoryToMap("/nobackup/afodor_research/topeOneAtATime/file3/mergedOut", map);
		addDirectoryToMap("/nobackup/afodor_research/topeOneAtATime/file4/mergedOut", map);
		
		List<Holder> list = new ArrayList<Holder>();
		for(String s : map.keySet())
		{
			Holder h = new Holder();
			h.sequence = s;
			h.count = map.get(s);
			list.add(h);
		}
		
		Collections.sort(list);
		
		System.out.println("Writing");
		BufferedWriter writer =new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/topeOneAtATime/mergedForSwarmIncludingSingletons.txt")));
		
		int index=0;
		for(Holder h : list)
			//if( h.count > 1)
			{
				index++;
				writer.write(">seq" + index + "_" + h.count + "\n");
				writer.write(h.sequence + "\n");
			}
			
		writer.flush();  writer.close();
		System.out.println("Finished");
	}
	
	private static class Holder implements Comparable<Holder>
	{
		String sequence;
		int count;
		
		@Override
		public int compareTo(Holder o)
		{
			return o.count - this.count;
		}
	}
	
	private static void addDirectoryToMap(String directoryPath, HashMap<String, Integer> map) 
				throws Exception
	{
		File topDir = new File(directoryPath);
		
		String[] files =topDir.list();
		
		for(String s : files)
		{
			System.out.println(s);
			
			FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(
				topDir.getAbsolutePath() + File.separator + s	);	
			
			for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs = fsoat.getNextSequence())
			{
				Integer count = map.get(fs.getSequence());
				
				if( count == null)
					count =0;
				
				count++;
				map.put(fs.getSequence(), count);
			}
			
			fsoat.close();
		}
		
	}
}
