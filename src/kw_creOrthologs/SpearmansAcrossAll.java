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
import java.util.StringTokenizer;
import creOrthologs.kmers.*;

import utils.Spearman;

public class SpearmansAcrossAll
{
	private static final File GATHERED_DIR =
			new File("/nobackup/afodor_research/af_broad/orthologs/gatheredKmerDistanceMatrices");

	public static final int EXPECTED_NUM_LINES = 340;
	
	public static void main(String[] args) throws Exception
	{
		HashMap<File, List<Float>> map = getFilesToDo();
		
		System.out.println("Got " + map.size() + " files ");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			"/nobackup/afodor_research/af_broad/orthologs/symmetric1.txt"	)));
		
		List<File> list = new ArrayList<File>(map.keySet());
		Collections.sort(list);

		HashMap<String, Float> cache = new HashMap<String,Float>();
		
		for(int x=0; x < list.size(); x++)
		{
			writer.write(getIdentifier(list.get(x)));
			
			List<Float> xVals = map.get(list.get(x));
			
			for( int y=0; y < list.size(); y++)
			{
				if( x== y)
				{
					writer.write("\t" + 1);
				}
				else if ( x < y )
				{
					List<Float> yVals = getVals(list.get(y));
					
					if( xVals.size() != yVals.size())
						throw new Exception("No " + xVals.size() + " " + yVals.size());
					
					float val = (float) Spearman.getSpear(yVals, xVals).getRs();
					writer.write("\t" + val);
					
					String key = y + "@" + x;
					
					if( cache.containsKey(key))
						throw new Exception("No " + key);
					
					cache.put(key, val);
				}
				else
				{
					String key = x + "@" + y;
					Float val = cache.get(key);
					
					if( val == null)
						throw new Exception("No " + key);
					
					writer.write("\t" + val);
					
					map.remove(key);
				}
				
			}
			
			writer.write("\n");
			writer.flush();
		}
		
		writer.flush();  writer.close();
	}
	
	private static int getFileLines(File f ) throws Exception
	{
		int count =0;
		
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		for(String s=  reader.readLine(); s != null; s= reader.readLine())
			count++;
		
		reader.close();
		
		return count;
	}
	
	private static List<Float> getVals(File f) throws Exception
	{
		List<Float> vals = new ArrayList<Float>();
		
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		for(String s=  reader.readLine(); s != null; s= reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);

			sToken.nextToken();
			
			while(sToken.hasMoreTokens())
				vals.add(Float.parseFloat(sToken.nextToken()));
		}
		
		reader.close();
		
		return vals;
	}
	
	private static HashMap<File, List<Float>> getFilesToDo() throws Exception
	{
		HashMap<File, List<Float>> map = new HashMap<File, List<Float>>();
		HashSet<Integer> includeSet = getIncludeSet();
		String[] files = GATHERED_DIR.list();
		
		
		for(String s : files)
			if( s.startsWith("orthogroups_carolina_klebsiella_pneumoniae_chs_11.0") && s.endsWith("_dist.txt"))
			{
				File f= new File(GATHERED_DIR.getAbsolutePath() + File.separator +s);
				
				if( getFileLines(f) == EXPECTED_NUM_LINES ) 
				{
					if( includeSet.contains(Integer.parseInt(f.getName().split("_")[5])))
						map.put(f, getVals(f));
				}
					
			}
		
		return map;
	}
	
	private static HashSet<Integer> getIncludeSet() throws Exception
	{
		HashSet<Integer> set = new HashSet<Integer>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				"/nobackup/afodor_research/af_broad/initialConstrainedMap.txt")));
		
		reader.readLine();
		
		for(String s= reader.readLine(); s != null; s= reader.readLine())
		{
			String[] splits = s.split("\t");
			
			if( splits.length !=3 )
				throw new Exception("No " + s);
			
			if( Double.parseDouble(splits[2]) < 0.98)
				set.add(Integer.parseInt(splits[0]));
		}
		
		reader.close();
		
		return set;
	}
	
	private static String getIdentifier(File f) throws Exception
	{
		String[] splits = f.getName().split("_");
		
		if( splits.length != 8)
			throw new Exception( f.getName() + " " + splits.length );
		
		return splits[4] + "_" + splits[5] + "_" + splits[6];
	}
}
