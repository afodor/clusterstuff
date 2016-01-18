package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

import utils.Spearman;

public class SpearmansAcrossAll
{
	public static final int EXPECTED_NUM_LINES = 340;
	
	public static void main(String[] args) throws Exception
	{
		List<File> list = getFilesToDo();
		
		System.out.println("Got " + list.size() + " files ");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			"/nobackup/afodor_research/af_broad/symmetric1.txt"	)));
		
		for(int x=0; x < list.size(); x++)
		{
			writer.write(getIdentifier(list.get(x)));
			
			List<Double> xVals = getVals(list.get(x));
			
			for( int y=0; y < list.size(); y++)
			{
				List<Double> yVals = getVals(list.get(y));
				
				if( xVals.size() != yVals.size())
					throw new Exception("No " + xVals.size() + " " + yVals.size());
				
				writer.write("\t" + Spearman.getSpearFromDouble(yVals, xVals).getRs());
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
	
	private static List<Double> getVals(File f) throws Exception
	{
		List<Double> vals = new ArrayList<Double>();
		
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		for(String s=  reader.readLine(); s != null; s= reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			
			while(sToken.hasMoreTokens())
				vals.add(Double.parseDouble(sToken.nextToken()));
		}
		
		reader.close();
		
		return vals;
	}
	
	private static List<File> getFilesToDo() throws Exception
	{
		HashSet<Integer> includeSet = getIncludeSet();
		String[] files = GatherDistanceMatrix.GATHERED_DIR.list();
		
		List<File> list = new ArrayList<File>();
		
		for(String s : files)
			if( s.startsWith("klebsiella_pneumoniae_chs_11.0_") && s.endsWith("_dist.txt"))
			{
				File f= new File(GatherDistanceMatrix.GATHERED_DIR.getAbsolutePath() + File.separator +s);
				
				if( getFileLines(f) == EXPECTED_NUM_LINES ) 
				{
					if( includeSet.contains(Integer.parseInt(f.getName().split("_")[5])))
						list.add(f);
				}
					
			}
		
		Collections.sort(list);
		return list;
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
