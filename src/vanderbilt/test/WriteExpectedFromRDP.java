package vanderbilt.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.StringTokenizer;

public class WriteExpectedFromRDP
{
	public static final double THRESHOLD =0.49999;
	public static final String TAXA_LEVEL = "family";
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 2)
		{
			System.out.println("Usage rdpFile summaryFile");
			System.exit(1);
		}
		
		HashMap<String, Integer> map =getCountMap(args[0]);
		
		File outFile = new File(args[1]);
		
		if( outFile.exists())
			throw new Exception("out file already exists");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		
		for(String s : map.keySet())
			writer.write(s + "\t" + map.get(s) + "\n");
		
		writer.flush();  writer.close();
	}
	
	private static HashMap<String, Integer> getCountMap(String filePath) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(new File(filePath)));
		
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		
		for( String s= reader.readLine(); s != null; s = reader.readLine())
		{
			String taxa = getAtThreshold(s, TAXA_LEVEL, THRESHOLD);
			
			if( taxa != null)
			{
				Integer val = map.get(taxa);
				
				if( val == null)
					val = 0;
				
				val++;
				
				map.put(taxa, val);
			}
		}
		
		return map;
	}
	
	private static String getAtThreshold(String s, String level, double threshold) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer(s, "\t");
		
		System.out.println("First " +  sToken.nextToken());
		
		String taxaNameToken = null;
		String taxaLevelToken = null;
		Double score = null;
		
		while(sToken.hasMoreTokens())
		{
			taxaNameToken = sToken.nextToken();
			taxaLevelToken = sToken.nextToken();
			score = Double.parseDouble(sToken.nextToken());
			
			if( taxaLevelToken.equals(level) &&  score >= THRESHOLD)
				return taxaNameToken;
		}
		
		return null;
	}
}
