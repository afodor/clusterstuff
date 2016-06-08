package cmcDistances.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.StringTokenizer;

import utils.Translate;

/*
 * Example SNP file:
 * Num match = 7376725.0
Fraction match = 0.8962499465413206
sequence        context1        context2        distance
440069386814719587      [0,0,6,0]       [0,5,0,0]       7.810249675906654
 */

public class TestASnpFile
{
	private static final char[] NUCLEOTIDES = { 'A', 'C', 'G', 'T' };
	
	private static final int BOTH_CUTOFF = 10;
	
	private static HashMap<String, String> getMostMap(File contextFile) throws Exception
	{
		HashMap<String, String> map = new HashMap<String,String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(contextFile));
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String kmer = new StringTokenizer(s).nextToken();
			
			String val =  map.get(kmer);
			
			if( val != null)
				throw new Exception("Duplicate entry " + kmer);
				
			String reverse = Translate.reverseTranscribe(kmer);
				
			if( ! kmer.equals(reverse) && map.containsKey(reverse))
					throw new Exception("Duplicate entry " + kmer + " " + reverse);
			
			map.put(kmer, getMostLine(s));
		}
		
		reader.readLine();
		
		return map;
	}
	
	private static String getMostLine(String fileLine) throws Exception
	{
		Integer highest = null;
		String collected = null;
		StringTokenizer sToken = new StringTokenizer(fileLine);
		sToken.nextToken();
		
		for( int x=0; x < NUCLEOTIDES.length; x++ )
		{
			Integer count = Integer.parseInt(sToken.nextToken());
			
			if( highest == null || count >= highest)
			{
				if( highest != null && highest.equals(count))
				{
					collected += NUCLEOTIDES[x];
				}
				else if(  highest == null || count > highest )
				{
					collected = "" + NUCLEOTIDES[x];
					highest = count;
				}
				else throw new Exception("Logic error");
			}
		}
		
		if( collected == null)
			throw new Exception("Logic error");
		
		if( sToken.hasMoreTokens())
			throw new Exception("Parsing error");
		
		return collected;
	}
	
	
	
	private static void verifyAFile(File file) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer(file.getName(), "@");
		
		File aFile = new File(WriteContextScripts.OUTPUT_DIRECTORY.getAbsolutePath() + 
				File.separator + sToken.nextToken() + ".context");
		
		if( ! aFile.exists())
			throw new Exception("Could not find " + aFile.getAbsolutePath());
		
		File bFile = new File(WriteContextScripts.OUTPUT_DIRECTORY.getAbsolutePath() + 
				File.separator + sToken.nextToken() + ".context");
		
		if( ! bFile.exists())
			throw new Exception("Could not find " + aFile.getAbsolutePath());
		
		if( sToken.hasMoreTokens())
			throw new Exception("Unexpected file name " +bFile.getAbsolutePath() );
		
		HashMap<String, String> aMap = getMostMap(aFile);
		HashMap<String, String> bMap =getMostMap(bFile);
		
		
		
	}
}
