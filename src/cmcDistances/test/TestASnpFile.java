package cmcDistances.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;

import cmcDistances.Encode;
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
	
	private static int getMax(String s) throws Exception
	{
		int max =0;
		
		s =s.replace("[", "").replace("]", "");
		StringTokenizer sToken = new StringTokenizer(s, ",");
		
		for( int x=0; x < 4; x++ )
		{
			max = Math.max( Integer.parseInt(sToken.nextToken()), max);
		}
		
		if( sToken.hasMoreTokens())
			throw new Exception("No");
		
		return max;
		
	}
	
	private static List<String> getExpectedDifferences(File file)
		throws Exception
	{
		List<String> list = new ArrayList<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		reader.readLine(); 
		reader.readLine(); 
		reader.readLine(); 
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			if( splits.length != 4)
				throw new Exception("Parsing error " + s);
			
			if( getMax(splits[1]) >= BOTH_CUTOFF && getMax(splits[2]) >= BOTH_CUTOFF )
			{
				String key = Encode.getKmer(Long.parseLong(splits[0]), 30);
				list.add(key);
			}
		}
		
		reader.close();
		
		return list;
	}
	
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
		int max = 0;
		
		for( int x=0; x < NUCLEOTIDES.length; x++ )
		{
			Integer count = Integer.parseInt(sToken.nextToken());
			max = Math.max(count, max);
			
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
		
		if( max < BOTH_CUTOFF)
			return "";
		
		return collected;
	}
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("Usage fileToVerify");
			System.exit(1);
		}
		
		verifyAFile(new File(args[0]));
	}
	
	private static void verifyAFile(File file) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer(file.getName(), "@");
		
		File aFile = new File(cmcDistances.WriteContextScripts.OUTPUT_DIRECTORY.getAbsolutePath() + 
				File.separator + sToken.nextToken() + ".context");
		
		if( ! aFile.exists())
			throw new Exception("Could not find " + aFile.getAbsolutePath());
		
		File bFile = new File(cmcDistances.WriteContextScripts.OUTPUT_DIRECTORY.getAbsolutePath() + 
				File.separator + sToken.nextToken().replace(".txt", "") + ".context");
		
		if( ! bFile.exists())
			throw new Exception("Could not find " + bFile.getAbsolutePath());
		
		if( sToken.hasMoreTokens())
			throw new Exception("Unexpected file name " +file.getName());
		
		HashMap<String, String> aMap = getMostMap(aFile);
		HashMap<String, String> bMap =getMostMap(bFile);
		
		List<String> expectedDiffs = getExpectedDifferences(file);
		
		int numMatching = 0;
		int numTooSmall =0;
		
		for(String s : aMap.keySet())
		{
			if( bMap.containsKey(s))
			{
				String aMost = aMap.get(s);
				String bMost = bMap.get(s);
				 
				if( aMost.length() == 0 || bMost.length() ==0 )
				{
					numTooSmall++;
				}
				else if( aMost.equals(bMost))
				{
					numMatching++;
				}
				else
				{
					boolean gotOne = false;
					
					for(String s2 : expectedDiffs)
					{
						if( ! gotOne)
						{
							if( s.equals(s2) || s.equals(Translate.reverseTranscribe(s2)))
							{
								System.out.println("MATCH " + s + " " + s2);
								gotOne = true;
							}
						}
						
						if( ! gotOne)
						{
							System.out.println("Did not match " + s + " " + s2 );
						}
					}
				}
			}
		}
		
		System.out.println("Matched "  + numMatching + " skipped " + numTooSmall );
		
	}
}
