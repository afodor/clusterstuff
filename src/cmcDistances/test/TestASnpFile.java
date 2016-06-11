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
	
	private static class Holder
	{
		int numA=0;
		int numC=0;
		int numG=0;
		int numT=0;
		
		private void flip()
		{
			int temp = numA;
			this.numA = this.numT;
			this.numT = temp;
		
			temp = numC;
			this.numC = this.numG;
			this.numG = temp;
		}
		
		private boolean charIsZero(char c) throws Exception
		{
			if( c == 'A' )
			{
				return numA ==0;
			}
			else if ( c =='C')
			{
				return numC ==0;
			} else if ( c== 'G')
			{
				return numG ==0;
			} else if ( c == 'T')
			{
				return numT ==0;
			}
			else throw new Exception("Unknown character");
				
		}
		
		private char mostChar() throws Exception
		{
			Character c = null;
			int max=-1;
			
			if( numA > max)
			{
				max = numA;
				c = 'A';
			}
			
			if( numC > max)
			{
				max = numC;
				c = 'C';
			}
			
			if( numG > max )
			{
				max = numG;
				c = 'G';
			}
			
			if( numT > max )
			{
				max = numT;
				c = 'T';
			}
			
			if( c== null)
				throw new Exception("Logic error");
			
			return c;
		}
	}
	
	private static HashMap<String, Holder> getMostMap(File contextFile) throws Exception
	{
		HashMap<String, Holder> map = new HashMap<String,Holder>();
		
		BufferedReader reader = new BufferedReader(new FileReader(contextFile));
		reader.readLine();
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String kmer = new StringTokenizer(s).nextToken();
			
			Holder h =  map.get(kmer);
			
			if( h != null)
				throw new Exception("Duplicate entry " + kmer);
				
			String reverse = Translate.reverseTranscribe(kmer);
				
			if( ! kmer.equals(reverse) && map.containsKey(reverse))
					throw new Exception("Duplicate entry " + kmer + " " + reverse);
			
			h = getHolder(s);
			
			if( h != null)
				map.put(kmer, h);
		}
		
		reader.close();
		return map;
	}
	
	private static Holder getHolder(String fileLine) throws Exception
	{
		Holder h= new Holder();
		
		StringTokenizer sToken = new StringTokenizer(fileLine);
		sToken.nextToken();
		int max = 0;
		
		h.numA = Integer.parseInt(sToken.nextToken());
		max = Math.max(h.numA, max);
		
		h.numC = Integer.parseInt(sToken.nextToken());
		max = Math.max(h.numC, max);
		
		h.numG = Integer.parseInt(sToken.nextToken());
		max = Math.max(h.numG, max);
		
		h.numT = Integer.parseInt(sToken.nextToken());
		max = Math.max(h.numT, max);
		
		if( sToken.hasMoreTokens())
			throw new Exception("Parsing error");
		
		if( max < BOTH_CUTOFF)
			return null;
		
		return h;
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
		
		HashMap<String, Holder> aMap = getMostMap(aFile);
		HashMap<String, Holder> bMap =getMostMap(bFile);
		
		List<String> expectedDiffs = getExpectedDifferences(file);
		
		int numMatching = 0;
		int numNotMatching =0;
		int numOffSnp =0;
		
		for(String s : aMap.keySet())
		{
			Holder bHolder = bMap.get(s);
			
			if( bHolder != null)
			{
				String flip = Translate.reverseTranscribe(s);
				bHolder = bMap.get(flip);
				
				if( bHolder != null)
					bHolder.flip();
			}
			
			if( bHolder != null)
			{
				numMatching++;
				Holder aHolder = aMap.get(s);
				
				char aChar = aHolder.mostChar();
				char bChar = bHolder.mostChar();
				
				if( aChar==bChar)
				{
					if( aHolder.charIsZero(bChar) && bHolder.charIsZero(aChar) )
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
						}
						
						if( gotOne)
							System.out.println("Matched " + s);
						else
							System.out.println("Failed to match " + s);
					}
					else
					{
						numOffSnp++;
					}
				}
				
			}
			else
			{
				numNotMatching++;
			}
		}
		
		System.out.println("Matched "  + numMatching + " skipped " + numNotMatching + " off SNP " + 
						numOffSnp);
		
	}
}
