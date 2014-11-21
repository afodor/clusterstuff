package hmp;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.StringTokenizer;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;
import utils.ConfigReader;

public class PivotOutBySite
{
	private static String getSampleId(FastaSequence fseq) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer(fseq.getHeader());
		
		while(sToken.hasMoreTokens())
		{
			StringTokenizer innerToken = new StringTokenizer(sToken.nextToken(), "=");
			
			if( innerToken.nextToken().equals("sample"))
				return innerToken.nextToken();
		}
		
		throw new Exception("Could not find sample id");
	}
	
	private static String getPrimer(FastaSequence fs) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer(fs.getHeader());
		
		while(sToken.hasMoreTokens())
		{
			StringTokenizer innerToken = new StringTokenizer(sToken.nextToken(), "=");
			
			if( innerToken.nextToken().equals("primer"))
				return innerToken.nextToken();
		}
		
		throw new Exception("Could not find primer");
	}
	
	private static boolean isStool(FastaSequence fseq) throws Exception
	{
		StringTokenizer sToken = new StringTokenizer(fseq.getHeader());
		
		while(sToken.hasMoreTokens())
		{
			StringTokenizer innerToken = new StringTokenizer(sToken.nextToken(), "=");
			
			if( innerToken.nextToken().equals("body_site"))
			{
				if( innerToken.nextToken().equals("Stool"))
					return true;
				else
					return false;
			}
		}
		
		throw new Exception("Could not find body site");
	}
	
	public static void main(String[] args) throws Exception
	{
		FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(ConfigReader.getHMPDir() + 
				File.separator + "SRP002395-7514-cs-nbp-rc.fsa.gz");
		
		File parentDir = new File(ConfigReader.getHMPDir() + File.separator + 
										"stoolBySample" + File.separator  );
		
		HashMap<String, BufferedWriter> fileMap = new HashMap<String, BufferedWriter>();
		
		int numDone =0;
		for( FastaSequence fs = fsoat.getNextSequence(); fs != null; 
							fs = fsoat.getNextSequence())
		{
			if( isStool(fs) )
			{
				String sampleID = getSampleId(fs) + "_" + getPrimer(fs);
				BufferedWriter writer =  fileMap.get(sampleID);
				
				if( writer == null)
				{
					writer = new BufferedWriter(new FileWriter(new File(
							parentDir.getAbsolutePath() + File.separator + sampleID)));
					fileMap.put(sampleID, writer);
					
				}
				
				writer.write(fs.getHeader() + "\n");
				writer.write(fs.getSequence() + "\n");
				writer.flush();
				
				if( ++numDone % 10000 == 0)
					System.out.println(numDone);
				
			}
		}
		
		for(BufferedWriter writer : fileMap.values())
		{
			writer.flush();  writer.close();
		}
		System.out.println("finished");
	}
}
