package creOrthologs.kmers.chunks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.StringTokenizer;

import utils.Spearman;
import creOrthologs.kmers.GatherDistanceMatrix;
import creOrthologs.kmers.WriteSpearmanFromRandom;

public class ChunkSpearman
{
	public static class Holder
	{
		int startPos;
		int endPos;
		boolean isBaseline;
		File file;
	}
	
	static void addFileOrNull(Holder h) throws Exception
	{
		String genomePath = WriteScriptsForChunkDistanceMatrices.GENOME_PATH;
		String outFileBase =  genomePath.substring(genomePath.lastIndexOf("/")+1)
				.replace(".scaffolds.fasta", "") + "_" + WriteScriptsForChunkDistanceMatrices.CONTIG 
				+ "_" + h.startPos + "_" + h.endPos;

		File file = new File(GatherDistanceMatrix.GATHERED_DIR.getAbsolutePath() 
				+ File.separator + outFileBase + "_dist.txt");
		
		if( ! file.exists())
		{
			System.out.println("Could not find " + file.getAbsolutePath());
			h.file = null;
			return;
		}
		
		BufferedReader reader =new BufferedReader(new FileReader(file));
		
		int expected =340;
		int found =1;
		
		reader.readLine();
		
		for(String s = reader.readLine(); s != null; s= reader.readLine())
		{
			int numTokens =0;
			
			StringTokenizer sToken = new StringTokenizer(s);
			
			while( sToken.hasMoreTokens())
			{
				sToken.nextToken();
				numTokens++;
			}
			
			if( numTokens != expected)
			{
				System.out.println("Wrong number of tokens");
				h.file = null;
				reader.close();
				return;
			}
			
			found++;
		}
		
		reader.close();
		
		if( found != expected)
		{
			System.out.println("Wrong number of lines");
			h.file = null;
			return;
		}
		
		h.file = file;

	}
	
	static List<Holder> getComparisonMatrices() throws Exception
	{
		List<Holder> list = new ArrayList<Holder>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				"/nobackup/afodor_research/af_broad/chunks/pneuOnlyChunks_0.85_0.9.txt")));
		
		reader.readLine();
		
		String[] lastSplits= null;
		for(String s = reader.readLine() ; s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");

			if( splits.length != 4)
				throw new Exception("No");
			
			if( lastSplits != null )
			{
				Holder h = new Holder();
				h.startPos = Integer.parseInt(lastSplits[1]) + 6000;
				h.endPos = Integer.parseInt(splits[0]) -1000;
				h.isBaseline = true;
				
				
				if( h.endPos - h.startPos > 5000)
				{
					addFileOrNull(h);
					
					if( h.file != null)
						list.add( h);
				}
			}
			
			Holder h= new Holder();
			h.startPos = Integer.parseInt(splits[0]);
			h.endPos = Integer.parseInt(splits[1])+5000;
			h.isBaseline = false;
			addFileOrNull(h);
			
			if( h.file != null)
				list.add(h);
			
			lastSplits = splits;
		}
		
		
		return list;
	}
	
	public static void main(String[] args) throws Exception
	{
		List<Holder> list = getComparisonMatrices();
		System.out.println("Got " + list.size());
		List<String> genomeNames =  WriteSpearmanFromRandom.getGenomeNames();
		
		HashSet<Integer> includeIndex = WriteSpearmanFromRandom.getIncludeIndexes(genomeNames);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/chunkedSpearman.txt"));
		
		writer.write("aFile\tbFile\taStart\taEnd\tbStart\tbEnd\t" + 
							"aIsBaseline\tbIsBaseline\tbothBaseline\tneitherBaseline\tcrossed\t" + 
				"distancePneuOnly\n");
		
		for(int x=0; x < list.size()-1; x++)
		{
			Holder aHolder = list.get(x);
			
			List<Float> aVals = WriteSpearmanFromRandom.getValsOrNull(aHolder.file,includeIndex);
				
			if( aVals != null)
			{		
				for( int y=x+1; y < list.size(); y++)
				{
					Holder bHolder = list.get(y);
					
					List<Float> bVals = 
							WriteSpearmanFromRandom.getValsOrNull(bHolder.file, includeIndex);
									
					if(bVals != null)
					{
						writer.write(aHolder.file.getName()+ "\t");
						writer.write(bHolder.file.getName()+ "\t");
						writer.write(aHolder.startPos + "\t");
						writer.write(aHolder.endPos + "\t");
						writer.write(bHolder.startPos + "\t");
						writer.write(bHolder.endPos + "\t");
						writer.write(aHolder.isBaseline + "\t");
						writer.write(bHolder.isBaseline + "\t");
						writer.write( (aHolder.isBaseline && bHolder.isBaseline )+ "\t" );
						boolean neither = aHolder.isBaseline == false && bHolder.isBaseline == false;
						writer.write(  neither + "\t" );
						
						writer.write( (aHolder.isBaseline != bHolder.isBaseline) + "\t" );

						writer.write(Spearman.getSpear(aVals, bVals).getRs() + "\n");
						writer.flush();
						}
					}
			}
		}
		
		writer.flush(); writer.close();
	}
}
