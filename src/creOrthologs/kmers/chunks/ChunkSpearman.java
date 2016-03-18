package creOrthologs.kmers.chunks;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import creOrthologs.kmers.GatherDistanceMatrix;

public class ChunkSpearman
{
	private static class Holder
	{
		int startPos;
		int endPos;
		boolean isBaseline;
		File file;
	}
	
	private static void addFileOrNull(Holder h) throws Exception
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
	
	private static List<Holder> getComparisonMatrices() throws Exception
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
		}
		
		
		return list;
	}
	
	public static void main(String[] args) throws Exception
	{
		List<Holder> list = getComparisonMatrices();
		System.out.println("Got " + list.size());
	}
	
	/*
	public static void main(String[] args) throws Exception
	{
		List<String> genomeNames = getGenomeNames();
		HashSet<Integer> indexes = getIncludeIndexes(genomeNames);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/randomSpearman.txt"));
		
		writer.write("aFile\tbFile\tdistanceAll\tdistancePneuOnly\n");
		
		String[] list = KMER_DIST_DIRECTORY.list();
		
		for(int x=0; x < list.length-1; x++)
		{
			String xName = list[x];
			
			if( xName.endsWith("dist.txt"))
			{
				List<Float> aVals = getValsOrNull(new File(
						KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName),null);
				
				if( aVals != null)
				{
					List<Float> aValsPneuOnly =
							getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName),indexes);
					
					for( int y=x+1; y < list.length; y++)
					{
						String yName = list[y];
						
						if( yName.endsWith("dist.txt"))
						{
							List<Float> bVals = getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),null);
									
							if(bVals != null)
							{
								List<Float> bValsPneuOnly = 
										getValsOrNull(new File(
												KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),indexes);
								
								writer.write(xName + "\t");
								writer.write(yName + "\t");
								writer.write(Spearman.getSpear(aVals, bVals).getRs() + "\t");
								writer.write(Spearman.getSpear(aValsPneuOnly,bValsPneuOnly).getRs() + "\n");
								writer.flush();
							}
						}
					}
				}
			}
		}
		writer.flush(); writer.close();
	}	*/
}
