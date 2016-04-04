package creOrthologs.kmers.chunks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import utils.Spearman;
import creOrthologs.kmers.GatherDistanceMatrix;
import creOrthologs.kmers.WriteSpearmanFromRandom;

public class ChunkSpearmanToSquareAllContigs
{
	private static List<File> getFiles() throws Exception
	{
		List<File> list = new ArrayList<File>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				WriteScriptsForChunkDistanceMatricesWithAllContigs.LOG_PATH	)));
			
		reader.readLine();
		
		for(String s= reader.readLine(); s != null; s= reader.readLine())
		{
			String[] splits = s.split("\t");
			
			File f= new File(
					GatherDistanceMatrix.GATHERED_DIR.getAbsolutePath() 
					+ File.separator + splits[0]+ "_dist.txt");
			
			if( WriteSpearmanFromRandom.getValsOrNull(f, null) != null)
				list.add(f );
		}
			
		return list;
	}
	
	public static void main(String[] args) throws Exception
	{
		List<String> genomeNames = WriteSpearmanFromRandom.getGenomeNames();
		HashSet<Integer> indexes = WriteSpearmanFromRandom.getIncludeIndexes(genomeNames);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/chunkedSpearmanSquareAllContigs.txt"));
		
		List<File> list= getFiles();
		
		for(int x=0; x < list.size(); x++)
		{
			List<Float> aVals = WriteSpearmanFromRandom.getValsOrNull(list.get(x), indexes);
				
			writer.write(list.get(x).getName());
			for( int y=0; y < list.size(); y++)
			{
				List<Float> bVals = 
							WriteSpearmanFromRandom.getValsOrNull(list.get(y), indexes);
									
				writer.write("\t" + Spearman.getSpear(aVals, bVals).getRs());
			}
			writer.write("\n"); writer.flush();
		}
		
		writer.flush(); writer.close();
	}
}
	