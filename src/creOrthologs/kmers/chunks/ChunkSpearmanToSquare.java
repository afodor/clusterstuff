package creOrthologs.kmers.chunks;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;

import utils.Spearman;
import creOrthologs.kmers.WriteSpearmanFromRandom;
import creOrthologs.kmers.chunks.ChunkSpearman.Holder;

public class ChunkSpearmanToSquare
{
	public static void main(String[] args) throws Exception
	{
		List<Holder> list = ChunkSpearman.getComparisonMatrices();
		System.out.println("Got " + list.size());
		List<String> genomeNames =  WriteSpearmanFromRandom.getGenomeNames();
		
		HashSet<Integer> includeIndex = WriteSpearmanFromRandom.getIncludeIndexes(genomeNames);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/chunkedSpearmanSquare.txt"));
		
		for(int x=0; x < list.size(); x++)
		{
			Holder aHolder = list.get(x);
			
			List<Float> aVals = WriteSpearmanFromRandom.getValsOrNull(aHolder.file,includeIndex);
				
			if( aVals != null)
			{	
				writer.write(aHolder.file.getName());
				for( int y=0; y < list.size(); y++)
				{
					Holder bHolder = list.get(y);
					
					List<Float> bVals = 
							WriteSpearmanFromRandom.getValsOrNull(bHolder.file, includeIndex);
									
					if(bVals != null)
					{
						writer.write("\t" + Spearman.getSpear(aVals, bVals).getRs());
					}
				}
				writer.write("\n"); writer.flush();
			}
		}
		
		writer.flush(); writer.close();
	}
}
	