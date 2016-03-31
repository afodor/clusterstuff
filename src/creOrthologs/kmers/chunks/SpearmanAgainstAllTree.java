package creOrthologs.kmers.chunks;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;

import creOrthologs.kmers.WriteSpearmanFromRandom;
import utils.Spearman;

public class SpearmanAgainstAllTree
{
	private static final 
		File KMER_DIST_DIRECTORY = new File( "/nobackup/afodor_research/af_broad/gatheredKmerMatrices");
	
	public static List<String> getGenomeNames() throws Exception
	{
		String[] list = KMER_DIST_DIRECTORY.list();
		
		List<String> returnVals = null;
		
		for( String s : list)
		{
			if( s.endsWith("key.txt"))
			{
				File f = new File(KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + 
											s);
				
				if(returnVals == null)
				{
					returnVals = WriteSpearmanFromRandom.getNames(f);
				}
				else
				{
					if( ! returnVals.equals(WriteSpearmanFromRandom.getNames(f)))
						throw new Exception("No " + f.getAbsolutePath());
				}
			}
		}
		
		return returnVals;
	}
	
	public static void main(String[] args) throws Exception
	{
		List<String> genomeNames = WriteSpearmanFromRandom.getGenomeNames();
		HashSet<Integer> indexes = WriteSpearmanFromRandom.getIncludeIndexes(genomeNames);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/allContigsChs11PneuPlusNoPneu.txt"));
		
		writer.write("aFile\tbFile\tdistanceAll\tdistancePneuOnly\n");
		
		String[] list = KMER_DIST_DIRECTORY.list();
		
		for(int x=0; x < list.length-1; x++)
		{
			String xName = list[x];
			
			if( xName.endsWith("dist.txt") )
			{
				List<Float> aVals = WriteSpearmanFromRandom.getValsOrNull(new File(
						KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName),null);
				
				if( aVals != null)
				{
					List<Float> aValsPneuOnly =
							WriteSpearmanFromRandom.getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + xName),indexes);
					
					for( int y=x+1; y < list.length; y++)
					{
						String yName = list[y];
						
						if( yName.endsWith("dist.txt"))
						{
							List<Float> bVals = WriteSpearmanFromRandom.getValsOrNull(new File(
									KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),null);
									
							if(bVals != null)
							{
								List<Float> bValsPneuOnly = 
										WriteSpearmanFromRandom.getValsOrNull(new File(
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
	}	
}
