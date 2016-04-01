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
	
	private static final File ALL_TREE_MATRIX
		= new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allDist.txt");
	
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
		
		writer.write("distFile\tdistanceAll\tdistancePneuOnly\n");
		
		String[] list = KMER_DIST_DIRECTORY.list();
		
		List<Float> aVals = WriteSpearmanFromRandom.getValsOrNull(ALL_TREE_MATRIX,null);
		List<Float> aValsPneuOnly =WriteSpearmanFromRandom.getValsOrNull(ALL_TREE_MATRIX,indexes);
					
		for( int x=0; x < list.length; x++)
		{
			String yName = list[x];
						
			if( yName.endsWith("dist.txt"))
			{
				List<Float> bVals = WriteSpearmanFromRandom.getValsOrNull(new File(
					KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),null);
									
				if(bVals != null)
				{
					List<Float> bValsPneuOnly = 
							WriteSpearmanFromRandom.getValsOrNull(new File(
							KMER_DIST_DIRECTORY.getAbsolutePath() + File.separator + yName),indexes);
								
					writer.write(yName + "\t");
					writer.write(Spearman.getSpear(aVals, bVals).getRs() + "\t");
					writer.write(Spearman.getSpear(aValsPneuOnly,bValsPneuOnly).getRs() + "\n");
					writer.flush();
				}
			}
		}
		writer.flush(); writer.close();
	}	
}
