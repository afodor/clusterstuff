package creOrthologs.kmers.chunks;


public class ChunkDistanceMatrix
{
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
