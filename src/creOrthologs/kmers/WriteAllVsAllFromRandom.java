package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;

public class WriteAllVsAllFromRandom
{
	private static final File KMER_DIST_DIRECTORY = 
			new File("/nobackup/afodor_research/af_broad/randomKMerMatrices");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(
				new FileWriter("/nobackup/afodor_research/af_broad/randomDist.txt"));
		
		writer.write("aFile\tbFile\tdistance\n");
		
		String[] list = KMER_DIST_DIRECTORY.list();
		
		for( int x=0; x < list.length-1; x++)
		{
			String xFile = list[x];
			
			if( xFile.endsWith("dist.txt"))
			{
				HashMap<String, Integer> aMap = WriteDistanceMatrixOneVsAll.getCounts(
					new File(KMER_DIST_DIRECTORY + File.separator + xFile));
				long aSumSquared = WriteDistanceMatrixOneVsAll.getSumSquare(aMap);
				
				for( int y=x+1; y < list.length; y++ )
				{
					String yFile = list[y];
					
					if( yFile.endsWith("dist.txt"))
					{
						HashMap<String, Integer> bMap = WriteDistanceMatrixOneVsAll.getCounts(
							new File(KMER_DIST_DIRECTORY + File.separator + yFile));
						
						writer.write(xFile + "\t");
						writer.write(yFile  + "\t");
						double distance = WriteDistanceMatrixOneVsAll.getDistance(aMap, bMap, aSumSquared);
						writer.write(distance + "\n");
						writer.flush();
					}
				}
			}
		}
		
		writer.flush();  writer.close();
	}
}