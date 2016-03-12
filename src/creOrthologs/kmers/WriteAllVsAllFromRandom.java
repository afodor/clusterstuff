package creOrthologs.kmers;

import java.io.File;

public class WriteAllVsAllFromRandom
{
	private static final File KMER_DIST_DIRECTORY = 
			new File("/nobackup/afodor_research/af_broad/randomKMerMatrices");
	
	public static void main(String[] args) throws Exception
	{
		String[] list = KMER_DIST_DIRECTORY.list();
		
		for( int x=0; x < list.length-1; x++)
		{
			String xFile = list[x];
			
			if( xFile.endsWith("dist.txt"))
			{
				for( int y=x+1; y < list.length; y++ )
				{
					String yFile = list[y];
					
					if( yFile.endsWith("dist.txt"))
					{
						
					}
				}
			}
		}
	}
}