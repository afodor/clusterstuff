package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteScriptsForAllDistance
{
	public static final File KMER_RUN_DIR = new File( "/nobackup/afodor_research/af_broad/scripts/kmerRunDir");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				KMER_RUN_DIR.getAbsolutePath() + File.separator + "runAll.sh")));
		
		
		int index = 1;
		for(String s : MakeKmers.KMER_DIR.list())
		{
			if( s.endsWith("_kmers.txt"))
			{
				File aFile = new File(KMER_RUN_DIR.getAbsolutePath() + File.separator + 
											"run" + index + ".sh");
				
				BufferedWriter aWriter = new BufferedWriter(new FileWriter(aFile));
				
				aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
					+ "creOrthologs.kmers.WriteDistanceMatrixOneVsAll "  + s + "\n");
				
				writer.write("qsub -q \"viper_batch\" " + aFile.getAbsolutePath() + "\n");
				
				aWriter.flush();  aWriter.close();
				index++;
			}
		}
		
		writer.flush();  writer.close();
	}
}
