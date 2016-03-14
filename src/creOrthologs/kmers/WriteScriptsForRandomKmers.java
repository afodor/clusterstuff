package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteScriptsForRandomKmers
{
	public static final File KMER_RUN_DIR = 
			new File( "/nobackup/afodor_research/af_broad/scripts/runRandomKmers");
	
	public static final File QUERY_FILE = 
			new File("/nobackup/afodor_research/af_broad/kmers/klebsiella_pneumoniae_chs_11.0_kmers.txt");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				KMER_RUN_DIR.getAbsolutePath() + File.separator + "runAll.sh")));
		
		int index = 1;
		for(int x= 5000; x <= 30000; x =x + 100)
		{
			File aFile = new File(KMER_RUN_DIR.getAbsolutePath() + File.separator + 
											"run" + index + ".sh");
				
			BufferedWriter aWriter = new BufferedWriter(new FileWriter(aFile));
				
			aWriter.write("java -mx20000m  -cp " + 
			"/users/afodor/gitInstall/clusterstuff/bin creOrthologs.kmers.RandomKMers " + 
				x + " " + QUERY_FILE.getAbsolutePath() + 
				" /nobackup/afodor_research/af_broad/randomKMerMatrices/chs11_" + x +"_dist.txt " + 
				"/nobackup/afodor_research/af_broad/randomKMerMatrices/chs11_" + x + "_key.txt \n");
				
			writer.write("qsub -q \"viper_batch\" " + aFile.getAbsolutePath() + "\n");
				
			aWriter.flush();  aWriter.close();
			index++;
		}
		
		writer.flush();  writer.close();
	}
}
