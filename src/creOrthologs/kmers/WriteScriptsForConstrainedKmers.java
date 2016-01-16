package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteScriptsForConstrainedKmers
{
	private static final File CONSTRAINED_KMER_SCRIPT_DIR = 
				new File( "/nobackup/afodor_research/af_broad/scripts/constrainedKmerRunDir");
	
	private static BufferedWriter makeNewWriter( BufferedWriter allWriter, int fileNum ) throws Exception
	{
		File aFile = new File(
				CONSTRAINED_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + 
					"run_" + fileNum + ".sh");
			
		BufferedWriter aWriter = new BufferedWriter(new FileWriter(aFile));
		
		allWriter.write("qsub -q \"viper_batch\" " + aFile.getAbsolutePath() + "\n");
		
		return aWriter;
			
	}
	
	public static void main(String[] args) throws Exception
	{
		String genomePath = 
				"/nobackup/afodor_research/af_broad/carolina/klebsiella_pneumoniae_chs_11.0.scaffolds.fasta";
		
		String contig = "7000000220927538";
		
		BufferedWriter allWriter  = new BufferedWriter(new FileWriter(new File(
			CONSTRAINED_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + 
			"runAll.sh")));
		
		int fileNum =0;
		int index = 0;
		
		BufferedWriter aWriter = makeNewWriter(allWriter, fileNum);
		
		
		for( int x=0; x < 3985959 - 6000; x = x + 1000)
		{
			aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
					+ "creOrthologs.kmers.ConstrainKMersToRegion "  + 
					genomePath + " " + contig + " " +  x + " " + (x + 5001)+ "\n");
				
			aWriter.flush();
			index++;
			
			if( index == 20)
			{
				index=0;
				fileNum++;
				
				aWriter.flush(); aWriter.close();
				
				aWriter = makeNewWriter(allWriter, fileNum);
				
			}
		}
		
		allWriter.flush();  allWriter.close();
	}
}
