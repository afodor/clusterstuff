package creOrthologs.kmers;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

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
		
		List<FastaSequence> list = FastaSequence.readFastaFile(genomePath);
		
		int fileNum =0;
		int index = 0;
				
		BufferedWriter allWriter  = new BufferedWriter(new FileWriter(new File(
			CONSTRAINED_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + 
			"runAll.sh")));
		
		BufferedWriter aWriter = makeNewWriter(allWriter, fileNum);
		
		for( FastaSequence fs : list)
		{
			String contig = fs.getFirstTokenOfHeader();
			String seq = fs.getSequence();
			int length = seq.length();
			int slice = Math.min(5000, length-1);
			
			for( int x=0; x < Math.max(length-6000, 1); x = x + 1000)
			{
				aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
						+ "creOrthologs.kmers.ConstrainKMersToRegion "  + 
						genomePath + " " + contig + " " +  x + " " + (x + slice)+ "\n");
					
				aWriter.flush();
				index++;
				
				if( index == 25)
				{
					index=0;
					fileNum++;
					
					aWriter.flush(); aWriter.close();
					
					aWriter = makeNewWriter(allWriter, fileNum);
					
				}
			}
			
		}
		
		
		allWriter.flush();  allWriter.close();
	}
}
