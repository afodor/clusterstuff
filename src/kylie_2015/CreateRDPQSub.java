package kylie_2015;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class CreateRDPQSub
{
	private static int countNum =0 ;
	public static final File FASTQ_DIR =new File("/projects/afodor_research/kylie_2015/fastq/");
	public static final File FASTA_DIR = new File("/projects/afodor_research/kylie_2015/fasta/");
	public static final File RUN_RDP_DIR = new File("/projects/afodor_research/kylie_2015/runRDPDir/");
	public static final File RDP_OUT_DIR = new File("/projects/afodor_research/kylie_2015/rdpOutDir/");
	
	public static void main(String[] args) throws Exception
	{
			
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				 RUN_RDP_DIR + File.separator +  "runAll.sh")));
		
		int x=0;
		for(String s : FASTQ_DIR.list())
		{
			File seqFile = new File(FASTQ_DIR.getAbsolutePath() + File.separator + s);
			File shFile = writeCommandsForAFile(seqFile);
			
			x++;
			writer.write("qsub -q \"viper\" -N \"CountJob" 
					+ x + "\" " + shFile.getAbsolutePath() +  "\n"  );
		}
		
		writer.flush();  writer.close();
		
	}

	
	private static File writeCommandsForAFile( File aFile  ) 
			throws Exception
	{
		countNum++;
		File outFile =  new File( RUN_RDP_DIR + File.separator +  "qsubTarget" + countNum + ".sh");
		
		BufferedWriter writer = new BufferedWriter( 
			new FileWriter(outFile ));
		
		String basePath = FASTA_DIR.getAbsolutePath() + File.separator +  
							aFile.getName().replace(".fastq.gz", "") + 
								".fasta";
		
		writer.write("java -cp /users/afodor/gitInstall/clusterstuff/bin " + 
			"parsers.FastQToFastA" + " " + aFile.getAbsolutePath() +
			" " + basePath + "\n"
						);
		writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
				"-o " + RDP_OUT_DIR.getAbsolutePath() + File.separator + aFile.getName().replace(".fastq.gz", "") 
									+ "_TO_RDP.txt" + " -q " + basePath+ "\n" );
				
		writer.flush();  writer.close();
		
		return outFile;
	}
}
