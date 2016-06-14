package adenomasRelease;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class CreateRDPQSub
{
	private static int countNum =0 ;
	
	public static final File SEQ_DIR = new File("/nobackup/afodor_research/adenomasRelease/seqs");
	public static final File RDP_OUT_DIR = new File( "/nobackup/afodor_research/adenomasRelease/rdpOut");
	public static final File SCRIPT_DIR = new File("/nobackup/afodor_research/adenomasRelease/scripts/runRDP");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR.getAbsolutePath() + File.separator + "runAll.sh")));
			
		for(String s : SEQ_DIR.list())
		{
			countNum++;
			File fastaFile = new File(SEQ_DIR + File.separator + 
						s);
			
			File rdpOutFile = new File(RDP_OUT_DIR.getAbsolutePath() + File.separator + 
					s  + "toRDP.txt");
			
			File runFile = new File(SCRIPT_DIR.getAbsoluteFile() + File.separator + "run_" + 
						countNum + ".sh");
			
			BufferedWriter writer = new BufferedWriter( new FileWriter(runFile));
			
			writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
					"-o " + rdpOutFile.getAbsolutePath()  + " -q " + fastaFile+ "\n" );
			
			allWriter.write("qsub -q \"viper\" " + runFile.getAbsolutePath() +  "\n"  );
					
			writer.flush();  writer.close();
			allWriter.flush();
		}
		
		allWriter.flush();  allWriter.close();
	}
}
