package jobinCardio;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class CreateRDPQSub
{
	private static int countNum =0 ;
	
	public static final File RDP_OUT_DIR = new File( "/projects/afodor_research/JobinCardio/rdpOut");
	public static final File RDP_RUN_DIR = new File("/projects/afodor_research/JobinCardio/rdpRunDir");
	
	public static void main(String[] args) throws Exception
	{
		List<File> allShFiles = new ArrayList<File>();
		
		for(String s : Demultiplex.fastaOutDir.list())
		{
			countNum++;
			File fastaFile = new File(Demultiplex.fastaOutDir.getAbsolutePath() + File.separator + 
						s);
			
			File rdpOutFile = new File(RDP_OUT_DIR.getAbsolutePath() + File.separator + 
					s  + "toRDP.txt");
			
			File runFile = new File(RDP_RUN_DIR.getAbsoluteFile() + File.separator + "run_" + 
						countNum + ".sh");
			
			BufferedWriter writer = new BufferedWriter( new FileWriter(runFile));
			
			writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
					"-o " + rdpOutFile.getAbsolutePath()  + " -q " + fastaFile+ "\n" );
					
			writer.flush();  writer.close();
			
			allShFiles.add(runFile);
			
		}
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			RDP_RUN_DIR.getAbsoluteFile() + File.separator + "runAll.sh")));
		
		int x=0;
		for(File f : allShFiles)
		{
			x++;
			writer.write("qsub -q \"viper\" -N \"CountJob" 
					+ x + "\" " + f.getAbsolutePath() +  "\n"  );
		}
		
		writer.flush();  writer.close();
		
	}
}
