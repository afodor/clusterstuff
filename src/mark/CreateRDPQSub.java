package mark;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class CreateRDPQSub
{
	private static int countNum =0 ;
	
	public static final File RDP_OUT_DIR = new File( "/projects/afodor_research/mark/rdpOut");
	public static final File RUN_RDP_DIR = new File("/projects/afodor_research/mark/runRDP");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			RUN_RDP_DIR + File.separator + 	"runAll.sh")));
		
		int x=0;
		for(String s : BreakOutBySample.OUT_SEQUENCE_DIR.list())
		{
			if( s.startsWith("sample"))
			{
				File f = new File(BreakOutBySample.OUT_SEQUENCE_DIR + File.separator + s);
				File shFile = writeCommandsForAFile(f);
				x++;
				writer.write("qsub -q \"viper\" -N \"CountJob" 
						+ x + "\" " + shFile.getAbsolutePath() +  "\n"  );
			}
			
		}
		
		writer.flush();  writer.close();
		
	}
	
	private static File writeCommandsForAFile( File aFile  ) 
			throws Exception
	{
		countNum++;
		File outFile =  new File( RUN_RDP_DIR.getAbsolutePath() + File.separator +  "run_" +  countNum + ".sh");
		
		BufferedWriter writer = new BufferedWriter( 
			new FileWriter(outFile ));
		
		writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
				"-o " + RDP_OUT_DIR.getAbsolutePath() + File.separator + 
							aFile.getName() + "toRdp.txt"
				+ " -q " + aFile.getAbsolutePath()+ "\n" );
				
		writer.flush();  writer.close();
		
		return outFile;
	}
}
