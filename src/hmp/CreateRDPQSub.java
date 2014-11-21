package hmp;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class CreateRDPQSub
{
	private static int countNum =0 ;
	
	public static void main(String[] args) throws Exception
	{
		List<File> toRun = new ArrayList<File>();
		File file = new File("/projects/afodor/hmp/stoolBySample/");
		
		for(String s : file.list())
			toRun.add(writeCommandsForAFile(s));
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/users/afodor/hmp/runAll.sh")));
		
		int x=0;
		for(File f : toRun)
		{
			x++;
			writer.write("qsub -q \"viper\" -N \"CountJob" 
					+ x + "\" " + f.getAbsolutePath() +  "\n"  );
		}
		
		writer.flush();  writer.close();
	}
	
	private static File writeCommandsForAFile( String filename  ) 
			throws Exception
	{
		countNum++;
		File outFile =  new File("/users/afodor/runRDP/run_" + countNum + ".sh");
		
		BufferedWriter writer = new BufferedWriter( 
			new FileWriter(outFile ));
	
		writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
				"-o " + "/projects/afodor/rdpResults/" + filename + 
				"_TO_RDP.txt" + " -q " + "/projects/afodor/rdpResults/stoolBySample" + 
						filename+ "\n" );
				
		writer.flush();  writer.close();
		
		return outFile;
	}
}
