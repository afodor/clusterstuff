package urbanVsRural;

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
		File file = new File("/projects/afodor/hmp/stoolBySample");
		
		for(String s : file.list())
			writeCommandsForAFile(s);
	}
	
	private static File writeCommandsForAFile( String filename ) 
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
