package vanderbilt.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteTestScripts
{
	private static String RDP_DIR = "/projects/afodor/vanderbilt/rdpResults";
	private static String CLASSPATH_DIR = "/users/afodor/gitInstall/clusterstuff/src";
	private static String TEST_SCRIPT_DIR = "/projects/afodor/vanderbilt/test/scripts";
	private static String TEST_RECOUNT_DIR = "/projects/afodor/vanderbilt/test/recountDir";
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter= new BufferedWriter(new FileWriter(new File(TEST_SCRIPT_DIR + File.separator + "runAll.sh")));
		
		File rdpDir =new File(RDP_DIR);
		
		for( String s : rdpDir.list())
		{
			File rdpFile =  new File( rdpDir.getAbsolutePath() + File.separator + s);
			
			File scriptFile = new File(TEST_SCRIPT_DIR + File.separator + rdpFile.getName()+ ".sh");
			
			File outputFile = new File(TEST_RECOUNT_DIR + File.separator + rdpFile.getName() + ".counts");
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(scriptFile));
			
			writer.write("java -cp " + CLASSPATH_DIR + "vanderbilt.test.WriteExpectedFromRDP " +  rdpFile.getAbsolutePath() + " "+ 
							outputFile.getAbsolutePath());
			
			allWriter.write("qsub -q \"viper\""  + scriptFile.getAbsolutePath() +  "\n"  );
			
			writer.flush();  writer.close();
		}
		
		allWriter.flush();  allWriter.close();
	}
}
