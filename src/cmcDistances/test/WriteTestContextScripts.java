package cmcDistances.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteTestContextScripts
{
	private static final File CONTEXT_DIR = new File(
			"/nobackup/afodor_research/fromNury52016/contexts"
			);
	
	private static final File TEST_CONTEXT_SCRIPT_DIR = new File(
				"/nobackup/afodor_research/fromNury52016/scripts/compareContexts");

	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
			TEST_CONTEXT_SCRIPT_DIR + File.separator + "runAll.sh"	)));
		
		String[] files = WriteContextScripts.OUTPUT_DIRECTORY.list();
		
		int index = 1;
		
		for(String s : files)
		{
			if ( s.endsWith(".context"))
			{
				File file1 = new File(WriteContextScripts.OUTPUT_DIRECTORY.getAbsolutePath()
							+ File.separator + s);
				
				File file2 = new File(CONTEXT_DIR.getAbsolutePath() + File.separator + 
						s);
						
				File scriptFile = new File(TEST_CONTEXT_SCRIPT_DIR.getAbsolutePath() 
						+ File.separator 
						+ "run_" + index + ".sh");
				
				BufferedWriter writer = new BufferedWriter(new FileWriter(scriptFile));
				
				writer.write("java  -Xmx20g -cp /users/afodor/gitInstall/clusterstuff/bin cmcDistances.test.CompareTwoContexts " + 
								 file1.getAbsolutePath() + " " + file2.getAbsolutePath() 
								+"\n"
						);
				
				allWriter.write("qsub -q\"viper\" " + scriptFile.getAbsolutePath() + "\n");
				
				writer.flush();  writer.close();
				index++;
			}
		}
		
		allWriter.flush();  allWriter.close();
	}
}
