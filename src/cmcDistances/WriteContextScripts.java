package cmcDistances;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteContextScripts
{
	private static File INPUT_DIRECTORY = new File(
			"/projects/afodor_chs/From Nury 5-2016");
	
	private static File OUTPUT_DIRECTORY = new File(
			"/nobackup/afodor_research/fromNury52016/contexts");
	
	private static File SCRIPT_DIRECTORY = new File(
			"/nobackup/afodor_research/fromNury52016/scripts/makeContexts");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
			SCRIPT_DIRECTORY + File.separator + "runAll.sh"	)));
		
		String[] files = INPUT_DIRECTORY.list();
		
		int index = 1;
		
		for(String s : files)
		{
			if ( s.endsWith(".fastq"))
			{
				File inFile = new File(INPUT_DIRECTORY.getAbsolutePath() + File.separator + s);
				
				File outFile = new File(OUTPUT_DIRECTORY.getAbsolutePath() + File.separator + 
						s.replaceAll(".fastq", "") + ".context");
						
				File scriptFile = new File(SCRIPT_DIRECTORY.getAbsolutePath() + File.separator 
						+ "run_" + index + ".sh");
				
				BufferedWriter writer = new BufferedWriter(new FileWriter(scriptFile));
				
				writer.write("java -cp /users/afodor/gitInstall/clusterstuff/bin -Xmx20g " + 
								inFile.getAbsolutePath() + " " + outFile.getAbsolutePath() + " 15 " 
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
