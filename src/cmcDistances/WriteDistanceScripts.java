package cmcDistances;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteDistanceScripts
{
	private static File SCRIPT_DIRECTORY = new File(
			"/nobackup/afodor_research/fromNury52016/scripts/snpScripts");
	
	private static File SNP_DIRECTORY = new File(
			"/nobackup/afodor_research/fromNury52016/snpFiles");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
			SCRIPT_DIRECTORY + File.separator + "runAll.sh"	)));
		
		String[] files = WriteContextScripts.OUTPUT_DIRECTORY.list();
		
		int index = 1;
		
		for(int x=0; x < files.length-1; x++)
		{
			if ( files[x].endsWith(".context"))
			{
				File xFile = new File(WriteContextScripts.OUTPUT_DIRECTORY
								+ File.separator + files[x]);
				
				for( int y=x+1; y < files.length; y++)
					if( files[y].endsWith(".context"))
					{
						File yFile = new File(WriteContextScripts.OUTPUT_DIRECTORY
										+ File.separator + files[y]);
						
						File outFile = new File(SNP_DIRECTORY.getAbsolutePath() + File.separator + 
								xFile.getName().replaceAll(".context", "") + "@" + 
									yFile.getName().replaceAll(".context", "") + ".txt");
						
						File scriptFile = new File(SCRIPT_DIRECTORY.getAbsolutePath() + File.separator 
								+ "run_" + index + ".sh");
						
						BufferedWriter writer = new BufferedWriter(new FileWriter(scriptFile));
						
						writer.write("java  -Xmx20g -cp /users/afodor/gitInstall/clusterstuff/bin cmcDistances.WriteSNPFile " + 
										 xFile.getAbsolutePath() + " " + yFile.getAbsolutePath() + " "
										 		+ outFile.getAbsolutePath() +"\n");
						
						allWriter.write("qsub -q\"viper_batch\" " + scriptFile.getAbsolutePath() + "\n");
						
						writer.flush();  writer.close();
						index++;
					}
			}
		}
		
		allWriter.flush();  allWriter.close();
	}
}
