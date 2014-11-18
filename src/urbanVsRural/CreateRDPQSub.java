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
		List<File> allShFiles = new ArrayList<File>();
		
		writeCommandsForAllSubDirectories(allShFiles,
		"/projects/afodor/ChinaSequences/first/microbiome/F14FTSUSAT0494_HUMmaxM/Clean");

		writeCommandsForAllSubDirectories(allShFiles,
		"/projects/afodor/ChinaSequences/first/microbiome/F14FTSUSAT0494_HUMmaxM_1022/Clean");
		
		writeCommandsForAllSubDirectories(allShFiles,
		"/projects/afodor/ChinaSequences/second/microbiome/F14FTSUSAT0494_HUMssuM/Clean");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/users/afodor/runChina/runAll.sh")));
		
		int x=0;
		for(File f : allShFiles)
		{
			x++;
			writer.write("qsub -q \"viper\" -N \"CountJob" 
					+ x + "\" " + f.getAbsolutePath() +  "\n"  );
		}
		
		writer.flush();  writer.close();
		
	}

	private static void 
	writeCommandsForAllSubDirectories( List<File> allShFiles , String parentPath )
		throws Exception
	{
		File parentdir = new File(parentPath);
		for(String s : parentdir.list() )
		{
			File childDir = new File( parentdir + File.separator + s );
			
			if( childDir.isDirectory())
			{
				for(String s2 : childDir.list())
				{
					allShFiles.add(writeCommandsForAFile(
							new File(
						childDir.getAbsoluteFile() + File.separator + s2 )));
				}
			}
		}
	}
	
	private static File writeCommandsForAFile( File aFile  ) 
			throws Exception
	{
		countNum++;
		File outFile =  new File("/users/afodor/runChina/qsubTarget" + countNum + ".sh");
		
		BufferedWriter writer = new BufferedWriter( 
			new FileWriter(outFile ));
		
		String basePath = "/projects/afodor/ChinaSequences/rdpResults" + 
							aFile.getName() + "_" + aFile.getParentFile().getName() +
								".fasta";
				
		
		writer.write("java -cp /users/afodor/gitInstall/clusterstuff/bin " + 
			"parsers.FastQToFastA" + " " + aFile.getAbsolutePath() +
			" " + basePath + "\n"
						);
		writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
				"-o " + basePath + "_TO_RDP.txt" + " -q " + basePath+ "\n" );
				
		writer.flush();  writer.close();
		
		return outFile;
	}
}
