package vanderbilt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class CreateRDPQSub
{
	private static int countNum =0 ;
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/projects/afodor/vanderbilt/runAll/runAll.sh")));
		
		
		int x=0;
		for(String s : SplitIntoSeparateFiles.outDir.list())
		{
			if( s.endsWith(".fasta"))
			{
				File f = new File(SplitIntoSeparateFiles.outDir.getAbsolutePath() + File.separator + s);
				writeCommandsForAFile(f);
				x++;
				writer.write("qsub -q \"viper\" -N \"CountJob" 
						+ x + "\" " + f.getAbsolutePath() +  "\n"  );
			}
			
		}
		
		writer.flush();  writer.close();
		
	}

	
	
	private static File writeCommandsForAFile( File aFile  ) 
			throws Exception
	{
		countNum++;
		File outFile =  new File("/projects/afodor/vanderbilt/runAll/" + countNum + ".sh");
		
		BufferedWriter writer = new BufferedWriter( 
			new FileWriter(outFile ));
		
		writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
				"-o " + "/projects/afodor/vanderbilt/rdpResults/" + aFile.getName().replace(".fasta", "toRdp.txt") 
				+ " -q " + aFile.getAbsolutePath()+ "\n" );
				
		writer.flush();  writer.close();
		
		return outFile;
	}
}
