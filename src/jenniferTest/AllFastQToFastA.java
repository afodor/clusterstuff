package jenniferTest;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class AllFastQToFastA
{
	private static void recurse(File dir, List<File> list) throws Exception
	{
		for(String s : dir.list())
		{
			File f = new File(dir.getAbsolutePath() + File.separator + s);
			
			if( f.isDirectory())
				recurse(f, list);
			else
				list.add(f);
		}
	}

	public static void main(String[] args) throws Exception
	{
		File startDir = new File("/projects/afodor_research/jwellerTestRunRawData");
		File rdpDir = new File("/projects/afodor_research/jwellerTest/rdpOutDir");
		File runDir = new File("/projects/afodor_research/jwellerTest/runDir");
		File fastaDir = new File("/projects/afodor_research/jwellerTest/fastaDir");
		
		BufferedWriter allSHWriter = new BufferedWriter(new FileWriter(new File(
				runDir.getAbsolutePath() + File.separator + "runAll.sh")));
		
		List<File> list = new ArrayList<File>();
		recurse(startDir, list);
		
		int index=1;
		for( File f : list)
		{
			File shFile = new File(runDir.getAbsolutePath() + "run_" + index + ".sh");
			index++;
			allSHWriter.write("qsub -q \"viper\" " + shFile.getAbsolutePath() +  "\n"  );
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(shFile));
			
			if (f.getName().endsWith("fastq.gz"))
			{
				System.out.println(f.getAbsolutePath());
				
				String newName = f.getName().replace("fastq.gz", "fasta.txt");
				
				File fastAFile = new File(fastaDir.getAbsolutePath() + File.separator + 
						newName	);

				writer.write("java -cp /users/afodor/gitInstall/clusterstuff/bin " + 
						"parsers.FastQToFastA" + " " + f.getAbsolutePath() +
							" " + fastAFile + "\n");

				File rdpFile = new File(rdpDir.getAbsolutePath() + File.separator + 
								newName.replace("fasta.txt", "")+ "_rdpOut.txt"	);
				
				writer.write("java -jar /users/afodor/rdp/rdp_classifier_2.10.1/dist/classifier.jar " + 
						"-o "+ rdpFile.getAbsolutePath() + " -q " + fastAFile+ "\n" );
				
				writer.flush();  writer.close();

			}
		}
		
		allSHWriter.flush();  allSHWriter.close();

	}
}