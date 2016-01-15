package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class MakeScriptsToReRunErrors
{
	public static void main(String[] args) throws Exception
	{
		File errorFile = new File("/nobackup/afodor_research/af_broad/scripts/kmerRunDir/rerun/errors.txt");
		
		BufferedReader reader = new BufferedReader(new FileReader(errorFile));
		
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/scripts/kmerRunDir/rerun/runAll.txt")));
		
		int count=1;
		for(String s= reader.readLine(); s != null; s= reader.readLine())
		{
			File shFile = new File("/nobackup/afodor_research/af_broad/scripts/kmerRunDir/rerun/" + count
							+".sh");
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(shFile));
			
			s = s.substring(s.indexOf("java") -1, s.length());
			
			writer.write(s + "\n");
			writer.flush();  writer.close();
			
			allWriter.write("qsub -q \"viper\" " + shFile.getAbsolutePath());
			
		}
		
		allWriter.flush();  allWriter.close();
		reader.close();
	}
}
