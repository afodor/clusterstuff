package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class RerunIfNotAtTarget
{
	private static int getNumLines(File file) throws Exception
	{
		int count=0;
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		for(String s=  reader.readLine(); s != null; s= reader.readLine())
			count++;
		
		reader.close();
		
		return count;
	}
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			WriteScriptsForAllDistance.KMER_RUN_DIR.getAbsolutePath() + File.separator + "runAll.sh")));
		
		int index = 1;
		for(String s : WriteDistanceMatrixOneVsAll.KMER_DISTANCE_MATRIX_DIR.list())
		{
			if( s.endsWith("_toAll.txt"))
			{
				if( getNumLines(new File(WriteDistanceMatrixOneVsAll.KMER_DISTANCE_MATRIX_DIR.getAbsolutePath() + 
						File.separator + s)) != 339)
				{
					File aFile = new File(WriteScriptsForAllDistance.KMER_RUN_DIR.getAbsolutePath() + 
							File.separator + "run" + index + ".sh");

					BufferedWriter aWriter = new BufferedWriter(new FileWriter(aFile));
					
					aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
						+ "creOrthologs.kmers.WriteDistanceMatrixOneVsAll "  
							+ s.replace("_toAll.txt", "_kmers.txt") + "\n");
					
					writer.write("qsub -q \"viper_batch\" " + aFile.getAbsolutePath() + "\n");
					
					aWriter.flush();  aWriter.close();
					index++;
				}
			
			}
		}
		
		writer.flush();  writer.close();
	}
}
