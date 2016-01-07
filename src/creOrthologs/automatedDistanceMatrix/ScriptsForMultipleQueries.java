package creOrthologs.automatedDistanceMatrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class ScriptsForMultipleQueries
{
	private static File SCRIPT_DIR = 
		new File("/nobackup/afodor_research/af_broad/scripts/runQueriesToDistance");
	
	private static File OUTPUT_DIR = 
			new File("/nobackup/afodor_research/af_broad/distanceMatrices");
	
	private static File INPUT_GENOME = 
			new File("/nobackup/afodor_research/af_broad/carolina/klebsiella_pneumoniae_chs_11.0.scaffolds.fasta");
	
	public static void main(String[] args) throws Exception
	{
		writeOneFile(INPUT_GENOME, "7000000220927533", 729729, 749719);
	}
	
	private static File writeOneFile( File genomePath, String contig, 
										int startPos, int endPos)  throws Exception
	{
		File shFile = new File(SCRIPT_DIR + File.separator + 
				genomePath.getName() + "_" + contig + "_" + startPos + "_"+ endPos + ".sh");
		
		File outFile = new File(OUTPUT_DIR + File.separator +
				genomePath.getName() + "_" + contig + "_" + startPos + "_"+ endPos + "_dist.txt");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(shFile));
		
		writer.write(
				"java -cp /users/afodor/gitInstall/clusterstuff/bin/ " + 
		"creOrthologs.automatedDistanceMatrix.GenerateDistanceMatrix " 
		+ genomePath.getAbsolutePath() + " " + contig + " " + startPos + " " + endPos + " " + 
					outFile.getAbsolutePath() +  "\n" );
				
		writer.flush(); writer.close();
		
		return shFile;
	}
	
}
