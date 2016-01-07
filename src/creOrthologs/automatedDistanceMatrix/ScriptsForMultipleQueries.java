package creOrthologs.automatedDistanceMatrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import creOrthologs.RunBlastAll;

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
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File(
			SCRIPT_DIR.getAbsolutePath() + File.separator + "runAll.sh"	)));
		
		File queryFile = writeOneExtractionFile(INPUT_GENOME, "7000000220927533", 729729, 749719);
		
		allWriter.write("qsub -q \"viper\" " + queryFile.getAbsolutePath());
		
		allWriter.flush();  allWriter.close();
	}
	
	private static void addBlastRuns(BufferedWriter writer, File queryFile,
			File genomePath, String contig, 
			int startPos, int endPos
			) throws Exception
	{
		File topDir = new File(  ExtractOne.WORKING_DIR.getAbsolutePath() + File.separator + 
				genomePath.getName() +"_" + contig + "_" + startPos + "_" + endPos );
		
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/nobackup/afodor_research/af_broad" + File.separator + d);
				
			String[] list = genomeDir.list();

			for( String s : list)
			{	
				if( s.endsWith("fasta"))
				{
					File inSeqs= new File( genomeDir.getAbsolutePath() + File.separator + s);
					File outFile = new File(topDir.getAbsolutePath()+ 
						"contig_" + contig+ File.separator + startPos + "_to_" + endPos + "_" + 
								s.replaceAll(".scaffolds.fasta", ".txt"));
						
					writer.write("module load blast\n");
					writer.write("/apps/pkg/ncbi-blast-2.2.29+/rhel6_u5-x86_64/gnu/bin/blastn -db " + 
								inSeqs.getAbsolutePath() + " -out " + 
								outFile.getAbsolutePath() +  
								" -query " + queryFile.getAbsolutePath() +
								" -outfmt 6\n");
				
						writer.flush();
				}
			}
		}
	}
	
	private static File writeOneExtractionFile( File genomePath, String contig, 
										int startPos, int endPos)  throws Exception
	{
		File shFile = new File(SCRIPT_DIR + File.separator + 
				genomePath.getName() + "_" + contig + "_" + startPos + "_"+ endPos + ".sh");
		
		File outFile = new File(OUTPUT_DIR + File.separator +
				genomePath.getName() + "_" + contig + "_" + startPos + "_"+ endPos + "_dist.txt");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(shFile));
		
		writer.write(
				"java -cp /users/afodor/gitInstall/clusterstuff/bin/ " + 
		"creOrthologs.automatedDistanceMatrix.ExtractOne " 
		+ genomePath.getAbsolutePath() + " " + contig + " " + startPos + " " + endPos + " " + 
					outFile.getAbsolutePath() +  "\n" );
				
		addBlastRuns(writer, outFile.getAbsoluteFile(), genomePath, contig, startPos, endPos);
		
		writer.flush(); writer.close();
		
		return shFile;
	}
	
}
