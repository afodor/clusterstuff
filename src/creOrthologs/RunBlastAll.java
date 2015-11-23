package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class RunBlastAll
{
	public static final String BLAST_RESULTS_PATH = "/projects/afodor_research/af_broad/blastResultsByDir";
	
	public static final String[] DIRECTORIES = { "carolina", "resistant", "susceptible" };
	
	private static final String SCRIPT_DIR = 
			"/projects/afodor_research/af_broad/scripts/runBlastAllScripts";
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter runAllWriter = new BufferedWriter(new FileWriter(new File( 
			SCRIPT_DIR + File.separator +  "runAll.sh"	)));
		String[] innerList = MakeBlastDB.ORTHOLOG_DIR.list();
		
		int index = 1;
		for(String d : DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + d);
			
			String[] list = genomeDir.list();

			for( String s : list)
			{	
				if( s.endsWith("fasta"))
				{
					File outSubDir = new File(BLAST_RESULTS_PATH + File.separator + s.replaceAll(".scaffolds.fasta",""));
					
					if( ! outSubDir.exists())
						outSubDir.mkdir();
					
					File shellFile = new File( SCRIPT_DIR + File.separator + "run_" + index + ".sh" );
					runAllWriter.write("qsub -q \"viper_batch\" -N \"CountJob" 
							+ index + "\" " + shellFile.getAbsolutePath() +  "\n" );
					runAllWriter.flush();
					
					index = index + 1;

					BufferedWriter writer = new BufferedWriter(new FileWriter(shellFile));
					writer.write("module load blast\n");
					
					for(String s2 : innerList)
					{
						if( s2.endsWith("dnaseq"))
						{
							File databaseFile = new File( MakeBlastDB.ORTHOLOG_DIR.getAbsoluteFile() + File.separator + 
													s2);
							
							
							File outFile = new File( outSubDir.getAbsolutePath()+ File.separator + 
									s.replaceAll(".scaffolds.fasta","") + "_" + d + "_to_" + 
									s2.replaceAll("dnaseq", "") + "txt");
							
							writer.write("/apps/pkg/ncbi-blast-2.2.29+/rhel6_u5-x86_64/gnu/bin/blastn -db " + 
							databaseFile.getAbsolutePath() + " -out " + 
								outFile.getAbsolutePath() +  
								" -query " + genomeDir.getAbsolutePath() + File.separator + s + 
								" -outfmt 6\n");
							
							writer.write("gzip " + outFile.getAbsolutePath() + "\n");
						}
					}
					
					
					writer.flush();  writer.close();
				}
			}
			
		}
		
		runAllWriter.flush(); runAllWriter.close();
	}
}
