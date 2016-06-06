package seerPipeline;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import creOrthologs.RunBlastAll;

public class MakePhenoFile
{
	public static final File SEER_DIR= new File( "/nobackup/afodor_research/seerStuff");
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				SEER_DIR.getAbsolutePath() + File.separator + "All3.pheno")));
		
		for( String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + d);
			
			String[] list = genomeDir.list();

			for( String s : list)
			{	
				if( s.endsWith("fasta"))
				{
					String cleanS = s.replace(".scaffolds.fasta", "");
					writer.write(cleanS + "\t" + cleanS + "\t" + d + "\n");
				}
			}
		}
		
		writer.flush();  writer.close();
	}
}
