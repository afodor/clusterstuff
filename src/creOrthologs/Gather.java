package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class Gather
{
	public static void main(String[] args) throws Exception
	{
	
		long numDone =0;
		long numFound = 0;
		
		BufferedWriter logWriter =new BufferedWriter(new FileWriter(new File( 
				"/projects/afodor_research/af_broad/log.txt")));
		
		String[] innerList = MakeBlastDB.ORTHOLOG_DIR.list();
		
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + d);
			
			String[] list = genomeDir.list();

			for( String s : list)
			{	
				if( s.endsWith("fasta"))
				{
					for(String s2 : innerList)
					{
						if( s2.endsWith("dnaseq"))
						{
							
							File outFile = new File( RunBlastAll.BLAST_RESULTS_PATH + File.separator + 
									s.replaceAll(".scaffolds.fasta","") + "_" + d + "_to_" + 
									s2.replaceAll("dnaseq", "") + "txt");
							
							if( outFile.exists())
								numFound++;
							
							numDone++;
							
							if( numDone % 1000 == 0 )
							{
								logWriter.write(numDone + " " + numFound + " " +((double)numDone / numFound));
								logWriter.flush();
							}
							
						}
					}
				}
			}
			
		}
		
		logWriter.flush();  logWriter.close();
	}
}
