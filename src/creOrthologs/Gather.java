package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import parsers.HitScores;

public class Gather
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(
				new File("/projects/afodor_research/af_broad/blastResults.txt")));
	
		writer.write( "query\ttarget\tgenomeFileName\tmaxBitScore\n");
		
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
					File outSubDir = new File( RunBlastAll.BLAST_RESULTS_PATH + File.separator + s.replaceAll(".scaffolds.fasta",""));
					
					for(String s2 : innerList)
					{
						if( s2.endsWith("dnaseq"))
						{
							File outFile = new File( outSubDir.getAbsolutePath()+ File.separator + 
									s.replaceAll(".scaffolds.fasta","") + "_" + d + "_to_" + 
									s2.replaceAll("dnaseq", "") + "txt.gz");
							
							if( outFile.exists())
							{
								numFound++;
								
								HitScores hs = HitScores.getTopHitByBitScore(outFile);
								
								//writer.write( "query\ttarget\tparentDir\tmaxBitScore\n");
								
								writer.write(hs.getQueryId() + "\t" + hs.getTargetId() + "\t" + 
										s.replaceAll(".scaffolds.fasta","") + "\t" + hs.getBitScore() + "\n"
											);
							}
							else
							{
								//logWriter.write("Could not find " + outFile.getAbsoluteFile() + "\n");
							}
								
							numDone++;
							
							if( numDone % 1000 == 0 )
							{
								logWriter.write(numDone + " " + numFound + " " +((double)numFound / numDone) + "\n");
								logWriter.flush();
							}	
						}
					}
				}
			}
			
		}
		
		logWriter.flush();  logWriter.close();
		writer.flush();  writer.close();
	}
}
