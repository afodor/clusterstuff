package creOrthologs;

import java.io.File;

public class ExtractAlignment
{
	public static void main(String[] args) throws Exception
	{
		File topDir = new File( "/projects/afodor_research/af_broad/individualBlastRuns/contig_7000000220927531");
		
		String[] list = topDir.list();
		
		for( String s : list) 
		{
			if( ! s.endsWith("fasta"))
			{
				System.out.println(s);
				String[] splits = s.split("_");
				String aName = splits[2] + "_" + splits[3] + "_" + splits[4] + "_" + splits[5] + ".scaffolds.fasta";
				File aFile = findFile(aName);
				System.out.println(aFile.getAbsolutePath());	
			}
			
		}
		
		
	}
	
	private static File findFile(String genome) throws Exception
	{
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + genome);
			
			if(genomeDir.exists())
				return genomeDir;
		}	
		
		
		throw new Exception("Could not find " + genome);
	}
}
