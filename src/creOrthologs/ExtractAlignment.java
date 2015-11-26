package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;
import parsers.HitScores;

public class ExtractAlignment
{
	public static void main(String[] args) throws Exception
	{
		File topDir = new File( "/projects/afodor_research/af_broad/individualBlastRuns/contig_7000000220927531");
		
		String[] list = topDir.list();
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter( "/projects/afodor_research/af_broad/contig_7000000220927531_forAlign.txt" ));
		
		for( String s : list) 
		{
			if( ! s.endsWith("fasta"))
			{
				System.out.println(s);
				String[] splits = s.split("_"
						+ "");
				
				StringBuffer b = new StringBuffer();
				
				for( int x=2; x < splits.length; x++)
					b.append(splits[x] + (x < splits.length - 1 ? "_" : ""));
				
				String aName = 
						(b.toString() + ".scaffolds.fasta").replace(".txt", "");
				File aFile = findFile(aName);
				System.out.println(aFile.getAbsolutePath());
				
				List<HitScores> hitList = HitScores.getTopHits(topDir.getAbsolutePath() + File.separator + s);
				
				if( hitList.size() != 1)
					throw new Exception("No");
			}
			else
			{
				List<FastaSequence> fastaList = FastaSequence.readFastaFile(topDir.getAbsolutePath() + File.separator + s);
				
				if( fastaList.size() != 1)
					throw new Exception("No");
				
				writer.write(">" + fastaList.get(0).getHeader()  + "\n");
				writer.write(fastaList.get(1).getSequence() + "\n");
				
			}
			
		}
		
		writer.flush();  writer.close();
		
	}
	
	private static File findFile(String genome) throws Exception
	{
		for(String d : RunBlastAll.DIRECTORIES)
		{
			File genomeDir = new File("/projects/afodor_research/af_broad" + File.separator + d + File.separator +  genome);
			
			if(genomeDir.exists())
				return genomeDir;
		}	
		
		
		throw new Exception("Could not find " + genome);
	}
}
