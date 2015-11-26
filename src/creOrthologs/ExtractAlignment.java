package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.List;

import parsers.FastaSequence;
import parsers.HitScores;
import utils.OrderedSequenceRange;

public class ExtractAlignment
{
	private static int getQuerySequenceLength() throws Exception
	{
		List<FastaSequence> list = 
				FastaSequence.readFastaFile(
						"/projects/afodor_research/af_broad/individualBlastRuns/contig_7000000220927531/7000000220927531.fasta");
		
		if( list.size() != 1)
			throw new Exception("No");
		
		return list.get(0).getSequence().length();
	}
	
	public static void main(String[] args) throws Exception
	{
		int queryLength = getQuerySequenceLength();
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
				
				if( hitList.size() > 1)
					throw new Exception("No");
				
				System.out.println(hitList.size());
				
				if( hitList.size() >0 )
				{
					HitScores hs = hitList.get(0);
					
					if( ((float)hs.getQueryAlignmentLength())/ queryLength > 0.9f)
					{
						HashMap<String, FastaSequence> fastaMap = FastaSequence.getFirstTokenSequenceMap(aFile);
						FastaSequence aSeq = fastaMap.get(hs.getTargetId());
						
						if( aSeq == null)
							throw new Exception("Could not find " + hs.getTargetId() + " in " + aFile.getAbsolutePath() );
						
						OrderedSequenceRange osr = hs.getTargetRange();
						
						writer.write(">" + hs.getTargetId() + " " + aFile.getName() + " " 
								+ osr.getStartPosition() + " " + osr.getEndPosition() + " bit=" + hs.getBitScore() + " escore=" + 
											hs.getEScore() + "\n");
						
						writer.write(aSeq.getSequence().substring(osr.getStartPosition() - 1, osr.getEndPosition()) + "\n");
				
					}
				}
			}
			else
			{
				/*  will be queried to itself so will be in the file that way
				List<FastaSequence> fastaList = FastaSequence.readFastaFile(topDir.getAbsolutePath() + File.separator + s);
				
				if( fastaList.size() != 1)
					throw new Exception("No");
				
				writer.write(">" + fastaList.get(0).getHeader()  + "\n");
				writer.write(fastaList.get(0).getSequence() + "\n");
				*/
				
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
