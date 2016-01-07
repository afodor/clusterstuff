package creOrthologs.automatedDistanceMatrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.List;

import creOrthologs.RunBlastAll;
import parsers.FastaSequence;
import parsers.HitScores;
import utils.OrderedSequenceRange;

public class ExtractForAlignment
{
	private static int getQuerySequenceLength(File querySequence) throws Exception
	{
		List<FastaSequence> list = 
				FastaSequence.readFastaFile(querySequence);
		
		if( list.size() != 1)
			throw new Exception("No");
		
		return list.get(0).getSequence().length();
	}
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 1)
		{
			System.out.println("usage workingDirectory");
		}
		
		File workingDir = new File(args[0]);
		
		if(!  workingDir.exists()  || workingDir.isDirectory())
			throw new Exception(workingDir.getAbsolutePath() + " does not appear to be a directory");
		
		File queryFile = null;
		
		for(String s : workingDir.list())
		{
			if( s.endsWith("fasta"))
			{
				if( queryFile != null)
					throw new Exception("Duplicate " + s);
				
				queryFile = new File(workingDir.getAbsolutePath() + File.separator + s);
			}
		}
		
		if( queryFile == null)
			throw new Exception("Could not find any fasta in " + workingDir.getAbsolutePath());
		
		int queryLength = getQuerySequenceLength(queryFile);
		
		BufferedWriter writer = new BufferedWriter(
				new FileWriter( workingDir.getAbsolutePath() + File.separator + "forAlign.align"));
		
		for( String s : workingDir.list()) 
		{
			if(  s.endsWith("_blastOut.txt"))
			{
				String aName = s.replace("_blastOut.txt", "");
				File aFile = findFile(aName);
				
				List<HitScores> hitList = 
						HitScores.getTopHits(workingDir.getAbsolutePath() + File.separator + s);
				
				if( hitList.size() > 1)
					throw new Exception("No");
				
				if( hitList.size() >0 )
				{
					HitScores hs = hitList.get(0);
					
					System.out.println( hs.getTargetId() + " " +  hs.getQueryAlignmentLength() + " "+ queryLength + " " +
							((float)hs.getQueryAlignmentLength()) + " " + hs.getBitScore() + " " + hs.getEScore());
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
