package creOrthologs.automatedDistanceMatrix;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

public class GenerateDistanceMatrix
{
	private static final File WORKING_DIR = new File("/nobackup/afodor_research/af_broad/workingDir");
	
	public static void main(String[] args) throws Exception
	{
		if(args.length != 5)
		{
			System.out.println("Usage inputGenomeFilePath contigName startPos endPos outDistanceMatrix");
			System.exit(1);
		}
		
		File genomePath = new File(args[0]);
		
		File topDir = new File(WORKING_DIR.getAbsolutePath() + File.separator + 
						genomePath.getName() +"_" + args[1] + "_" + args[2] + "_" + args[3] );
		
		topDir.mkdir();
		writeQuerySequence(args, topDir);
	}
	
	private static File writeQuerySequence( String[] args , File topDir )
		throws Exception
	{
		List<FastaSequence> list = 
				FastaSequence.readFastaFile(args[0]);
		
		File genomePath = new File(args[0]);
		String toFind = args[1];
		
		File outFile = 
				new File(
						topDir.getAbsolutePath() + File.separator + 
						 genomePath.getName() +"_" + args[1] + "_" + args[2] + "_" + args[3] + ".fasta");
		
		boolean foundOne = false;
			
		for(FastaSequence fs : list)
		{
			if( !foundOne && fs.getFirstTokenOfHeader().equals(toFind))
			{
				BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
				
				writer.write(fs.getHeader() + "\n");
				writer.write(fs.getSequence().substring(Integer.parseInt(args[2]), 
						Integer.parseInt(args[3])) + "\n");
				
				foundOne = true;

				writer.flush();  writer.close();
			}
		}
			
		if(!foundOne )
			throw new Exception("Could not find query sequence " + 
					args[0] +"_" + args[1] + "_" + args[2] + "_" + args[3] );
		
		return outFile;
	}		
}
