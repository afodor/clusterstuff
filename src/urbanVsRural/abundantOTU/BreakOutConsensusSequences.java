package urbanVsRural.abundantOTU;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

public class BreakOutConsensusSequences
{
	public static void main(String[] args) throws Exception
	{
		List<FastaSequence> list = 
				FastaSequence.readFastaFile("/projects/afodor/ChinaSequences/abundantOtuForwardResults");
		
		for(FastaSequence fs : list)
		{
			File outFile = new File("/projects/afodor/ChinaSequences/abundantOtuForwardResults/" +
						fs.getFirstTokenOfHeader() + ".fasta");
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
			
			
			writer.write(">" + fs.getFirstTokenOfHeader() + "\n");
			writer.write(fs.getSequence() + "\n");
			writer.flush();  writer.close();
		}
	}
}
