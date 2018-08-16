package farnazAbundantOTU;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import parsers.FastaSequence;
import parsers.FastaSequenceOneAtATime;

public class MergeSequences
{
	private static final String INPUT_DIRECTORY = "/scratch/afodor_research/datasets/GastricBypass/mbp_ready";
	private static final String OUTPUT_DIRECTORY = "/scratch/afodor_research/datasets/GastricBypass/abundantOTU";
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			OUTPUT_DIRECTORY + File.separator + "combinedFasta.txt"	)));
		
		String[] files = new File(INPUT_DIRECTORY).list();
		
		for(String s : files)
			if( s.endsWith(".fasta.gz"))
			{
				String sampleName = s.replace(".fastq.gz", "");
				System.out.println(sampleName);
				
				int index = 1;
				File inFile = new File(INPUT_DIRECTORY + File.separator + s);
				
				FastaSequenceOneAtATime fsoat = new FastaSequenceOneAtATime(inFile,true);
				
				for(FastaSequence fs = fsoat.getNextSequence(); fs != null; fs =fsoat.getNextSequence())
				{
					writer.write(">" + sampleName + "_" + index  + "\n");
					writer.write( fs.getSequence() + "\n");
					index++;
				}
			}
		
		writer.flush();  writer.close();
	}
}
