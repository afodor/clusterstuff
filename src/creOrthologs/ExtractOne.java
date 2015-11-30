package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

public class ExtractOne
{
	public static void main(String[] args) throws Exception
	{
		List<FastaSequence> list = 
				FastaSequence.readFastaFile(
						"/projects/afodor_research/af_broad/carolina/klebsiella_pneumoniae_chs_11.0.scaffolds.fasta");
		
		String toFind = "7000000220927531";
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
			"/projects/afodor_research/af_broad/individualBlastRuns/contig_" + toFind + 
			File.separator + toFind + ".fasta")));
		
		
		for(FastaSequence fs : list)
		{
			if(fs.getFirstTokenOfHeader().equals(toFind))
			{
				writer.write(fs.getHeader() + "\n");
				writer.write(fs.getSequence().substring(2482, 4032) + "\n");
			}
		}
			
		writer.flush();  writer.close();
	}
}
