package jamesDada2Stream;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteBlastSH
{
	/*
	 * 
#this is blast 2.2.30
module load blast
/apps/pkg/ncbi-blast-2.3.0+/rhel7_u2-x86_64/gnu/bin/blastn/blastn 
 * -db /scratch/afodor_research/silva/SILVA_132_SSURef_tax_silva.fasta -out /scratch/afodor_research/datasets/uegp/dada2ToBlast/blastOut
  -outfmt 6
	 */
	public static File blastOutDir = new File("/scratch/afodor_research/datasets/uegp/dada2ToBlast/blastOut");
	public static File scriptsDir = new File("/scratch/afodor_research/datasets/uegp/dada2ToBlast/scripts");
	
	public static void main(String[] args) throws Exception
	{
		String [] list = new File(BreakUpFastaFile.OUT_FILE_DIR).list();
		
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(scriptsDir.getAbsolutePath() +
						File.separator + "all.sh"));
		
		for(String s  : list)
		{
			if ( s.endsWith(".fasta"))
			{
				File inFile = new File(BreakUpFastaFile.OUT_FILE_DIR + File.separator + s);
				
				File outFile = new File(blastOutDir + File.separator + 
						s.replace(".fasta", "toSilva.txt"));
				
				allWriter.write("qsub -q \"copperhead\" " + outFile.getAbsolutePath() +  "\n"  );
				

				BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
				
				writer.write("#!/bin/sh\n");
				writer.write("#PBS -q copperhead\n");
				writer.write("PBS -l nodes=1:ppn=1\n");
				
				writer.write("module load blast\n");
				writer.write("/apps/pkg/ncbi-blast-2.3.0+/rhel7_u2-x86_64/gnu/bin/blastn ");
				writer.write(" -db /scratch/afodor_research/silva/SILVA_132_SSURef_tax_silva.fasta ");
				writer.write(" -out " + outFile.getAbsolutePath() );
				writer.write(" -query " + inFile.getAbsolutePath()  );
				writer.write(" -outfmt 6 \n" );
				writer.flush();  writer.close();
				
			
			}
		}
		
		allWriter.flush();  allWriter.close();
	}
	
}
