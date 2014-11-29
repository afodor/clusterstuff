package urbanVsRural.abundantOTU;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

public class WriteBlastSH
{
	/*
	 * 
#this is blast 2.2.29+
module load blast
/apps/pkg/ncbi-blast-2.2.29+/rhel6_u5-x86_64/gnu/bin/blastn -db /projects/afodor
/silva/SILVA_119_SSURef_tax_silva.fasta -out /projects/afodor/ChinaSequences/abu
ndantOtuForwardResults/runBlast/cons1ToSilva.txt -query /projects/afodor/ChinaSe
quences/abundantOtuForwardResults/runBlast/Consensus27.fasta -outfmt 6
	 */
	
	public static void main(String[] args) throws Exception
	{
		List<FastaSequence> list = 
				FastaSequence.readFastaFile("/projects/afodor/ChinaSequences/abundantOtuForwardResults/chinaForward.cons");
		
		BufferedWriter allWriter= new BufferedWriter(new FileWriter(new File( 
				"/projects/afodor/ChinaSequences/abundantOtuForwardResults/runBlast/runAll.sh")));
		
		for(FastaSequence fs : list)
		{
			File outFile = new File("/projects/afodor/ChinaSequences/abundantOtuForwardResults/runBlast/" +
					"run_" + fs.getFirstTokenOfHeader() + ".sh");
			
			allWriter.write("qsub -q \"viper\" " + outFile.getAbsolutePath() +  "\n"  );
			
			BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		
			writer.write("module load blast\n");
			writer.write("/apps/pkg/ncbi-blast-2.2.29+/rhel6_u5-x86_64/gnu/bin/blastn ");
			writer.write("-db /projects/afodor/silva/SILVA_119_SSURef_tax_silva.fasta ");
			writer.write(" -out /projects/afodor/ChinaSequences/abundantOtuForwardResults/runBlast/" +
				fs.getFirstTokenOfHeader() + 	"ToSilva.txt ");
			writer.write(" -query /projects/afodor/ChinaSequences/abundantOtuForwardResults/runBlast/" + 
					fs.getFirstTokenOfHeader() + ".fasta " );
			writer.write(" -outfmt 6 \n" );
			writer.flush();  writer.close();
		}
		
		allWriter.flush();  allWriter.close();
	}
	
}
