package chs_snps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;


public class WriteScriptsForCreatingBinaryFiles
{
	public static void main(String[] args) throws Exception
	{
		int index=1;
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("/projects/afodor_chs/run/runAll.sh")));
		
		for( int x=1; x <=9; x++)
		{
			File outFile = new File("/projects/afodor_chs/run/run_" + index);
			
			BufferedWriter aWriter = new BufferedWriter(new FileWriter(outFile));
			
			writer.write("qsub -q \"viper\" run_" + index);
			
			aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin -mx2800m chs_snps.WriteBinaryContexts "
					+ "/projects/afodor_chs/fasta/chs_241_" + index + " /projects/afodor_chs/context/chs_241_" + 
						index + "_context.gz");
			
			aWriter.flush();  aWriter.close();
			index++;
		}
		
		for( int x=1; x <=8; x++)
		{
			File outFile = new File("/projects/afodor_chs/run/run_" + index);
			
			
			BufferedWriter aWriter = new BufferedWriter(new FileWriter(outFile));
			
			writer.write("qsub -q \"viper\" run_" + index);
			
			aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin -mx2800m chs_snps.WriteBinaryContexts "
					+ "/projects/afodor_chs/fasta/chs_242_" + index + " /projects/afodor_chs/context/chs_242_" + 
						index + "_context.gz");
			
			aWriter.flush();  aWriter.close();
			index++;
		}
		
		writer.flush();  writer.close();
	}
}
