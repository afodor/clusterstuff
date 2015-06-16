/**
 * generate new run files for the fasta files that ran out of memory
 */
package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import utils.ConfigReader;


public class WriteScriptsForCreatingBinaryFilesBigCHSFiles
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/cophylog_all80chs/run/runAllBig.sh")));
		
		File folder = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/run/");
		File[] files = folder.listFiles();
		
		for( int i=0; i < files.length; i++)
		{
			String inFile = files[i].getName();
			if(inFile.contains(".e")) {
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				String line = br.readLine();
				if(line != null && line.length() > 1) {//had an error

				inFile = inFile.replaceAll("\\.e[0-9]+", "");
				File outFile = new File("/projects/afodor_research/kwinglee/cophylog_all80chs/run/run_big_" + inFile);
				
				BufferedWriter aWriter = new BufferedWriter(new FileWriter(outFile));
				
				//aWriter.write("#PBS -l nodes=1:ppn=8\n");
				//aWriter.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_ddr\n");
				aWriter.write("#PBS -l nodes=1:ppn=12\n");
				
				writer.write("qsub -q \"Cobra_batch\" run_big_" + inFile  + "\n");
				
				
				aWriter.write("java -cp /users/kwinglee/git/clusterstuff/bin -mx30000m chs_snps.WriteBinaryContextsFromFastQ "
						+ "/projects/afodor_research/mjzapata/CRE/CHS_raw/" + inFile + ".fastq.gz /projects/afodor_research/kwinglee/cophylog_all80chs/context/context" + 
							inFile + "_context.gz");
				
				aWriter.flush();  aWriter.close();
				}
			}
		}
		
		writer.flush();  writer.close();
	}
}
