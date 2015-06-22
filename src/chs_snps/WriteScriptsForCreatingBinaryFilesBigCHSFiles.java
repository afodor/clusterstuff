/**
 * generate new run files for the fasta files that ran out of memory
 */
package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import utils.ConfigReader;


public class WriteScriptsForCreatingBinaryFilesBigCHSFiles
{
	//based on whether the context file exists
	public static void main(String[] args) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/cophylog_all80chs/run/runAllBig.sh")));
		
		String seqDir = "/projects/afodor_research/mjzapata/CRE/CHS_raw";
		String contextDir = "/projects/afodor_research/kwinglee/cophylog_all80chs/context/";
		File seqFolder = new File(seqDir);
		File[] seqs = seqFolder.listFiles();
		for(int i=0; i < seqs.length; i++) {
			String seq = seqs[i].getName();
			if(seq.endsWith(".fastq.gz")) {
				seq = seq.replaceAll(".fastq.gz", "");
				File con = new File(contextDir + "context"+seq+"_context.gz");
				System.out.println(con.getName());
				System.out.println(con.exists());
				if(!con.exists()) {
					BufferedWriter aWriter = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/cophylog_all80chs/run/run_big_" + seq)));
					
					aWriter.write("#PBS -l nodes=1:ppn=16\n");
					aWriter.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_qdr2\n");
					aWriter.write("java -cp /users/kwinglee/git/clusterstuff/bin -Xmx124000m chs_snps.WriteBinaryContextsFromFastQ "
							+ "/projects/afodor_research/mjzapata/CRE/CHS_raw/" + seq + ".fastq.gz /projects/afodor_research/kwinglee/cophylog_all80chs/context/context" + 
								seq + "_context.gz");
					
					aWriter.flush();  aWriter.close();
					
					writer.write("qsub -q \"viper_batch\" run_big_" + seq  + "\n");
					
				}
			}
		}
		writer.close();
	}
	
	//based on whether there were any errors
	/*public static void main(String[] args) throws Exception
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

				inFile = inFile.replaceAll("\\.e[0-9]+", "").replace("run_", "");
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
	}*/
}
