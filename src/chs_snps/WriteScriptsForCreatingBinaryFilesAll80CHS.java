package chs_snps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import utils.ConfigReader;


public class WriteScriptsForCreatingBinaryFilesAll80CHS
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File("/projects/afodor_chs/kwinglee/cophylog_all80chs/run/runAll.sh")));
		
		File folder = new File("/projects/afodor_research/mjzapata/CRE/CHS_raw/");
		File[] files = folder.listFiles();
		
		for( int i=0; i < files.length; i++)
		{
			String inFile = files[i].getName();
			if(inFile.endsWith(".fastq.gz")) {
				inFile = inFile.replace(".fastq.gz", "");
				File outFile = new File("/projects/afodor_chs/kwinglee/cophylog_all80chs/run/run_" + inFile);
				
				BufferedWriter aWriter = new BufferedWriter(new FileWriter(outFile));
				
				//aWriter.write("#PBS -l nodes=1:ppn=8\n");
				//aWriter.write("#PBS -W x=NODESET:ONEOF:FEATURE:ib_ddr\n");
				aWriter.write("#PBS -l nodes=1:ppn=12\n");
				
				writer.write("qsub -q \"Cobra_batch\" run_" + inFile  + "\n");
				
				
				aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin -mx20000m chs_snps.WriteBinaryContextsFromFastQ "
						+ "/projects/afodor_research/mjzapata/CRE/CHS_raw/" + inFile + ".fastq.gz /projects/afodor_chs/kwinglee/cophylog_all80chs/context" + 
							inFile + "_context.gz");
				
				aWriter.flush();  aWriter.close();
			}
		}
		
		writer.flush();  writer.close();
	}
}
