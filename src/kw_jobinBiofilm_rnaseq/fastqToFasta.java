/*
 * convert fastq files to fasta, getting number of reads in process
 * and look at number of reads containing the given index
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class fastqToFasta {
	public static String BASE_DIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/";
	public static String FQ_DIR = BASE_DIR + "RNAseqTestRunFastqs/";
	public static String FA_DIR = BASE_DIR + "RNAseqTestRunFastas/";
	
	public static void main(String[] args) throws IOException {
		convert("RNABFplusTnumber2M115_bc1_L001-39898859/RNABFplusTnumber2M115-bc1-L001_S1_L001_",
				"ATCACG");
		convert("RNABFplusTnumber2M317_bc4_L001-39898862/RNABFplusTnumber2M317-bc4-L001_S2_L001_",
				"TGACCA");
		convert("RNABFplusTnumber2M426_bc5_L001-39898863/RNABFplusTnumber2M426-bc5-L001_S3_L001_",
				"ACAGTG");
		convert("RNABFplusTnumber2M527_bc6_L001-39898861/RNABFplusTnumber2M527-bc6-L001_S4_L001_",
				"GCCAAT");
		convert("RNABFplusTnumber2M635_bc7_L001-39898864/RNABFplusTnumber2M635-bc7-L001_S5_L001_",
				"CAGATC");
		convert("RNABFplusTnumber2M937_bc8_L001-39898860/RNABFplusTnumber2M937-bc8-L001_S6_L001_",
				"ACTTGA");

	}
	
	//for the given filePrefix convert forward and reverse read to fasta
	//look at whether the index is present and the number of reads
	public static void convert(String filePrefix, String index) throws IOException {
		String[] reads = new String[]{"R1", "R2"};
		String sID = filePrefix.split("-")[0];
		for(String r : reads) {
			int numReads = 0;
			int numIndex = 0;
			BufferedReader fq = new BufferedReader(
					new InputStreamReader(new GZIPInputStream(
							new FileInputStream(new File(
					FQ_DIR + filePrefix + r + "_001.fastq.gz")))));
			BufferedWriter fa = new BufferedWriter(new FileWriter(new File(
					FA_DIR + sID + "_" + r + ".fasta")));
			for(String line1 = fq.readLine(); line1 != null; line1 = fq.readLine()) {
				numReads++;
				String line2 = fq.readLine();
				fq.readLine();
				fq.readLine();
				
				fa.write(line1.replace("@", ">") + "\n" + line2 + "\n");
				
				if(line2.contains(index)) {
					numIndex++;
				}
			}
			fq.close();
			fa.close();
			System.out.println(sID + "\t" + r + "\t" + numReads + "\t" + numIndex);
		}
	}
}
