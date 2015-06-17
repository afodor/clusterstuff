/**
 * determine the number of reads in each fasta for each strain
 * also write a file that indicates the biggest file for each strain
 */

package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class ReadCounts {
	public static String conversionFile = "/projects/afodor_research/mjzapata/CRE/CHS_raw/chs_batch_download_results.csv";//file containing the conversion
	public static String outFile = "/projects/afodor_research/kwinglee/cophylog_all80chs/ReadCounts.txt";//Name of file to write results to
	public static String seqDir = "/projects/afodor_research/mjzapata/CRE/CHS_raw/";//directory with sequences
	public static String outBig = "/projects/afodor_research/kwinglee/cophylog_all80chs/BiggestSRRPerStrain.txt";//Name of file to write results to
	
	public static void main(String[] args) throws IOException {
		//get conversion file
		BufferedReader convert = new BufferedReader(new FileReader (new File(conversionFile)));
		String line = convert.readLine();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(outFile)));
		out.write("Strain\tSRRFile\tNum_reads\n");
		BufferedWriter big = new BufferedWriter(new FileWriter(new File(outBig)));
		while(line != null) {
			String[] csp = line.split("\t");
			String chs;
			if(csp[0].length() == 1) {
				chs = "CHS0" + csp[0];
			} else {
				chs = "CHS" + csp[0];
			}
			int max = 0;//most number of reads for this chs
			String maxSRR = "";//srr file with max reads
				
			//read counts for each file
			String[] srrlist = csp[1].replace("[", "").replace("]", "").split(",");
			for(int i = 0; i < srrlist.length; i++) {//for each file
				for(int j = 1; j < 2; j++) {//forward and reverse reads
					int num_reads = 0;
					String srr = srrlist[i].trim()+"_"+j;
					
					//read file
					BufferedReader file = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(seqDir+srr+".fastq.gz"))));
					String s = file.readLine();
					while(s != null) {
						num_reads++;
						s = file.readLine();
					}
					file.close();
					
					num_reads/=4;
					if(num_reads > max) {
						max = num_reads;
						maxSRR = srr;
					}
					
					//write results
					out.write(chs+"\t"+srr+"\t"+num_reads+"\n");
				}
			}
			big.write(chs+"\t"+maxSRR+"\n");
			line = convert.readLine();
		}
		out.close();
		big.close();
		convert.close();
	}
}
				
