/*
 * Convert run2 4066 P1 joined reads from fastq to fasta
 */
package kw_meyer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ConvertRun2JoinToFasta_4066P1 {
	private static String DIR = "/nobackup/afodor_research/kwinglee/meyer/run2joinedReads/";
	
	public static void main(String[] args) throws IOException {
		File[] fqs = new File(DIR).listFiles();
		for(File fq : fqs) {
			if(fq.getName().endsWith("join.fastq") && fq.getName().startsWith("P1")) {
				BufferedReader br = new BufferedReader(new FileReader(fq));
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(
						fq.getAbsolutePath().replace(".fastq", ".fasta"))));
				String head = br.readLine();
				while(head != null) {
					String seq = br.readLine();
					br.readLine();//plus
					br.readLine();//qual
					out.write(head.replace("@", ">") + "\n" + seq + "\n");
					head = br.readLine();
				}
				br.close();
				out.close();
			}
		}
	}

}
