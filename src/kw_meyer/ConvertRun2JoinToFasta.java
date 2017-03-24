/*
 * Convert run2 joined reads from fastq to fasta
 */
package kw_meyer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ConvertRun2JoinToFasta {
	private static String DIR = "/nobackup/afodor_research/kwinglee/meyer/run2joinedReads/";
	
	public static void main(String[] args) throws IOException {
		File[] fqs = new File(DIR).listFiles();
		for(File fq : fqs) {
			if(fq.getName().endsWith("join.fastq")) {
				BufferedReader br = new BufferedReader(new FileReader(fq));
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(
						fq.getAbsolutePath().replace(".fastq", "fasta"))));
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
