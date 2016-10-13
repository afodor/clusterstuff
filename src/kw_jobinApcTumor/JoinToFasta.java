/*
 * convert output of the join to fasta (is currently fastq)
 */
package kw_jobinApcTumor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class JoinToFasta {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/apcTumor/joinedReads/";
	
	public static void main(String[] args) throws IOException {
		String[] fqs = new File(DIR).list();
		for(String f : fqs) {
			if(f.endsWith("join.fastq")) {
				BufferedReader fq = new BufferedReader(new FileReader(new File(DIR + f)));
				BufferedWriter fa = new BufferedWriter(new FileWriter(new File(
						DIR + f.replace(".fastq", ".fasta"))));
				String header = fq.readLine();
				while(header != null) {
					String seq = fq.readLine();
					fq.readLine();
					fq.readLine();
					fa.write(header.replace("@", ">") + "\n"
							+ seq + "\n");
					header = fq.readLine();
				}
				fq.close();
				fa.close();
			}
		}
	}

}
