/*
 * Convert fastq.gz files to fasta for RDP
 */
package kw_jobinDolphin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class FastqToFasta {
	public static final String FQ_DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/fastqs/";
	public static final String FA_DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/fastas/";
	
	public static void main(String[] args) throws IOException {
		String[] samples = new File(FQ_DIR).list();
		
		for(String s : samples) {
			String path = FQ_DIR + s + "/Data/Intensities/BaseCalls/";
			String[] reads = new File(path).list();
			for(String r : reads) {
				BufferedReader fq = new BufferedReader(
						new InputStreamReader(
								new GZIPInputStream(
										new FileInputStream(
												path + r))));
				String faName = r.replace(".fastq.gz", ".fasta");
				BufferedWriter fa = new BufferedWriter(
						new FileWriter(
								new File(FA_DIR + faName)));
				String line1 = fq.readLine();
				while(line1 != null) {
					String line2 = fq.readLine();
					fq.readLine();//+
					fq.readLine();//quality score
					String header = line1.replaceFirst("@", ">");
					fa.write(header + "\n" + line2 + "\n");					
					line1 = fq.readLine();
				}
				fa.close();
				fq.close();
			}
		}
	}

}
