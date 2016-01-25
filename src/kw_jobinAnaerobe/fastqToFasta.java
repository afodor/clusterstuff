/**
 * convert joined reads from fastq to fasta
 */
package kw_jobinAnaerobe;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class fastqToFasta {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/anaerobe/join/";
	
	public static void main(String[] args) throws IOException {
		File[] files = new File(DIR).listFiles();
		for(File f : files) {
			String name = f.getName();
			if(name.endsWith("join.fastq")) {
				BufferedReader in = new BufferedReader(new FileReader(f));
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(
						DIR + name.replace("fastq", "fasta"))));
				String header = in.readLine();
				while(header != null) {
					String seq = in.readLine();
					in.readLine();
					in.readLine();
					out.write(header.replaceFirst("@", ">") + "\n" + seq + "\n");
					header = in.readLine();
				}
				in.close();
				out.close();
			}
		}
	}

}
