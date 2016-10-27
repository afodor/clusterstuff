/*
 * Get the number of 16S reads for each sample
 */

package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class GetNumber16sReads {
	public static String DIR = "/nobackup/afodor_research/kwinglee/china/";
	public static String FQDIR = DIR + "fastqs_16s/";
	
	public static void main(String[] args) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "number16Sreads.txt")));
		out.write("file\tnumberReads\n");
		File[] fqs = new File(FQDIR).listFiles();
		for(File fq : fqs) {
			if(fq.getName().endsWith(".fq.gz")) {
				BufferedReader reads = new BufferedReader(
						new InputStreamReader(
								new GZIPInputStream(
										new FileInputStream(fq))));
				int numReads = 0;
				for(String line = reads.readLine(); line != null; line = reads.readLine()) {
					numReads++;
				}
				reads.close();
				out.write(fq.getName() + "\t" + (numReads / 4) + "\n");
			}
		}
		out.close();
	}

}
