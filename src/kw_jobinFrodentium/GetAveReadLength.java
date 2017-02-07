/*
 * Calculate average read lengths for each group
 */
package kw_jobinFrodentium;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GetAveReadLength {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/F_rodentium/";
	private static String OUTDIR = DIR + "readLengths";
	
	public static void main(String[] args) throws IOException {
		File temp = new File(OUTDIR);
		if(!temp.exists()) {
			temp.mkdirs();
		}
		
		//run scripts for each set of data
		String[] fastaDirs = new String[]{"/nobackup/afodor_research/kwinglee/jobin/apcTumor/fastas/",
				"/nobackup/afodor_research/kwinglee/jobin/biofilm/fastas/",
				"/nobackup/afodor_research/kwinglee/jobin/gemcitabine/demultiplexedReads/"};
		String[] names = new String[]{"apc", "biofilm", "biofilmReassoc"};
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				OUTDIR + "datasetAverage.txt")));
		all.write("dataset\tnumSamples\taveLength\n");
		for(int i = 0; i < fastaDirs.length; i++) {
			String fd = fastaDirs[i];
			File[] fastas = new File(fd).listFiles();
			double setAve = 0;
			int setNumSamps = 0;
			BufferedWriter set = new BufferedWriter(new FileWriter(new File(
					OUTDIR + names[i] +"_average.txt")));
			set.write("file\tnumReads\ttotLength\taveLength\n");
			for(File f : fastas) {
				String n = f.getName();
				if(n.endsWith(".fasta") && ! n.contains("PCR") &&
						!n.contains("NC101") && !n.contains("H20") &&
						!n.contains("water") && !n.contains("other") &&
						(!fd.contains("gemcitabine") || n.startsWith("B"))) {
					double len = 0;
					int numReads = 0;
					BufferedReader fa = new BufferedReader(new FileReader(f));
					for(String line = fa.readLine(); line != null; line = fa.readLine()) {
						if(line.startsWith(">")) {
							numReads++;
						} else {
							len += line.length();
						}
					}
					fa.close();
					double ave = len/numReads;
					set.write(n.replace(".fasta", "") + "\t" + numReads + "\t" +
							len + "\t" + ave + "\n");
					setNumSamps++;
					setAve += ave;
				}
			}
			set.close();
			all.write(names[i] + "\t" + setNumSamps + "\t"
					+ (setAve/setNumSamps) + "\n");
		}
		all.close();
	}
}
