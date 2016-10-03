/*
 * filter mouse reads from fastq files
 */

package kw_jobinBiofilm_rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

public class FilterMouse {
	public static String BASE = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	public static String ALIGNDIR = BASE + "alignToMouse/";
	public static String OUTDIR = BASE + "mouseFilteredFastq/";

	public static void main(String[] args) throws Exception {
		File[] files = new File(ALIGNDIR).listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".mapped.sam")) {
				int numMapped = 0;
				int numReads = 0;
				String fq = "";
				String name = f.getName().replace(".mapped.sam", "").replace("mm10aln.", "");
				
				//get set of mapped reads
				HashSet<String> mapped = new HashSet<String>();
				BufferedReader mouse = new BufferedReader(new FileReader(f));
				for(String line = mouse.readLine(); line != null; line = mouse.readLine()) {
					if(!line.startsWith("@")) {
						mapped.add(line.split("\t")[0]);
						numMapped++;
					} else if(line.startsWith("@PG")) {//get file
						String[] sp = line.split(" ");
						fq = sp[sp.length-1];
					}
				}
				mouse.close();
				if(numMapped != mapped.size()) {
					throw new Exception("Different number to filter: "
							+ name + "\t" + numMapped + "\t" + mapped.size());
				}
				
				//filter fastq
				int numLeft = 0;
				BufferedReader fastq = new BufferedReader(
						new InputStreamReader(new GZIPInputStream(
								new FileInputStream(new File(fq)))));
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(
						OUTDIR + name + ".mouseFiltered.fastq")));
				for(String line1 = fastq.readLine(); line1 != null; line1 = fastq.readLine()) {
					String line2 = fastq.readLine();
					String line3 = fastq.readLine();
					String line4 = fastq.readLine();
					
					String[] sp = line1.split(" ");
					if(!mapped.contains(sp[0].replace("@", ""))) {
						out.write(line1 + "\n" + line2 + "\n" + 
								line3 + "\n" + line4 + "\n");
						numLeft++;
					}
					numReads++;
				}
				
				fastq.close();
				out.close();
				if(numLeft != (numReads - numMapped)) {
					throw new Exception("Filtering failed: " + name + "\t"
							+ numReads + "\t" + numMapped + "\t" + numLeft);
				}
				
				System.out.println(name + "\t" + numReads + "\t" + numMapped + "\t" 
						+ (100.0*numMapped/numReads) + "\t" + numLeft);
			}
		}
	}
}
