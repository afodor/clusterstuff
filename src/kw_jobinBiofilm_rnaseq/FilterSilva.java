/*
 * filter silva reads from already mouse filtered fastq files
 */

package kw_jobinBiofilm_rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashSet;
import java.util.zip.GZIPInputStream;

public class FilterSilva {
	public static String BASE = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	public static String ALIGNDIR = BASE + "alignToSilva/";//location of alignment results
	public static String OUTDIR = BASE + "mouseAndSilvaFiltered/";//where to write the filtered files
	public static String INDIR = BASE + "mouseFilteredFastq/";//location of original fastq files

	public static void main(String[] args) throws Exception {
		File[] files = new File(ALIGNDIR).listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".mapped.sam") && f.getName().startsWith("lsu")) {
				int numReads = 0;
				String name = f.getName().replace(".mapped.sam", "").replace("lsu.aln.", "");
				
				//get set of mapped reads
				HashSet<String> mapped = getMappedReads(f);//large subunit
				HashSet<String> ssu = getMappedReads(new File(ALIGNDIR 
						+ f.getName().replace("lsu", "ssu")));//small subunit
				mapped.addAll(ssu);
				int numMapped = mapped.size();
				
				//filter fastq
				int numLeft = 0;
				BufferedReader fastq = new BufferedReader(
						new FileReader(new File(INDIR + name + ".fastq")));
				BufferedWriter fqout = new BufferedWriter(new FileWriter(new File(
						OUTDIR + name + ".silvaFiltered.fastq")));
				BufferedWriter faout = new BufferedWriter(new FileWriter(new File(
						OUTDIR + name + ".silvaFiltered.fasta")));
				for(String line1 = fastq.readLine(); line1 != null; line1 = fastq.readLine()) {
					String line2 = fastq.readLine();
					String line3 = fastq.readLine();
					String line4 = fastq.readLine();
					
					String[] sp = line1.split(" ");
					if(!mapped.contains(sp[0].replace("@", ""))) {
						fqout.write(line1 + "\n" + line2 + "\n" + 
								line3 + "\n" + line4 + "\n");
						faout.write(line1.replace("@", ">") + "\n"
								+ line2 + "\n");
						numLeft++;
					}
					numReads++;
				}
				
				fastq.close();
				fqout.close();
				faout.close();
				if(numLeft != (numReads - numMapped)) {
					throw new Exception("Filtering failed: " + name + "\t"
							+ numReads + "\t" + numMapped + "\t" + numLeft);
				}
				
				System.out.println(name + "\t" + numReads + "\t" + numMapped + "\t" 
						+ (100.0*numMapped/numReads) + "\t" + numLeft);
			}
		}
	}
	
	private static HashSet<String> getMappedReads(File f) throws IOException {
		HashSet<String> mapped = new HashSet<String>();
		BufferedReader mouse = new BufferedReader(new FileReader(f));
		for(String line = mouse.readLine(); line != null; line = mouse.readLine()) {
			if(!line.startsWith("@")) {
				mapped.add(line.split("\t")[0]);
			}
		}
		mouse.close();
		return(mapped);
	}
}
