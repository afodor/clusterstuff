/*
 * align 100 reads from two 16S samples to E. coli to figure out which variable region
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

public class Determine16SVariableRegion {
	private static String DIR = "/nobackup/afodor_research/kwinglee/china/";
	private static String VDIR = DIR + "16SvariableRegion/";
	private static String FQDIR = DIR + "fastqs_16s/";
	private static String REF = VDIR + "E_coli_J0169516S_rRNA.fasta";
	private static int MAXREADS = 100;

	public static void main(String[] args) throws IOException {
		//set up BWA script
		BufferedWriter bwa = new BufferedWriter(new FileWriter(new File(
				VDIR + "bwaScript")));
		String bwaDB = VDIR + "rRNAbwaDB";
		bwa.write("#PBS -l procs=1\n");
		bwa.write("module load bwa\n");
		bwa.write("bwa index -p " + bwaDB + " " + REF + "\n");
		
		//set up BLAST script
		BufferedWriter blast = new BufferedWriter(new FileWriter(new File(
				VDIR + "blastScript")));
		String blastDB = VDIR + "rRNAblastDB";
		blast.write("#PBS -l procs=1\n");
		blast.write("module load blast\n");
		blast.write("makeblastdb -in " + REF + " -dbtype nucl -out "
				+ blastDB + "\n");
		
		//make 100 read fasta files and add to scripts
		String[] files = new String[]{"first_15A_1.fq.gz", "second_B120A_2.fq.gz"};
		for(String f : files) {
			String faName = VDIR + f.replace(".fq.gz", ".fasta");
			bwa.write("bwa mem " + bwaDB + " " + faName
					+ " > " + faName.replace(".fasta", ".bwa.sam") + "\n");
			blast.write("blastn -outfmt 6 -db " + blastDB + 
					" -query " + faName
					+ " -out " + faName.replace(".fasta", ".blast.txt") + "\n");
			BufferedReader fq = new BufferedReader(new InputStreamReader(
					new GZIPInputStream(new FileInputStream(FQDIR + f))));
			BufferedWriter fa = new BufferedWriter(new FileWriter(new File(faName)));
			String line1 = fq.readLine();
			int numReads = 0;
			while(numReads < MAXREADS) {
				numReads++;
				String line2 = fq.readLine();
				/*String line3 =*/ fq.readLine();
				/*String line4 =*/ fq.readLine();
				fa.write(line1.replaceFirst("@", ">") + "\n" + line2 + "\n");
				line1 = fq.readLine();
			}
			fq.close();
			fa.close();
		}
		
		bwa.close();
		blast.close();
	}
}
