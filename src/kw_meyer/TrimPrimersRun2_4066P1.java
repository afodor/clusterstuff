/*
 * Trim primers from run2 4066-P1 reads and generate both fastq and fasta
 */
package kw_meyer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class TrimPrimersRun2_4066P1 {
	private static String DIR = "/nobackup/afodor_research/kwinglee/meyer/";
	private static String INDIR = DIR + "cardia_seq2_fastqs/4066-KAM-P1-34087111/";
	public static String OUTDIR = DIR + "run2filteredSeqs/";

	public static void main(String[] args) throws IOException {
		File[] samples = new File(INDIR).listFiles();
		System.out.println(samples.length);
		for(File s : samples) { // for each sample in each run, filter adapter
			filter(s);
		}
	}

	private static void filter(File file) throws IOException {
		File[] list = file.listFiles();
		File read1 = null;
		File read2 = null;
		if(list[0].getName().contains("_R1_")) {
			read1 = list[0];
			read2 = list[1];
		} else if(list[0].getName().contains("_R2_")) {
			read1 = list[1];
			read2 = list[0];
		} else {
			System.out.println("Bad set of reads: " + file.getName());
			System.exit(1);
		}

		BufferedReader br1 = new BufferedReader(
				new InputStreamReader(
						new GZIPInputStream(
								new FileInputStream(read1))));
		BufferedReader br2 = new BufferedReader(
				new InputStreamReader(
						new GZIPInputStream(
								new FileInputStream(read2))));
		BufferedWriter fq1 = new BufferedWriter(
				new FileWriter(
						new File(TrimPrimersRun2_4066P1.OUTDIR + "P1" + read1.getName().replace(".gz", ""))));
		BufferedWriter fq2 = new BufferedWriter(
				new FileWriter(
						new File(TrimPrimersRun2_4066P1.OUTDIR + "P1" + read2.getName().replace(".gz", ""))));
		BufferedWriter fa1 = new BufferedWriter(
				new FileWriter(
						new File(TrimPrimersRun2_4066P1.OUTDIR + "P1" + read1.getName().replace("q.gz", "a"))));
		BufferedWriter fa2 = new BufferedWriter(
				new FileWriter(
						new File(TrimPrimersRun2_4066P1.OUTDIR + "P1" + read2.getName().replace("q.gz", "a"))));
		String head1 = br1.readLine();
		String head2 = br2.readLine();
		int numRemR1 = 0;
		int numRemR2 = 0;
		while(head1 != null && head2 != null) {
			String seq1 = br1.readLine();
			String seq2 = br2.readLine();
			String plus1 = br1.readLine();
			String plus2 = br2.readLine();
			String qual1 = br1.readLine();
			String qual2 = br2.readLine();

			//remove primers
			int len1 = seq1.length();
			int len2 = seq2.length();
			seq1 = seq1.replaceAll("^CCTACGGG[AGTC]GGC[AT]GCAG", "");
			seq2 = seq2.replaceAll("^GACTAC[ACT][ACG]GGGTATCTAATCC", "");

			if(seq1.length() != len1 - 17) {
				numRemR1++;
			} else if(seq2.length() != len2 - 21) {
				numRemR2++;
			} else {
				qual1 = qual1.substring(17);
				qual2 = qual2.substring(21);
				//write fastq
				fq1.write(head1 + "\n" + seq1 + "\n" + plus1 + "\n" + qual1 + "\n");
				fq2.write(head2 + "\n" + seq2 + "\n" + plus2 + "\n" + qual2 + "\n");

				//write fasta
				fa1.write(head1.replace("@", ">") + "\n" + seq1 + "\n");
				fa2.write(head2.replace("@", ">") + "\n" + seq2 + "\n");
			}

			head1 = br1.readLine();
			head2 = br2.readLine();
		}


		br1.close();
		br2.close();
		fq1.close();
		fq2.close();
		fa1.close();
		fa2.close();

		System.out.println(file.getName() + " removedR1: " + numRemR1
				+ " removedR2: " + numRemR2 + " totalRemoved: " + (numRemR1+numRemR2));
	}

}