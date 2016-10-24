/*
 * Merge the blast outputs and write a script to generate the linked files for humann
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class MergeBlastResults {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/rnaseq/";
	private static String RESDIR = DIR + "kegg_split_blastx_mouseSilvaFiltered/";
	private static String OUTDIR = DIR + "kegg_blastx_mouseSilvaFiltered_merge/";
	private static String FASTA_DIR = DIR + "mouseAndSilvaFiltered/";
	private static String HUMINPUTDIR = "/nobackup/afodor_research/kwinglee/humann-0.99/input/";
	private static int NUM_SAMPLES = 6;
	
	public static void main(String[] args) throws Exception {
		//get list of genomes
		String[] fas = new File(FASTA_DIR).list();
		String[] samples = new String[NUM_SAMPLES];
		int numSamps = 0;
		for(String fa : fas) {
			if(fa.endsWith(".fasta") && fa.contains("_R1")) {
				samples[numSamps] = fa.replace(".mouseFiltered.silvaFiltered.fasta", "");
				numSamps++;
			}
		}
		if(NUM_SAMPLES != numSamps) {
			throw new Exception("Wrong number of samples: " + NUM_SAMPLES + " " + numSamps);
		}
		
		//set up script to write links
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(OUTDIR + "humannLink.sh")));
		
		//combine genomes and add to script
		for(String samp : samples) {
			String fileName = "kegg_" + samp + ".txt";
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(
					OUTDIR + fileName)));
			script.write("ln -s " + OUTDIR + fileName + HUMINPUTDIR + ".\n");
			String[] reads = new String[]{"_R1", "_R2"};
			for(String read : reads) {
				int split = 0;
				File file = new File(RESDIR + "kegg_" + samp + read + "_split" + split + ".txt");
				while(file.exists()) {
					BufferedReader in = new BufferedReader(new FileReader(file));
					for(String line = in.readLine(); line != null; line = in.readLine()) {
						out.write(line + "\n");
					}
					in.close();
					split++;
					file = new File(RESDIR + "kegg_" + samp + read + "_split" + split + ".txt");
				}
				System.out.println(samp + read + " numSplits: " + (split-1));
			}
			out.close();
		}
		script.close();
	}
}
