/*
 * generates the fasta files and scripts to align beta-lactamase genes
 * Klebsiella pneumonia only
 */
package kw_cre;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class AlignSelectCardsHitsSHV {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static String FASTA_DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/geneFastas/";

	public static void main(String[] args) throws Exception {
		String[] files = new String[]{"groupPos4A", "groupPos4C", 
				"groupPos4other", "groupPos4otherPlusLen", "groupAllPlusLen",
				"groupAllNoTem", "groupAllNoTemPlusLen"};//files to analyze

		String outDir = DIR + "betaLactamaseAlignments/";

		//generate fasta file for each file (representing a group of similar genes)
		for(String f : files) {
			BufferedReader shvFile = new BufferedReader(new FileReader(new File(
					DIR + "bestBlastHit_SHV_" + f + ".txt")));
			BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
					outDir + "SHV_" + f + ".fasta")));
			String line = shvFile.readLine();//header
			for(line = shvFile.readLine(); line != null; line = shvFile.readLine()) {
				String gene = line.split("\t")[0];

				//get fasta sequence
				String[] name = gene.split("_A[E|F]");
				if(name.length != 2) {
					fasta.close();
					throw new Exception("bad gene split " + name.length + " " + gene);
				}
				BufferedReader seq = new BufferedReader(new FileReader(new File(
						FASTA_DIR + name[0] + File.separator + gene + ".fasta")));
				for(String seqline = seq.readLine(); seqline != null; seqline = seq.readLine()) {
					fasta.write(seqline + "\n");
				}
				seq.close();
			}
			shvFile.close();
			fasta.close();
		}
		
		//generate scripts
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				outDir + "alignSHVCommands")));
		script.write("module load mafft\n");
		for(String f : files) {
			script.write("echo \"\" >&2\n");
			script.write("echo \"Starting " + f + "\" >&2\n");//write gene to standard error so can tell differences
			script.write("echo \"\" >&2\n");
			script.write("mafft --auto --adjustdirection " + outDir + "SHV_" + f + ".fasta" 
					+ " > " + outDir + "SHV_" + f + ".mafftAlign.fasta\n");
		}
		script.close();
	}
}

