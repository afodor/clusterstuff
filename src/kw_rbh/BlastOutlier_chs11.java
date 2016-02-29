/*
 * For orthologGroup outliers identified from all v. all MDS plot, blast
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class BlastOutlier_chs11 {
	public static final String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";
	public static final String BLAST_DIR = DIR + "blastOutliers/";
	public static final String DB = "/users/kwinglee/non-redundant-db-nt/nt";

	public static void main(String[] args) throws IOException {
		//get list of outliers
		HashSet<String> outliers = new HashSet<String>();//set of outliers
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + "chs11_kmer_allvall_mds1_v_p_outliers.txt")));
		String line = br.readLine();
		while(line != null) {
			outliers.add(line);
			line = br.readLine();
		}
		br.close();

		//set up blast -> use the sequence of the first gene in the list
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				BLAST_DIR + "runAll_orthogroup.sh")));
		for(String gene : outliers) {
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					BLAST_DIR + "nrBLAST_" + gene)));
			script.write("module load blast\n");
			script.write("blastn -query " + DIR + "geneFastas/carolina_klebsiella_pneumoniae_chs_11.0/" + gene + ".fasta"
					+ " -db " + DB + " -outfmt 7 -out " +
					BLAST_DIR + "/nrBLAST_" + gene + ".txt\n");
			script.close();

			runAll.write("qsub -q \"viper_batch\" nrBLAST_" + gene + "\n");
		}
		runAll.close();

	}
}
