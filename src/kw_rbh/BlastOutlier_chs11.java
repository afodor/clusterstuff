/*
 * For chs11 outliers identified from all v. all MDS plot, blast
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
	public static final String DB_NT = "/nobackup/afodor_research/kwinglee/ncbi-non-redundant-db-nt/nt";
	public static final String DB_NR = "/nobackup/afodor_research/kwinglee/ncbi-non-redundant-db-nr/nr";

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
		/*BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				BLAST_DIR + "runAll_chs11.sh")));*/
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				BLAST_DIR + "ntBLAST_chs11")));
		BufferedWriter scriptx = new BufferedWriter(new FileWriter(new File(
				BLAST_DIR + "nrBLAST_chs11")));
		script.write("module load blast\n");
		scriptx.write("module load blast\n");
		for(String gene : outliers) {
			//String[] g = gene.split("_");
			//String name = g[g.length-2] + "_" + g[g.length-1];
			/*BufferedWriter script = new BufferedWriter(new FileWriter(new File(
					BLAST_DIR + "nrBLAST_" + name)));
			script.write("module load blast\n");*/
			script.write("blastn -query " + DIR + "geneFastas/carolina_klebsiella_pneumoniae_chs_11.0/" + gene + ".fasta"
					+ " -db " + DB_NT + " -outfmt 7 -out " +
					BLAST_DIR + "ntBLASTresults_" + gene + ".txt\n");
			scriptx.write("blastx -query " + DIR + "geneFastas/carolina_klebsiella_pneumoniae_chs_11.0/" + gene + ".fasta"
					+ " -db " + DB_NR + " -outfmt 7 -out " +
					BLAST_DIR + "nrBLASTresults_" + gene + ".txt\n");
			//script.close();
			
			//get the fasta file
			BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
					BLAST_DIR + gene + ".fasta")));
			BufferedReader fa = new BufferedReader(new FileReader(new File(
					DIR + "geneFastas/carolina_klebsiella_pneumoniae_chs_11.0/" + gene + ".fasta")));
			String s = fa.readLine();
			while(s != null) {
				fasta.write(s + "\n");
				s = fa.readLine();
			}
			fasta.close();
			fa.close();

			//runAll.write("qsub -q \"viper_batch\" nrBLAST_" + name + "\n");
		}
		//runAll.close();
		script.close();
		scriptx.close();

	}
}
