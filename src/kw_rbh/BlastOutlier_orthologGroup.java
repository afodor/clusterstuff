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

public class BlastOutlier_orthologGroup {
	public static final String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";
	public static final String BLAST_DIR = DIR + "blastOutliers/";
	public static final String DB = "/nobackup/afodor_research/kwinglee/ncbi-non-redundant-db-nt/nt";

	public static void main(String[] args) throws IOException {
		//get list of outliers
		HashSet<String> outliers = new HashSet<String>();//set of outliers
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + "orthologGroup_kmer_allvall_mds1_v_p_outliers.txt")));
		String line = br.readLine();
		while(line != null) {
			outliers.add(line);
			line = br.readLine();
		}
		br.close();
		
		//set up blast -> use the sequence of the first gene in the list
		BufferedReader orthogroups = new BufferedReader(new FileReader(new File(
				DIR + "rbhOrthologs/orthologGroups150.txt")));
		/*BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				BLAST_DIR + "runAll_orthogroup.sh")));*/
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				BLAST_DIR + "ntBLAST_orthogroup.sh")));
		script.write("module load blast\n");
		line = orthogroups.readLine();//header
		line = orthogroups.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String orth = sp[0];
			if(outliers.contains(orth)) {//orthogroup is an outlier
				//get names for sequence
				String gene = sp[2].split(";")[0];
				String[] g = gene.split("_");
				String genome = g[0];
				for(int i = 1; i < g.length-2; i++) {
					genome += "_" + g[i];
				}
				
				/*BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						BLAST_DIR + "nrBLAST_" + orth)));
				script.write("module load blast\n");*/
				script.write("blastn -query " + DIR + "geneFastas/" + genome + "/" + gene + ".fasta"
						+ " -db " + DB + " -outfmt 7 -out " +
						BLAST_DIR + "ntBLASTresults_" + orth + ".txt\n");
				/*script.close();
				
				runAll.write("qsub -q \"viper_batch\" nrBLAST_" + orth + "\n");*/
				
			}
			line = orthogroups.readLine();
		}
		orthogroups.close();
		//runAll.close();
		script.close();
	}
}
