/*
 * Generates a table of genome name to class (carolina, resistant, susceptible)
 */
package kw_rbh;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class GenomeToClass {
	public static String GEN_DIR = "/nobackup/afodor_research/af_broad/";
	public static String OUT_DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";
	
	public static void main(String[] args) throws IOException {
		String[] classes = {"carolina", "resistant", "susceptible"};
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				OUT_DIR + "GenomeToClass.txt")));
		out.write("Genome\tClass\n");
		for(String c : classes) {
			String[] files = new File(GEN_DIR + c).list();
			for(String f : files) {
				if(f.endsWith("genes.gtf")) {
					out.write(f.replace(".genes.gtf", "") +  "\t" + c  + "\n");
				}
			}
		}
		
		out.close();
	}

}
