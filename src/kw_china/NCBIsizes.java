/*
 * get a table of genome to number of genes and total gene length
 */
package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class NCBIsizes {
	public static void main(String[] args) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/kwinglee/china/ncbiSizes.txt")));
		File ncbiDir = new File("/nobackup/afodor_research/ncbi");
		
		out.write("genome\tnumberGenes\ttotalGeneLength\n");
		File[] flist1 = ncbiDir.listFiles();
		for(File f1 : flist1) {
			if(f1.isDirectory()) {
				File[] flist2 = f1.listFiles();
				for(File f2 : flist2) {
					int numGenes = 0;
					int totGeneLength = 0;
					if(f2.getName().endsWith(".frn")) {
						BufferedReader fa = new BufferedReader(new FileReader(f2));
						for(String line = fa.readLine(); line != null; line = fa.readLine()) {
							if(line.startsWith(">")) {
								numGenes++;
							} else {
								totGeneLength += line.length();
							}
						}
						fa.close();
					} else {
						System.out.println("Skipped file: " + f2.getAbsolutePath());
					}
					out.write(f1.getName() + "\t" + numGenes + 
							"\t" + totGeneLength + "\n");
				}
			}
		}
		
		out.close();
	}

}
