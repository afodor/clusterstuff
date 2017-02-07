/*
 * Convert U to T in Silva sequences
 */
package kw_jobinFrodentium;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ConvertToDNA {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/F_rodentium/";
	
	public static void main(String[] args) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + "arb-silva.de_2017-02-07_id405408_tax_silva.fasta")));
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "Faecalibaculum_rodentium_SilvaSSU_DNA.fasta")));
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			if(line.startsWith(">")) {
				out.write(line + "\n");
			} else {
				out.write(line.replaceAll("U", "T").replaceAll("u", "t") + "\n");
			}
		}
		br.close();
		out.close();
	}

}
