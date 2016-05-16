/*
 * parse RDP classification of abundant OTU results
 * for QIIME blast_fragments chimera detector
 */
package kw_jobinDolphin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ParseRDPforQIIME {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/abunOTU/";
	
	public static void main(String[] args) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(new File(
				DIR + "dolphinAbundantOTU.cons.rdpTaxonomy")));
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "dolphinAbundantOTU.cons.rdpTaxonomy.forQiime")));
		for(String line = in.readLine(); line != null; line=in.readLine()) {
			String[] sp = line.split("\t");
			out.write(sp[0] + "\t");//the id
			String taxa = "";
			int index = 5;
			//domain
			while(!sp[index].equals("domain")) {
				index += 3;
			}
			if(Double.parseDouble(sp[index+1]) > 0.8) {
				taxa += sp[index-1].replace("\"", "") + ";";
			} else {
				taxa += ";";
			}
			index += 3;
			
			//phylum
			while(!sp[index].equals("phylum")) {
				index += 3;
			}
			if(Double.parseDouble(sp[index+1]) > 0.8) {
				taxa += sp[index-1].replace("\"", "") + ";";
			} else {
				taxa += ";";
			}
			index += 3;
			
			//class
			while(!sp[index].equals("class")) {
				index += 3;
			}
			if(Double.parseDouble(sp[index+1]) > 0.8) {
				taxa += sp[index-1].replace("\"", "") + ";";
			} else {
				taxa += ";";
			}
			index += 3;
			
			//order
			while(!sp[index].equals("order")) {
				index += 3;
			}
			if(Double.parseDouble(sp[index+1]) > 0.8) {
				taxa += sp[index-1].replace("\"", "") + ";";
			} else {
				taxa += ";";
			}
			index += 3;
			
			//family
			while(!sp[index].equals("family")) {
				index += 3;
			}
			if(Double.parseDouble(sp[index+1]) > 0.8) {
				taxa += sp[index-1].replace("\"", "") + ";";
			} else {
				taxa += ";";
			}
			index += 3;
			
			//genus
			while(!sp[index].equals("genus")) {
				index += 3;
			}
			if(Double.parseDouble(sp[index+1]) > 0.8) {
				taxa += sp[index-1].replace("\"", "");
			} 
			index += 3;
			
			out.write(taxa + "\n");
		}
		
		in.close();
		out.close();
	}

}
