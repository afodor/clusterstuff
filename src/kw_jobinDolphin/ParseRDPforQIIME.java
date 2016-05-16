/*
 * parse RDP classification of abundant OTU results
 * for QIIME blast_fragments chimera detector
 * and in format like metaphlan mpa output
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
		BufferedWriter outmpa = new BufferedWriter(new FileWriter(new File(
				DIR + "dolphinAbundantOTU.cons.rdpTaxonomy.mpa")));
		for(String line = in.readLine(); line != null; line=in.readLine()) {
			String[] sp = line.split("\t");
			out.write(sp[0] + "\t");//the id
			String taxa = "";
			String mpa = "";
			//domain
			taxa += getTaxonomyName(sp, "domain");
			mpa += getTaxonomyNameMpa(sp, "domain");
			
			//phylum
			taxa += getTaxonomyName(sp, "phylum");
			mpa += getTaxonomyNameMpa(sp, "phylum");
			
			//class
			taxa += getTaxonomyName(sp, "class");
			mpa += getTaxonomyNameMpa(sp, "class");
			
			//order
			taxa += getTaxonomyName(sp, "order");
			mpa += getTaxonomyNameMpa(sp, "order");
			
			//family
			taxa += getTaxonomyName(sp, "family");
			mpa += getTaxonomyNameMpa(sp, "family");
			
			//genus
			taxa += getTaxonomyName(sp, "genus").replace(";", "");
			mpa += getTaxonomyNameMpa(sp, "genus").replace("\\|", "");
			
			out.write(taxa + "\n");
			outmpa.write(mpa + "\n");
		}
		
		in.close();
		out.close();
		outmpa.close();
	}

	public static String getTaxonomyName(String[] splits, String level) {
		int index = 6;
		while(index < splits.length && !splits[index].equals(level)) {
			index += 3;
		}
		if(index >= splits.length) {
			System.out.println("Missing " + level + ":");
			for(int i = 0; i < splits.length; i++) {
				System.out.print(splits[i] + "\t");
			}
			System.out.println();
			return(";");
		} else if(Double.parseDouble(splits[index+1]) > 0.8) {
			return(splits[index-1].replace("\"", "") + ";");
		} else {
			return(";");
		}
	}
	
	public static String getTaxonomyNameMpa(String[] splits, String level) {
		int index = 6;
		String prefix = level.charAt(0) + "__";
		while(index < splits.length && !splits[index].equals(level)) {
			index += 3;
		}
		if(index >= splits.length) {
			System.out.println();
			return(prefix + "|");
		} else if(Double.parseDouble(splits[index+1]) > 0.8) {
			return(prefix + splits[index-1].replace("\"", "") + "|");
		} else {
			return(prefix + "|");
		}
	}
}
