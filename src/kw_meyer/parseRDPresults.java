/**
 * Parse hierarchical file output from RDP into separate files for each taxonomic level
 * remove the taxid, lineage, and rank columns
 */

package kw_meyer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class parseRDPresults {
	//public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\Meyer\\";
	public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\Meyer\\RDP analysis trimmed\\";
	
	public static void main(String[] args) throws IOException {
		//get hierarchical file results
		//BufferedReader hier = new BufferedReader(new FileReader(new File(DIR + "hierarch_merge_all.txt")));
		BufferedReader hier = new BufferedReader(new FileReader(new File(DIR + "hierarch_merge_all_trimmed.txt")));
		
		//set up output files
		BufferedWriter genus = new BufferedWriter(new FileWriter(new File(DIR + "rdp_genus.txt")));
		BufferedWriter fam = new BufferedWriter(new FileWriter(new File(DIR + "rdp_family.txt")));
		BufferedWriter order = new BufferedWriter(new FileWriter(new File(DIR + "rdp_order.txt")));
		BufferedWriter cl = new BufferedWriter(new FileWriter(new File(DIR + "rdp_class.txt")));
		BufferedWriter phy = new BufferedWriter(new FileWriter(new File(DIR + "rdp_phylum.txt")));
		HashMap<String, BufferedWriter> out = new HashMap<String, BufferedWriter>();//map of level to output file
		out.put("genus", genus);
		out.put("family", fam);
		out.put("order", order);
		out.put("class", cl);
		out.put("phylum", phy);
		
		//set up header and metadata in each file
		//rows: file name, sample ID, run, read, sample type
		String line = hier.readLine();
		String[] sp = line.split("\t");
		for(String key: out.keySet()) {
			BufferedWriter file = out.get(key);
			//write file names
			file.write("fileName");
			for(int i = 4; i < sp.length; i++) {
				file.write("\t" + sp[i]);
			}
			file.write("\n");
		}
		
		//split counts into separate files based on rank (and remove taxid, lineage and rank columns)
		line = hier.readLine();
		while(line != null) {
			sp = line.split("\t");
			String rank = sp[3];
			if(out.containsKey(rank)) {
				writeLine(out.get(rank), sp);				
			} 
			line = hier.readLine();
		}
		
		
		//close files
		hier.close();
		genus.close();
		fam.close();
		order.close();
		cl.close();
		phy.close();
	}
	
	/**
	 * Writes the given line of counts (split by tabs) to the given file
	 * @param file	output file
	 * @param sp	line of counts for a taxa split by tabs
	 * @throws IOException
	 */
	public static void writeLine(BufferedWriter file, String[] sp) throws IOException {
		file.write(sp[2]);//name
		//write counts
		for(int i=4; i < sp.length; i++) {
			file.write("\t" + sp[i]);
		}
		file.write("\n");
	}

}
