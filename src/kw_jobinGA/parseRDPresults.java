/**
 * Parse hierarchical file output from RDP into separate files for each taxonomic level
 * Add metadata on samples; remove the taxid, lineage, and rank columns
 * @author kwinglee
 * @date 10/22/15
 */

package kw_jobinGA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class parseRDPresults {
	public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\JobinCollaboration\\GA-stools\\rdpResults\\";
	
	public static void main(String[] args) throws IOException {
		//get hierarchical file results
		BufferedReader hier = new BufferedReader(new FileReader(new File(DIR + "hier_merge.txt")));
		
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
			//write sample id
			file.write("\nsampleID");
			for(int i = 4; i < sp.length; i++) {
				String[] name = sp[i].split("_");
				file.write("\t" + name[2].replace(".fasta", ""));
			}
			//write run
			file.write("\nrun");
			for(int i = 4; i < sp.length; i++) {
				String[] name = sp[i].split("_");
				file.write("\t" + name[0].replace("Run", ""));
			}
			//write read
			file.write("\nread");
			for(int i = 4; i < sp.length; i++) {
				String[] name = sp[i].split("_");
				file.write("\t" + name[1].replace("R", ""));
			}
			//write sample type
			file.write("\nsampleType");
			for(int i = 4; i < sp.length; i++) {
				String[] name = sp[i].split("_");
				if(name[2].startsWith("G")) {
					file.write("\tgastric_aspirate");
				} else if(name[2].startsWith("S")) {
					file.write("\tstool");
				} else if(name[2].startsWith("H20") || name[2].startsWith("NC101") ||
						name[2].startsWith("NC101") || name[2].startsWith("neg")) {
					file.write("\tcontrol");
				} else {
					System.out.println("BAD NAME: " + name[2]);
				}
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
			} /*else if(rank.length()==0) {//blank rank; the taxa is unclassified but write based on the lineage (ex if unclassified family, write as genus)
				String lin = sp[1];//lineage
				if(lin.contains("family")) {
					writeLine(genus, sp);
				} else if(lin.contains("order")) {//lines that are unclassified but end in suborder are included on the total order count and so should be added to family
					writeLine(fam, sp);
				} else if(lin.contains(";class;")) {
					writeLine(order, sp);
				} else if(lin.contains("phylum")) {
					writeLine(cl, sp);
				} else {
					writeLine(phy, sp);
				}
			} else {//skip subclass and suborder; these are included in the class and order counts
				System.out.println("Skipped rank: " + rank);
			}(*/ /////currently ignoring unclassified sequences
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
