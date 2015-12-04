/**
 * code to make list of all taxa seen in RDP files in one place
 * output is table of each taxa from RDP database to be used in writing out taxonomies
 * 12/4/15
 */
package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class getRDPtaxonomy {
	
	public static String getTaxonomy(String[] sp) {
		String phy = "";
		String cl = "";
		String ord = "";
		String fam = "";
		String gen = "";
		for(int i = 0; i < sp.length; i++) {
			if(sp[i].equals("phylum")) {
				phy = sp[i-1];
			} else if(sp[i].equals("class")) {
				cl = sp[i-1];
			} else if(sp[i].equals("order")) {
				ord = sp[i-1];
			} else if(sp[i].equals("family")) {
				fam = sp[i-1];
			} else if(sp[i].equals("genus")) {
				gen = sp[i-1];
			}
		}
		return(phy + "\t" + cl + "\t" + ord + "\t" + fam + "\t" + gen);
	}
	
	public static void main(String[] args) throws IOException {
		//folders with RDP results
		String[] folders = {"/projects/afodor_research/kwinglee/jobin/anaerobe/rdpResults",
				"/projects/afodor_research/kwinglee/jobin/biofilm/rdpResults", 
				"/projects/afodor_research/kwinglee/jobin/ga-stool/rdpResults"};
		//files to analyze
		HashSet<File> files = new HashSet<File>();
		//files.add(new File("C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\China\\abundantOTU.chinaForward.taxonomy.txt"));
		files.add(new File("/projects/afodor_research/kwinglee/china/abundantOTU.chinaForward.taxonomy.txt"));
		for(String fname: folders) {
			File f = new File(fname);
			File[] flist = f.listFiles();
			for(File file : flist) {
				if(file.getName().startsWith("rdp_")) {
					files.add(file);
				}
			}
		}
		
		//analyze files
		HashMap<String, String> taxonomy = new HashMap<String, String>();//map of genus to taxonomy
		for(File file : files) {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line != null) {
				String[] sp = line.split("\t");
				int s = sp.length - 3;
				if(!sp[s+1].equals("genus")) {//check this is genus
					System.out.println("Need to update s: " + file.getName() + "\t" + sp[s+1] + "\t" + line);
				} else {
					String gen = sp[s];
					String tax = getTaxonomy(sp);
					if(taxonomy.containsKey(gen)) {
						//check taxonomies are same
						String tax2 = taxonomy.get(gen);
						if(!tax.equals(tax2)) {
							System.out.println("MIXED TAXONOMIES: ");
							System.out.println(tax);
							System.out.println(tax2 + "\n");
						}
					} else {
						taxonomy.put(gen,  tax);
					}
				}
			}
			br.close();
		}
		
		//write results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/china/RDPtaxonomy.txt")));
		Set<String> gen = taxonomy.keySet();
		for(String g : gen) {
			System.out.println(taxonomy.get(g));
		}
		out.close();
	}

}
