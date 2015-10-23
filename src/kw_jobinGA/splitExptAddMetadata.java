/**
 * split the log normalized results into gastric aspirate vs stool
 * add the relevant metadata
 * just focus on run2 read1
 * @author kwinglee
 */

package kw_jobinGA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;

import utils.ConfigReader;

public class splitExptAddMetadata {

	public static void main(String[] args) throws FileNotFoundException, Exception {
		//get stool sample metadata
		HashMap<String, String> stool = new HashMap<String, String>();//hashmap of stool id to check group
		BufferedReader metadata = new BufferedReader(new InputStreamReader (new FileInputStream(new File(
				ConfigReader.getJobinGAStoolDir() + File.separator + "Stool Microbiome for Miseq no-DOB 10-5-2015.txt"))));
		for(int i = 0; i < 3; i++) {
			metadata.readLine(); //read extra lines
		}
		String line = metadata.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String id = sp[7];
			if(line.contains("no check")) {
				stool.put(id, "no_check");
			} else {
				stool.put(id, "check");
			}
			line = metadata.readLine();
		}
		metadata.close();
		/*for(String key: stool.keySet()) {
			System.out.println(key + "\t" + stool.get(key));
		}*/		
		
		//get gastric aspirate metadata
		HashMap<String, String> ga = new HashMap<String, String>();//hashmap of ga id to analysis group
		metadata = new BufferedReader(new InputStreamReader (new FileInputStream(new File(
				ConfigReader.getJobinGAStoolDir() + File.separator + "Gast.Asp. Microbiome for Miseq 10-5-2015.txt"))));
		for(int i = 0; i < 2; i++) {
			metadata.readLine(); //read extra header lines
		}
		line = metadata.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			if(sp.length == 8) {//skip blank lines; end up with extra blank key from the extra groups listed at the end
				ga.put(sp[6], sp[7]);
			}
			line = metadata.readLine();
		}
		metadata.close();
		/*for(String key: ga.keySet()) {
			System.out.println(key + "\t" + ga.get(key));
		}*/
		
		//for each taxa level, split into ga or stool, only write run2 read1, add metadata
		String[] taxaLevels = {"phylum", "class", "order", "family", "genus"};
		for(int i = 0; i < taxaLevels.length; i++) {
			String taxa = taxaLevels[i];
			BufferedReader table = new BufferedReader(new InputStreamReader (new FileInputStream(new File(
					ConfigReader.getJobinGAStoolRDPDir() + File.separator + "rdp_taxaAsCol_logNorm_" + taxa + ".txt"))));
			
			//set up output files
			BufferedWriter stoolOut = new BufferedWriter(new FileWriter (new File(
					ConfigReader.getJobinGAStoolRDPDir() + File.separator + "stool" + File.separator + "rdp_stool_" + taxa + ".txt")));
			BufferedWriter gaOut = new BufferedWriter(new FileWriter (new File(
					ConfigReader.getJobinGAStoolRDPDir() + File.separator + "gastricAspirate" + File.separator + "rdp_ga_" + taxa + ".txt")));
			
			//set up headers
			line = table.readLine();
			String[] sp = line.split("\t");
			stoolOut.write("sampleID\tgroup");
			gaOut.write("sampleID\tgroup");
			for(int j = 1; j < sp.length; j++) {
				stoolOut.write("\t" + sp[j]);
				gaOut.write("\t" + sp[j]);
			}
			stoolOut.write("\n");
			gaOut.write("\n");
			
			//get sample counts, write if run2 read1, add metadata
			line = table.readLine();
			while(line != null) {
				if(line.startsWith("Run2_R1")) {
					sp = line.split("\t");
					String id = sp[0].replace("Run2_R1_", "").replace(".fasta", "");
					if(id.startsWith("S")) {
						stoolOut.write(id + "\t" + stool.get(id));
						for(int j = 1; j < sp.length; j++) {
							stoolOut.write("\t" + sp[j]);
						}
						stoolOut.write("\n");
					} else if(id.startsWith("G")) {
						gaOut.write(id + "\t" + ga.get(id));
						for(int j = 1; j < sp.length; j++) {
							gaOut.write("\t" + sp[j]);
						}
						gaOut.write("\n");
					} /*else {
						System.out.println(id);
					}*/
				}
				line = table.readLine();
			}
			
			//close files
			table.close();
			stoolOut.close();
			gaOut.close();
		}
	}
}
