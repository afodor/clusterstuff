/**
 * For each Broad genome, generate a fasta file of all genes from the gtf file
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

public class geneFastas {
	public static String GenomeTopDir = "/nobackup/afodor_research/af_broad/";

	public static void main(String[] args) throws IOException {
		analyzeFolder("carolina");
		/*analyzeFolder("susceptible");
		analyzeFolder("resistant");*/
	}
	
	//analyze the genomes in the given folder
	public static void analyzeFolder(String folder) throws IOException {
		String outDir = "/nobackup/afodor_research/kwinglee/cre/rbh/" + folder + "_";
		String dirName = GenomeTopDir + folder;
		File dir = new File(dirName);
		File[] allFiles = dir.listFiles();
		for(File gtf : allFiles) {
			if(gtf.getName().endsWith(".genes.gtf")) {
				String name = gtf.getName().replace(".genes.gtf", "");//genome to analyze
				//make directory for output files
				File gtfDir = new File(outDir + name);
				gtfDir.mkdirs();
				
				//get corresponding fasta file as hash map of scaffold name to sequence
				HashMap<String, String> scaff = getScaffolds(new File(dirName + "/" + name + ".scaffolds.fasta"));
				Iterator<String> itSet = scaff.keySet().iterator();
				while(itSet.hasNext()) {
					String p = itSet.next();
					System.out.println(scaff.get(p) + "\t" + p);
				}
				break;
				//make gene files
			}
		}
	}
	
	//takes the given fasta file, reads it and returns a hash map of scaffold name to sequence
	public static HashMap<String, String> getScaffolds(File fasta) throws IOException {
		HashMap<String, String> map = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fasta)));
		String line = br.readLine();
		String header = line.replaceFirst(">", "");
		String seq = "";
		while(line != null) {
			if(line.startsWith(">")) { //header/start of new sequence
				if(seq.length() > 1) { //add to map unless its the beginning of the file
					map.put(header, seq);
				}
				header = line.replaceFirst(">", "");
				seq = "";
			} else {
				seq += line;
			}
		}
		br.close();
		return(map);
	}
}
