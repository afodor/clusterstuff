/**
 * merge lists output by ex, the intermediate steps of orthologGroups into one list
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class MergeOrthologLists {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/";
	public static int MIN = 10;//minimum numbers of members in set to keep orthogroup

	public static void main(String[] args) throws IOException {
		BufferedWriter log = new BufferedWriter(new FileWriter(new File("/nobackup/afodor_research/kwinglee/cre/rbh/orthologMergeLog")));//log to track progress
		File[] tables = new File(DIR + "orthogroupTablesByGenome").listFiles();

		//merge
		List<Set<String>> orthogroups = new ArrayList<Set<String>>();
		int numMerge = 1;
		for(File f : tables) {
			List<Set<String>> list = fileToGeneSets(f);
			log.write("merging orthogroups " + numMerge + "\n");
			numMerge++;
			log.flush();

			for(Set<String> set : list) {
				boolean merged = false;
				for(Set<String> orth : orthogroups) {
					Set<String> copy = new HashSet<String>();
					copy.addAll(orth);
					copy.retainAll(set);//see if sets intersect
					if(!copy.isEmpty()) {//sets intersect
						orth.retainAll(set);
						if(merged) {
							System.out.println("set merged twice");
						}
						merged = true;
					}
				}
				if(!merged) {//not already in orthogroup list
					orthogroups.add(set);
				}
			}
			log.write("   number orthogroups = " + orthogroups.size() + "\n");
			log.flush();
		}

		//write orthogroups and get fasta files
		log.write("writing...\n");
		log.flush();
		String fastaDir = DIR + "orthologGroupFastas/";
		new File(fastaDir).mkdirs();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "orthologGroups.txt")));
		out.write("orthogroupNumber\tnumberOfGenes\tgeneList\n");//header
		for(int i = 0; i < orthogroups.size(); i++) {
			Set<String> set = orthogroups.get(i);
			if(!set.isEmpty() && set.size() > MIN) {
				//write list
				out.write("orthogroup" + i + "\t" + set.size() + "\t");
				for(String s : set) {
					out.write(s + ";");
				}
				out.write("\n");

				//write fastas
				BufferedWriter fasta = new BufferedWriter(new FileWriter(new File(
						fastaDir + "orthogroup" + i + ".fasta")));
				for(String s : set) {
					//String genome = s.replaceFirst("_[A-Z]*[0-9]*_[0-9]*$", "");//works on carolina, haven't checked other groups
					String[] sp = s.split("_");
					String genome = "";
					for(int j = 0; j < sp.length-2; j++) {
						genome += sp[j] + "_";
					}
					genome = genome.replaceFirst("_$", "");

					try {
						BufferedReader br = new BufferedReader(new FileReader(new File(
								"/nobackup/afodor_research/kwinglee/cre/rbh/geneFastas/" + 
										genome + "/" + s + ".fasta")));
						String line = br.readLine();
						while(line != null) {
							fasta.write(line + "\n");
							line = br.readLine();
						}
						br.close();
					} catch(Exception e) {
						System.err.println(e.getMessage());
					}
				}

				fasta.close();
			}
		}
		out.close();

		log.close();
	}
	
	/*
	 * for the given file, where each row is a set of genes in an orthogroup, 
	 * convert to a list of set of strings
	 */
	private static List<Set<String>> fileToGeneSets(File file) throws IOException {
		List<Set<String>> list = new ArrayList<Set<String>>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = br.readLine();
		while(line != null) {
			String[] sp = line.split(";");
			Set<String> set = new HashSet<String>(Arrays.asList(sp));
			list.add(set);
			line = br.readLine();
		}
		br.close();
		return(list);
	}
}
