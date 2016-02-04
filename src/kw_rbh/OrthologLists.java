/**
 * from the ortholog tables, generate list of orthogroups, basing membership of set intersection
 * also generate fasta file for each orthogroup
 * very slow to run
 * Speed up improvements over orthologLists: check set size before doing analysis
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class OrthologLists {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/";
	public static int MIN = 10;//minimum numbers of members in set to keep orthogroup

	public static void main(String[] args) throws IOException {
		BufferedWriter log = new BufferedWriter(new FileWriter(new File("/nobackup/afodor_research/kwinglee/cre/rbh/orthologListLog")));//log to track progress
		File[] tables = new File(DIR + "orthologTables").listFiles();
		List<Set<String>> intersect = new ArrayList<Set<String>>();//intersection of everything
		List<Set<String>> union = new ArrayList<Set<String>>();//once an intersection is found, union is everything in that set so that don't later add sets that were already removed
		for(File t : tables) {
			if(t.getName().startsWith("orthologNameTable_")) {
				log.write(t.getName() + "\n");
				log.flush();
				List<Set<String>> tab = tableToSet(t);
				//add sets to current lists
				for(int i = 0; i < tab.size(); i++) {
					if(tab.get(i).size() > MIN) {
						boolean matched = false; //if true, this has been matched to what is in list
						for(int j = 0; j < union.size(); j++) {
							Set<String> temp = new HashSet<String>();
							temp.addAll(union.get(j));
							temp.retainAll(tab.get(i));//intersection
							if(!temp.isEmpty()) {//if intersection > 0, add to union and update intersection
								union.get(j).addAll(tab.get(i));
								intersect.get(j).retainAll(tab.get(i));
								
								//check if previously had union
								if(matched) {
									System.err.println("overlapping orthogroups: " + t.getName() + "\t" + i + "\t" + j);
								}
								matched = true;
							}
						}
						if(!matched) {//set not previously seen, add to lists
							union.add(tab.get(i));
							intersect.add(tab.get(i));
						}
					}
				}
			}
		}
		log.write("writing...\n");
		log.flush();
		//write lists and fastas for orthogroups
		String fastaDir = DIR + "orthologListFastas/";
		new File(fastaDir).mkdirs();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "orthogroupList.txt")));
		out.write("orthogroupNumber\tnumberOfGenes\tgeneList\n");//header
		for(int i = 0; i < intersect.size(); i++) {
			Set<String> set = intersect.get(i);
			if(!set.isEmpty() && set.size() > MIN) {
				//write list
				out.write("Orthogroup" + i + "\t" + set.size() + "\t");
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
	 * for the given table in f, return a list of the sets for each gene
	 */
	private static List<Set<String>> tableToSet(File f) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = br.readLine();
		//initialize the sets with the genes
		List<Set<String>> sets = new ArrayList<Set<String>>();
		String[] sp = line.split("\t");
		for(int i = 1; i < sp.length; i++) {//skip first element, aka the genome ids
			Set<String> s = new HashSet<String>();
			s.add(sp[i]);
			sets.add(s);
		}
		line = br.readLine();
		//add genes to set
		while(line != null) {
			sp = line.split("\t");
			for(int i = 1; i < sp.length; i++) {
				sets.get(i-1).add(sp[i]);
			}
			line = br.readLine();
		}
		br.close();
		//remove NA
		for(int i = 0; i < sets.size(); i++) {
			sets.get(i).remove("NA");
		}
		
		return(sets);
	}
}
