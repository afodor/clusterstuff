/**
 * from the ortholog tables, generate list of orthogroups, basing membership of set intersection
 * also generate fasta file for each orthogroup
 * alternative to orthologLists
 */

package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class orthologGroups {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/";
	public static int MIN = 5;//minimum numbers of members in set to keep orthogroup

	public static void main(String[] args) throws IOException {
		BufferedWriter log = new BufferedWriter(new FileWriter(new File("/nobackup/afodor_research/kwinglee/cre/rbh/orthologGroupLog")));//log to track progress
		File[] tables = new File(DIR + "orthologTables").listFiles();
		
		//get hash of genome to hash of gene in that genome to orthologs
		HashMap<String, HashMap<String, Set<String>>> genomeToOrth = new HashMap<String, HashMap<String, Set<String>>>();
		for(File t : tables) {
			if(t.getName().startsWith("orthologNameTable_")) {
				log.write("Getting map for " + t.getName() + "\n");
				log.flush();
				String name = t.getName().replace("orthologNameTable_", "").replace(".txt", "");
				HashMap<String, Set<String>> map = tableToHash(t);
				genomeToOrth.put(name, map);
			}
		}
		
		//for each table, look at each set: for each member, get intersection of that member with its set in other tables
		String intDir = DIR + "orthogroupTablesByGenome/";//where to write these intermediate results
		new File(intDir).mkdirs();
		List<List<Set<String>>> orthogroupLists = new ArrayList<List<Set<String>>>();//list of list of orthogroups
		for(String genome : genomeToOrth.keySet()) {
			log.write("Looking at sets for " + genome + "\n");
			log.flush();
			HashMap<String, Set<String>> geneMap = genomeToOrth.get(genome);//copy so don't change table for next round
			List<Set<String>> orthogroups = new ArrayList<Set<String>>();//list of orthogroups
			
			//for each set of orthologs, for each member get intersection with that member's set then add to orthogroup list
			for(String gene : geneMap.keySet()) {
				Set<String> geneSet = geneMap.get(gene);
				if(!geneSet.isEmpty()) {
					Set<String> orth = new HashSet<String>();//the orthogroup (make a copy so don't modify the original table)
					orth.addAll(geneSet);				
					//get intersection of each member
					for(String gene2 : geneSet) {
						//String genome2 = gene2.replaceFirst("_[A-Z]*[0-9]*_[0-9]*$", "");//works on carolina, haven't checked other groups
						String[] sp = gene2.split("_");
						String genome2 = "";
						for(int i = 0; i < sp.length-2; i++) {
							genome2 += sp[i] + "_";
						}
						genome2 = genome2.replaceFirst("_$", "");
						HashMap<String, Set<String>> map2 = genomeToOrth.get(genome2);
						Set<String> set2 = map2.get(gene2);
						orth.retainAll(set2);
					}
					if(!orth.isEmpty() && orth.size() > MIN) {
						orthogroups.add(orth);
					}
				}
			}
			log.write("  number orthogroups: " + orthogroups.size() + "\n");
			log.flush();
			orthogroupLists.add(orthogroups);
			
			//write intermediate tables
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(intDir + "orthogroups_" +  genome + ".txt")));
			for(Set<String> set : orthogroups) {
				for(String gene : set) {
					out.write(gene + ";");
				}
				out.write("\n");
			}
			out.close();			
		}
		
		//try to clear memory
		genomeToOrth = null;
		//System.gc();
		
		//merge
		List<Set<String>> orthogroups = new ArrayList<Set<String>>();
		int numMerge = 1;
		for(List<Set<String>> list : orthogroupLists) {
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
					String genome = s.replaceFirst("_.*_.*$", "");//genome name is everything except the last two strings when separating on _
					
					BufferedReader br = new BufferedReader(new FileReader(new File(
							"/nobackup/afodor_research/kwinglee/cre/rbh/geneFastas/" + 
							genome + "/" + s + ".fasta")));
					String line = br.readLine();
					while(line != null) {
						fasta.write(line + "\n");
						line = br.readLine();
					}
					br.close();
				}
				
				fasta.close();
			}
		}
		out.close();
		
		log.close();
	}
	
	/*
	 * returns map of gene to set of genes that are an ortholog for the given table
	 */
	public static HashMap<String, Set<String>> tableToHash(File table) throws IOException {
		HashMap<String, Set<String>> map = new HashMap<String, Set<String>>();
		
		//read file
		BufferedReader br = new BufferedReader(new FileReader(table));
		String line = br.readLine();
		
		//initialize the sets with the genes
		List<Set<String>> sets = new ArrayList<Set<String>>();
		String[] genes = line.split("\t");
		for(int i = 1; i < genes.length; i++) {//skip first element, aka the genome ids
			Set<String> s = new HashSet<String>();
			s.add(genes[i]);
			sets.add(s);
		}
		line = br.readLine();
		
		//add genes to set
		while(line != null) {
			String[] sp = line.split("\t");
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
		
		//put in map
		for(int i = 1; i < genes.length; i++) {
			map.put(genes[i], sets.get(i-1));
		}
		
		//write results to system.out (to check)
		/*for(String k : map.keySet()) {
			Set<String> set = map.get(k);
			System.out.println(k + ":");
			for(String s : set) {
				System.out.print(s + "; ");
			}
			System.out.println("\n");
		}*/
		
		return(map);
	}
}
