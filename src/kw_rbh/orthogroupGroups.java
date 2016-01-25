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

public class orthogroupGroups {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/rbhOrthologs/";

	public static void main(String[] args) throws IOException {
		BufferedWriter log = new BufferedWriter(new FileWriter(new File("/nobackup/afodor_research/kwinglee/cre/rbh/orthologGroupLog")));//log to track progress
		File[] tables = new File(DIR + "orthologTables").listFiles();
		
		//get hash of genome to hash of gene in that genome to orthologs
		HashMap<String, HashMap<String, Set<String>>> genomeToOrth = new HashMap<String, HashMap<String, Set<String>>>();
		for(File t : tables) {
			if(t.getName().startsWith("orthologNameTable_")) {
				log.write("Getting map for " + t.getName());
				log.flush();
				String name = t.getName().replace("orthologNameTable_", "").replace(".txt", "");
				HashMap<String, Set<String>> map = tableToHash(t);
				genomeToOrth.put(name, map);
				break;
			}
		}
		
		//for each table, look at each set: for each member, get intersection of that member with its set in other tables
		
		//merge
		
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
		for(String k : map.keySet()) {
			Set<String> set = map.get(k);
			System.out.println(k + ":");
			for(String s : set) {
				System.out.print(s + "; ");
			}
			System.out.println("\n");
		}
		
		return(map);
	}
}
