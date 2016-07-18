/*
 * Uses the results from BWAbetaLactamaseScripts (idxstats and depth) to get
 * the proportion of reads that mapped to each reference, and the average depth.
 * Also combines individual SRR file results into one number for each isolate.
 */
package kw_cre;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class BWAbetaLactamaseStats {
	//public static final String DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/bwaAlignToBetaLactamases/";
	public static final String DIR = "/nobackup/test_afodor_research/kwinglee/cre/chs_v_cards/bwaAlignToBetaLactamases/";
	public static final String CONVERT = "/nobackup/afodor_research/mjzapata/CRE/CHS_raw/chs_batch_download_results.csv";//file used to convert to from SRR to CHS
	//public static final String REF = "/users/kwinglee/card/beta_lactamase.protein_homolog.fasta";//reference used for alignment
	public static final String REF = "/users/kwinglee/card/beta_lactamase2.protein_homolog.fasta";//reference used for alignment
	
	public static void main(String[] args) throws Exception {
		//get SRR to CHS conversion (map of CHS number to list of corresponding SRR numbers
		//note: SRR numbers are all duplicated in the conversion file
		HashMap<String, HashSet<String>> chs = new HashMap<String, HashSet<String>>();
		BufferedReader conv = new BufferedReader(new FileReader(new File(CONVERT)));
		ArrayList<Integer> strains = new ArrayList<Integer>();
		for(String line = conv.readLine(); line != null; line = conv.readLine()) {
			String[] sp = line.replace("[", "").replace("]", "").split("\t");
			HashSet<String> set = new HashSet<String>(Arrays.asList(sp[1].split(", ")));
			chs.put(sp[0], set);
			strains.add(Integer.parseInt(sp[0]));
		}
		conv.close();
		Collections.sort(strains);
		
		//list of references and their lengths
		HashMap<String, Integer> refLengths = new HashMap<String, Integer>();
		BufferedReader brRef = new BufferedReader(new FileReader(new File(REF)));
		String name = "";
		int len = 0;
		for(String line = brRef.readLine(); line != null; line = brRef.readLine()) {
			if(line.startsWith(">")) {
				if(name.length() > 0) {
					refLengths.put(name, len);
				}
				name = line.replace(">", "").split(" ")[0];
				len = 0;
			} else {
				len += line.length();
			}
		}
		refLengths.put(name, len);
		brRef.close();
		List<String> refs = new ArrayList<String>();
		refs.addAll(refLengths.keySet());
		Collections.sort(refs);
		
		//list of files
		String[] results = new File(DIR).list();
		
		//sum of depths (from depth command)
		HashMap<String, HashMap<String, Integer>> sumDepth = new HashMap<String, HashMap<String, Integer>>();
		//map of SRR to map of reference to sum of number reads that aligned at each position in that reference
		for(String file : results) {
			if(file.endsWith("2.depth.txt")) {
				String srr = file.split("\\.")[0];
				HashMap<String, Integer> map = new HashMap<String, Integer>();
				BufferedReader br = new BufferedReader(new FileReader(new File(
						DIR + file)));
				String line = br.readLine();
				if(line == null || line.length() < 1) {//file is empty
					br.close();
					throw new Exception("File " + file + "is empty");
				}
				while(line != null) {
					String[] sp = line.split("\t");
					String key = sp[0];
					int count = Integer.parseInt(sp[2]);
					if(!map.containsKey(key)) {
						map.put(key, count);
					} else {
						map.put(key, map.get(key) + count);
					}
					line = br.readLine();
				}
				br.close();
				//if no reads mapped reference, ref is not in depth output file,
				//so add missing refs in
				if(refs.size() > map.size()) {
					Set<String> keys = map.keySet();
					Set<String> refset = new HashSet<String>();
					refset.addAll(refs);
					refset.removeAll(keys);
					for(String r : refset) {
						map.put(r, 0);
					}
				} else if(map.size() > refs.size()) {
					System.out.println("File " + file + 
							" has incorrect number of refs: " + map.size());
				}
				sumDepth.put(srr, map);
			}
		}
		
		//number of reads mapped (from idxstats) and total number reads
		HashMap<String, HashMap<String, Integer>> numMapped = new HashMap<String, HashMap<String, Integer>>();
		     //map of SRR to map of reference to number of reads that mapped that reference
		HashMap<String, Integer> numReads = new HashMap<String, Integer>();
		    //map of SRR to total number of reads
		for(String file : results){
			if(file.endsWith("2.idxstats.txt")) {
				String srr = file.split("\\.")[0];
				HashMap<String, Integer> map = new HashMap<String, Integer>();
				BufferedReader br = new BufferedReader(new FileReader(new File(
						DIR + file)));
				String line = br.readLine();
				if(line == null || line.length() < 1 || line.startsWith("*")) {//file is empty or no reads mapped
					br.close();
					throw new Exception("File " + file + "is empty");
				}
				int totReads = 0;
				while(line != null) {
					String[] sp = line.split("\t");
					String key = sp[0];
					//check lengths
					if(!refLengths.containsKey(key) && !key.equals("*")) {
						br.close();
						throw new Exception("Extra key " + key);
					}
					if(!key.equals("*") && refLengths.get(key) != Integer.parseInt(sp[1])) {
						br.close();
						throw new Exception("Mismatching lengths: " + file + " "
								+ key + " " + Integer.parseInt(sp[1]) + " " 
										+ refLengths.get(key));
					}
					int numMap = Integer.parseInt(sp[2]);
					totReads += numMap + Integer.parseInt(sp[3]);
					if(!key.equals("*")) {
						if(!map.containsKey(key)) {
							map.put(key, numMap);
						} else {
							map.put(key, map.get(key) + numMap);
						}
					}
					line = br.readLine();
				}
				br.close();
				numReads.put(srr, totReads);
				numMapped.put(srr, map);
				if(refs.size() != map.size()) {
					System.out.println("File " + file + 
							" has incorrect number of refs: " + map.size());
				}
			}
		}
		
		//get list of srr files
		HashSet<String> srrSet = new HashSet<String>();
		srrSet.addAll(sumDepth.keySet());
		srrSet.addAll(numMapped.keySet());
		List<String> srr = new ArrayList<String>();
		srr.addAll(srrSet);
		Collections.sort(srr);
		
		//write intermediate table to check numbers
		BufferedWriter srrOut = new BufferedWriter(new FileWriter(new File(
				DIR + "betaLactamaseAlignmentStatsBySRR_ref2.txt")));
		srrOut.write("id\ttotReads\treference\tnumMappedReads\tsumDepth\n");
		for(String s : srr) {
			for(String r : refs) {
				srrOut.write(s + "\t" + numReads.get(s) + "\t" + r + 
						"\t" + numMapped.get(s).get(r) + "\t" +
						sumDepth.get(s).get(r) + "\n");
			}
		}
		srrOut.close();
		
		//combine individual read files into one
		BufferedWriter chsOut = new BufferedWriter(new FileWriter(new File(
				DIR + "betaLactamaseAlignmentStatsByCHS_ref2.txt")));
		chsOut.write("chs\ttotReads\treference\tnumMappedReads\tpropMappedRead\tsumDepth\taveDepth\n");
		for(int s : strains) {
			String c = Integer.toString(s);
			String[] files = chs.get(c).toArray(new String[chs.get(c).size()]);
			int totReads = 0;
			for(String f : files) {
				totReads += numReads.get(f);
			}
			for(String r : refs) {
				int mapReads = 0;
				int dep = 0;
				for(String f : files) {
					mapReads += numMapped.get(f).get(r);
					dep += sumDepth.get(f).get(r);
				}
				chsOut.write(c + "\t" + totReads + "\t" + r + "\t" +
						mapReads + "\t" + (1.0 * mapReads / totReads) + "\t" +
						dep + "\t" + (1.0 * dep / refLengths.get(r)) + "\n");
			}
		}
		chsOut.close();
	}

}
