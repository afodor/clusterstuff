/*
 * get count table from BWA alignments
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;

public class ParseBWA {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/miRBaseBWA/";
	
	public static void main(String[] args) throws IOException {
		HashMap<String, ArrayList<Integer>> hairpin = new HashMap<String, ArrayList<Integer>>();
		HashMap<String, ArrayList<Integer>> mature = new HashMap<String, ArrayList<Integer>>();
		String[] sams = new File(DIR).list();
		for(String sam : sams) {
			if(sam.endsWith(".mature.sam")) {
				getCounts(sam, mature);
			} else if(sam.endsWith(".hairpin.sam")) {
				getCounts(sam, hairpin);
			}
		}
		
		writeCounts(mature, "BWAMature.txt");
		writeCounts(hairpin, "BWAHairpin.txt");
	}

	private static void writeCounts(HashMap<String, ArrayList<Integer>> map,
			String name) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + name)));
		out.write("sampleID\tnumReads\tnumMapped\tnumUnique\tnumPerfect\n");
		ArrayList<String> samples = new ArrayList<String>(map.keySet());
		Collections.sort(samples);
		for(String s : samples) {
			ArrayList<Integer> counts = map.get(s);
			out.write(s + "\t" + counts.get(0) + "\t" + counts.get(1)
					+ "\t" + counts.get(2) + "\t" + counts.get(3) + "\n");
		}
		out.close();		
	}

	private static void getCounts(String sam,
			HashMap<String, ArrayList<Integer>> map) throws IOException {
		String perfect = "";
		for(int i = 0; i < 25; i++) {
			perfect += "M";
		}//perfect match is 25 matches
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + sam)));
		String id = sam.split("\\.")[0];
		int numReads = 0;
		int numMapped = 0;
		int numUnique = 0;
		int numPerfect = 0;
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			if(!line.startsWith("@")) {
				numReads++;
				String[] sp = line.split("\\s+");
				if(!sp[1].equals("4")) {//read is mapped
					numMapped++;
					System.err.println(line);
					if(!sp[1].equals("256")) {//not secondary alignment
						numUnique++;
					}
					if(sp[5].equals(perfect)) {
						numPerfect++;
						System.out.println(id + "\t" + sp[3] + "\tperfect");//print which ref aligned
					} else {
						System.out.println(id + "\t" + sp[3]);//print which ref aligned
					}
				}
			}
		}
		br.close();
		ArrayList<Integer> results = new ArrayList<Integer>(
				Arrays.asList(numReads, numMapped, numUnique, numPerfect));
		map.put(id, results);
	}
}
