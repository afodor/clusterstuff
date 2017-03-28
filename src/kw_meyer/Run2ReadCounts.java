/*
 * Look at the number of reads for each step of run2
 */

package kw_meyer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class Run2ReadCounts {
	private static String DIR = "/nobackup/afodor_research/kwinglee/meyer/";

	public static void main(String[] args) throws IOException {
		////raw reads
		HashMap<String, Integer> raw = new HashMap<String, Integer>();//sample to number raw reads
		File[] readDirs = new File(DIR + "cardia_seq2_fastqs/").listFiles();
		for(File rd : readDirs) {
			File[] samples = rd.listFiles();
			for(File s : samples) {
				if(s.getName().endsWith(".fastq.gz") && s.getName().contains("_R1_")) {
					int count = 0;
					String id = s.getName().replace("_L001_R1_001.fastq", "");
					BufferedReader br = new BufferedReader(new FileReader(s));
					for(String line = br.readLine(); line != null; line = br.readLine()) {
						count++;
					}
					br.close();
					raw.put(id, count/4);
				}
			}
		}

		////filtered reads
		HashMap<String, Integer> filter = new HashMap<String, Integer>();//sample to number of reads after filtering primers
		File[] filtFiles = new File(DIR + "run2filteredSeqs/").listFiles();
		for(File f : filtFiles) {
			if(f.getName().endsWith(".fasta") && f.getName().contains("_R1_")) {
				int count = 0;
				String id = f.getName().replace("_L001_R1_001.fasta", "");
				BufferedReader br = new BufferedReader(new FileReader(f));
				for(String line = br.readLine(); line != null; line = br.readLine()) {
					if(line.startsWith(">")) {
						count++;						
					}
				}
				br.close();
				filter.put(id, count);
			}
		}

		////joined reads
		HashMap<String, Integer> join = new HashMap<String, Integer>();//sample to number of stitched reads
		File[] jFiles = new File(DIR + "run2joinedReads/").listFiles();
		for(File f : jFiles) {
			if(f.getName().endsWith("join.fasta")) {
				int count = 0;
				String id = f.getName().replace("join.fasta", "");
				BufferedReader br = new BufferedReader(new FileReader(f));
				for(String line = br.readLine(); line != null; line = br.readLine()) {
					if(line.startsWith(">")) {
						count++;						
					}
				}
				br.close();
				join.put(id, count);
			}
		}

		////write
		Set<String> sampSet = new HashSet<String>();
		sampSet.addAll(raw.keySet());
		sampSet.addAll(filter.keySet());
		sampSet.addAll(join.keySet());
		ArrayList<String> samples = new ArrayList<String>();
		samples.addAll(sampSet);
		Collections.sort(samples);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "Run2ReadCounts.txt")));
		out.write("sampleID\tnumRawReads\tnumFilteredReads\tnumJoinedReads\n");
		for(String key : samples) {
			out.write(key + "\t");
			if(raw.containsKey(key)) {
				out.write(raw.get(key) + "\t");
			} else {
				out.write("0\t");
			}
			if(filter.containsKey(key)) {
				out.write(filter.get(key) + "\t");
			} else {
				out.write("0\t");
			}
			if(join.containsKey(key)) {
				out.write(join.get(key) + "\n");
			} else {
				out.write("0\n");
			}
		}
		out.close();
	}
}
