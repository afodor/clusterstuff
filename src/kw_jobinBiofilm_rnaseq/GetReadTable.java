/*
 * Make a table of the number of reads for each sample after each step
 */
package kw_jobinBiofilm_rnaseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class GetReadTable {
	public static String BASEDIR = "/nobackup/afodor_research/kwinglee/jobin/biofilm/";

	public static void main(String[] args) throws IOException {
		HashMap<String, ArrayList<Integer>> counts = new HashMap<String, ArrayList<Integer>>();
		//get number of reads in SILVA filtered
		String silvaDir = BASEDIR + "rnaseq/mouseAndSilvaFiltered/";
		File[] files = new File(silvaDir).listFiles();
		for(File f : files) {
			if(f.getName().endsWith("_R1.mouseFiltered.silvaFiltered.fasta")) {
				String name = f.getName().replace("_R1.mouseFiltered.silvaFiltered.fasta", "");
				int r1 = getNumFastaReads(f);
				int r2 = getNumFastaReads(new File(f.getAbsolutePath().replace("_R1", "_R2")));

				ArrayList<Integer> reads = new ArrayList<Integer>();
				reads.add(r1+r2);
				counts.put(name, reads);	
			}
		}

		//get number of reads in mouse filtered
		String mouseDir = BASEDIR + "rnaseq/mouseFilteredFastq/";
		files = new File(mouseDir).listFiles();
		for(File f : files) {
			if(f.getName().endsWith("_R1.mouseFiltered.fastq")) {
				String name = f.getName().replace("_R1.mouseFiltered.fastq", "");
				int r1 = getNumFastqReads(f);
				int r2 = getNumFastqReads(new File(f.getAbsolutePath().replace("_R1", "_R2")));

				if(counts.containsKey(name)) {
					counts.get(name).add(r1+r2);
				} else {
					System.out.println("Bad mouse key: " + name);
					ArrayList<Integer> reads = new ArrayList<Integer>();
					reads.add(r1+r2);
					counts.put(name, reads);
				}
			}
		}

		//get number of reads in initial file 
		String allDir = BASEDIR + "RNAseqTestRunFastas/";
		files = new File(allDir).listFiles();
		for(File f : files) {
			if(f.getName().endsWith("_R1.fasta")) {
				String name = f.getName().replace("_R1.fasta", "");
				int r1 = getNumFastaReads(f);
				int r2 = getNumFastaReads(new File(f.getAbsolutePath().replace("_R1", "_R2")));

				if(counts.containsKey(name)) {
					counts.get(name).add(r1+r2);
				} else {
					System.out.println("Bad all key: " + name);
					ArrayList<Integer> reads = new ArrayList<Integer>();
					reads.add(r1+r2);
					counts.put(name, reads);
				}
			}
		}

		//write results
		ArrayList<String> keys = new ArrayList<String>(counts.keySet());
		Collections.sort(keys);
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				BASEDIR + "rnaseq/microbiomeReadTable.txt")));
		out.write("sample\tnumInitialReads\tnumReadsAfterFilterMouse\tnumReadsAfterFilterMouseThenrRNA\n");
		for(String k : keys) {
			ArrayList<Integer> reads = counts.get(k);
			out.write(k + "\t" + reads.get(2) + "\t" + reads.get(1) + "\t" + reads.get(0) + "\n");
		}
		out.close();
	}

	private static int getNumFastaReads(File file) throws IOException {
		int numReads = 0;
		BufferedReader br = new BufferedReader(new FileReader(file));
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			if(line.startsWith(">")) {
				numReads++;
			}
		}
		br.close();
		return(numReads);
	}

	private static int getNumFastqReads(File file) throws IOException {
		int numReads = 0;
		BufferedReader br = new BufferedReader(new FileReader(file));
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			numReads++;
			line = br.readLine();
			line = br.readLine();
			line = br.readLine();
		}
		br.close();
		return(numReads);
	}
}

