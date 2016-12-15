/*
 * generate table of number of reads aligned to each reference
 * one table per blast, bowtie, and using Java to get exact matches
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class ResultsTables {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	private static HashMap<String, Integer> NUMREADS = new HashMap<String, Integer>();//map of id to number of reads

	public static double PID = 99;//minimum percent ID allowed
	public static int MISMATCH = 1;//maximum number of mismatches
	
	public static void main(String[] args) throws IOException {
		//bowtie (and get number of reads)
		analyzeBowtie();

		//blast
		analyzeBlast();

		//java
		//analyzeJava();
	}

	private static void analyzeBowtie() throws IOException {
		HashMap<String, Integer> miRmature = new HashMap<String, Integer>();
		HashMap<String, Integer> miRhairpin = new HashMap<String, Integer>();
		HashMap<String, Integer> mouse = new HashMap<String, Integer>();
		HashMap<String, Integer> piR = new HashMap<String, Integer>();
	
		//miR
		File[] files = new File(DIR + "miRBaseBowtie").listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".sam")) {
				String id = f.getName().split("\\.")[0];
				int count = getBowtieCounts(f);
				if(f.getName().contains("mature")) {
					miRmature.put(id, count);
				} else if(f.getName().contains("hairpin")) {
					miRhairpin.put(id, count);
				} else {
					System.out.println("Extra sam: " + f.getName() + ' ' + count);
				}
			}
		}
	
		//mouse
		files = new File(DIR + "mouseBowtie/").listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".sam")) {
				mouse.put(f.getName().split("\\.")[0], getBowtieCounts(f));
			}
		}
	
		//piR
		/*files = new File(DIR + "piRBaseBowtie/").listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".sam")) {
				piR.put(f.getName().split("\\.")[0], getBowtieCounts(f));
			}
		}*/
		for(int i = 1; i <= 10; i++) {
			String key = "Sample" + i;
			piR.put(key, -1);
		}
	
		writeTable("bowtie", miRmature, miRhairpin, mouse, piR);
	}

	private static void analyzeBlast() throws IOException {
		HashMap<String, Integer> miRmature = new HashMap<String, Integer>();
		HashMap<String, Integer> miRhairpin = new HashMap<String, Integer>();
		HashMap<String, Integer> mouse = new HashMap<String, Integer>();
		HashMap<String, Integer> piR = new HashMap<String, Integer>();

		//miR
		File[] files = new File(DIR + "miRBaseBlast").listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".txt")) {
				String id = f.getName().split("\\.")[0];
				int count = getBlastCounts(f);
				if(f.getName().contains("mature")) {
					miRmature.put(id, count);
				} else if(f.getName().contains("hairpin")) {
					miRhairpin.put(id, count);
				} else {
					System.out.println("Extra sam: " + f.getName() + ' ' + count);
				}
			}
		}

		//mouse
		files = new File(DIR + "mouseBlast").listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".mouse.blast.txt")) {
				mouse.put(f.getName().split("\\.")[0], getBlastCounts(f));
			}
		}

		//piR
		files = new File(DIR + "piRBaseBlast").listFiles();
		for(File f : files) {
			if(f.getName().endsWith(".piR.blast.txt")) {
				piR.put(f.getName().split("\\.")[0], getBlastCounts(f));
			}
		}

		writeTable("blast", miRmature, miRhairpin, mouse, piR);
	}

	private static void analyzeJava() throws IOException {
		HashMap<String, Integer> miRmature = new HashMap<String, Integer>();
		HashMap<String, Integer> miRhairpin = new HashMap<String, Integer>();
		HashMap<String, Integer> mouse = new HashMap<String, Integer>();
		HashMap<String, Integer> piR = new HashMap<String, Integer>();
	
		//mouse
		/*for(int i = 1; i <= 10; i++) {
			String key = "Sample" + i;
			mouse.put(key, NA);
		}*/
	
		writeTable("java", miRmature, miRhairpin, mouse, piR);
	}

	private static int getBowtieCounts(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String id = file.getName().split("\\.")[0];
		int numReads = 0;
		HashSet<String> mapped = new HashSet<String>();//set of reads mapped (in case read mapped multiple times)
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			if(!line.startsWith("@")) {
				numReads++;
				String[] sp = line.split("\\s+");
				if(!sp[1].equals("4")) {//read is mapped
					mapped.add(sp[0]);
				}
			}
		}
		br.close();
		if(NUMREADS.containsKey(id)) {
			if(NUMREADS.get(id) != numReads) {
				System.err.println("Discordant read number: " + id + " " + NUMREADS.get(id)
						+ " " + file.getName() + " " + numReads);
				if(numReads < NUMREADS.get(id)) {
					NUMREADS.put(id, numReads);
				}
			}
		} else {
			NUMREADS.put(id, numReads);
		}
		return(mapped.size());
	}

	private static Integer getBlastCounts(File file) throws IOException {
		HashSet<String> reads = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			String[] sp = line.split("\t");
			//qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
			if(Double.parseDouble(sp[2]) > PID &&
					Integer.parseInt(sp[4]) <= MISMATCH) {
				reads.add(sp[0]);
			}
		}
		br.close();		
		return(reads.size());
	}

	private static void writeTable(String name,
			HashMap<String, Integer> miRmature,
			HashMap<String, Integer> miRhairpin,
			HashMap<String, Integer> mouse, HashMap<String, Integer> piR) throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "Results_" + name + ".txt")));
		out.write("sampleID\tgroup\ttotalNumReads\tnumMapMiRmature\t"
				+ "numMapMiRhairpin\tnumMapMouse\tnumMapPiR\n");
		for(int i = 1; i <= 10; i++) {
			String key = "Sample" + i;
			out.write(key + "\t");
			if(i <= 5) {
				out.write("GF\t");
			} else {
				out.write("SPF\t");
			}
			out.write(NUMREADS.get(key) + "\t" + miRmature.get(key) + "\t" +
					miRhairpin.get(key) + "\t" + mouse.get(key) + "\t" +
					piR.get(key) + "\n");
		}
		out.close();		
	}
}
