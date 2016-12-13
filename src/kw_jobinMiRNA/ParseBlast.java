/*
 * Get number of reads that aligned using blast
 */
package kw_jobinMiRNA;

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

public class ParseBlast {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static double PID = 99;//minimum percent ID allowed
	public static int MISMATCH = 1;//maximum number of mismatches
	
	public static void main(String[] args) throws IOException {
		String[] dirs = new String[]{"miRBaseBlast/", "miRBaseBlastReadsAsDB/"};
		int[] qids = new int[]{0,1};
		for(int i = 0; i < dirs.length; i++) {
			String d = dirs[i];
			HashMap<String, Integer> mature = new HashMap<String, Integer>();//id -> number of reads mapped
			HashMap<String, Integer> hairpin = new HashMap<String, Integer>();//id -> number of reads mapped
			File[] files = new File(DIR + d).listFiles();
			for(File f : files) {
				if(f.getName().endsWith(".hairpin.txt")) {
					getCounts(f, hairpin, qids[i]);
				} else if(f.getName().endsWith(".mature.txt")) {
					getCounts(f, mature, qids[i]);
				}
			}
			
			writeCounts(mature, new File(DIR + d + "BlastMature.txt"));
			writeCounts(hairpin, new File(DIR + d + "BlastHairpin.txt"));
		}
	}

	/*
	 * map = map to write
	 * file = output file
	 */
	private static void writeCounts(HashMap<String, Integer> map, File file) throws IOException {
		ArrayList<String> samples = new ArrayList<String>(map.keySet());
		Collections.sort(samples);
		BufferedWriter out = new BufferedWriter(new FileWriter(file));
		out.write("sampleID\tnumReadsMapped\n");
		for(String s : samples) {
			out.write(s + "\t" + map.get(s) + "\n");
		}
		out.close();
	}

	/*
	 * file = file to read
	 * map = map to add results to
	 * qseqid = which column contains query sequence id
	 */
	private static void getCounts(File file, HashMap<String, Integer> map,
			int qseqid) throws IOException {
		HashSet<String> reads = new HashSet<String>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			String[] sp = line.split("\t");
			//qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
			if(Double.parseDouble(sp[2]) > PID &&
					Integer.parseInt(sp[4]) <= MISMATCH) {
				reads.add(sp[qseqid]);
			}
		}
		br.close();		
		map.put(file.getName().split("\\.")[0], reads.size());
	}
}
