/*
 * parse F. rodentium blast alignments
 */
package kw_jobinFrodentium;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class parseBlastResults {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/F_rodentium/blastResults/";
	//public static int MINLEN = 260;//average reads are all around 275
	private static int MINLEN1 = 250;
	private static double MINPID1 = 90;
	private static int MINLEN2 = 260;
	private static double MINPID2 = 95;
	
	public static void main(String[] args) throws IOException {
		File[] files = new File(DIR).listFiles();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "blastResults_Frodentium.txt")));
		out.write("sample\tnoFilter\tPID" + MINPID1 + "len" + MINLEN1 + "\tPID" + 
				MINPID2 + "len" + MINLEN2 + "\n");
		for(File f : files) {
			String sid = f.getName();
			if(sid.startsWith("Frod_") && sid.endsWith(".txt")) {
				sid = sid.replace("Frod_", "").replace(".txt", "");
				HashSet<String> noFilt = new HashSet<String>();
				HashSet<String> filt1 = new HashSet<String>();
				HashSet<String> filt2 = new HashSet<String>();
				BufferedReader res = new BufferedReader(new FileReader(f));
				for(String line = res.readLine(); line != null; line = res.readLine()) {
					if(!line.startsWith("#")) {
						String[] sp = line.split("\t");
						//Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
						String read = sp[0];
						double pid = Double.parseDouble(sp[2]);
						int length = Integer.parseInt(sp[3]);
						noFilt.add(read);
						if(length > MINLEN1 && pid > MINPID1) {
							filt1.add(read);
						}
						if(length > MINLEN2 && pid > MINPID2) {
							filt2.add(read);
						}
					}
				}
				out.write(sid + "\t" + noFilt.size() + "\t" + filt1.size() + "\t"
						+ filt2.size() + "\n");
				res.close();
			}
		}
		out.close();
	}

}
