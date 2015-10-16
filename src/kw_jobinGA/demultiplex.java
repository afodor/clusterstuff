/**
 * Demultiplex GA-stools sample from Jobin group, writing results as fasta
 */

package kw_jobinGA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.zip.GZIPInputStream;

public class demultiplex {
	//directory containing the needed files; also where will write results
	//public static final String DIR = "/projects/afodor_research/kwinglee/jobin/ga-stool/";
	public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\JobinCollaboration\\GA-stools\\V1_V3_16S_GA+stools_2-25611692\\Sample_1-29344834\\Data\\Intensities\\BaseCalls\\";

	public static void main(String[] args) throws IOException {
		analyze("Sample-Name-1_S1_L001_R1_001.fastq.gz", "Run2_R1_");
	}
	
	/**
	 * Reverse complements the sequence seq
	 * @param seq	the sequence to reverse complement
	 * @return the reverse complement
	 */
	public static String revComp(String seq) {
		char[] fwd = seq.toUpperCase().toCharArray();
		String rev = "";
		for(int i = fwd.length-1; i >= 0; i--) {
			if(fwd[i] == 'A') {
				rev += "T";
			} else if(fwd[i] == 'T') {
				rev += "A";
			} else if(fwd[i] == 'G') {
				rev += "C";
			} else if(fwd[i] == 'C') {
				rev += "G";
			} else {
				throw new IllegalArgumentException("Invalid base: " + fwd[i]);
			}
		}
		return(rev);
	}
	
	public static void analyze(String fastqFile, String outPrefix) throws IOException {
		//set up dictionary of primer name to primer sequence
		HashMap<String, String> pToSeq = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "primer_to_sequence.txt"))));
		br.readLine();//header
		String line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			pToSeq.put(sp[0].split("_")[2], sp[1]);
			line = br.readLine();
		}
		br.close();
		//check
		/*Iterator<String> it = pToSeq.keySet().iterator();
		while(it.hasNext()) {
			String k = it.next();
			System.out.println(k + "\t" + pToSeq.get(k));
		}
		System.out.println(pToSeq.get("F12"));
		System.out.println(revComp(pToSeq.get("F12")));*/
		
		//set up dictionary of forward-reverse primer pair to sample (ex F12R1 -> S1)
		//and set up dictionary of sample id to writer
		HashMap<String, String> pToSamp = new HashMap<String, String>();//primer pair to sample
		HashMap<String, BufferedWriter> sToFile = new HashMap<String, BufferedWriter>();//sample to sample's output file
		br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "primer_to_sample.txt"))));
		br.readLine();//header
		line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String samp = sp[1]; //sample name
			pToSamp.put(sp[2]+sp[3], samp);
			sToFile.put(samp, 
					new BufferedWriter(new FileWriter(new File(DIR + outPrefix + samp))));
			line = br.readLine();
		}
		br.close();
		//check
		/*Iterator<String> it = pToSamp.keySet().iterator();
		while(it.hasNext()) {
			String k = it.next();
			System.out.println(k + "\t" + pToSamp.get(k));
		}*/
		
		//add extra "other" file for unmatched reads
		sToFile.put("other",
				new BufferedWriter(new FileWriter(new File(DIR + outPrefix + "other"))));
		
		//list of primers used
		String[] forward = {"F12", "F13", "F14", "F15"};
		String[] reverse = new String[16];
		for(int i = 1; i < 17; i++) {
			reverse[i-1] = "R" + Integer.toString(i);
		}
		
		
		//read fastq file
		//convert each read to fast and determine the correct read
		BufferedReader fastq = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(DIR + fastqFile))));
		String l1 = fastq.readLine();//first line of read (@Seqid)
		while(l1 != null) {
			//get rest of read info
			String l2 = fastq.readLine(); //second line of read (actual sequence)
			fastq.readLine();//third line (+)
			fastq.readLine();//4th line (quality scores)
			
			String header = l1.replace("@", ">");//fasta header
			
			//figure out what the sample is
			/**
			 * currently just assuming if fwd and rev is present is correct
			 * not checking that both are reverse complemented or both forward
			 */
			String samp = "other";
			String fwd = "";
			String rev = "";
			for(int f = 0; f < forward.length; f++) {//iterate over forward primers
				String fseq = pToSeq.get(forward[f]);
				if(l1.startsWith(fseq) ||
						l1.startsWith(revComp(fseq)) ||
						l1.endsWith(fseq) ||
						l1.endsWith(revComp(fseq))){ //f is the beginning of the sequence
					fwd = forward[f];
					//remove forward primer
					l1 = l1.replace(fseq, "");
					l1 = l1.replace(revComp(fseq), "");
					for(int r = 0; r < reverse.length; r++) {
						String rseq = pToSeq.get(reverse[r]);
						if(l1.startsWith(rseq) ||
								l1.startsWith(revComp(rseq)) ||
								l1.endsWith(rseq) ||
								l1.endsWith(revComp(rseq))) {
							rev = reverse[r];
							//remove reverse primer
							l1 = l1.replace(rseq, "");
							l1 = l1.replace(revComp(rseq), "");
							break;
						}
					}
					break;
				}
			}
			String key = fwd+rev;
			if(pToSamp.containsKey(key)) {
				samp = pToSamp.get(key);
			}
			
			//write to file
			BufferedWriter out = sToFile.get(samp);
			out.write(header + "\n" + l1 + "\n");
			
			l1 = fastq.readLine();
		}
		
		//close all files
		fastq.close();
		Iterator<String> it = sToFile.keySet().iterator();
		while(it.hasNext()) {
			String k = it.next();
			sToFile.get(k).close();
		}
	}
}
