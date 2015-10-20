/**
 * Demultiplex GA-stools sample from Jobin group, writing results as fasta
 */

package kw_jobinGA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

public class demultiplex {
	//directory containing the needed files; also where will write results
	//public static final String DIR = "/projects/afodor_research/kwinglee/jobin/ga-stool/";
	public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\JobinCollaboration\\GA-stools\\V1_V3_16S_GA+stools_2-25611692\\Sample_1-29344834\\Data\\Intensities\\BaseCalls\\";

	public static final String[] FWD_PRIMERS = {"F12", "F13", "F14", "F15"};//forward primers used
	public static String[] REV_PRIMERS = {};//reverse primers; filled in in main
	private static HashMap<String, String> pToSeq;//hash of primer to primer sequence
	private static int numMultiple = 0;
	
	public static void main(String[] args) throws IOException {
		//make list of reverse primers used
		REV_PRIMERS = new String[16];
		for(int i = 1; i < 17; i++) {
			REV_PRIMERS[i-1] = "R" + Integer.toString(i);
		}
				
		analyze("Sample-Name-1_S1_L001_R1_001.fastq.gz", "Sample-Name-1_S1_L001_R2_001.fastq.gz", "Run2_");
		//analyze("Run2-Sample-Name-1_S1_L001_R1_001.fastq.gz", "Run2-Sample-Name-1_S1_L001_R2_001.fastq.gz", "Run2_");
		//analyze("Run1-Sample-Name-1_S1_L001_R1_001.fastq.gz", "Run1-Sample-Name-1_S1_L001_R2_001.fastq.gz", "Run1_");
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
	
	/**
	 * Return the primer sequence for a given read (seq)
	 * @param seq
	 * @return
	 */
	public static String getPrimer(String seq) {
		String primer = null;
		//iterate over primers
		Iterator<String> it = pToSeq.keySet().iterator();
		while(it.hasNext()) {
			String p = it.next();
			String pseq = pToSeq.get(p);
			if(seq.startsWith(pseq) ||
					seq.startsWith(revComp(pseq)) ||
					seq.endsWith(pseq) ||
					seq.endsWith(revComp(pseq))){
				if(primer != null) {//multiple primers
					numMultiple++;
					return(null);
				}
				primer = p;
			}
		}
		return(primer);
		
		
		/*for(int f = 0; f < FWD_PRIMERS.length; f++) {//iterate over forward primers
			String fseq = pToSeq.get(FWD_PRIMERS[f]);
			if(seq.startsWith(fseq) ||
					seq.startsWith(revComp(fseq)) ||
					seq.endsWith(fseq) ||
					seq.endsWith(revComp(fseq))){ //f is the beginning of the sequence
				//remove forward primer
				seq = seq.replace(fseq, "");
				seq = seq.replace(revComp(fseq), "");
				for(int r = 0; r < REV_PRIMERS.length; r++) {
					String rseq = pToSeq.get(REV_PRIMERS[r]);
					if(seq.startsWith(rseq) ||
							seq.startsWith(revComp(rseq)) ||
							seq.endsWith(rseq) ||
							seq.endsWith(revComp(rseq))) {
						//remove reverse primer
						seq = seq.replace(rseq, "");
						seq = seq.replace(revComp(rseq), "");
						numRevEnd++;
						break;
					}
				}
				break;
			}
		}*/
	}
	
	public static void analyze(String fastqFileF, String fastqFileR, String outPrefix) throws IOException {
		//set up dictionary of primer name to primer sequence
		pToSeq = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "primer_to_sequence.txt"))));
		br.readLine();//header
		String line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String pseq = sp[1].replace("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "");//this initial part of the forward primer seems to have been already removed
			pseq = pseq.replace("CAAGCAGAAGACGGCATACGAGATCGGCATTCCTGCTGAACCGCTCTTCCGATCT", "");//5' reverse read
			//pseq = pseq.replace("ATTACCGCGGCTGCTGG", "");//3' of reverse
			pToSeq.put(sp[0].split("_")[2], pseq);
			line = br.readLine();
		}
		br.close();
		//check
		/*Iterator<String> test = pToSeq.keySet().iterator();
		while(test.hasNext()) {
			String k = test.next();
			System.out.println(k + "\t" + pToSeq.get(k));
		}
		System.out.println(pToSeq.get("F12"));
		System.out.println(revComp(pToSeq.get("F12")));*/
		
		//set up dictionary of forward-reverse primer pair to sample (ex F12R1 -> S1)
		//and set up dictionary of sample id to writer
		HashMap<String, String> pToSamp = new HashMap<String, String>();//primer pair to sample
		HashMap<String, BufferedWriter[]> sToFile = new HashMap<String, BufferedWriter[]>();//sample to sample's output file
		br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "primer_to_sample.txt"))));
		br.readLine();//header
		line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String samp = sp[1]; //sample name
			pToSamp.put(sp[2]+sp[3], samp);
			BufferedWriter[] files = {new BufferedWriter(new FileWriter(new File(DIR + outPrefix + "R1_" + samp + ".fasta"))),
					new BufferedWriter(new FileWriter(new File(DIR + outPrefix + "R2_" + samp + ".fasta")))};
			sToFile.put(samp, files);
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
		BufferedWriter[] files = {new BufferedWriter(new FileWriter(new File(DIR + outPrefix + "R1_other.fasta"))),
				new BufferedWriter(new FileWriter(new File(DIR + outPrefix + "R2_other.fasta")))};
		sToFile.put("other", files);
		
		
		//read fastq file
		//convert each read to fast and determine the correct read
		BufferedReader fastqF = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(DIR + fastqFileF))));
		BufferedReader fastqR = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(DIR + fastqFileR))));
		String headF = fastqF.readLine();//first line of read (@Seqid) for forward reads
		String headR = fastqR.readLine();
		int numRead = 0;
		int numMatch = 0;
		while(headF != null && headR != null) {
			//get rest of read info
			String readF = fastqF.readLine(); //second line of read (actual sequence) for forward reads
			fastqF.readLine();//third line (+)
			fastqF.readLine();//4th line (quality scores)
			String readR = fastqR.readLine();
			fastqR.readLine();
			fastqR.readLine();
			numRead++;
			//System.out.println(seq);
			
			String headerF = headF.replaceFirst("@", ">");//fasta header
			String headerR = headR.replaceFirst("@", ">");
			
			//figure out what the sample is
			/**
			 * currently just assuming if fwd and rev is present is correct
			 * not checking that both are reverse complemented or both forward
			 */
			String samp = "other";
			String pF = getPrimer(readF);
			String pR = getPrimer(readR);
			
			
			//trim primers from sequence
			String key = "";
			if(pF != null && pR != null) {				
				//check that one is a fwd primer and one is rev
				String fwd = "";
				String rev = "";
				if(pF.startsWith("F") && pR.startsWith("R")) {
					fwd = pF;
					rev = pR;
					numMatch++;
				} else if(pF.startsWith("R") && pR.startsWith("F")) {
					fwd = pR;
					rev = pF;
					numMatch++;
				}
				
				key = fwd+rev;
			}
			
			if(pToSamp.containsKey(key)) {
				samp = pToSamp.get(key);
				//only remove primers if not going into other category
				readF = readF.replace(pToSeq.get(pF), "");
				readF = readF.replace(revComp(pToSeq.get(pF)), "");
				readR = readR.replace(pToSeq.get(pR), "");
				readR = readR.replace(revComp(pToSeq.get(pR)), "");
			}
			
			if(numRead % 1000000 == 0) {
				System.out.println("numread = " + numRead + " num sorted = " + numMatch + " num multiple " + numMultiple);
			}
			
			//write to file
			//Forward read
			BufferedWriter out = sToFile.get(samp)[0];
			out.write(headerF + "\n" + readF + "\n");
			//Reverse read
			out = sToFile.get(samp)[1];
			out.write(headerR + "\n" + readR + "\n");
			
			headF = fastqF.readLine();
			headR = fastqR.readLine();
		}
		
		if((headF != null && headR == null) ||
				headF == null && headR != null) {
			System.out.println("Uneven number of reads");
		}
		
		//close all files
		System.out.println("Total reads = " + numRead);//13896021
		System.out.println("Reads matched = " + numMatch);//8081744
		System.out.println("Reads with multiple primers = " + numMultiple);//1166833
		fastqF.close();
		fastqR.close();
		Iterator<String> it = sToFile.keySet().iterator();
		while(it.hasNext()) {
			String k = it.next();
			sToFile.get(k)[0].close();
			sToFile.get(k)[1].close();
		}
	}
}
