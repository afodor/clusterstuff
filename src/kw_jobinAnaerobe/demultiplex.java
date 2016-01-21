/**
 * Demultiplex C jejuni aerobe vs anaerobe samples from Jobin group, writing results as fasta
 * 12/4/15
 */

package kw_jobinAnaerobe;

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
	public static final String DIR = "/nobackup/afodor_research/kwinglee/jobin/anaerobe/";
	//public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\JobinCollaboration\\C jejuni anaerobe v aerobe\\JG_V1_V3_16S_Aerobe_vs_Anaerobe-26697712\\Sample_1-30857322\\Data\\Intensities\\BaseCalls\\";

	private static HashMap<String, String> P_TO_SEQ;//hash of primer to primer sequence
	private static int numMultiple = 0;
	
	public static void main(String[] args) throws IOException {
		analyze("Sample-Name-1_S1_L001_R1_001.fastq.gz", "Sample-Name-1_S1_L001_R2_001.fastq.gz");
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
				System.out.println(seq);
				System.out.println(fwd[i]);
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
		Iterator<String> it = P_TO_SEQ.keySet().iterator();
		while(it.hasNext()) {
			String p = it.next();
			String pseq = P_TO_SEQ.get(p);
			if(seq.startsWith(pseq) ||
					seq.startsWith(revComp(pseq)) ||
					seq.endsWith(pseq) ||
					seq.endsWith(revComp(pseq))){
				if(primer != null) {//multiple primers
					/*System.out.println(p);
					System.out.println(primer);
					System.out.println(seq);
					System.out.println();*/
					numMultiple++;
					return(null);
				}
				primer = p;
			}
		}
		return(primer);
	}
	
	public static void analyze(String fastqFileF, String fastqFileR) throws IOException {
		//set up dictionary of primer name to primer sequence
		P_TO_SEQ = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "AnaerobePrimers.txt"))));
		br.readLine();//header
		String line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String pseq = sp[1].replace("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "");//this initial part of the forward primer seems to have been already removed
			pseq = pseq.replace("CAAGCAGAAGACGGCATACGAGATCGGCATTCCTGCTGAACCGCTCTTCCGATCT", "");//5' reverse read
			String pname = sp[0].replace("PE1_27", "").replace("PE2_534", "").replace("-", "");//primer name -> reduce to ex F1/R1
			P_TO_SEQ.put(pname, pseq);
			line = br.readLine();
		}
		br.close();
		/*Iterator<String> itSet = P_TO_SEQ.keySet().iterator();
		while(itSet.hasNext()) {
			String p = itSet.next();
			System.out.println(P_TO_SEQ.get(p) + "\t" + p);
		}*/
		
		//set up dictionary of forward-reverse primer pair to sample (ex F12R1 -> S1)
		//and set up dictionary of sample id to writers
		HashMap<String, String> pToSamp = new HashMap<String, String>();//primer pair to sample
		HashMap<String, BufferedWriter[]> sToFile = new HashMap<String, BufferedWriter[]>();//sample to sample's output file
		br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "AnaerobeMetadata.txt"))));
		br.readLine();//header
		line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String samp = sp[0]; //sample name
			pToSamp.put(sp[sp.length-2]+sp[sp.length-1], samp);
			BufferedWriter[] files = {new BufferedWriter(new FileWriter(new File(DIR + File.separator + "fastas" + File.separator + samp + "_R1.fasta"))),
					new BufferedWriter(new FileWriter(new File(DIR + File.separator + "fastas" + File.separator + samp + "_R2.fasta")))};
			sToFile.put(samp, files);
			line = br.readLine();
		}
		br.close();
		/*Iterator<String> itSet = pToSamp.keySet().iterator();
		while(itSet.hasNext()) {
			String p = itSet.next();
			System.out.println(p + "\t" + pToSamp.get(p));
		}*/
		
		//add extra "other" file for unmatched reads
		BufferedWriter[] files = {new BufferedWriter(new FileWriter(new File(DIR + File.separator + "fastas" + File.separator + "other_R1.fasta"))),
				new BufferedWriter(new FileWriter(new File(DIR + File.separator + "fastas" + File.separator + "other_R2.fasta")))};
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
				} /*else {
					System.out.println(pF + "\t" + pR);
					System.out.println(readF);
					System.out.println(readR);
				}*/
				
				key = fwd+rev;
			}
			
			if(pToSamp.containsKey(key)) {
				samp = pToSamp.get(key);
				/*if(samp.equals("183")) {
					System.out.println("Full sequence:");
					System.out.println(readF);
					System.out.println(readR);
				}*/
				//only remove primers if not going into other category
				readF = readF.replace(P_TO_SEQ.get(pF), "");
				readF = readF.replace(revComp(P_TO_SEQ.get(pF)), "");
				readR = readR.replace(P_TO_SEQ.get(pR), "");
				readR = readR.replace(revComp(P_TO_SEQ.get(pR)), "");
				/*if(samp.equals("183")) {
					System.out.println("Primers removed:");
					System.out.println(readF);
					System.out.println(readR);
					System.out.println();
				}*/
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
		System.out.println("Total reads = " + numRead);
		System.out.println("Reads matched = " + numMatch);
		System.out.println("Reads with multiple primers = " + numMultiple);
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
