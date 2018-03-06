/**
 * Demultiplex GA-stools sample from Jobin group, writing results as fasta
 */

package kw_jobinGA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

public class demultiplexFastQ {
	//directory containing the needed files; also where will write results
	public static final String DIR = "/nobackup/afodor_research/archive/kwinglee/jobin/ga-stool/";
	//public static final String DIR = "C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\JobinCollaboration\\GA-stools\\V1_V3_16S_GA+stools_2-25611692\\Sample_1-29344834\\Data\\Intensities\\BaseCalls\\";

	private static HashMap<String, String> P_TO_SEQ;//hash of primer to primer sequence
	private static int numMultiple = 0;
	
	public static void main(String[] args) throws Exception {
		//analyze("Sample-Name-1_S1_L001_R1_001.fastq.gz", "Sample-Name-1_S1_L001_R2_001.fastq.gz", "Run2_");
		analyze("Run2-Sample-Name-1_S1_L001_R1_001.fastq.gz", "Run2-Sample-Name-1_S1_L001_R2_001.fastq.gz", "Run2_");
		analyze("Run1-Sample-Name-1_S1_L001_R1_001.fastq.gz", "Run1-Sample-Name-1_S1_L001_R2_001.fastq.gz", "Run1_");
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
		Iterator<String> it = P_TO_SEQ.keySet().iterator();
		while(it.hasNext()) {
			String p = it.next();
			String pseq = P_TO_SEQ.get(p);
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
		
	}
	
	public static void analyze(String fastqFileF, String fastqFileR, String outPrefix) throws Exception {
		//set up dictionary of primer name to primer sequence
		P_TO_SEQ = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "primer_to_sequence.txt"))));
		br.readLine();//header
		String line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String pseq = sp[1].replace("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT", "");//this initial part of the forward primer seems to have been already removed
			pseq = pseq.replace("CAAGCAGAAGACGGCATACGAGATCGGCATTCCTGCTGAACCGCTCTTCCGATCT", "");//5' reverse read
			//pseq = pseq.replace("ATTACCGCGGCTGCTGG", "");//3' of reverse
			P_TO_SEQ.put(sp[0].split("_")[2], pseq);
			line = br.readLine();
		}
		br.close();
		
		//set up dictionary of forward-reverse primer pair to sample (ex F12R1 -> S1)
		//and set up dictionary of sample id to writers
		HashMap<String, String> pToSamp = new HashMap<String, String>();//primer pair to sample
		HashMap<String, BufferedWriter[]> sToFile = new HashMap<String, BufferedWriter[]>();//sample to sample's output file
		br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(DIR + "primer_to_sample.txt"))));
		br.readLine();//header
		line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String samp = sp[1]; //sample name
			pToSamp.put(sp[2]+sp[3], samp);
			System.out.println("CREATING " + DIR + File.separator + 
					"fastqSeqs" + File.separator + outPrefix + "R1_" + samp + ".fastq" );
			BufferedWriter[] files = {new BufferedWriter(new FileWriter(new File(DIR + File.separator + 
					"fastqSeqs" + File.separator + outPrefix + "R1_" + samp + ".fastq"))),
					new BufferedWriter(new FileWriter(new File(DIR + File.separator + "fastqSeqs" + 
							File.separator + outPrefix + "R2_" + samp + ".fastq")))};
			sToFile.put(samp, files);
			line = br.readLine();
		}
		br.close();
		
		//add extra "other" file for unmatched reads
		BufferedWriter[] files = {new BufferedWriter(new FileWriter(new File(DIR + File.separator + 
				"fastqSeqs" + outPrefix + "R1_other.fastq"))),
				new BufferedWriter(new FileWriter(new File(DIR + File.separator + 
						"fastqSeqs" + outPrefix + "R2_other.fastq")))};
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
			String thirdF =fastqF.readLine();//third line (+)
			String fourthF = fastqF.readLine();//4th line (quality scores)
			String readR = fastqR.readLine();
			String thirdR = fastqR.readLine();
			String fourthR = fastqR.readLine();
			numRead++;
			
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
				
				Holder h1 = trim( readF, fourthF,P_TO_SEQ.get(pF));
				readF = h1.seqString;
				fourthF = h1.qualScoreString;
				
				Holder h2 = trim( readF, fourthF,revComp(P_TO_SEQ.get(pF)));
				readF = h2.seqString;
				fourthF = h2.qualScoreString;
				
				Holder h3 = trim( readR, fourthR,P_TO_SEQ.get(pR));
				readR = h3.seqString;
				fourthR = h3.qualScoreString;
				
				Holder h4 = trim( readR, fourthR,revComp(P_TO_SEQ.get(pR)));
				readR = h4.seqString;
				fourthR = h4.qualScoreString;
			}
			
			if(numRead % 1000000 == 0) {
				System.out.println("numread = " + numRead + " num sorted = " + numMatch + " num multiple " + numMultiple);
			}
			
			//write to file
			//Forward read
			BufferedWriter out = sToFile.get(samp)[0];
			out.write(headF + "\n" + readF + "\n"+ thirdF + "\n" + fourthF + "\n");
			//Reverse read
			out = sToFile.get(samp)[1];
			out.write(headR + "\n" + readR + "\n" + thirdR + "\n" + fourthR +"\n");
			
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
	
	private static class Holder 
	{
		String seqString;
		String qualScoreString;
	}
	
	private static Holder trim(String seqString, String qualScoreString, String primer)
		throws Exception
	{
		int index = seqString.indexOf(primer);
		
		if( seqString.length() != qualScoreString.length())
			throw new Exception("Parsing error");
		
		if( index != -1)
		{
			int endIndex = index + primer.length();
			endIndex = Math.min(endIndex, seqString.length() -1);
			seqString = seqString.substring(0, index) + seqString.substring(endIndex);
			qualScoreString = qualScoreString.substring(0,index) + qualScoreString.substring(endIndex);
		}
		
		if( seqString.length() != qualScoreString.length())
			throw new Exception("logic error");
		
		Holder h =new Holder();
		h.seqString = seqString;
		h.qualScoreString = qualScoreString;
		return h;
	}
}
