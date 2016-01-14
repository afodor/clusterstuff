/**
 * For each Broad genome, generate a fasta file of all genes from the gtf file
 */
package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;

public class geneFastas {
	public static String GenomeTopDir = "/nobackup/afodor_research/af_broad/";

	public static void main(String[] args) throws IOException {
		analyzeFolder("carolina");
		analyzeFolder("susceptible");
		analyzeFolder("resistant");
	}
	
	//analyze the genomes in the given folder
	public static void analyzeFolder(String folder) throws IOException {
		System.out.println("ANALYZING: " + folder);
		String outDir = "/nobackup/afodor_research/kwinglee/cre/rbh/" + folder + "/";
		String dirName = GenomeTopDir + folder;
		File dir = new File(dirName);
		File[] allFiles = dir.listFiles();
		for(File gtf : allFiles) {
			if(gtf.getName().endsWith(".genes.gtf")) {
				String name = gtf.getName().replace(".genes.gtf", "");//genome to analyze
				System.out.println(name);
				
				//make directory for output files
				File gtfDir = new File(outDir + name);
				gtfDir.mkdirs();
				
				//get corresponding fasta file as hash map of scaffold name to sequence
				HashMap<String, String> scaff = getScaffolds(new File(dirName + "/" + name + ".scaffolds.fasta"));
				
				//make gene files
				gtfToGene(gtf, folder + "_" + name, gtfDir.getAbsolutePath() + "/", scaff);
			}
		}
	}
	
	//takes the given fasta file, reads it and returns a hash map of scaffold name to sequence
	public static HashMap<String, String> getScaffolds(File fasta) throws IOException {
		HashMap<String, String> map = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fasta)));
		String line = br.readLine();
		String header = line.replaceFirst(">", "");
		String seq = "";
		while(line != null) {
			if(line.startsWith(">")) { //header/start of new sequence
				if(seq.length() > 1) { //add to map unless its the beginning of the file
					map.put(header, seq);
				}
				header = line.replaceFirst(">", "");
				seq = "";
			} else {
				seq += line;
			}
			line = br.readLine();
		}
		br.close();
		map.put(header, seq);//last sequence
		return(map);
	}
	
	//takes the given gtf file, determines the gene sequence from scaff 
	//and writes the result as a fasta in gtfDir
	public static void gtfToGene(File gtf, String name, String gtfDir, HashMap<String, String> scaff) throws IOException {
		BufferedWriter allGenes = new BufferedWriter(new FileWriter(
				new File(gtfDir + name + "_allGenes.fasta")));
		BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(gtf)));
		String line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			if(sp[2].equals("exon")) { //only include exons (which include start and stop codon)
				//get sequence
				String scaffSeq = scaff.get(sp[0]);
				String seq = scaffSeq.substring(
						Integer.parseInt(sp[3])-1, Integer.parseInt(sp[4]));//minus 1 to go from ones to zero based counting; end should include the last base so don't substract
				if(sp[6].equals("-")) {//reverse orientation
					seq = revComp(seq);
				}
				
				//use gene_id
				String id = sp[8].split(";")[0].split("\"")[1];

				//write as in file with all genes
				allGenes.write(">" + id + "\n" + seq + "\n");
				
				//write as separate fasta file
				BufferedWriter gene = new BufferedWriter(new FileWriter(
						new File(gtfDir + name + "_" + id +".fasta")));
				gene.write(">" + id + "\n" + seq + "\n");
				gene.close();
			}
			line = br.readLine();
		}
		allGenes.close();
		br.close();
	}
	
	//returns the reverse complement of the given sequence seq
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
			} else if(fwd[i] == 'N') {
				rev += "N";
			} else {
				System.out.println(seq);
				System.out.println(fwd[i]);
				throw new IllegalArgumentException("Invalid base: " + fwd[i]);
			}
		}
		return(rev);
	}
}
