/*
 * get table of genome and genome stats from NCBI blast
 */
package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ParseNCBIresults {
	private static String BLAST_RESULTS = "/nobackup/afodor_research/ChinaSequences/ncbi/allForwardToNCBI_Oct15.txt";
	private static String GENE_TO_GENOME = "/nobackup/afodor_research/ncbi/geneIDToGenome.txt";
	private static String NCBI_STATS = "/nobackup/afodor_research/kwinglee/china/ncbiSizes.txt";
	private static String OUTFILE = "/nobackup/afodor_research/kwinglee/china/ncbiBlastResultsWithGeneStats.txt";
	
	public static void main(String[] args) throws IOException {
		//get best blast hit, based on best bit score
		HashMap<String, String[]> blastHits = new HashMap<String, String[]>();//map of consensus to reference + bit score
		BufferedReader blast = new BufferedReader(new FileReader(new File(BLAST_RESULTS)));
		for(String line = blast.readLine(); line != null; line = blast.readLine()) {
			String[] sp = line.split("\t");
			String con = sp[0];
			String ref = sp[1];
			String bit = sp[11];
			
			if(blastHits.containsKey(con)) {
				String[] oldHit = blastHits.get(con);
				Double newBit = Double.parseDouble(bit);
				Double oldBit = Double.parseDouble(oldHit[1]);
				if(newBit > oldBit) {
					blastHits.put(con, new String[]{ref, bit});
				}
			} else {
				blastHits.put(con, new String[]{ref, bit});
			}
		}
		blast.close();
		
		//get map of reference gene ID to genome name
		HashMap<String, String> gene2genome = new HashMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(new File(GENE_TO_GENOME)));
		br.readLine();
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			String[] sp = line.split("\t");
			gene2genome.put(sp[0].replace(">>", "").split(" ")[0], sp[1]);
		}
		br.close();
		
		//get map of genome to stats
		HashMap<String, String> stats = new HashMap<String, String>();//genome to number genes + gene length
		br = new BufferedReader(new FileReader(new File(NCBI_STATS)));
		br.readLine();//header
		for(String line = br.readLine(); line != null; line = br.readLine()) {
			String[] sp = line.split("\t");
			stats.put(sp[0], sp[1] + "\t" + sp[2]);
		}
		br.close();
		
		//write results
		ArrayList<String> otus = new ArrayList<String>();
		otus.addAll(blastHits.keySet());
		Collections.sort(otus);
		
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(OUTFILE)));
		out.write("OTU\t16Sgene\tgenome\tnumberGenes\ttotalGeneLength\n");
		for(String otu : otus) {
			String gene = blastHits.get(otu)[0];
			String genome = gene2genome.get(gene);
			out.write(otu + "\t" + gene + "\t" + genome 
					+ "\t" + stats.get(genome) + "\n");
		}
		out.close();
	}

}
