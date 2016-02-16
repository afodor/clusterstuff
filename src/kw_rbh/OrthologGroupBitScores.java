/*
 * Get the bit score (and gene name) tables for the ortholog groups
 * (use the list of groups with 150 or more genomes)
 * Bit score is the average bit score for all pairs
 */

package kw_rbh;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class OrthologGroupBitScores {
	public static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/";
	public static String RBHDIR = DIR + "rbhOrthologs/";
	public static int NUM_GENOME = 339;
	public static String[] GENOMES = new String[NUM_GENOME];
	
	public static void main(String[] args) throws Exception {
		//log file
		BufferedWriter log = new BufferedWriter(new FileWriter(new File(
				DIR + "orthologGroupBitScoreLog.txt")));
		log.write("get list of genomes\n");
		log.flush();
		
		//get list of genomes
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + "GenomeToClass.txt")));
		String line = br.readLine();//header
		line = br.readLine();
		int pos = 0;
		while(line != null) {
			String[] sp = line.split("\t");
			GENOMES[pos] = sp[1] + "_" + sp[0];
			pos++;
			line = br.readLine();
		}
		br.close();
		if(pos != NUM_GENOME) {
			log.close();
			throw new Exception("Incorrect number genomes " + pos + " " + NUM_GENOME);
		}
		
		//get list of groups (represented as set of members)
		log.write("Get ortholog group list\n");
		log.flush();
		HashMap<String, HashSet<String>> groupMap = new HashMap<String, HashSet<String>>();
		br = new BufferedReader(new FileReader(new File(
				RBHDIR + "orthologGroups150.txt")));
		line = br.readLine();//header
		line = br.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String[] group = sp[2].split(";");
			HashSet<String> set = new HashSet<String>(Arrays.asList(group));
			groupMap.put(sp[0], set);
			line = br.readLine();
		}
		br.close();
		
		//for each genome, get bit scores and gene name for each orthogroup
		String[] orthogroups = groupMap.keySet().toArray(new String[groupMap.keySet().size()]);
		double[][] bitScore = new double[NUM_GENOME][orthogroups.length];//for each genome in genome (first index), the sum of the bit score for that orthogroup (second index)
		//bitScore varies for each genome comparison -> take average
		int[][] numScores = new int[NUM_GENOME][orthogroups.length];//number of times have added to this value in bitScore
		String[][] gene = new String[NUM_GENOME][orthogroups.length];//for each genome in genome (first index), the gene name for that orthogroup (second index)
		String[] folders = {"carolina_v_carolina", "carolina_v_resistant",
				"carolina_v_susceptible", "resistant_v_carolina", "resistant_v_resistant",
				"resistant_v_susceptible", "susceptible_v_carolina", 
				"susceptible_v_resistant", "susceptible_v_susceptible"};
		for(String f : folders) {
			log.write("Starting folder " + f + "\n");
			log.flush();
			String path = RBHDIR + f;
			String[] genList = new File(path).list();
			for(String g : genList) {//each genome folder in folder
				String genPath = path + "/" + g;
				//System.out.println(path);
				String[] tableList = new File(genPath).list();
				log.write("Starting genome " + genPath + " " + tableList.length + "\n");
				log.flush();
				for(String table : tableList) {//each table
					log.write("table " + table + "\n");
					log.flush();
					br = new BufferedReader(new FileReader(new File(
							path + "/" + table)));
					line = br.readLine();//header
					line = br.readLine();
					while(line != null) {
						String[] sp = line.split("\t");
						String gene1 = sp[0];
						String gene2 = sp[1];
						double bit1 = Double.parseDouble(sp[2]);
						double bit2 = Double.parseDouble(sp[3]);
						if(bit1 != bit2) {
							br.close();
							log.close();
							throw new Exception("Unequal bit scores: " + line);
						}
						
						for(int i = 0; i < orthogroups.length; i++) { //for each line, see if those genes are in a particular orthogroup
							HashSet<String> set = groupMap.get(orthogroups[i]);
							if(set.contains(gene1) && set.contains(gene2)) {
								String[] genomeNames = table.split("_v_");//table name is genome1_v_genome2.txt
								int gen1 = getGenomeNumber(genomeNames[0].replace("rbhResults_", ""));
								int gen2 = getGenomeNumber(genomeNames[1].replace(".txt", ""));
								/*if(bitScore[gen1][i] != 0 && bitScore[gen1][i] != bit1) {
									br.close();
									throw new Exception("Inconsistent bit score " + bitScore[gen1][i] +
											" " + bit1 + " " + gene1 + " " + orthogroups[i]);
								}
								if(bitScore[gen2][i] != 0 && bitScore[gen2][i] != bit1) {
									br.close();
									throw new Exception("Inconsistent bit score " + bitScore[gen2][i] +
											" " + bit2 + " " + gene2 + " " + orthogroups[i]);
								}*/
								if(gene[gen1][i] != null && !gene[gen1][i].equals(gene1)) {
									log.close();
									br.close();
									throw new Exception("Inconsistent gene name " + gene[gen1][i] +
											" " + gene1 + " " + orthogroups[i]);
								}
								if(gene[gen2][i] != null && !gene[gen2][i].equals(gene2)) {
									log.close();
									br.close();
									throw new Exception("Inconsistent gene name " + gene[gen2][i] +
											" " + gene2 + " " + orthogroups[i]);
								}
								bitScore[gen1][i] += bit1;
								bitScore[gen2][i] += bit1;
								numScores[gen1][i]++;
								numScores[gen2][i]++;
								gene[gen1][i] = gene1;
								gene[gen2][i] = gene2;
							}
						}
						
						line = br.readLine();
					}
					br.close();
				}
			}
		}
		
		//write results
		log.write("writing results\n");
		log.flush();
		BufferedWriter bitOut = new BufferedWriter(new FileWriter(new File(
				RBHDIR + "bitScoreTable_orthologGroups150.txt")));
		BufferedWriter geneOut = new BufferedWriter(new FileWriter(new File(
				RBHDIR + "orthologNameTable_orthologGroups150.txt")));
		//write header
		bitOut.write("genomeID");
		geneOut.write("genomeID");
		for(String o : orthogroups) {
			bitOut.write("\t" + o);
			geneOut.write("\t" + o);
		}
		bitOut.write("\n");
		geneOut.write("\n");
		//write values for each genome
		for(int i = 0; i < NUM_GENOME; i++) {
			bitOut.write(GENOMES[i]);
			geneOut.write(GENOMES[i]);
			double[] bit = bitScore[i];
			int[] count = numScores[i];
			String[] geneName = gene[i];
			for(int j = 0; j < bit.length; j++) {
				bitOut.write("\t" + (bit[j]/count[j]));
				if(geneName[j] == null) {
					geneOut.write("\tNA");
				} else {
					geneOut.write("\t" + geneName[j]);
				}
			}
			bitOut.write("\n");
			geneOut.write("\n");
		}
		
		bitOut.close();
		geneOut.close();
		log.close();
	}
	
	//return the position of the genome gen in GENOMES
	public static int getGenomeNumber(String gen) throws Exception {
		int pos = -1;
		for(int i = 0; i < GENOMES.length; i++) {
			if(gen.equals(GENOMES[i])) {
				pos = i;
				break;
			}
		}
		if(pos < 0) {
			throw new Exception("Missing genome " + gen);
		}
		return(pos);
	}
}
