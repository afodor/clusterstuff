/*
 * Calculate the average distance between all genomes, or within carolina,
 * for each gathered kmer 
 * 
 * Calculate for all genomes and Kleb only
 */
package kw_creOrthologs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class AveDist {
	public static String DIR = "/nobackup/afodor_research/af_broad/orthologs/gatheredKmerDistanceMatrices";
	public static HashMap<String, String> GenToClass = new HashMap<String, String>();

	public static void main(String[] args) throws Exception {
		//get map of genome to class
		BufferedReader convert = new BufferedReader(new FileReader(new File(
				"/nobackup/afodor_research/kwinglee/cre/rbh/GenomeToClass.txt")));
		String line = convert.readLine();//header
		line = convert.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			GenToClass.put(sp[0], sp[1]);
			line = convert.readLine();
		}
		convert.close();
		
		//set up output
		BufferedWriter outAll = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/orthologs/aveDist.txt")));
		outAll.write("orthogroup\taveAll\taveCarolina\taveResistant\taveSusceptible\n");
		BufferedWriter outKleb = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/orthologs/aveKlebDist.txt")));
		outKleb.write("orthogroup\taveAllKleb\taveCarolinaKleb\taveResistantKleb\taveSusceptibleKleb\n");
		
		//include table of all kmers
		getAve(new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allDist.txt"), 
				new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allKey.txt"),
				outAll, false);
		getAve(new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allDist.txt"), 
				new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allKey.txt"),
				outKleb, true);
		
		File[] tables = new File(DIR).listFiles();
		for(File t : tables) {
			if(t.getName().endsWith("dist.txt")) {
				String name = t.getAbsolutePath();
				getAve(t, new File(name.replace("dist.txt", "Key.txt")), outAll, false);
				getAve(t, new File(name.replace("dist.txt", "Key.txt")), outKleb, true);
			}
		}
		
		outAll.close();	
		outKleb.close();
	}
	
	//for given distance file dist and its corresponding key, calculate the 
	//averages and write the results
	//if kleb is true, only look at kleb genomes
	private static void getAve(File dist, File key, BufferedWriter out, boolean kleb) throws Exception {
		int numGenom = 339;
		if(kleb) {
			numGenom = 206;
		}
		
		//get key
		String[] genomes = new String[numGenom];
		HashSet<Integer> klebs = new HashSet<Integer>();//set of indices of which genomes are kleb
		BufferedReader k = new BufferedReader(new FileReader(key));
		String line;
		if(!key.getName().equals("allKey.txt")) {
			line = k.readLine(); //header
		}
		line = k.readLine();
		int count = 0;
		//int numCar = 0;
		while(line != null) {
			String[] sp = line.split("\\s");
			if(sp.length == 3 || 
					(key.getName().equals("allKey.txt") && sp.length==2)) {
				String name = sp[1];
				if(!kleb || name.contains("kleb")) {
					if(name.contains("kleb")) {
						klebs.add(count);
					}
					genomes[count] = name;
					count++;
				}
				/*if(name.contains("chs")) {
					numCar++;
				}*/
			} else {
				k.close();
				throw new Exception("Incorrect split length: " + key.getName() + 
						" " + sp.length + " " + line);
			}
			line = k.readLine();
		}
		k.close();
		if(count != numGenom) {
			throw new Exception("Incorrect number of genomes in: " + key.getName());
		}
		/*if(numCar != 76) {
			throw new Exception("Missing Carolina: " + key.getName());
		}*/
		
		//get table
		BufferedReader d = new BufferedReader(new FileReader(dist));
		line = d.readLine();
		line = d.readLine();
		String[][] table = new String[numGenom][];
		count = 0;
		while(line != null) {
			String[] sp = line.split(" ");
			if(sp.length != numGenom + 1) {
				d.close();
				throw new Exception("Incorrect length in distance: " + dist.getName() +
						" " + sp.length + " " + sp[0]);
			}
			table[count] = Arrays.copyOfRange(sp, 1, sp.length);
			count++;
			line = d.readLine();
		}
		d.close();
		if(count != numGenom) {
			throw new Exception("Incorrect dist: " + dist.getName());
		}
		
		//get sums
		double allSum = 0;
		double carSum = 0;
		double resSum = 0;
		double susSum = 0;
		int countAll = 0;
		int countCar = 0;
		int countRes = 0;
		int countSus = 0;
		for(int r = 0; r < numGenom; r++) {//row
			for(int c = 0; c < r; c++) {//column
				if(!kleb || (klebs.contains(r) && klebs.contains(c))) {
					allSum += Double.parseDouble(table[r][c]);
					countAll++;
					String cl1 = GenToClass.get(genomes[r]);
					String cl2 = GenToClass.get(genomes[c]);
					if(cl1 == null) {
						System.err.println(genomes[r]);
					}
					if(cl2 == null) {
						System.err.println(genomes[c]);
					}
					if(cl1.equals("carolina") && cl2.equals("carolina")) {
					//if(genomes[r].contains("chs") && genomes[c].contains("chs")) {
						carSum += Double.parseDouble(table[r][c]);
						countCar++;
					} else if(cl1.equals("resistant") && cl2.equals("resistant")) {
						resSum += Double.parseDouble(table[r][c]);
						countRes++;
					} else if(cl1.equals("susceptible") && cl2.equals("susceptible")) {
						susSum += Double.parseDouble(table[r][c]);
						countSus++;
					}
				}
			}
		}
		String name = dist.getName();
		System.out.println(name + "\t" + allSum + "\t" + carSum);
		out.write(name.replace(".fasta_dist.txt", "") + "\t" 
				+ (allSum/countAll) + "\t" + (carSum/countCar) + "\t"
				+ (resSum/countRes) + "\t" + (susSum/countSus) + "\n");
		out.flush();
	}
}
