/*
 * Calculate the mean distance between all genomes, or within carolina,
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class MedianDist {
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
				"/nobackup/afodor_research/af_broad/orthologs/medianDist.txt")));
		outAll.write("orthogroup\tmeanAll\tmeanCarolina\tmeanResistant\tmeanSusceptible\n");
		BufferedWriter outKleb = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/orthologs/medianKlebDist.txt")));
		outKleb.write("orthogroup\tmeanAllKleb\tmeanCarolinaKleb\tmeanResistantKleb\tmeanSusceptibleKleb\n");
		
		//include table of all kmers
		getMedian(new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allDist.txt"), 
				new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allKey.txt"),
				outAll, false);
		getMedian(new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allDist.txt"), 
				new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allKey.txt"),
				outKleb, true);
		
		File[] tables = new File(DIR).listFiles();
		for(File t : tables) {
			if(t.getName().endsWith("dist.txt")) {
				String name = t.getAbsolutePath();
				getMedian(t, new File(name.replace("dist.txt", "Key.txt")), outAll, false);
				getMedian(t, new File(name.replace("dist.txt", "Key.txt")), outKleb, true);
			}
		}
		
		outAll.close();	
		outKleb.close();
	}
	
	//for given distance file dist and its corresponding key, calculate the 
	//means and write the results
	//if kleb is true, only look at kleb genomes
	private static void getMedian(File dist, File key, BufferedWriter out, boolean kleb) throws Exception {
		int numAllGenom = 339;
		
		//get key
		String[] genomes = new String[numAllGenom];
		HashSet<Integer> klebs = new HashSet<Integer>();//set of indices of which genomes are kleb
		BufferedReader k = new BufferedReader(new FileReader(key));
		String line;
		if(!key.getName().equals("allKey.txt")) {
			line = k.readLine(); //header
		}
		line = k.readLine();
		int count = 0;
		while(line != null) {
			String[] sp = line.split("\\s");
			if(sp.length == 3 || 
					(key.getName().equals("allKey.txt") && sp.length==2)) {
				String name = sp[1];
				if(!kleb || name.contains("kleb")) {
					if(name.contains("kleb")) {
						klebs.add(count);
					}
				}
				genomes[count] = name;
				count++;
			} else {
				k.close();
				throw new Exception("Incorrect split length: " + key.getName() + 
						" " + sp.length + " " + line);
			}
			line = k.readLine();
		}
		k.close();
		if(count != numAllGenom) {
			throw new Exception("Incorrect number of genomes in: " + key.getName());
		}
		/*if(numCar != 76) {
			throw new Exception("Missing Carolina: " + key.getName());
		}*/
		
		//get table
		BufferedReader d = new BufferedReader(new FileReader(dist));
		line = d.readLine();
		if(Integer.parseInt(line) != numAllGenom) {
			d.close();
			throw new Exception("Missing genomes in distance: " + dist.getName());
		}
		line = d.readLine();
		String[][] table = new String[numAllGenom][];
		count = 0;
		while(line != null) {
			String[] sp = line.split(" ");
			if(sp.length != numAllGenom + 1) {
				d.close();
				throw new Exception("Incorrect length in distance: " + dist.getName() +
						" " + sp.length + " " + sp[0]);
			}
			table[count] = Arrays.copyOfRange(sp, 1, sp.length);
			count++;
			line = d.readLine();
		}
		d.close();
		if(count != numAllGenom) {
			throw new Exception("Incorrect dist: " + dist.getName());
		}
		
		//get list of all distances
		List<Double> allSum = new ArrayList<Double>();
		List<Double>  carSum = new ArrayList<Double>();
		List<Double>  resSum = new ArrayList<Double>();
		List<Double>  susSum = new ArrayList<Double>();
		for(int r = 1; r < numAllGenom; r++) {//row
			if(!kleb || klebs.contains(r)) {
				for(int c = 0; c < r; c++) {//column
					if(!kleb || klebs.contains(c)) {
						allSum.add(Double.parseDouble(table[r][c]));
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
							carSum.add(Double.parseDouble(table[r][c]));
						} else if(cl1.equals("resistant") && cl2.equals("resistant")) {
							resSum.add(Double.parseDouble(table[r][c]));
						} else if(cl1.equals("susceptible") && cl2.equals("susceptible")) {
							susSum.add(Double.parseDouble(table[r][c]));
						}
					}
				}
			}
		}
		
		//write results
		String name = dist.getName();
		System.out.println(name + "\t" + allSum + "\t" + carSum);
		out.write(name.replace(".fasta_dist.txt", "") + "\t" 
				+ calculateMedian(allSum) + "\t" + calculateMedian(carSum) + "\t"
				+ calculateMedian(resSum) + "\t" + calculateMedian(susSum) + "\n");
		out.flush();
	}
	
	public static double calculateMedian(List<Double> list) {
		//sort
		Collections.sort(list);
		
		//return median
		double median;
		int mid = list.size()/2;
		if(list.size() % 2 == 0) {
			median = (list.get(mid) + list.get(mid-1))/2;
		} else {
			median = list.get(mid);
		}
		return(median);
	}
}
