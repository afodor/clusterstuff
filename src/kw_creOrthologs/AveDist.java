/*
 * Calculate the average distance between all genomes, or within carolina,
 * for each gathered kmer 
 */
package kw_creOrthologs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;

public class AveDist {
	public static String DIR = "/nobackup/afodor_research/af_broad/orthologs/gatheredKmerDistanceMatrices";

	public static void main(String[] args) throws Exception {
		//set up output
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				"/nobackup/afodor_research/af_broad/orthologs/aveDist.txt")));
		out.write("orthogroup\taveAll\taveCarolina\n");
		
		//include table of all kmers
		getAve(new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allDist.txt"), 
				new File("/nobackup/afodor_research/af_broad/gatheredKmerMatrices/allKey.txt"),
				out);
		
		File[] tables = new File(DIR).listFiles();
		for(File t : tables) {
			if(t.getName().endsWith("dist.txt")) {
				String name = t.getAbsolutePath();
				getAve(t, new File(name.replace("dist.txt", "Key.txt")), out);
			}
		}
		
		out.close();	
	}
	
	//for given distance file dist and its corresponding key, calculate the 
	//averages and write the results
	private static void getAve(File dist, File key, BufferedWriter out) throws Exception {
		int numGenom = 339;
		
		//get key
		String[] genomes = new String[numGenom];
		BufferedReader k = new BufferedReader(new FileReader(key));
		String line = k.readLine(); //header
		line = k.readLine();
		int count = 0;
		int numCar = 0;
		while(line != null) {
			String[] sp = line.split(" ");
			if(sp.length == 3) {
				String name = sp[1];
				genomes[count] = name;
				count++;
				if(name.contains("chs")) {
					numCar++;
				}
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
		if(numCar != 76) {
			throw new Exception("Missing Carolina: " + key.getName());
		}
		
		//get table
		BufferedReader d = new BufferedReader(new FileReader(dist));
		line = d.readLine();
		if(Integer.parseInt(line) != numGenom) {
			d.close();
			throw new Exception("Missing genomes in distance: " + dist.getName());
		}
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
		for(int r = 0; r < numGenom; r++) {//row
			for(int c = 0; c < r; c++) {//column
				allSum += Double.parseDouble(table[r][c]);
				if(genomes[r].contains("chs") && genomes[c].contains("chs")) {
					carSum += Double.parseDouble(table[r][c]);
				}
			}
		}
		out.write((allSum/numGenom) + "\t" + (carSum/numCar) + "\n");
	}
}
