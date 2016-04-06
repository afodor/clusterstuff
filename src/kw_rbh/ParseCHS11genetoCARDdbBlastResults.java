/*
 * parse results of blasting card database to chs11 genes
 * only include hits that are greater than 50% of the length of both the CHS11 and 
 * CARD gene
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
import java.util.Iterator;
import java.util.Set;


public class ParseCHS11genetoCARDdbBlastResults {
	private static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/chs11_v_card/";
	private static String ARO_KEY = "/users/kwinglee/card/aro.csv";
	private static String CHS11_GENE = "/nobackup/afodor_research/kwinglee/cre/rbh/carolina_klebsiella_pneumoniae_chs_11.0_genePositions.txt";
	private static HashMap<String, String[]> ARO;//ARO to [name, description]
	private static HashMap<String, String[]> CHS11;//gene name to [length, scaffold, start, stop]
	
	public static void main(String[] args) throws Exception {
		//get map of gene to gene length for CHS11
		CHS11 = new HashMap<String, String[]>();
		//map will be gene name -> [length, scaffold, start, stop]
		BufferedReader chs11gene = new BufferedReader(new FileReader(new File(
				CHS11_GENE)));
		String line = chs11gene.readLine();//header
		line = chs11gene.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			String scaff = sp[1];
			int start = Integer.parseInt(sp[2]);
			int stop = Integer.parseInt(sp[3]);
			CHS11.put(sp[0], new String[]{Integer.toString(stop-start), 
					scaff, Integer.toString(start), Integer.toString(stop)});
			line = chs11gene.readLine();
		}
		chs11gene.close();
		
		//get map of aro to annotation
		//map will be aro number to [name, description]
		ARO = new HashMap<String, String[]>();
		BufferedReader arocsv = new BufferedReader(new FileReader(new File(
				ARO_KEY)));
		line = arocsv.readLine();//header
		line = arocsv.readLine();
		while(line != null) {
			String[] sp = line.split(",");
			ARO.put(sp[0], new String[]{sp[1], sp[2]});
			line = arocsv.readLine();
		}
		arocsv.close();
		
		
		analyzeResults("chs11genes_v_cardHomolog.txt", 
				"nucleotide_fasta.protein_homolog.fasta");
		analyzeResults("chs11genes_v_cardVariant.txt", 
				"nucleotide_fasta.protein_variant.fasta");
	}
	
	/*
	 * for the given blast results table, for each gene give the best hit
	 * and any other hits, with annotation for the best
	 * genes must match > 50% gene length
	 * database is the name of the fasta database in card folder
	 */
	public static void analyzeResults(String table, String database) throws Exception {
		//get database gene length
		HashMap<String, Integer> dbLen = new HashMap<String, Integer>();//map of aro number to length
		BufferedReader dbbr = new BufferedReader(new FileReader(new File(
				"/users/kwinglee/card/database")));
		String line = dbbr.readLine();
		int len = 0;
		String name = null;
		while(line != null) {
			if(line.startsWith(">")) {
				if(name != null) {//put in map if not the first line
					dbLen.put(name, len);
				}
				//get aro as name
				name = line.split("\\|")[3];
				if(!name.startsWith("ARO")) {
					dbbr.close();
					throw new Exception("database ARO number wrong: " + name + "\n" + line);
				}
				len = 0;
			} else {
				len += line.length();
			}
			line = dbbr.readLine();
		}
		dbbr.close();
		
		//get gene to hit
		//map is chs11 gene name to set of card hits (represented as aro;name;bitscore)
		HashMap<String, Set<String>> hits = new HashMap<String, Set<String>>();
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + table)));
		line = br.readLine();
		while(line != null) {
			if(!line.startsWith("#")) {
				String[] sp = line.split("\t");
				if(sp.length != 11) {
					br.close();
					throw new Exception("results wrong length: " + sp.length + "\n" + line);
				}
				String chs = sp[0];
				String[] sp2 = sp[1].split("\\|");
				String aro = sp2[3];//aro
				String shortName = sp2[4];
				String bit = sp[11];
				if(!aro.startsWith("ARO")) {
					br.close();
					throw new Exception("results ARO number wrong: " + name + "\n" + line);
				}
				int alLen = Integer.parseInt(sp[2]);//alignment length
				if(alLen > Integer.parseInt(CHS11.get(chs)[0]) / 2 &&
						alLen > dbLen.get(aro)) {//alignment length is more than 50% of gene lengths
					String value = aro + ";" + shortName + ";" + bit;
					if(hits.containsKey(chs)) {
						hits.get(chs).add(value);
					} else {
						Set<String> set = new HashSet<String>();
						set.add(value);
						hits.put(chs, set);
					}
				}
			}
			line = br.readLine();
		}
		br.close();

		//write
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + table.replace(".txt", "_summary.txt"))));
		//header
		out.write("CHS11 gene\tscaffold\tgene start\tgene stop\t"
				+ "best CARD hit ARO\tARO name\tARO description\t"
				+ "best hit bit score\tother hits (ARO;ARO name;bit score)\n");
		String[] genes = (String[]) CHS11.keySet().toArray();
		Arrays.sort(genes);
		for(String key : genes) {
			if(hits.containsKey(key)) {
				Set<String> set = hits.get(key);
				Iterator<String> it = set.iterator();
				if(set.size() == 1) {
					out.write(parseOutput(key, it.next()) + "NA\n");
				} else {
					//identify best hit
					String best = "";
					int max = 0;//highest bit score
					while(it.hasNext()) {
						String h = it.next();
						String[] sp = h.split(";");
						if(Integer.parseInt(sp[2]) > max) {
							max = Integer.parseInt(sp[2]);
							best = h;
						}
					}
					//write best hit
					out.write(parseOutput(key, best));
					//write other hits
					it = set.iterator();
					while(it.hasNext()) {
						String h = it.next();
						if(!h.equals(best)) {
							out.write(h + "|");
						}
					}
					out.write("\n");
				}
			}
		}
		out.close();
	}
	
	/*
	 * for the given CHS gene key and its best CARD match hit,
	 * return the output for all columns except the other hits
	 */
	public static String parseOutput(String key, String hit) {
		String[] chs = CHS11.get(key);
		String[] sp = hit.split(";");
		String a = sp[0];
		String[] aro = ARO.get(a);
		if(!sp[1].equals(aro[0])) {//see if name is same from fasta and csv file
			System.out.println("Different names for " + a + ": " + 
					sp[1] + "\t" + aro[0]);
		}
		return(key + "\t" + chs[1] + "\t" + chs[2] + "\t" + chs[3] + "\t" + 
				a + "\t" + aro[0] + "\t" + aro[1] + "\t" + sp[2] + "\t");
	}

}
