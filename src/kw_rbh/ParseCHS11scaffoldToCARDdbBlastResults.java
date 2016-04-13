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
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;


public class ParseCHS11scaffoldToCARDdbBlastResults {
	private static String DIR = "/nobackup/afodor_research/kwinglee/cre/rbh/chs11_v_card/";
	private static String ARO_KEY = "/users/kwinglee/card/aro.csv";
	private static HashMap<String, String[]> ARO;//ARO to [name, description]
	private static Double LENDIFF = 0.5;//percent of gene that must align
	
	public static void main(String[] args) throws Exception {
		//get map of aro to annotation
		//map will be aro number to [name, description]
		ARO = new HashMap<String, String[]>();
		BufferedReader arocsv = new BufferedReader(new FileReader(new File(
				ARO_KEY)));
		String line = arocsv.readLine();//header
		line = arocsv.readLine();
		while(line != null) {
			if(line.startsWith("ARO")) {//there are empty rows and rows start with Gene orientation that throw errors
				String[] sp = line.split(",");
				if(sp.length == 3) {
					ARO.put(sp[0], new String[]{sp[1], sp[2]});		
				} else {
					line = line.replace(sp[0] + ",", "");
					String[] sp2 = line.split("\"");
					if(sp2.length == 2) {//extra commas in description
						line = line.replace(sp2[0], "");
						ARO.put(sp[0], new String[]{sp2[0].replace(",", ""), line});
					} else {//some lines don't seem to have a name
						ARO.put(sp[0], new String[]{"", line});	
					}
				}
			}
			line = arocsv.readLine();
		}
		arocsv.close();
		
		
		analyzeResults("chs11scaff_v_cardHomolog.txt", 
				"nucleotide_fasta.protein_homolog.fasta");
		analyzeResults("chs11scaff_v_cardVariant.txt", 
				"nucleotide_fasta.protein_variant.fasta");
	}
	
	/*
	 * for the given blast results table, get position of hit
	 * genes must match > 50% gene length
	 * database is the name of the fasta database in card folder
	 */
	public static void analyzeResults(String table, String database) throws Exception {
		//get database gene length
		HashMap<String, Integer> dbLen = new HashMap<String, Integer>();//map of aro number to length
		BufferedReader dbbr = new BufferedReader(new FileReader(new File(
				"/users/kwinglee/card/" + database)));
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
		//map is scaffold to set of card hits (represented as aro;name;bitscore;start;stop)
		HashMap<String, Set<String>> hits = new HashMap<String, Set<String>>();
		BufferedReader br = new BufferedReader(new FileReader(new File(
				DIR + table)));
		line = br.readLine();
		while(line != null) {
			if(!line.startsWith("#")) {
				String[] sp = line.split("\t");
				if(sp.length != 12) {
					br.close();
					throw new Exception("results wrong length: " + sp.length + "\n" + line);
				}
				String scaff = sp[0];
				String[] sp2 = sp[1].split("\\|");
				String aro = sp2[3];//aro
				String shortName = sp2[4];
				int start = Integer.parseInt(sp[6]);
				int stop = Integer.parseInt(sp[7]);
				if(stop < start) {//reversed match
					int temp = stop;
					start = stop;
					stop = temp;
				}
				String bit = sp[11].replaceAll("\\s+","");//remove extra white space ex for carolina_klebsiella_pneumoniae_chs_11.0_AE67_00312      ARO:3001190
				if(!aro.startsWith("ARO")) {
					br.close();
					throw new Exception("results ARO number wrong: " + name + "\n" + line);
				}
				double alLen = Double.parseDouble(sp[3]);//alignment length
				if(alLen / dbLen.get(aro) > LENDIFF) {//alignment length is more than 50% of gene length
					String value = aro + ";" + shortName + ";" + bit + ";" + start
							+ ";" + stop;
					if(hits.containsKey(scaff)) {
						hits.get(scaff).add(value);
					} else {
						Set<String> set = new HashSet<String>();
						set.add(value);
						hits.put(scaff, set);
					}
				}
			}
			line = br.readLine();
		}
		br.close();
		System.out.println("1000 bp postions with hits: " + hits.size() + " for " + table);

		//write
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + table.replace(".txt", "_summary.txt"))));
		//header
		out.write("scaffold\thit start\thit stop\t"
				+ "best CARD hit ARO\tARO name\tARO description\t"
				+ "best hit bit score\tother hits (ARO;ARO name;bit score)\n");
		String[] hitSet = hits.keySet().toArray(new String[hits.keySet().size()]);
		Arrays.sort(hitSet);
		for(String key : hitSet) {
			//write best hit if among overlapping hits
			Set<String> set = hits.get(key);
			checkOverlappingSet(set, key, out);
		}
		out.close();
	}
	
	/*
	 * for the given set, write best hit and any overlapping genes
	 * repeat for genes not overlapping the best hit
	 */
	public static void checkOverlappingSet(Set<String> set, 
			String key, BufferedWriter out) throws IOException {
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
			
			//write other hits if they are overlapping; otherwise make new set
			it = set.iterator();
			Set<String> nonOverlap = new HashSet<String>();
			String[] b = best.split(";");
			int bstart = Integer.parseInt(b[3]);
			int bstop = Integer.parseInt(b[4]);
			while(it.hasNext()) {
				String h = it.next();
				if(!h.equals(best)) {
					String[] sp = h.split(";");
					int hstart = Integer.parseInt(sp[3]);
					int hstop = Integer.parseInt(sp[4]);
					if((hstart <= bstart & hstop >= bstart) ||
							(hstart <= bstop & hstop >= bstop)) {//overlapping
						out.write(h + "|");
					} else {
						nonOverlap.add(h);
					}
				}
			}
			out.write("\n");
			if(!nonOverlap.isEmpty()) {
				checkOverlappingSet(nonOverlap, key, out);
			}
		} 
	}
	
	/*
	 * for the given CHS gene key and its best CARD match hit,
	 * return the output for all columns except the other hits
	 */
	public static String parseOutput(String key, String hit) {
		String[] sp = hit.split(";");
		String a = sp[0];
		String[] aro = ARO.get(a);
		if(!sp[1].equals(aro[0])) {//see if name is same from fasta and csv file
			System.out.println("Different names for " + a + ": " + 
					sp[1] + "\t" + aro[0]);
		}
		return(key + "\t" + sp[3] + "\t" + sp[4] + "\t" + 
				a + "\t" + aro[0] + "\t" + aro[1] + "\t" + sp[2] + "\t");
	}

}
