/*
 * make table of best (highest bit score) hit for each consensus sequence
 * from blasting abundant otu consensus sequences against SILVA database
 * table contains hit, bit score and percent identity
 */
package kw_jobinDolphin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class ParseSilvaBlastResults {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/abunOTU/";
	
	public static void main(String[] args) throws IOException {
		//get best hit
		HashMap<String, Integer> bitScore = new HashMap<String, Integer>();//map of consensus to bit score
		HashMap<String, String> hits = new HashMap<String, String>();//map of consensus to hit and percent identity
		BufferedReader blast = new BufferedReader(new FileReader(new File(
				DIR + "blastSilva_v_dolphinAbundantOTUcons.txt")));
		for(String line = blast.readLine(); line != null; line = blast.readLine()) {
			if(!line.startsWith("#")) {
				String[] sp = line.split("\t");
				String cons = sp[0];
				Integer bit = Integer.parseInt(sp[11]);
				String hit = sp[1] + "\t" + sp[2];
				//add to map only if consensus sequence not in map or this bit score is higher than previously seen
				if(!bitScore.containsKey(cons) || bit > bitScore.get(cons)) {
					bitScore.put(cons, bit);
					hits.put(cons, hit);
				} 
			}
		}
		blast.close();
		
		//print results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				DIR + "blastSilva_v_dolphinAbundantOTUcons.table.txt")));
		out.write("consensusSequence\tSilvaHit\tpercentIdentity\tbitScore\n");
		ArrayList<String> cons = new ArrayList<String>(hits.keySet());
		Collections.sort(cons);
		for(String c : cons) {
			out.write(c + "\t" + hits.get(c) + "\t" + bitScore.get(c) + "\n");
		}
		out.close();
	}

}
