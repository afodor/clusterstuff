/*
 * virulence database annotations were lost during analysis -> produce list of genes
 * in each database
 */

package kw_china_wgs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;

public class GetVirulenceAnnotation {
	public static final String DIR = "/nobackup/afodor_research/kwinglee/software/virulence/";
	
	public static void main(String[] args) throws IOException {
		//mvirdb
		BufferedReader mvirdb = new BufferedReader(new FileReader(new File(
				DIR + "MvirDB/virulenceDB.nucleic.fasta")));
		BufferedWriter mout = new BufferedWriter(new FileWriter(new File(
				DIR + "MvirDB_headers.txt")));
		for(String line = mvirdb.readLine(); line != null; line = mvirdb.readLine()) {
			if(line.startsWith(">")) {
				mout.write(line.replace(">", "") + "\n");
			}
		}
		mvirdb.close();
		mout.close();
		
		////VFDB
		//full
		BufferedReader full = new BufferedReader(new FileReader(new File(
				DIR + "VFDB/VFDB_full.fas")));
		HashSet<String> vfdb = new HashSet<String>();
		for(String line = full.readLine(); line != null; line = full.readLine()) {
			if(line.startsWith(">")) {
				vfdb.add(line.replace(">", ""));
			}
		}
		System.out.println("Full size: " + vfdb.size());
		full.close();
		//check core is a subset of full
		int numMissing = 0;
		BufferedReader core = new BufferedReader(new FileReader(new File(
				DIR + "VFDB/VFDB_core.fas")));
		for(String line = core.readLine(); line != null; line = core.readLine()) {
			if(line.startsWith(">")) {
				line = line.replace(">", "");
				if(!vfdb.contains(line)) {
					numMissing++;
					vfdb.add(line);
				}
			}
		}
		core.close();
		System.out.println("Number unique to core: " + numMissing);
		//write
		BufferedWriter vout = new BufferedWriter(new FileWriter(new File(
				DIR + "VFDB_headers.txt")));
		for(String v : vfdb) {
			vout.write(v + "\n");
		}
		vout.close();
	}

}
