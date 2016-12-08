/*
 * Check the stool MD5 for HMP
 */
package kw_machineLearning;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class CheckHMPmd5 {
	public static String DIR = "/nobackup/afodor_research/kwinglee/machineLearning/hmp/";

	public static void main(String[] args) throws IOException {
		//get map of ID to HMP MD5
		HashMap<String, String> hmp = new HashMap<String, String>();
		BufferedReader list = new BufferedReader(new FileReader(new File
				(DIR + "HMIWGS_healthy.csv")));
		for(String line = list.readLine(); line != null; line = list.readLine()) {
			line = line.replaceAll("\"", "");
			String[] sp = line.split(",");
			if(sp[1].equals("stool")) {
				hmp.put(sp[0], sp[3]);
			}
		}
		list.close();
		
		//get map of ID to download MD5
		HashMap<String, String> download = new HashMap<String, String>();
		BufferedReader file = new BufferedReader(new FileReader(new File
				(DIR + "fastqs/stoolMD5")));
		for(String line = file.readLine(); line != null; line = file.readLine()) {
			String[] sp = line.split("  ");
			line.contains(" ");
			String id = sp[1].replace(".tar.bz2", "");
			download.put(id, sp[0]);
		}
		file.close();
		
		//check same set of IDs
		if(hmp.size() != download.size()) {
			System.err.println("Different sizes: " + hmp.size() + ' ' + download.size());
		} else {
			Set<String> hkeys = new HashSet<String>();
			hkeys.addAll(hmp.keySet());
			Set<String> dkeys = new HashSet<String>();
			dkeys.addAll(download.keySet());
			hkeys.removeAll(dkeys);
			if(hkeys.size() != 0) {
				System.err.println("Extra HMP keys: " + hkeys.size());
			}
			hkeys.addAll(hmp.keySet());
			dkeys.addAll(download.keySet());
			dkeys.removeAll(hkeys);
			if(dkeys.size() != 0) {
				System.err.println("Extra download keys: " + dkeys.size());
			}
		}
		
		//check MD5s
		Set<String> keys = hmp.keySet();
		for(String k : keys) {
			if(!hmp.get(k).equals(download.get(k))) {
				System.err.println("Bad download: " + k);
			}
		}
	}
}
