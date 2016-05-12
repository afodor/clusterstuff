/*
 * compare metaphlan runs with and without filtering human reads
 */
package kw_china_wgs_metaphlan;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class CompareMetaphlanRuns {
	public static String DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	
	public static void main(String[] args) throws IOException {
		String filterDir = DIR + "metaphlanResultsFilterHuman/";
		String otherDir = DIR + "metaphlanResults/";
		String[] files = new File(filterDir).list();
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				DIR + "metScripts/compare.sh")));
		for(String f : files) {
			if(f.startsWith("metaphlan_table_filterHuman")) {
				String sample = f.replace("metaphlan_table_filterHuman_", "");
				script.write("echo " + sample + "\n");
				script.write("diff " + filterDir + "metaphlan_table_filterHuman_"
						+ sample + " " + otherDir + "metaphlan_table_" + sample + "\n");
			}
		}
		script.close();
	}

}
