/**
 * shotmap gives an error about broken gzip pipe
 * -> check number of lines in raw folder is correct
 * 11/18/15
 */
package kw_china_wgs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class testShotmapSplit {
	//test 081A
	private static String dir = "/projects/afodor_research/kwinglee/china/wgs/shotmapSfamResults/081A_1/081A_1/raw";
	
	public static void main(String[] args) throws IOException {
		File[] fastas = new File(dir).listFiles();
		int total = 0;
		for(File f: fastas) {
			int fileTotal = 0;
			BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
			String line = br.readLine();
			while(line != null) {
				fileTotal++;
				total++;
				line = br.readLine();
			}
			br.close();
			System.out.println(f.getName() + "\t" + fileTotal);
		}
		System.out.println("total = " + total);
	}

}
