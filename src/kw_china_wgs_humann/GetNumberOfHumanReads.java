/*
 * outputs table of number of reads that aligned to Hg38 (number of human reads)
 * for each China WGS sample
 */
package kw_china_wgs_humann;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class GetNumberOfHumanReads {
	public static String HGDIR = "/nobackup/afodor_research/kwinglee/china/wgs/alignToHG38/";
	
	public static void main(String[] args) throws IOException {
		String[] files = new File(HGDIR).list();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				HGDIR + "numberHumanReadsPerSample.txt")));
		out.write("sample\tnumberReadsAlignedToHg38\n");
		int numSamps = 0;
		int totReads = 0;
		for(String f : files) {
			if(f.endsWith(".hg38.mapped.sam")) {
				numSamps++;
				int count = 0;
				BufferedReader br = new BufferedReader(new FileReader(new File(
						HGDIR + f)));
				String line = br.readLine();
				while(line != null) {
					count++;
					line = br.readLine();
				}
				br.close();
				out.write(f.replace("_1.hg38.mapped.sam", "") + 
						"\t" + count + "\n");
				totReads+=count;
			}
		}
		out.write("average\t" + ((double)totReads / numSamps) + "\n");
		out.close();
	}

}
