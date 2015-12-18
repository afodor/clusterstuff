/**
 * Generate table of number of reads and average read length in each fastq file for supplemental
 * 12/18/15
 */


package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class getNumberOfReads {
	public static void main(String[] args) throws FileNotFoundException, IOException {
		File fastaDir = new File("/nobackup/afodor_research/ChinaSequences/rdpResults");
		File[] fastas = fastaDir.listFiles();
		BufferedWriter out = new BufferedWriter(new FileWriter(new File("/nobackup/afodor_research/kwinglee/china/readStats.txt")));
		out.write("FileName\tNumberReads\tAverageReadLength\n");
		for(File f : fastas) {
			if(f.getName().endsWith(".fasta.gz")) {
				BufferedReader in =new BufferedReader( 
						new InputStreamReader(new GZIPInputStream(new FileInputStream(f))));
				String line = in.readLine();
				int numReads = 0;
				double lenReads = 0;
				while(line != null) {
					numReads++;
					String read = in.readLine();
					lenReads += read.length();
					line = in.readLine();
				}
				out.write(f.getName() + "\t" + numReads + "\t" + (lenReads/numReads) + "\n");
				in.close();
			}
		}
		out.close();
	}

}
