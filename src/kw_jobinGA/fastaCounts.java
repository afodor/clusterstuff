/**
 * generate table of fasta file to number of reads and average read length
 */

package kw_jobinGA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class fastaCounts {
	public static final String DIR = "/projects/afodor_research/kwinglee/jobin/ga-stool/";
	
	public static void main(String[] args) throws IOException {
		//output file
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(DIR + "readCounts.txt")));
		out.write("fileName\tnumberReads\taveReadLength\n");
		
		//list of fasta files
		File folder = new File(DIR + "sequences");
		File[] fileList = folder.listFiles();
		
		for(File f: fileList) {
			int numReads = 0;//number of reads
			int readLen = 0;//read length
			BufferedReader fasta = new BufferedReader(new InputStreamReader (new FileInputStream(f)));
			String head = fasta.readLine();//sequence header
			while(head != null) {
				String seq = fasta.readLine();//sequence
				numReads++;
				readLen += seq.length();
				head = fasta.readLine();
			}
			fasta.close();
			out.write(f.getName() + "\t" + numReads + "\t" + ((double)readLen/numReads) + "\n");
		}
		
		out.close();
	}

}
