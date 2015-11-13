/**
 * Test that got all reads in fastq to fasta conversion
 */
package kw_china_wgs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class testFastqToFasta {
	
	public static void main(String[] args) throws FileNotFoundException, IOException {
		String fqDir = "/projects/afodor_research/kwinglee/china/wgs/from_cpc/0.cleandata/";//directory containing fastq sequences from CPC
		String faDir = "/projects/afodor_research/kwinglee/china/wgs/fastas/";//directory to write fastas to
		
		//set up output
		BufferedWriter out = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/china/wgs/fastaCheck.txt")));
		out.write("file\tnumLines_fastq\tnumLines_fasta\tpass?\n");
		
		File inDir = new File(fqDir);
		File[] fqList = inDir.listFiles();
		for(int i = 0; i < fqList.length; i++) {
			if(fqList[i].getName().endsWith(".fq.gz")) {
				//fastq count
				String name = fqList[i].getName().replace(".fq.gz", "");
				BufferedReader fq = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(fqList[i]))));
				int numLinesFq = 0;
				String line = fq.readLine();
				while(line != null) {
					numLinesFq++;
					line = fq.readLine();
				}
				fq.close();
				
				//fasta count
				BufferedReader fa = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(new File(faDir + name + ".fa.gz")))));
				int numLinesFa = 0;
				line = fa.readLine();
				while(line != null) {
					numLinesFa++;
					line = fa.readLine();
				}
				fa.close();
				
				//write results
				out.write(name + "\t" + numLinesFq + "\t" + numLinesFa + "\t");
				if(numLinesFq == numLinesFa * 2) {
					out.write("pass\n");
				} else {
					out.write("FAIL\n");
					System.out.println(name);
				}
			}
		}
		out.close();
	}

}
