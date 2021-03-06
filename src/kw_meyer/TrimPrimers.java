/*
 * Remove primers from merged reads
 * Primers are:
 *  341F CCT ACG GGN GGC WGC AG
 *  785R GAC TAC HVG GGT ATC TAA TCC
 */
package kw_meyer;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class TrimPrimers {
	public static final String FASTA_FOLDER = "/nobackup/afodor_research/katieDataSep21/inFasta/";
	public static final String OUT_DIR = "/nobackup/afodor_research/kwinglee/meyer/fastas_merged_trimmed/";
	
	public static void main(String[] args) throws IOException {
		File[] fas = new File(FASTA_FOLDER).listFiles();
		for(File fa : fas) {
		//File fa = fas[0];
			String name = fa.getName().replace(".fna", "_trimmed.fna");
			BufferedReader in = new BufferedReader(new FileReader(fa));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(
					OUT_DIR + name)));
			int numSeqs = 0;
			int numRem = 0;
			String header = in.readLine();
			for(String line = in.readLine(); line != null; line = in.readLine()) {
				if(!line.startsWith(">")) {
					numSeqs++;
					String seq = line.replaceAll("^CCTACGGG[AGTC]GGC[AT]GCAG", "").replaceAll("GGATTAGATACCC[^A][^C]GTAGTC$", "");
					if(seq.length() == line.length() - 17 - 21) {
						out.write(header + "\n" + seq + "\n");
					} else {
						/*if(numRem < 10) {
							System.out.println(line);
							System.out.println(seq);
							System.out.println("Fwd: " + line.startsWith("CCTACGGG"));
							System.out.println("Rev: " + line.endsWith("GTAGTC"));
							System.out.println("");
						}*/
						numRem++;
					}
				} else {
					header = line;
				}
			}
			
			in.close();
			out.close();
			System.out.println(fa.getName() + " " + numSeqs + " " + numRem + " " + (100.0 * numRem / numSeqs) + " " + (numSeqs - numRem));
		}
	}
}
