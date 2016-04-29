/*
 * Convert fastq.gz files to fastq for stitching
 * Also remove 16S primers (barcode and adaptor has been removed)
 */
package kw_jobinDolphin;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class TrimFastqForStitching {
	public static final String FQ_DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/fastqs/";
	public static final String OUT_DIR = "/nobackup/afodor_research/kwinglee/jobin/dolphin/trimmed_fastqs/";
	/*public static final String F16S = "CCTACGGGNGGCWGCAG";//forward 16s primer
	public static final String R16S = "GACTACHVGGGTATCTAATCC";//reverse 16s primer*/
	public static final String F16S = "[NC]CTACGGG[ACTG]GGC[AT]GCAG";//forward 16s primer
	public static final String R16S = "[NG]ACTAC[ACT][ACG]GGGTATCTAATCC";//reverse 16s primer
	
	public static void main(String[] args) throws Exception {
		String[] samples = new File(FQ_DIR).list();
			
		for(String s : samples) {
			String path = FQ_DIR + s + "/Data/Intensities/BaseCalls/";
			String[] reads = new File(path).list();
			//fastq readers
			BufferedReader fwd_in = null;
			BufferedReader rev_in = null;
			//fasta writers
			BufferedWriter fwd_out = null;
			BufferedWriter rev_out = null;
			for(String r : reads) {
				String[] rName = r.split("_");
				if(r.contains("R1")) {
					fwd_in = new BufferedReader(
							new InputStreamReader(
									new GZIPInputStream(
											new FileInputStream(
													path + r))));
					fwd_out = new BufferedWriter(
							new FileWriter(
									new File(
											OUT_DIR + rName[0] + "_R1.fasta")));
				} else if(r.contains("R2")) {
					rev_in = new BufferedReader(
							new InputStreamReader(
									new GZIPInputStream(
											new FileInputStream(
													path + r))));
					rev_out = new BufferedWriter(
							new FileWriter(
									new File(
											OUT_DIR + rName[0] + "_R2.fasta")));
				}
			}
			
			String f1 = fwd_in.readLine();//line 1 = header
			String r1 = rev_in.readLine();
			while(f1 != null && r1 != null) {
				//get additional read info
				String f2 = fwd_in.readLine();//line 2 = sequence
				String f3 = fwd_in.readLine();//line 3 = +
				String f4 = fwd_in.readLine();//line 4 = quality score
				String r2 = rev_in.readLine();
				String r3 = rev_in.readLine();
				String r4 = rev_in.readLine();
				
				//if both primers match, trim and write
				if(f2.matches("^" + F16S) && r2.matches("^" + R16S)) {
					//trim sequences
					String fseq = f2.replaceFirst(F16S, "");
					String rseq = r2.replaceFirst(R16S, "");
					//trim quals
					String fqual = f4.substring(17);//fwd primer is 17 bases
					String rqual = r4.substring(21);//rev primer is 21 bases
					//write
					fwd_out.write(f1 + "\n"
							+ fseq + "\n"
							+ f3 + "\n"
							+ fqual + "\n");
					rev_out.write(r1 + "\n"
							+ rseq + "\n"
							+ r3 + "\n"
							+ rqual + "\n");
				}
				f1 = fwd_in.readLine();
				r1 = rev_in.readLine();
			}
			if(f1 != null || r1 != null) {
				System.err.println("Uneven number of lines in " + s);
			}
			fwd_in.close();
			rev_in.close();
			fwd_out.close();
			rev_out.close();
		}
	}

}
