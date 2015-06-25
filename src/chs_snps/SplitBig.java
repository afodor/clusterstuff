/*
 * Split the 12 biggest files, which keep running out of memory, into 2 files of reads, 
 * and generate scripts to run contexts
 */
package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class SplitBig {
	public static String outdir = "/projects/afodor_research/kwinglee/cophylog_all80chs/splitBig/";
	public static String indir = "/projects/afodor_research/mjzapata/CRE/CHS_raw/";
	
	public static void main(String[] args) throws FileNotFoundException, IOException {
		String[] big = {"SRR1159063_1", "SRR1159216_1", "SRR1159123_1", "SRR1159223_2",
				"SRR1159345_1", "SRR1159182_2", "SRR1159182_1", "SRR1159063_2", 
				"SRR1159123_2", "SRR1159216_2", "SRR1159223_1", "SRR1159345_2"};//the 12 genomes
		
		//log file
		BufferedWriter log = new BufferedWriter(new FileWriter(new File(outdir+"SplitBig_log")));
		log.write("file\tNumLines1\tNumLines2\tNumLines3\tNumLines4TotalLines\n");
		
		//command to run scripts
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(outdir+"runAllContext.sh")));
		
		for(int i = 0; i < big.length; i++) {
			String f = big[i];
			//set up scripts to run contexts
			all.write("qsub -q \"Cobra_batch\" run_" + f  + "A\n");
			all.write("qsub -q \"Cobra_batch\" run_" + f  + "B\n");
			all.write("qsub -q \"Cobra_batch\" run_" + f  + "C\n");
			all.write("qsub -q \"Cobra_batch\" run_" + f  + "D\n");
			
			BufferedWriter script = new BufferedWriter(new FileWriter(new File(outdir+"run_"+f+"A")));
			script.write("#PBS -l nodes=1:ppn=12\n");
			script.write("java -cp /users/kwinglee/git/clusterstuff/bin -mx30000m chs_snps.WriteBinaryContextsFromFastQ "
					+ outdir + f + "A.fastq.gz " + outdir + "context" + 
						f + "A_context.gz");
			script.close();
			
			script = new BufferedWriter(new FileWriter(new File(outdir+"run_"+f+"B")));
			script.write("#PBS -l nodes=1:ppn=12\n");
			script.write("java -cp /users/kwinglee/git/clusterstuff/bin -mx30000m chs_snps.WriteBinaryContextsFromFastQ "
					+ outdir + f + "B.fastq.gz " + outdir + "context" + 
						f + "B_context.gz");
			script.close();
			
			script = new BufferedWriter(new FileWriter(new File(outdir+"run_"+f+"C")));
			script.write("#PBS -l nodes=1:ppn=12\n");
			script.write("java -cp /users/kwinglee/git/clusterstuff/bin -mx30000m chs_snps.WriteBinaryContextsFromFastQ "
					+ outdir + f + "C.fastq.gz " + outdir + "context" + 
						f + "C_context.gz");
			script.close();
			
			script = new BufferedWriter(new FileWriter(new File(outdir+"run_"+f+"D")));
			script.write("#PBS -l nodes=1:ppn=12\n");
			script.write("java -cp /users/kwinglee/git/clusterstuff/bin -mx30000m chs_snps.WriteBinaryContextsFromFastQ "
					+ outdir + f + "D.fastq.gz " + outdir + "context" + 
						f + "D_context.gz");
			script.close();
			
			//split file
			BufferedReader reader = 
					new BufferedReader(new InputStreamReader( 
							new GZIPInputStream( new FileInputStream(new File(indir+f+".fastq.gz")))));
			BufferedWriter out1 = new BufferedWriter(new OutputStreamWriter(
					new GZIPOutputStream(new FileOutputStream(new File(outdir+f+"A.fastq.gz")))));
			BufferedWriter out2 = new BufferedWriter(new OutputStreamWriter(
					new GZIPOutputStream(new FileOutputStream(new File(outdir+f+"B.fastq.gz")))));
			BufferedWriter out3 = new BufferedWriter(new OutputStreamWriter(
					new GZIPOutputStream(new FileOutputStream(new File(outdir+f+"C.fastq.gz")))));
			BufferedWriter out4 = new BufferedWriter(new OutputStreamWriter(
					new GZIPOutputStream(new FileOutputStream(new File(outdir+f+"D.fastq.gz")))));
			String line = reader.readLine();
			int count = 1; //number of lines read
			int c1 = 0;//num reads in 1
			int c2 = 0;//num read in 2
			int c3 = 0;//num reads in 3
			int c4 = 0;//num reads in 4
			while(line != null) {
				if(count <= 4) {
					out1.write(line+"\n");
					c1++;
				} else if(count <= 8) {
					out2.write(line + "\n");
					c2++;
				} else if(count <=12) {
					out3.write(line + "\n");
				} else if(count <=16) {
					out4.write(line+"\n");
					if(count == 16) {
						count = 0;
					}
				} else {
					System.err.println("Problem with count: "+count);
				}
				count++;
				line = reader.readLine();
			}
			log.write(f + "\t" + c1 + "\t" + c2 + "\t" + c3 + "\t" + c4 + "\t" + (c1+c2+c3+c4)+"\n");
			log.flush();
			reader.close();
			out1.close();
			out2.close();
			out3.close();
			out4.close();
		}
		log.close();
		all.close();
	}

}
