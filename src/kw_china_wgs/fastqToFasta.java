/**
 * convert the fastq files from CPC to fasta
 * (input to shotmap must be in .fa or .fa.gz format)
 */
package kw_china_wgs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class fastqToFasta implements Runnable {
	private File fName;
	public fastqToFasta(File name) {
		this.fName = name;
	}
	public void run() {
		String outDir = "/projects/afodor_research/kwinglee/china/wgs/fastas/";//directory to write fastas to
		String name = this.fName.getName().replace(".fq.gz", "");
		BufferedReader fq;
		try {
			fq = new BufferedReader(new InputStreamReader(new GZIPInputStream( new FileInputStream(fName))));
			BufferedWriter fa = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(new File(outDir + name + ".fa.gz")))));
			String line1 = fq.readLine();//seq id
			while(line1 != null) {
				String line2 = fq.readLine();//seqeunce
				fq.readLine();//+
				fq.readLine();//quality scores
				String id = line1.replaceFirst("@", ">");
				fa.write(id + "\n");
				fa.write(line2 + "\n");
				line1 = fq.readLine();
			}
			fq.close();
			fa.close();
		} catch (Exception e) {
			System.out.println("error in " + name);
			e.printStackTrace();
		} 
	}

	public static void main(String[] args) throws InterruptedException {
		String fqDir = "/projects/afodor_research/kwinglee/china/wgs/from_cpc/0.cleandata/";//directory containing fastq sequences from CPC
		
		File inDir = new File(fqDir);
		File[] fqList = inDir.listFiles();
		Thread[] allThreads = new Thread[fqList.length-1];//extra file in this directory: Cleandata.stat
		int tPos = 0;//position in allThreads to add next thread to
		for(int i = 0; i < fqList.length; i++) {
			if(fqList[i].getName().endsWith(".fq.gz")) {
				Runnable r = new fastqToFasta(fqList[i]);
				Thread t = new Thread(r);
				t.start();	
				allThreads[tPos] = t;
				tPos++;
			}
		}
		
		//make sure all threads finish
		for(int i = 0; i < allThreads.length; i++) {
			allThreads[i].join();
		}
	}

}
