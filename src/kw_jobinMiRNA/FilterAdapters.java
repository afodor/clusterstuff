/*
 * Filter reads that contain the adapters anywhere
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Semaphore;

public class FilterAdapters {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "fastqs/";
	public static String OUTDIR = DIR + "adapterFiltered/";
	
	public static void main(String[] args) throws InterruptedException {
		//set up output directory
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		
		//set up multithreaded
		int numThreads = 8;
		Semaphore s = new Semaphore(numThreads);
		
		File[] fqs = new File(FQDIR).listFiles();
		for(File fq : fqs) {
			s.acquire();
			Filter f = new Filter(fq, s);
			new Thread(f).start();
		}
		
		for(int i = 0; i < numThreads; i++) {
			s.acquire();
		}
	}
}

class Filter implements Runnable {
	private final File file;
	private final Semaphore semaphore;
	
	public Filter(File file, Semaphore semaphore) {
		this.file = file;
		this.semaphore = semaphore;
	}

	
	@Override
	public void run() {
		try {
			String name = file.getName().split("-")[0];
			BufferedReader fq = new BufferedReader(new FileReader(file));
			BufferedWriter outFA = new BufferedWriter(new FileWriter(new File(
					FilterAdapters.OUTDIR + name + ".adapterfiltered.fasta")));//fasta file
			BufferedWriter outFQ = new BufferedWriter(new FileWriter(new File(
					FilterAdapters.OUTDIR + name + ".adapterfiltered.fastq")));//fasta file
			int numReads = 0;
			int numFiltered = 0;
			for(String line1 = fq.readLine(); line1 != null; line1 = fq.readLine()) {
				numReads++;
				String line2 = fq.readLine();
				String line3 = fq.readLine();
				String line4 = fq.readLine();
				if(line2.contains("AGATCGGAAGAGCACACGTCT") || //3'
						line2.contains("TCTAGCCTTCTCGTGTGCAGA") || //3' reverse complement
						line2.contains("GTTCAGAGTTCTACAGTCCGACGATC") || //5'
						line2.contains("CAAGTCTCAAGATGTCAGGCTGCTAG")) {//5' reverse complement
					numFiltered++;
				} else {
					outFQ.write(line1 + "\n" + line2 + "\n" + line3 + "\n" + line4 + "\n");
					outFA.write(line1.replaceFirst("@", ">") + "\n" + line2 + "\n");
				}
			}
			outFA.close();
			outFQ.close();
			fq.close();
			System.out.println(name + ": " + numReads + " " + numFiltered + " " + (100.0 * numFiltered / numReads));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		} finally {
			semaphore.release();
		}
		
	}
	
}