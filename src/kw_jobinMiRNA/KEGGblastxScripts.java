/*
 * Scripts to run blastx on the microRNA reads against the KEGG database
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class KEGGblastxScripts {
	public static final String DB = "/nobackup/afodor_research/kwinglee/china/wgs/allKeggPep/allKegg";
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FADIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "keggBlast/";
	public static String SCRIPT_DIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/kegg/";
	private static String SPLITFA = DIR + "splitFastas/"; //directory for split fastas
	private static final int NUM_READS = 100000;//max number of reads per file
	
	public static void main(String[] args) throws IOException {
		//script to run all files
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "runAll.sh")));
		
		File[] fastas = new File(FADIR).listFiles();
		for(File f : fastas) {
			if(f.getName().endsWith(".fasta")) {
				int numSplits = 0;
				String name = f.getName().replace(".adapterfiltered.fasta", "")+ "_split" + numSplits;
								
				writeScript(name, runAll);
				
				int numReads = 0;
				BufferedReader fa = new BufferedReader(new FileReader(f));
				BufferedWriter out = new BufferedWriter(new FileWriter(new File(
						SPLITFA + name + ".fasta")));
				for(String line = fa.readLine(); line != null; line = fa.readLine()) {
					if(line.startsWith(">")) {
						numReads++;
					}
					if(numReads == NUM_READS) {
						out.close();
						numSplits++;
						name = f.getName().replace(".adapterfiltered.fasta", "")+ "_split" + numSplits;
						writeScript(name, runAll);
						out = new BufferedWriter(new FileWriter(new File(
								SPLITFA + name + ".fasta")));
					}
					out.write(line + "\n");
				}
				
				fa.close();
				out.close();
			}
		}
		
		runAll.close();
	}

	private static void writeScript(String name, BufferedWriter runAll) throws IOException {
		String scriptName = "kegg_" + name;
		//add to runAll
		runAll.write("qsub -q \"copperhead\" " + scriptName + "\n");
		//runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");
		
		//write individual script
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + scriptName)));
		script.write("#PBS -l walltime=300:00:00,mem=10GB,procs=1\n");
		script.write("module load blast/2.5.0+\n");
		script.write("blastx -outfmt 6 -db " + DB + " -query " +
				SPLITFA + name + ".fasta -out " +
				OUTDIR + name + ".kegg.blast.txt\n");//blast command
		script.close();
	}

}
