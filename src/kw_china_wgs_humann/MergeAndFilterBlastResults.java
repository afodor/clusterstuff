/*
 * Merge split blast results and filter out human reads
 * Designed so can run before all blast jobs complete, so have list of genomes
 * 	previously finished, not finished, and finished now
 */
package kw_china_wgs_humann;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class MergeAndFilterBlastResults {
	public static String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static int NUM_GENOMES = 40;//number of genomes
	public static int NUM_SPLITS = 100;//number of files genome split into
	public static String OUT_DIR = BASE_DIR + "kegg_split_blastx_results_merge_filter/";//directory to write results to
	
	public static void main(String[] args) throws Exception {
		//get list of genomes
		String[] fastas = new File(BASE_DIR + "fastas/").list();
		String[] genomes = new String[NUM_GENOMES];
		int num = 0;
		for(String f : fastas) {
			if(f.endsWith(".fa")) {
				genomes[num] = f.replace("_1.fa", "");
				num++;
			}
		}
		if(num != NUM_GENOMES) {
			throw new Exception("Incorrect number of genomes");
		}
		
		//get set of previously completed genomes
		File prev = new File(OUT_DIR + "genomes_finished_blasting.txt");
		Set<String> done = new HashSet<String>();
		if(prev.exists()) {
			BufferedReader br = new BufferedReader(new FileReader(prev));
			String line = br.readLine();
			while(line != null) {
				done.add(line);
				line = br.readLine();
			}
			br.close();
		}
		
		//set up output files
		BufferedWriter writeDone = new BufferedWriter(new FileWriter(prev));
		BufferedWriter writeNewDone = new BufferedWriter(new FileWriter(new File(
				OUT_DIR + "genomes_newly_finished.txt")));
		BufferedWriter writeNotDone = new BufferedWriter(new FileWriter(new File(
				OUT_DIR + "genomes_not_done.txt")));
		BufferedWriter log = new BufferedWriter(new FileWriter(new File(
				BASE_DIR + "runMerge.log")));//log file
		
		//for each genome, if hasn't been done before, start merging and filtering
		//see if done
		for(String gen : genomes) {
			log.write(gen);
			if(!done.contains(gen)) {//skip if was previously finished
				//check all output files exist
				String exists = checkAllOutputExists(gen);
				if(!exists.equals("TRUE")) {
					writeNotDone.write(gen + "\t" + exists + "\tmissing\n");
					log.write("\tNot done\n");
					log.flush();
				} else {
					//for each file, check if done by whether the last read in the blast 
					//		table is close to the last in the fasta
					//in the process, if read did not map to human (not filtered), 
					//		write to merged file
					boolean complete = mergeFilterCheck(gen);
					if(complete) {
						writeDone.write(gen + "\n");
						writeNewDone.write(gen + "\n");
						log.write("\tDone\n");
						log.flush();
					} else {
						writeNotDone.write(gen + "\tnot complete\n");
						log.write("\tNot done\n");
						log.flush();
					}
				}
			} else {
				writeDone.write(gen + "\n");
				log.write("\tPreviously done\n");
				log.flush();
			}
		}
		
		writeDone.close();
		writeNewDone.close();
		writeNotDone.close();
		log.close();
	}
	
	//for each of the 100 blast results, check that each blast result exists
	//if it does, return "TRUE", else return the missing file
	public static String checkAllOutputExists(String genome) {
		String path = BASE_DIR + "kegg_split_blastx_results/kegg_split_" + genome + "_";
		for(int i = 0; i < NUM_SPLITS; i++) {
			File f = new File(path + i + ".txt");
			if(!f.exists()) {
				return(f.getName());
			}
		}
		return("TRUE");
	}
	
	//for the given genome, return the set of reads that mapped to the human genome
	public static Set<String> getHumanReads(String genome) throws IOException {
		Set<String> reads = new HashSet<String>();
		BufferedReader map = new BufferedReader(new FileReader(new File(
				BASE_DIR + "alignToHG38/" + genome + "_1.hg38.mapped.sam")));
		String line = map.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			reads.add(sp[0] + "/1");//kegg results and fasta have extra /1
			line = map.readLine();
		}
		map.close();
		return(reads);
	}
	
	//for each genome, for each file, 
	//		check if done by whether the last read in the blast 
	//		table is close to the last in the fasta
	//in the process, if read did not map to human (not filtered), 
	//		write to merged file
	public static boolean mergeFilterCheck(String genome) throws IOException {
		BufferedWriter merge = new BufferedWriter(new FileWriter(new File(
				OUT_DIR + "kegg_merge_filter_human_" + genome + ".txt")));
		//get set of reads mapped to human genome
		Set<String> human = getHumanReads(genome);
		//for each file merge, filter, test if done
		for(int i = 0; i < NUM_SPLITS; i++) {
			String lastRead = "";
			BufferedReader br = new BufferedReader(new FileReader(new File(
					BASE_DIR + "kegg_split_blastx_results/kegg_split_" 
							+ genome + "_" + i + ".txt")));
			String line = br.readLine();
			//filter and merge
			while(line != null) {
				String read = line.split("\t")[0];
				lastRead = read;
				if(!human.contains(read)) {
					merge.write(line + "\n");
				}
				line = br.readLine();
			}
			br.close();
			
			//check if file was done
			int numReads = 0;
			br = new BufferedReader(new FileReader(new File(
					BASE_DIR + "splitFastas/split_" + genome + "_" + i + ".fa")));
			line = br.readLine();
			boolean seen = false;//if last read has been seen
			while(line != null) {
				if(line.startsWith(">")) {
					numReads++;
					line = line.replace(">", "");
					if(line.equals(lastRead)) {
						seen = true;
						if(numReads/100000 < .9995) {
							//last read is less than 99.95% through the file -> within 50 reads of end
							br.close();
							merge.close();
							System.out.println(genome + "_" + i + "not complete: " + numReads + "reads");
							return(false);
						} else {
							break;
						}
					}
				}
				line = br.readLine();
			}
			br.close();
			if(!seen) {
				System.err.println(genome + "_" + i + " last read not seen: " + lastRead);
			}
		}
		merge.close();
		return(true);
	}
}
