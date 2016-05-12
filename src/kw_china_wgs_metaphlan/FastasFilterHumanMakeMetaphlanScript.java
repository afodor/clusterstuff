/*
 * Filter reads that align to human genome from fasta files
 * make scripts to run metaphlan on resulting file
 * 5/12/16
 */

package kw_china_wgs_metaphlan;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class FastasFilterHumanMakeMetaphlanScript {
	public static String DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static String FAOUTDIR = DIR + "fastasFilterHuman/";
	public static String HG38DIR = DIR + "alignToHG38/";
	public static String FADIR = DIR + "fastas/";
	public static String SCRIPTDIR = DIR + "metScripts/";
	public static String METOUTDIR = DIR + "metaphlanResultsFilterHuman/";
	public static String MET_DIR = "/nobackup/afodor_research/kwinglee/software/metaphlan2";//location of metaphlan
	public static String BOW_DIR = "/nobackup/afodor_research/kwinglee/software/bowtie2-2.2.8";//bowtie2 location
	
	public static void main(String[] args) throws IOException {
		String[] files = new File(FADIR).list();
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "runAllFilterMet.sh")));
		for(String f : files) {
			if(f.endsWith(".fa")) {
				//generate fasta file with reads mapped to hg38 removed
				String sample = f.replace("_1.fa", "");
				Set<String> human = getHumanReads(sample);
				BufferedWriter out = new BufferedWriter(new FileWriter(new File
						(FAOUTDIR + sample + "_filterHuman.fa")));
				BufferedReader fa = new BufferedReader(new FileReader(new File(
						FADIR + f)));
				String line = fa.readLine();
				boolean isHuman = false;
				while(line != null) {
					if(line.startsWith(">")) {
						String head = line.replace(">", "");
						isHuman = human.contains(head);
						if(!isHuman) {
							out.write(line + "\n");
						}
					} else {
						if(!isHuman) {
							out.write(line + "\n");
						}
					}
					line = fa.readLine();
				}
				fa.close();
				out.close();
				
				//set up metaphlan scripts
				String scriptName = "metFilter_" + sample;
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				//load python
				script.write("module load python/2.7.10\n");
				//set path variables
				script.write("export PATH=" + MET_DIR + ":" + BOW_DIR + ":$PATH\n");
				script.write("export mpa_dir=" + MET_DIR + "\n");
				//run metaphlan
				script.write("metaphlan2.py " + FAOUTDIR + f + 
						" --bowtie2out " + METOUTDIR + "met_bowtie2_filterHuman_" + sample + ".bz2" +
						" --input_type fasta --nproc 2 > " + 
						METOUTDIR + "metaphlan_table_filterHuman_" + sample + "\n");
				script.close();
				
				runAll.write("qsub -q \"Cobra\" " + scriptName + "\n");
			}
		}
		runAll.close();
	}

	//for the given genome, return the set of reads that mapped to the human genome
	public static Set<String> getHumanReads(String sample) throws IOException {
		Set<String> reads = new HashSet<String>();
		BufferedReader map = new BufferedReader(new FileReader(new File(
				HG38DIR + sample + "_1.hg38.mapped.sam")));
		String line = map.readLine();
		while(line != null) {
			String[] sp = line.split("\t");
			reads.add(sp[0] + "/1");//kegg results and fasta have extra /1
			line = map.readLine();
		}
		map.close();
		return(reads);
	}
}
