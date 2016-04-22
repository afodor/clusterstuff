/*
 * create scripts to run metaphlan for all China WGS
 */
package kw_china_wgs_metaphlan;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class AltMetaphlanOutputs {
	public static String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";//directory for outputs
	public static String MET_DIR = "/nobackup/afodor_research/kwinglee/software/metaphlan2";//location of metaphlan
	public static String BOW_DIR = "/nobackup/afodor_research/kwinglee/software/bowtie2-2.2.8";//bowtie2 location
	
	public static void main(String[] args) throws IOException {
		String scriptDir = BASE_DIR + "metScripts/";
		String outDir = BASE_DIR + "metaphlanResults/";
		String fastaDir = BASE_DIR + "fastas/";
		
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				scriptDir + "altMet")));
		//load python
		script.write("module load python/2.7.10\n");
		//set path variables
		script.write("export PATH=" + MET_DIR + ":" + BOW_DIR + ":$PATH\n");
		script.write("export mpa_dir=" + MET_DIR + "\n");
		String[] fastas = new File(fastaDir).list();
		for(String f : fastas) {
			if(f.endsWith(".fa")) {
				String sample = f.replace("_1.fa", "");
				
				//run metaphlan with marker abundance output
				/*script.write("metaphlan2.py " + outDir + "met_bowtie2_" + sample + ".bz2 " +
						" --input_type bowtie2out --nproc 2 -t marker_ab_table > " + 
						outDir + "metaphlan_marker_ab_table_" + sample + "\n");
				//run metaphlan with clade profiles
				script.write("metaphlan2.py " + outDir + "met_bowtie2_" + sample + ".bz2 " +
						" --input_type bowtie2out --nproc 2 -t clade_profiles > " + 
						outDir + "metaphlan_clade_table_" + sample + "\n");*/
				//run metaphlan with relative abundance with read stats
				script.write("metaphlan2.py " + outDir + "met_bowtie2_" + sample + ".bz2 " +
						" --input_type bowtie2out --nproc 2 -t rel_ab_w_read_stats > " + 
						outDir + "metaphlan_rel_ab_w_read_stats_table_" + sample + "\n");
				//run metaphlan with clade profiles
				script.write("metaphlan2.py " + outDir + "met_bowtie2_" + sample + ".bz2 " +
						" --input_type bowtie2out --nproc 2 -t reads_map > " + 
						outDir + "metaphlan_reads_map_" + sample + "\n");
				
			}
		}
		
		script.close();
	}

}
