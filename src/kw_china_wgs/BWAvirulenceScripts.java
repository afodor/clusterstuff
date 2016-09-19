/*
 * Scripts to use BWA to align China WGS reads against the virulence databases
 * MvirDB and VFDB (full and core)
 */
package kw_china_wgs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BWAvirulenceScripts {
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static final String REF_DIR = "/nobackup/afodor_research/kwinglee/software/virulence/";
	public static final String SCRIPT_DIR = "/projects/afodor_research/kwinglee/scripts/china/wgs/virulence/";
	public static final String OUT_DIR = BASE_DIR + "alignToVirulenceDB/";//folder to write results
	public static BufferedWriter RUNALL;
	
	public static void main(String[] args) throws IOException {
		new File(OUT_DIR).mkdirs();
		String fastaDir = BASE_DIR + "fastas/";//folder with fasta files
		
		//script to run all files
		RUNALL = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + "runAllAlignToVirulence.sh")));
		
		//align each set of reads to each reference 
		//filter for mapped reads
		File[] fastas = new File(fastaDir).listFiles();
		for(File f : fastas) {
			if(f.getName().endsWith(".fa")) {
				String name = f.getName().replace("_1.fa", "");
				
				writeScript("VFDBfull_v_" + name, "VFDB/VFDB_full.fas", f);
				writeScript("VFDBcore_v_" + name, "VFDB/VFDB_core.fas", f);
				writeScript("MvirDB_v_" + name, "MvirDB/virulenceDB.nucleic.fasta", f);
			}
		}
		
		RUNALL.close();
	}
	
	/*
	 * writes a script with the given  name that aligns the given
	 * fasta file fa to the given database db, with the output in name
	 * adds the script to runall
	 */
	public static void writeScript(String name, String db, File fa) throws IOException {
		String sname = "runAlignToVir_" + name;
		BufferedWriter script = new BufferedWriter(new FileWriter(new File(
				SCRIPT_DIR + sname)));				

		//load needed modules
		//script.write("#PBS -l walltime=400:00:00\n");
		script.write("module load bwa\n");
		script.write("module load samtools\n");
		//align
		script.write("bwa mem " + REF_DIR + db + 
				" " + fa.getAbsolutePath() + " > " +
				OUT_DIR + name + ".sam\n");//command to align file
		//get mapped reads only
		script.write("samtools view -h -S -F 4 " + OUT_DIR + name + ".sam > " 
				+ OUT_DIR + name + ".mapped.sam\n"); //get mapped reads; use -f 4 for unmapped
		script.close();
		
		//add to run all
		RUNALL.write("qsub -q \"copperhead\" " + sname + "\n");
		
	}

}
