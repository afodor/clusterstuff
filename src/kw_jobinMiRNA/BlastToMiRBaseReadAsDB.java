/*
 * use blast to align all samples to miRBase
 */
package kw_jobinMiRNA;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BlastToMiRBaseReadAsDB {
	public static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	public static String FQDIR = DIR + "adapterFiltered/";
	public static String OUTDIR = DIR + "miRBaseBlastReadsAsDB/";
	public static String DBDIR = FQDIR + "readBlastDBs/";
	public static String SCRIPTDIR = "/projects/afodor_research/kwinglee/scripts/jobin/microRNA/";
	public static String MATREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/mature.fa";
	public static String PINREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/hairpin.fa";

	public static void main(String[] args) throws IOException {
		File odir = new File(OUTDIR);
		if(!odir.exists()) {
			odir.mkdirs();
		}
		File ddir = new File(DBDIR);
		if(!ddir.exists()) {
			ddir.mkdirs();
		}
		
		File[] files = new File(FQDIR).listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(
				SCRIPTDIR + "blastAlignAll.sh")));
		for(File fa : files) {
			if(fa.getName().endsWith(".fasta")){
				String id = fa.getName().replace(".fasta", "");
				String scriptName = "blastAlignReadAsDB_" + id;
				String name = id + ".blast.readAsDB";

				//add to run all script
				all.write("qsub -q \"copperhead\" " + scriptName + "\n");

				//write script
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						SCRIPTDIR + scriptName)));
				script.write("#PBS -l procs=1\n");
				script.write("module load blast/2.5.0+\n");
				//set up database
				script.write("makeblastdb -dbtype nucl -in " + fa.getAbsolutePath()
						+ " -out " + DBDIR + id + "\n");
				//align to mature
				script.write("blastn -outfmt 6 -db " + DBDIR + id + 
						" -query " + MATREF
						+ " -out " + OUTDIR + name + ".mature.txt\n");
				//align to hairpin
				script.write("blastn -outfmt 6 -db " + DBDIR + id +  
						" -query " + PINREF
						+ " -out " + OUTDIR + name + ".hairpin.txt\n");
				script.close();
			}
		}

		all.close();
	}
}
