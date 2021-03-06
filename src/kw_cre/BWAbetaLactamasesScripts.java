package kw_cre;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

public class BWAbetaLactamasesScripts {
	//public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static final String BASE_DIR = "/nobackup/test_afodor_research/kwinglee/cre/chs_v_cards/";
	public static final String REF_DIR = "/users/kwinglee/card/";
	//public static final String FQ_DIR = "/nobackup/afodor_research/mjzapata/CRE/CHS_raw/";
	public static final String FQ_DIR = "/nobackup/test_afodor_research/mjzapata/CRE/CHS_raw/";

	public static void main(String[] args) throws IOException {
		String outDir = BASE_DIR + "bwaAlignToBetaLactamases/";//folder to write results
		File outFolder = new File(outDir);
		if(!outFolder.exists()) {
			outFolder.mkdirs();
		}
		String scriptDir = BASE_DIR + "bwaScripts/";
		File scriptFile = new File(scriptDir);
		if(!scriptFile.exists()) {
			scriptFile.mkdirs();
		}
		String ref = REF_DIR + "beta_lactamase2.protein_homolog.fasta";

		//make reference file containing one representative for each beta lactam
		Set<String> refBL = new HashSet<String>();
		//reference is lowest number of that family that is from kleb pneumoniae
		//below were references used in first round
		/*refBL.add(">gb|NC_022346.1|17139-18021|ARO:3002312|KPC-2");
		refBL.add(">gb|X04515|285-1125|ARO:3002454|LEN-1");
		//refBL.add(">gb|AJ635401|0-861|ARO:3002418|OKP-A-1");
		refBL.add(">gb|M55547|0-825|ARO:3001404|OXA-9");
		refBL.add(">gb|FJ668814|76-937|ARO:3001059|SHV-1");
		refBL.add(">gb|X64523|476-1337|ARO:3000875|TEM-3");*/
		//below are references used in second round
		refBL.add(">gb|KM379100|0-882|ARO:3003180|KPC-22");
		refBL.add(">gb|AM850922|24-885|ARO:3002453|OKP-B-20");
		refBL.add(">gb|KF151169|0-831|ARO:3001793|OXA-320");//from Proteus mirabilis but mapped CHS66
		refBL.add(">gb|KP050489|0-861|ARO:3001364|SHV-182");
		refBL.add(">gb|HQ877606|0-861|ARO:3001056|TEM-197");
		BufferedReader db = new BufferedReader(new FileReader(new File(
				REF_DIR + "nucleotide_fasta.protein_homolog.fasta")));
		BufferedWriter dbout = new BufferedWriter(new FileWriter(new File(
				ref)));
		boolean write = false;
		for(String line = db.readLine(); line != null; line = db.readLine()) {
			if(line.startsWith(">")) {
				String[] sp = line.split(" ");
				write = refBL.contains(sp[0]); 
				if(write) {
					dbout.write(line + "\n");
				} 
			} else if(write) {
				dbout.write(line + "\n");
			}
		}
		db.close();
		dbout.close();

		//align each set of reads to the beta-lactamases and get stats
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				scriptDir + "runAllAlignToBetaLactamases.sh")));
		File[] fastqs = new File(FQ_DIR).listFiles();
		for(File f : fastqs) {
			if(f.getName().endsWith("_1.fastq.gz")) {
				String name = f.getName().replace("_1.fastq.gz", "");

				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + "runAlignToBLs_" + name)));		

				//add to run all
				runAll.write("qsub -q \"viper_batch\" runAlignToBLs_" + name + "\n");
				
				//load needed modules
				script.write("module load bwa\n");
				script.write("module load samtools\n");

				String fname = outDir + name + ".blalign2";
				//align
				script.write("bwa mem " + ref + 
						" " + f.getAbsolutePath() + " " + 
						f.getAbsolutePath().replace("_1", "_2") + "> " +
						fname + ".sam\n");//command to align file
				//convert to bam
				script.write("samtools view -bSh " + fname +".sam > " 
						+ fname + ".bam\n");
				//sort
				script.write("samtools sort " + fname + ".bam "
						+ fname + ".sort\n");
				//remove duplicates
				script.write("samtools rmdup " + fname + ".sort.bam " + fname +
						".sort.rmdup.bam\n");
				//index
				script.write("samtools index " + fname + ".sort.rmdup.bam\n");
				//stats
				script.write("samtools idxstats " + fname + ".sort.rmdup.bam > "
						+ fname + ".idxstats.txt\n");
				//depth
				script.write("samtools depth " + fname + ".sort.rmdup.bam > " 
						+ fname + ".depth.txt\n");
				script.close();
			}
		}
		runAll.close();
	}
}
