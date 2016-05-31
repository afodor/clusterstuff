/*
 * Scripts to use BWA to align China WGS reads against the cards database
 * and call variants
 */
package kw_china_wgs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BWAcardsScripts {
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/china/wgs/";
	public static final String REF_DIR = "/users/kwinglee/card/";
	
	public static void main(String[] args) throws IOException {
		String outDir = BASE_DIR + "alignToCards/";//folder to write results
		new File(outDir).mkdirs();
		String fastaDir = BASE_DIR + "fastas/";//folder with fasta files
		String scriptDir = BASE_DIR + "bwaScripts/";
		
		//script to run all files
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				scriptDir + "runAllAlignToCards.sh")));
		
		//index the models
		String[] models = {"nucleotide_fasta.protein_homolog.fasta",
				"nucleotide_fasta.protein_variant.fasta",
				"nucleotide_fasta.protein_wild_type.fasta",
				"nucleotide_fasta.rrna.fasta"};
		String[] modelNames = {"pro_homolog", "pro_variant", "pro_wt", "rrna"};
		BufferedWriter index = new BufferedWriter(new FileWriter(new File(
				scriptDir + "indexCards")));
		index.write("module load bwa\n");
		for(String mod : models) {
			index.write("bwa index " + REF_DIR + mod + "\n");
		}
		index.close();
		
		//align each set of reads to each model and then call variants
		//filter for mapped reads
		File[] fastas = new File(fastaDir).listFiles();
		for(File f : fastas) {
			if(f.getName().endsWith(".fa")) {
				String name = f.getName().replace(".fa", "");
				
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + "runAlignToCards_" + name)));				

				//load needed modules
				script.write("#PBS -l walltime=400:00:00\n");
				script.write("module load bwa\n");
				script.write("module load samtools\n");
				
				//add to run all
				runAll.write("qsub -q \"Cobra\" runAlignToCards_" + name + "\n");
				
				for(int i = 0; i < models.length; i++) {
					String fname = outDir + modelNames[i] + "_v_" + name;
					//align
					script.write("bwa mem " + REF_DIR + models[i] + 
							" " + f.getAbsolutePath() + " > " +
							fname + ".sam\n");//command to align file
					//get mapped reads only
					script.write("samtools view -S -F 4 " + fname + ".sam > " 
							+ fname + ".mapped.sam\n"); //get mapped reads; use -f 4 for unmapped
					//convert to bam
					script.write("samtools view -bS " + fname +".mapped.sam > " 
							+ fname + ".mapped.bam\n");
					//index
					script.write("samtools index " + fname + ".mapped.bam\n");
					//call snps
					script.write("samtools mpileup -f " + REF_DIR + models[i]
							+ " " + fname + ".mapped.bam\n");
				}
				script.close();
			}
		}
		
		runAll.close();
	}

}
