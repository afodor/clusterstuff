/*
 * Scripts to use BWA to align CHS reads against the cards database
 * and call variants
 * 
 * only run on forward reads
 */
package kw_cre;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class BWAcardsScripts {
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static final String REF_DIR = "/users/kwinglee/card/";
	public static final String FQ_DIR = "/nobackup/afodor_research/mjzapata/CRE/CHS_raw/";
	
	public static void main(String[] args) throws IOException {
		String outDir = BASE_DIR + "bwaAlignToCards/";//folder to write results
		new File(outDir).mkdirs();
		String scriptDir = BASE_DIR + "bwaScripts/";
		new File(scriptDir).mkdirs();
		
		//script to run all files
		BufferedWriter runAll = new BufferedWriter(new FileWriter(new File(
				scriptDir + "runAllAlignToCards.sh")));
		
		//set up the models
		String[] models = {"nucleotide_fasta.protein_homolog.fasta",
				"nucleotide_fasta.protein_variant.fasta",
				"nucleotide_fasta.protein_wild_type.fasta",
				"nucleotide_fasta.rrna.fasta"};
		String[] modelNames = {"pro_homolog", "pro_variant", "pro_wt", "rrna"};
		
		//align each set of reads to each model and then call variants
		//filter for mapped reads
		File[] fastqs = new File(FQ_DIR).listFiles();
		for(File f : fastqs) {
			if(f.getName().endsWith("_1.fastq.gz")) {
				String name = f.getName().replace("_1.fastq.gz", "");
				
				BufferedWriter script = new BufferedWriter(new FileWriter(new File(
						scriptDir + "runAlignToCards_" + name)));				

				//load needed modules
				script.write("#PBS -l walltime=400:00:00\n");
				script.write("module load bwa\n");
				script.write("module load samtools\n");
				
				//add to run all
				runAll.write("qsub -q \"viper\" runAlignToCards_" + name + "\n");
				
				for(int i = 0; i < models.length; i++) {
					String fname = outDir + modelNames[i] + "_v_" + name;
					//align
					script.write("bwa mem " + REF_DIR + models[i] + 
							" " + f.getAbsolutePath() + " > " +
							fname + ".sam\n");//command to align file
					//get mapped reads only
					script.write("samtools view -h -S -F 4 " + fname + ".sam > " 
							+ fname + ".mapped.sam\n"); //get mapped reads; use -f 4 for unmapped
					//convert to bam
					script.write("samtools view -bS " + fname +".mapped.sam > " 
							+ fname + ".mapped.bam\n");
					//sort
					script.write("samotools sort -o " + fname + ".mapped.sort.bam "
							+ fname + ".mapped.bam\n");
					//index
					script.write("samtools index " + fname + ".mapped.sort.bam\n");
					//call snps
					script.write("samtools mpileup -uf " + REF_DIR + models[i]
							+ " " + fname + ".mapped.sort.bam | bcftools call -mv > " 
							+ fname + ".vcf\n");
				}
				script.close();
			}
		}
		
		runAll.close();
	}

}
