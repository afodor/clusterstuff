/**
 * code to run SnpEff on the reduced VCF files
 * @author kwinglee
 * @date 6/22/15
 */
package alignment_analysis;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class SnpEffScripts {
	public static String indir = "/projects/afodor_research/kwinglee/red_vcf/";
	public static String outdir = "/projects/afodor_research/kwinglee/annotate_vcfs/";
	
	public static void main(String[] args) throws IOException {
		File folder = new File(indir);
		File[] files = folder.listFiles();
		BufferedWriter all = new BufferedWriter(new FileWriter(new File(outdir+"runAll.sh")));
		for(int i = 0; i < files.length; i++) {
			String name = files[i].getName();
			BufferedWriter br = new BufferedWriter(new FileWriter(new File(outdir+"run_"+name)));
			br.write("java -cp /projects/afodor_research/kwinglee/snpeff/snpEff/snpEff.jar GCA_000598005.1.26 "+
					indir+name + " > " + outdir+name.replace("vcf", "annotate.vcf"));
			br.close();
			all.write("qsub -q \"Cobra\" run_" + name);
		}
		all.close();
		
	}

}
