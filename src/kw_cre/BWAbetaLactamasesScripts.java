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
	public static final String BASE_DIR = "/nobackup/afodor_research/kwinglee/cre/chs_v_cards/";
	public static final String REF_DIR = "/users/kwinglee/card/";
	public static final String FQ_DIR = "/nobackup/afodor_research/mjzapata/CRE/CHS_raw/";
	
	public static void main(String[] args) throws IOException {
		String outDir = BASE_DIR + "bwaAlignToBetaLactamases/";//folder to write results
		new File(outDir).mkdirs();
		String scriptDir = BASE_DIR + "bwaScripts/";
		
		//make reference file containing one representative for each beta lactam
		Set<String> refBL = new HashSet<String>();
		//reference is lowest number of that family that is from kleb pneumoniae
		refBL.add(">gb|NC_022346.1|17139-18021|ARO:3002312|KPC-2");
		refBL.add(">gb|X04515|285-1125|ARO:3002454|LEN-1");
		refBL.add(">gb|AJ635401|0-861|ARO:3002418|OKP-A-1");
		refBL.add(">gb|JX893517|0-786|ARO:3001791|OXA-247");
		refBL.add(">gb|FJ668814|76-937|ARO:3001059|SHV-1");
		refBL.add(">gb|X64523|476-1337|ARO:3000875|TEM-3");
		BufferedReader db = new BufferedReader(new FileReader(new File(
				REF_DIR + "nucleotide_fasta.protein_homolog.fasta")));
		BufferedWriter dbout = new BufferedWriter(new FileWriter(new File(
				REF_DIR + "beta_lactamase.protein_homolog.fasta")));
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
		
		//index new reference file
	}
}
