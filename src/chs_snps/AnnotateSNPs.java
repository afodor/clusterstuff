/**
 * takes results from WriteSNPFiles and converts context to sequence, writing results to new annotated snps file
 *
 */
package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import coPhylog.BitHolder;

public class AnnotateSNPs {
	
	/**
	 * 
	 * @param counts	counts for each base, in format [A,C,G,T]
	 * @return
	 */
	public static String getMid(String counts) {
		String[] sp = counts.replace("[", "").replace("]", "").split(",");
		int maxi = 0;
		for(int i = 1; i < sp.length; i++) {
			if(Integer.parseInt(sp[i]) > Integer.parseInt(sp[maxi])) {
				maxi = i;
			}
		}
		if(maxi == 0) {
			return("A");
		} else if(maxi == 1) {
			return("C");
		} else if(maxi == 2) {
			return("G");
		} else if(maxi == 3) {
			return("T");
		} else {//should not get here
			return("N");
		}
		
	}
	
	
	/*
	 * expected input = path to comparison files
	 */
	public static void main(String[] args) throws Exception {
		if(args.length != 1) {
			System.out.println("Expected input: path_to_folder");
			System.exit(1);
		}
		File folder = new File(args[0]);
		if(!args[0].endsWith("/")) {
			args[0]+="/";
		}
		File[] files = folder.listFiles();
		for(int i = 0; i < files.length; i++) {
			if(files[i].isFile()) {
				BufferedReader br = new BufferedReader(new FileReader(files[i]));
				String line = br.readLine();
				if(line != null && line.equals("longID\tcontext1\tcontext2\tdistance")) { //header must be exact match; way to confirm only annotating files in correct format
					BufferedWriter out = new BufferedWriter(new FileWriter(new File(args[0]+files[i].getName()+".annotate.txt")));
					out.write(line+"\tseq1\tseq2\n");
					line = br.readLine();
					while(line != null) {
						String[] sp = line.split("\t");
						
						//get left and right sequences
						long con = Long.valueOf(sp[0]);
						String left = BitHolder.getASequence(con, true, 13);
						String right = BitHolder.getASequence(con, false, 13);
						
						//get middle sequence for each strain
						String mid1 = getMid(sp[1]);
						String mid2 = getMid(sp[2]);
						
						//write results
						out.write(line+"\t"+left+mid1+right+"\t"+left+mid2+right+"\n");
						line = br.readLine();
					}
					out.close();
				}
				br.close();
			}
		}
	}

}
