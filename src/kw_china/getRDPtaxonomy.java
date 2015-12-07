/**
 * code to make list of all taxa seen in RDP files in one place
 * output is table of each taxa from RDP database to be used in writing out taxonomies
 * 12/4/15
 */
package kw_china;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.*;


public class getRDPtaxonomy {
	
	public static class AnalyzeFileCallable implements Callable {
		private File file;
		private HashMap<String, String> taxonomy;
		public AnalyzeFileCallable(File f) {
			this.file = f;
		}
		
		public HashMap<String, String> call() throws IOException {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line != null) {
				String[] sp = line.split("\t");
				int s = sp.length - 3;
				if(!sp[s+1].equals("genus")) {//check this is genus
					System.out.println("Need to update s: " + file.getName() + "\t" + sp[s+1] + "\t" + line);
				} else {
					String gen = sp[s];
					String tax = getTaxonomy(sp);
					if(taxonomy.containsKey(gen)) {
						//check taxonomies are same
						String tax2 = taxonomy.get(gen);
						if(!tax.equals(tax2)) {
							System.out.println("MIXED TAXONOMIES: ");
							System.out.println(tax);
							System.out.println(tax2 + "\n");
						}
					} else {
						taxonomy.put(gen,  tax);
					}
				}
			}
			br.close();
			return(taxonomy);
		}
	}
	
	public static String getTaxonomy(String[] sp) {
		String phy = "";
		String cl = "";
		String ord = "";
		String fam = "";
		String gen = "";
		for(int i = 0; i < sp.length; i++) {
			if(sp[i].equals("phylum")) {
				phy = sp[i-1];
			} else if(sp[i].equals("class")) {
				cl = sp[i-1];
			} else if(sp[i].equals("order")) {
				ord = sp[i-1];
			} else if(sp[i].equals("family")) {
				fam = sp[i-1];
			} else if(sp[i].equals("genus")) {
				gen = sp[i-1];
			}
		}
		return(phy + "\t" + cl + "\t" + ord + "\t" + fam + "\t" + gen);
	}
	
	public static void main(String[] args) throws InterruptedException, ExecutionException, IOException {
		//folders with RDP results
		String[] folders = {"/projects/afodor_research/kwinglee/jobin/anaerobe/rdpResults",
				"/projects/afodor_research/kwinglee/jobin/biofilm/rdpResults", 
				"/projects/afodor_research/kwinglee/jobin/ga-stool/rdpResults"};
		//files to analyze
		HashSet<File> files = new HashSet<File>();
		//files.add(new File("C:\\Users\\kwinglee.cb3614tscr32wlt\\Documents\\Fodor\\China\\abundantOTU.chinaForward.taxonomy.txt"));
		files.add(new File("/projects/afodor_research/kwinglee/china/abundantOTU.chinaForward.taxonomy.txt"));
		for(String fname: folders) {
			File f = new File(fname);
			File[] flist = f.listFiles();
			for(File file : flist) {
				if(file.getName().startsWith("rdp_")) {
					files.add(file);
				}
			}
		}
		
		//analyze files
		ExecutorService pool = Executors.newCachedThreadPool();
		HashSet<Future<HashMap<String, String>>> taxa = new HashSet<Future<HashMap<String, String>>>();
		for(File file : files) {
			Callable<HashMap<String, String>> call= new AnalyzeFileCallable(file);
			Future<HashMap<String, String>> future = pool.submit(call);
			taxa.add(future);
		}
		
		//get results (get will block until thread is finished) and combine
		HashMap<String, String> taxonomy = new HashMap<String, String>();//map of genus to taxonomy
		for(Future<HashMap<String, String>> future : taxa) {
			HashMap<String, String> map = future.get();
			Set<String> keys = map.keySet();
			for(String k : keys) {
				String tax = map.get(k);
				if(taxonomy.containsKey(k)) {
					//check taxonomies are same
					String tax2 = taxonomy.get(k);
					if(!tax.equals(tax2)) {
						System.out.println("MIXED TAXONOMIES: ");
						System.out.println(tax);
						System.out.println(tax2 + "\n");
					}
				} else {
					taxonomy.put(k,  tax);
				}
			}
		}
		
		//write results
		BufferedWriter out = new BufferedWriter(new FileWriter(new File("/projects/afodor_research/kwinglee/china/RDPtaxonomyMultithread.txt")));
		out.write("phylum\tclass\torder\tfamily\tgenus\n");
		Set<String> gen = taxonomy.keySet();
		for(String g : gen) {
			out.write(taxonomy.get(g) + "\n");
		}
		out.close();
	}

}
