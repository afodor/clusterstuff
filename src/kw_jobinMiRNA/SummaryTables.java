/*
 * Generate tables summarizing the number of reads mapped to each database,
 * and the number of reads uniquely mapping to each database
 */

package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

public class SummaryTables {
	private static String DIR = "/nobackup/afodor_research/kwinglee/jobin/microRNA/";
	private static String RAW_READ_DIR = DIR + "adapterFiltered/";
	private static String OUTDIR = DIR + "summaryTables/";
	private static HashMap<String, Integer> TOT_READS;//sample to number of reads
	private static HashMap<String, HashSet<String>> ALL_READS;//sample to set of all reads

	//BLAST parameters
	private static double PID = 99;//minimum percent ID allowed
	private static int MISMATCH = 1;//maximum number of mismatches
	private static int LENGTH = 26;//minimum length of alignment; read length is 36, but average database sizes are for piR and miR are 21-27

	public static void main(String[] args) throws Exception {
		//get number of reads for each sample
		TOT_READS = new HashMap<String, Integer>();
		ALL_READS = new HashMap<String, HashSet<String>>();
		File[] faList = new File(RAW_READ_DIR).listFiles();
		for(File fa : faList) {
			if(fa.getName().endsWith(".fasta")) {
				String name = fa.getName().replace(".adapterfiltered.fasta", "");
				HashSet<String> reads = new HashSet<String>();
				int count = 0;
				BufferedReader fasta = new BufferedReader(new FileReader(fa));
				for(String line = fasta.readLine(); line != null; line = fasta.readLine()) {
					if(line.startsWith(">")) {
						count++;
						reads.add(line.replace(">", "").split(" ")[0]);
					}
				}
				fasta.close();
				TOT_READS.put(name, count);
				ALL_READS.put(name, reads);
				if(count != reads.size()) {
					throw new Exception("Different sizes: " + fa.getName() +
							" " + count + " " + reads.size());
				}
			}
		}

		//analyze sets of databases
		analyzeDatabaseSet(new String[]{"mouseBowtie", "mouseBowtie", "miniKraken"},
				"test2Same");
		analyzeDatabaseSet(new String[]{"piRBaseBowtie", "mouseBowtie"},
				"test2Mouse");
		analyzeDatabaseSet(new String[]{"miniKraken", "mouseBowtie"},
				"classifyBacteriaOrMouse");
		analyzeDatabaseSet(new String[]{"keggBlast", "miRBaseBowtie", "piRBaseBowtie",
				"silvaBowtie", "sRNATarBaseBowtie"},
				"smallRNAdatabasesBowtie");
		analyzeDatabaseSet(new String[]{"keggBlast", "miRBaseBlast", "piRBaseBlast",
				"silvaBlast", "sRNATarBaseBlast"},
				"smallRNAdatabasesBlast");
	}

	//for given set of databases (folder containing results), run the analysis
	//prefix is the prefix of the output files
	private static void analyzeDatabaseSet(String[] databases, String prefix) 
			throws Exception {
		HashMap<String, HashMap<String, HashSet<String>>> dbMapped = 
				new HashMap<String, HashMap<String, HashSet<String>>>();//map database name to map of sample to set of reads mapped
		for(String db : databases) {
			HashMap<String, HashSet<String>> mapped;//map of sample to set of reads mapped
			if(db.contains("Bowtie")) {
				mapped = getBowtieMap(db);
			} else if(db.contains("Blast")) {
				mapped = getBlastMap(db);
			} else if(db.contains("Kraken")) {
				mapped = getKrakenMap(db);
			} else {
				throw new Exception("Invalid database type: " + db);
			}
			dbMapped.put(db, mapped);
		}	

		//write number of reads mapped (regardless of overlap between databases)
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(
				OUTDIR + prefix + "_ignoreOverlap.txt")));
		//write header
		out.write("sampleID\ttotalNumberReads");
		for(String db : databases) {
			out.write("\t" + db);
		}
		out.write("\n");
		//write counts
		for(int s = 1; s <= 10; s++) {
			String sid = "Sample" + s;
			out.write(sid + "\t" + TOT_READS.get(sid));
			for(String db : databases) {
				if(dbMapped.get(db).containsKey(sid)) {
					out.write("\t" + dbMapped.get(db).get(sid).size());
				} else {
					out.write("\t0");
					System.out.println("Missing sample: " + db + " " + sid);
				}
			}
			out.write("\n");
		}
		out.close();

		//write number of unique reads mapped to each database/number unmapped/number mapped to multiple databases
		out = new BufferedWriter(new FileWriter(new File(
				OUTDIR + prefix + "_uniqueReads.txt")));
		//write header
		out.write("sampleID\ttotalNumberReads");
		for(String db : databases) {
			out.write("\t" + db);
		}
		out.write("\tnumMultipleMappedReads\tnumUnmappedReads\n");
		//write counts
		for(int s = 1; s <= 10; s++) {
			String sid = "Sample" + s;
			HashSet<String> unmapped = new HashSet<String>();
			HashSet<String> multiple = new HashSet<String>();
			unmapped.addAll(ALL_READS.get(sid));
			out.write(sid + "\t" + TOT_READS.get(sid));
			for(int db1 = 0; db1 < databases.length; db1++) {
				HashSet<String> unique = new HashSet<String>();
				unique.addAll(dbMapped.get(databases[db1]).get(sid));
				unmapped.removeAll(unique);
				for(int db2 = 0; db2 < databases.length; db2++) {
					if(db1 != db2) {
						HashSet<String> reads = new HashSet<String>();
						reads.addAll(dbMapped.get(databases[db2]).get(sid));

						//multiple set is intersection of reads and unique
						HashSet<String> intsn = new HashSet<String>(reads);
						intsn.retainAll(unique);
						multiple.addAll(intsn);
						
						//unique reads (unique set) is all reads not in reads set
						unique.removeAll(reads);//unique.removeAll(intsn);
					}
				}
				out.write("\t" + unique.size());
			}
			out.write("\t" + multiple.size() + "\t" + unmapped.size() + "\n");
		}
		out.close();
	}

	//returns set of reads mapped to family level with kraken
	private static HashMap<String, HashSet<String>> getKrakenMap(String db) throws IOException {
		HashMap<String, HashSet<String>> mapped = new HashMap<String, HashSet<String>>();
		File[] tables = new File(DIR + db).listFiles();
		for(File f : tables) {
			if(f.getName().endsWith("_mpa")) {
				String name = f.getName().split("\\.")[0];
				HashSet<String> reads;
				if(mapped.containsKey(name)) {
					reads = mapped.get(name);
				} else {
					reads = new HashSet<String>();
				}
				BufferedReader br = new BufferedReader(new FileReader(f));
				for(String line = br.readLine(); line != null; line = br.readLine()) {
					if(line.contains("f__")) {
						reads.add(line.split("\t")[0]);
					}
				}
				br.close();
				mapped.put(name, reads);
			}
		}
		return(mapped);
	}

	//returns set of reads mapped using BLAST
	private static HashMap<String, HashSet<String>> getBlastMap(String db) throws IOException {
		HashMap<String, HashSet<String>> mapped = new HashMap<String, HashSet<String>>();
		File[] tables = new File(DIR + db).listFiles();
		for(File f : tables) {
			if(f.getName().endsWith(".txt") && f.getName().contains(".blast.") &&
					f.getName().startsWith("Sample")) {
				String name = f.getName().split("\\.")[0].split("_")[0];
				HashSet<String> reads;
				if(mapped.containsKey(name)) {
					reads = mapped.get(name);
				} else {
					reads = new HashSet<String>();
				}
				BufferedReader br = new BufferedReader(new FileReader(f));
				for(String line = br.readLine(); line != null; line = br.readLine()) {
					String[] sp = line.split("\t");
					if(sp.length > 5) {
						//qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
						if(Double.parseDouble(sp[2]) > PID &&
								Integer.parseInt(sp[4]) <= MISMATCH &&
								Integer.parseInt(sp[3]) >= LENGTH) {
							reads.add(sp[0]);
						}
					}
				}
				br.close();
				mapped.put(name, reads);
			}
		}
		return(mapped);
	}

	private static HashMap<String, HashSet<String>> getBowtieMap(String db) throws IOException {
		HashMap<String, HashSet<String>> mapped = new HashMap<String, HashSet<String>>();
		File[] tables = new File(DIR + db).listFiles();
		for(File f : tables) {
			if(f.getName().endsWith(".sam")) {
				String name = f.getName().split("\\.")[0];
				HashSet<String> reads;
				if(mapped.containsKey(name)) {
					reads = mapped.get(name);
				} else {
					reads = new HashSet<String>();
				}
				BufferedReader br = new BufferedReader(new FileReader(f));
				for(String line = br.readLine(); line != null; line = br.readLine()) {
					if(!line.startsWith("@")) {
						String[] sp = line.split("\\s+");
						if(!sp[1].equals("4")) {//read is mapped
							//TODO do we need to filter on number of exact matches or match length?
							reads.add(sp[0]);
						}
					}
				}
				br.close();
				mapped.put(name, reads);
			}
		}
		return(mapped);
	}

}
