package chs_snps;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import coPhylog.CoPhylogBinaryFileReader;
import coPhylog.ContextCount;

public class Test {
	public static String conversionFile = "/projects/afodor_research/mjzapata/CRE/CHS_raw/chs_batch_download_results.csv";//file containing the conversion
	public static String outDir = "/projects/afodor_research/kwinglee/cophylog_all80chs/test/";//Name of directory to write results to
	public static String contextDir = "/projects/afodor_research/kwinglee/cophylog_all80chs/context/";//path to all context files
	public static int MIN_READS = 1; //minimum number of reads to keep key
	
	public static void main(String[] args) throws Exception {
		if(args.length != 1) {
			System.err.println("Usage: 1=write gzip file, 2=write unzipped file, 3=read gzip file, 4=read unzipped file");
			System.exit(1);
		}
		String strain = "CHS01";
		if(args[0].equals("1")) {
			writeZip(strain);
		} else if (args[0].equals("2")) {
			writeUnzip(strain);
		} else if(args[0].equals("3")) {
			read(outDir + strain + "_context.gz");
		} else if(args[0].equals("4")) {
			read(outDir + strain + "_context_unzip.gz");
		} else {
			System.err.println("Usage: 1=write gzip file, 2=write unzipped file, 3=read gzip file, 4=read unzipped file");
			System.exit(1);
		}
	}
	
	public static void read(String file) throws Exception {
		HashMap<Long, ContextCount> map1 = 
				CoPhylogBinaryFileReader.readBinaryFileRequireMin(new File(file), MIN_READS);
		int tot = 0;
		for(Long k : map1.keySet()) {
			ContextCount cc = map1.get(k);
			tot += cc.getSum();
		}
		System.out.println(tot);
	}
	
	public static void writeZip(String strain) throws IOException {
		BufferedReader convert = new BufferedReader(new FileReader (new File(conversionFile)));
		String line = convert.readLine();
		while(line != null) {
			String[] csp = line.split("\t");
			String chs;
			if(csp[0].length() == 1) {
				chs = "CHS0" + csp[0];
			} else {
				chs = "CHS" + csp[0];
			}
			if(chs.equals(strain)) {
				//set up hash
				HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();

				//combine files into hash
				String[] srrlist = csp[1].replace("[", "").replace("]", "").split(",");
				for(int i = 0; i < srrlist.length; i++) {//for each file
					for(int j = 1; j <= 2; j++) {//forward and reverse reads
						try {
							//read file
							HashMap<Long, ContextCount> m = CoPhylogBinaryFileReader.readBinaryFileRequireMin(new File(contextDir+"context"+srrlist[i].trim()+"_"+j+"_context.gz"), 2);

							System.out.println(srrlist[i].trim()+"_"+j);
							
							//merge maps
							for(Long key : m.keySet()) {
								ContextCount con = m.get(key);
								if(map.containsKey(key)) {
									ContextCount mapcon = map.get(key);
									mapcon.increaseA(con.getNumA());
									mapcon.increaseT(con.getNumT());
									mapcon.increaseC(con.getNumC());
									mapcon.increaseG(con.getNumG());
									map.put(key, mapcon);
								} else {
									map.put(key, con);
								}
							}
						} catch (Exception e) {
							System.err.println(e);
						}


						//update log
						System.gc();
					}
				}

				//write results
				DataOutputStream out =new DataOutputStream( new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(outDir + chs + "_context.gz"))));

				out.writeInt(map.size());

				for( Long l : map.keySet() ) {

					ContextCount cc = map.get(l);
					int sum = cc.getSum();
					if(sum > MIN_READS) {
						out.writeLong(l);
						out.writeByte(cc.getAAsByte());
						out.writeByte(cc.getCAsByte());
						out.writeByte(cc.getGAsByte());
						out.writeByte(cc.getTAsByte());
					}
				} 
				out.close();
			}
			
			line = convert.readLine();
			
		}
		convert.close();
	}

	public static void writeUnzip(String strain) throws IOException {
		BufferedReader convert = new BufferedReader(new FileReader (new File(conversionFile)));
		String line = convert.readLine();
		while(line != null) {
			String[] csp = line.split("\t");
			String chs;
			if(csp[0].length() == 1) {
				chs = "CHS0" + csp[0];
			} else {
				chs = "CHS" + csp[0];
			}
			if(chs.equals(strain)) {
				//set up hash
				HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();

				//combine files into hash
				String[] srrlist = csp[1].replace("[", "").replace("]", "").split(",");
				for(int i = 0; i < srrlist.length; i++) {//for each file
					for(int j = 1; j <= 2; j++) {//forward and reverse reads
						try {
							//read file
							HashMap<Long, ContextCount> m = CoPhylogBinaryFileReader.readBinaryFileRequireMin(new File(contextDir+"context"+srrlist[i].trim()+"_"+j+"_context.gz"), 2);

							System.out.println(srrlist[i].trim()+"_"+j);
							
							//merge maps
							for(Long key : m.keySet()) {
								ContextCount con = m.get(key);
								if(map.containsKey(key)) {
									ContextCount mapcon = map.get(key);
									mapcon.increaseA(con.getNumA());
									mapcon.increaseT(con.getNumT());
									mapcon.increaseC(con.getNumC());
									mapcon.increaseG(con.getNumG());
									map.put(key, mapcon);
								} else {
									map.put(key, con);
								}
							}
						} catch (Exception e) {
							System.err.println(e);
						}


						//update log
						System.gc();
					}
				}

				//write results
				DataOutputStream out =new DataOutputStream(new FileOutputStream(outDir + chs + "_context_unzip"));

				out.writeInt(map.size());

				for( Long l : map.keySet() ) {

					ContextCount cc = map.get(l);
					int sum = cc.getSum();
					if(sum > MIN_READS) {
						out.writeLong(l);
						out.writeByte(cc.getAAsByte());
						out.writeByte(cc.getCAsByte());
						out.writeByte(cc.getGAsByte());
						out.writeByte(cc.getTAsByte());
					}
				} 
				out.close();
			}
			
			line = convert.readLine();
			
		}
		convert.close();
	}

}
