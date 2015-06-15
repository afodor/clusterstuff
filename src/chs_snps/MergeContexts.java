package chs_snps;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import coPhylog.ContextCount;
import coPhylog.CoPhylogBinaryFileReader;

public class MergeContexts {
	public static String conversionFile = "/projects/afodor_research/mjzapata/CRE/CHS_raw/chs_batch_download_results.csv";//file containing the conversion
	public static String outDir = "/projects/afodor_research/kwinglee/cophylog_all80chs/contextCombined/";//Name of directory to write results to
	public static String contextDir = "/projects/afodor_research/kwinglee/cophylog_all80chs/context/";//path to all context files
	
	public static void main(String[] args) throws Exception {
		BufferedReader convert = new BufferedReader(new FileReader (new File(conversionFile)));
		String cline = convert.readLine();
		while(cline != null) {
			String[] csp = cline.split("\t");
			String chs;
			if(csp[0].length() == 1) {
				chs = "CHS0" + csp[0];
			} else {
				chs = "CHS" + csp[0];
			}
			
			//set up hash
			HashMap<Long, ContextCount> map = new HashMap<Long, ContextCount>();
			
			//combine files into hash
			String[] srrlist = csp[1].replace("[", "").replace("]", "").split(",");
			for(int i = 0; i < srrlist.length; i++) {//for each file
				for(int j = 1; j < 2; j++) {//forward and reverse reads
					//read file
					HashMap<Long, ContextCount> m = CoPhylogBinaryFileReader.readBinaryFile(new File(contextDir+"context"+srrlist[i]+"_"+j+"_context.gz"));
					
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
				}
			}
			
			//write results
			DataOutputStream out =new DataOutputStream( new BufferedOutputStream(new GZIPOutputStream(new FileOutputStream(outDir + chs + "_context.gz"))));
			
			out.writeInt(map.size());
			
			for( Long l : map.keySet() )
			{
				out.writeLong(l);
				
				ContextCount cc = map.get(l);
				
				out.writeByte( cc.getAAsByte() );
				out.writeByte( cc.getCAsByte());
				out.writeByte(cc.getGAsByte());
				out.writeByte( cc.getTAsByte());
			}
			
			out.flush();  out.close();
			
			
			cline = convert.readLine();
		}
		convert.close();
	}

}
