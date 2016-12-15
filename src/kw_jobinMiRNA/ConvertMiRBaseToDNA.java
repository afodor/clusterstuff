/*
 * convert Us to Ts (from RNA to DNA)
 */
package kw_jobinMiRNA;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.Semaphore;

public class ConvertMiRBaseToDNA {
	public static String MATREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/mature.fa";
	public static String PINREF = "/nobackup/afodor_research/kwinglee/mirbase_v21/hairpin.fa";

	public static void main(String[] args) throws InterruptedException {
		String[] dbs = new String[]{MATREF, PINREF};
		Semaphore s = new Semaphore(dbs.length);
		for(String d : dbs) {
			s.acquire();
			Convert con = new Convert(d, s);
			new Thread(con).start();
		}
		for(int i = 0; i < dbs.length; i++) {
			s.acquire();
		}
	}

}

class Convert implements Runnable {
	private final String db;
	private final Semaphore semaphore;
	
	public Convert(String db, Semaphore s) {
		this.db = db;
		semaphore = s;
	}

	@Override
	public void run() {
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(db)));
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(db.replace(".fa", ".dna.fasta"))));
			for(String line = in.readLine(); line != null; line = in.readLine()) {
				if(line.startsWith(">")) {
					out.write(line + "\n");
				} else {
					out.write(line.toUpperCase().replaceAll("U", "T") + "\n");
				}
			}
			in.close();
			out.close();
		} catch(IOException e) {
			System.out.println(e.getMessage());
		} finally {
			semaphore.release();
		}
	}

}
