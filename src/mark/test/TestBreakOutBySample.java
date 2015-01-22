package mark.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.StringTokenizer;

import mark.BreakOutBySample;

public class TestBreakOutBySample
{
	public static void main(String[] args) throws Exception
	{
		HashMap<String, BufferedReader> map =new HashMap<String, BufferedReader>();
		
		BufferedReader inFasta = new BufferedReader(new BufferedReader(new FileReader(new File(
				BreakOutBySample.IN_SEQUENCE_DIR.getAbsolutePath() 
					+ File.separator + "mrg_j200p5_q20_take2_l25_r25.fna"))));
		
		String line = inFasta.readLine();
		BufferedReader lastReader = null;
		StringBuffer stringBuff= null;
		String lastHeader = null;
		
		while( line != null)
		{
			if( line.startsWith(">"))
			{
				if( stringBuff!= null)
				{
					String brokenheader = lastReader.readLine();
					
					if( ! brokenheader.startsWith(lastHeader))
						throw new Exception("Unexpected " + lastHeader + " " + brokenheader);
					
					String nextSeq = lastReader.readLine();
					
					if( ! stringBuff.toString().equals(nextSeq))
						throw new Exception("Mismatch " + stringBuff+ " " + nextSeq);
				}
				
				lastHeader = new StringTokenizer(line.substring(1), "_").nextToken();
				
				BufferedReader brokenOutFasta = map.get(lastHeader);
				
				if( brokenOutFasta == null)
				{
					brokenOutFasta = new BufferedReader(new FileReader(new File(BreakOutBySample.OUT_SEQUENCE_DIR
							+File.separator + lastHeader)));
					map.put(lastHeader, brokenOutFasta);
				}
				
				lastReader = brokenOutFasta;
				stringBuff = new StringBuffer();				
				line = inFasta.readLine();
			}
			else while( line != null && ! line.startsWith(">"))
			{
				stringBuff.append(line.trim());
			}
		}
		
		for(String s : map.keySet())
		{
			BufferedReader reader = map.get(s);
			
			if ( reader.readLine() != null)
				throw new Exception("Not at end of file " + s);
		}
		
		System.out.println("Global pass");

	}
	
	
}
