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
		
		while( line != null)
		{
			if( line.startsWith(">"))
			{
				String name = new StringTokenizer(line.substring(1), "_").nextToken();
				
				BufferedReader brokenOutFasta = map.get(name);
				
				if( brokenOutFasta == null)
				{
					brokenOutFasta = new BufferedReader(new FileReader(new File(BreakOutBySample.OUT_SEQUENCE_DIR
							+File.separator + name)));
					map.put(name, brokenOutFasta);
				}
			}
		}
	}
	
}
