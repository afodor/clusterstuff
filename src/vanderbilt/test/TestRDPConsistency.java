package vanderbilt.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.StringTokenizer;

public class TestRDPConsistency
{
	public static void main(String[] args) throws Exception
	{
		File topDir = new File("/projects/afodor/vanderbilt/rdpResults");
		int numDups =0 ;
		
		for(String s : topDir.list())
		{
			System.out.println(s);
			File rdpFile = new File(topDir.getAbsolutePath() + File.separator + s);
			
			String expected = new StringTokenizer(s, "_").nextToken();
			
			HashSet<String> set = new HashSet<String>();
			
			BufferedReader reader = new BufferedReader(new FileReader(rdpFile));
			
			for(String s2= reader.readLine(); s2 != null; s2 = reader.readLine())
			{
				String name = new StringTokenizer(s2).nextToken();
				
				if( set.contains(name))
					throw new Exception("Duplciate " + name);
				
				set.add(name);
				
				if( ! name.startsWith(expected))
				{
					System.out.println("Mismatch " + name + " " +expected);
					numDups++;
				}
					
			}
			
			reader.close();
			System.out.println("Passed " + s + " " + set.size());
		}
		
		System.out.println("Passed all with " + numDups + " duplicates ");
	}
}
