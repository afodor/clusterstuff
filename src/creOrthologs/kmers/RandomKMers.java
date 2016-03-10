package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.StringTokenizer;

public class RandomKMers
{
	public static void main(String[] args) throws Exception
	{
		if( args.length != 3)
		{
			System.out.println("Usage numKMers inKmerPath fileOut");
			System.exit(1);
		}
		
		List<String> kmerList = getKmerList(args[1]);
		
		System.out.println("Got " + kmerList.size());
		
		Collections.shuffle(kmerList);
		
		System.out.println("Shuffled");
		
	}
	
	private static List<String> getKmerList(String inKmerFilePath) throws Exception
	{
		List<String> list = new ArrayList<String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(
				inKmerFilePath)));

		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			StringTokenizer sToken = new StringTokenizer(s);
			
			String seq = sToken.nextToken();
			int aNum = Integer.parseInt(sToken.nextToken());
			
			if( sToken.hasMoreTokens())
				throw new Exception("No");
			
			for(int x=0; x < aNum; x++)
				list.add(seq);
		}
		
		reader.close();
		
		return list;
	}
}
