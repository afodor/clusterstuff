package creOrthologs.kmers;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

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
		
		HashMap<String, Integer> map = MakeKmers.breakIntoKmers(new File(inKmerFilePath));
		
		for(String s : map.keySet())
		{
			int aNum = map.get(s);
			
			for( int x=0; x < aNum; x++)
			{
				list.add(s);
			}
		}
		
		return list;
	}
}
