package cmcDistances.test;

import java.io.File;
import java.util.HashMap;

import cmcDistances.test.VerifyContif.Holder;

public class CompareTwoContexts
{
	public static void main(String[] args) throws Exception
	{
		if( args.length != 2)
			throw new Exception("Usage file1 file2");
		
		HashMap<String, Holder> map1 = VerifyContif.parseFile(
				new File(args[0]), 15, 15);
		
		HashMap<String, Holder> map2 = VerifyContif.parseFile(
				new File(args[1]), 15, 15);
		
		if( map1.size() != map2.size())
			throw new Exception("unequal sizes " + map1.size() + " " + map2.size());
		
		for(String s : map1.keySet())
		{
			if(  map1.get(s).equals(map2.get(s)) )
				throw new Exception("No " + map1.get(s) + " " + map2.get(s) );
		}
	}
	
	
}
