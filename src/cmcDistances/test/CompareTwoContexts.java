package cmcDistances.test;

import java.io.File;
import java.util.HashMap;

import cmcDistances.test.VerifyContif.Holder;
import utils.Translate;

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
			Holder h1 = map1.get(s);
			
			if( h1 == null)
				throw new Exception("No");
			
			Holder h2 = map2.get(2);
			
			if( h2 == null)
			{
				String rev = Translate.reverseTranscribe(s);
				
				h2 = map2.get(rev);
				
				if( h2 == null)
					throw new Exception("Could not find " + s + " " + rev);
				
				h2.flip();
				
			}
			
			if(  ! h1.equals(h2) )
				throw new Exception("No " + h1 + " " + h2 );
		}
		
		System.out.println("passed " + args[0] + " " + args[1] 
				+ " " + map1.size() + " " + map2.size());
	}
	
	
}
