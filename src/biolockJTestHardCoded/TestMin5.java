package biolockJTestHardCoded;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.StringTokenizer;

public class TestMin5
{
	private static void getMapForDirectory( File inFile, HashMap<String,Long> innerMap, String level ) throws Exception
	{
		System.out.println("check " + inFile.getAbsolutePath());
		BufferedReader reader = new BufferedReader(new FileReader(inFile));
		
		int checked =0;
		
		for(String s= reader.readLine(); s!= null; s= reader.readLine() )
		{
			StringTokenizer sToken = new StringTokenizer(s, "\t");
			
			if( sToken.countTokens() != 2)
				throw new Exception("Expecting two tokens");
			
			String name = sToken.nextToken();
			
			if( name.toLowerCase().indexOf("unclassified") == -1)
			{
				name = name.substring(name.indexOf(level + "__"), name.length());
				name = name.replace(level + "__", "");
				
				Long parsedVal = innerMap.get(name);
				
				if( parsedVal == null)
					throw new Exception("Could not find " + name);
				
				Long thisVal = Long.parseLong(sToken.nextToken());
				
				if( ! thisVal.equals(parsedVal))
					throw new Exception("Mismatched for " + name + " " + parsedVal + " " + thisVal);
				
				checked++;
				
				innerMap.remove(name);
			}
		}
		
		System.out.println("Ok checked " + checked);
		
		if( innerMap.size() != 0)
			System.out.println("COULD NOT FIND  " + innerMap);
	}
}
