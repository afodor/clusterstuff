package kw_topeAdenoma;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

import utils.Translate;

public class SampleMap
{
	/*
	 * The key is the reverse transcription of the 4th column + "@" + the 9th column
	 * 
	 */
	public static HashMap<String, String> getPrimersToSampleMap() throws Exception
	{
		String dir = "/projects/afodor_research/kwinglee/tope/adenoma/";
		HashMap<String, String> map= new HashMap<String, String>();

		BufferedReader reader = new BufferedReader(new FileReader(new File(
				dir + "libraryMetadataForDemultiplexing.txt")));
		
		reader.readLine();
		for(String s = reader.readLine() ; s != null; s = reader.readLine())
		{
			String[] splits = s.split("\t");
			
			String key = Translate.reverseTranscribe(splits[4]) + "@" + 
								splits[9];
			
			if( map.containsKey(key)) {
				reader.close();
				throw new Exception("Duplicate key");
			}
			
			map.put(key, new String( splits[14])); 
			
		}
		
		reader.close();
		return map;
	}
	
	public static void main(String[] args) throws Exception
	{
		HashMap<String, String> map = getPrimersToSampleMap();
		
		for(String s : map.keySet())
			System.out.println(s +  " " + map.get(s));
	}
	
}
