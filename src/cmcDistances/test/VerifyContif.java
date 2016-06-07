package cmcDistances.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;

public class VerifyContif
{
	static class Holder
	{
		int numA=0;
		int numC=0;
		int numG=0;
		int numT=0;
		
		void add(char c) throws Exception
		{
			if( c== 'A')
				numA++;
			else if ( c== 'C')
				numC++;
			else if ( c== 'G')
				numG++;
			else if ( c== 'T')
				numT++;
			else
				throw new Exception("No " + c);
		}
		
		@Override
		public boolean equals(Object obj)
		{
			Holder o = (Holder) obj;
			
			return this.numA == o.numA && this.numC == o.numC && this.numG == o.numG 
						&& this.numT == o.numT;
		}
		
		@Override
		public String toString()
		{
			return "[" + this.numA + "," + this.numC + "," + this.numG + "," + this.numT + "]";
		}
		
	}
	
	private static boolean isValid(String s) throws Exception
	{
		for( int x=0; x < s.length(); x++)
		{
			char c = s.charAt(x);
			
			if( c != 'A' && c != 'C' && c != 'T' && c != 'G')
				return false;
		}
		
		return true;
	}
	
	private static String reverseTranscribe( String s ) throws Exception
	{
		StringBuffer buff = new StringBuffer();
		
		for( int x=s.length()-1; x >= 0 ; x--)
			buff.append(flip(s.charAt(x)));
		
		return buff.toString();
	}
	
	private static char flip(char c) throws Exception
	{
		if( c == 'A')
			return 'T';
		
		if( c== 'T')
			return 'A';
		
		if( c == 'G')
			return 'C';
		
		if ( c== 'C')
			return 'G';
		
		throw new Exception("No " + c);
	}
	
	public static void main(String[] args) throws Exception
	{
		if( args.length != 4)
		{
			System.out.println("Usage inFile outFile leftIndexLength rightIndexLength");
			System.exit(1);
		}
		
		HashMap<String, Holder> map = parseFile(new File(args[0]), Integer.parseInt(args[2]), 
					Integer.parseInt(args[3]));
		
		writeResults(new File(args[1]), map);
	}
	
	private static void writeResults(File outFile, HashMap<String, Holder> map )
		throws Exception
	{
		BufferedWriter writer =new BufferedWriter(new FileWriter(outFile));
		
		for(String s : map.keySet())
		{
			writer.write(s + "\t" );
			Holder h = map.get(s);
			
			writer.write(h.numA + "\t" + h.numC + "\t" + h.numG + "\t" + h.numT + "\n");
		}
		
		writer.flush(); writer.close();
	}
	
	
	static HashMap<String, Holder> parseFile(File f, int leftIndex, int rightIndex) 
				throws Exception
	{
		HashMap<String, Holder> map = new HashMap<String,Holder>();
		
		BufferedReader reader = new BufferedReader(new FileReader(f));
		
		for(String s = reader.readLine(); s != null; s = reader.readLine())
		{
			String seq = reader.readLine();
			reader.readLine();  reader.readLine();
			
			for( int x=0; x < seq.length(); x++)
			{
				if( x + leftIndex + rightIndex + 1 <= seq.length() )
				{
					String leftKey = seq.substring(x, leftIndex+ x);
					char middle= seq.charAt(leftIndex + x);
					String rightKey = seq.substring(leftIndex + x + 1,
							x + leftIndex + rightIndex + 1 );
					
					String fullKey = leftKey + rightKey;
					
					if( isValid(fullKey))
					{
						Holder h = map.get(fullKey);
						
						if( h == null)
						{
							fullKey = reverseTranscribe(fullKey);
							middle = flip(middle);
							h = map.get(fullKey);
							
							if( h==null)
							{
								h = new Holder();
								map.put(fullKey, h );
							}
						}
						
						h.add(middle);
					}
				}
			}
		}
		
		reader.close();
		
		return map;
	}
}
