package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class WriteSummaryFile
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("D:\\MelHospital\\snpSummary.txt"));
		
		writer.write("file\tcategory\tnumberOfSnps\n");
		
		File topDir = new File("D:\\MelHospital\\comparison");
		
		for(String s : topDir.list())
		{
			File toParse = new File(topDir.getAbsolutePath() + File.separator + s);
			
			BufferedReader reader = new BufferedReader(new FileReader(toParse));
			
			int numLines =0;
			
			reader.readLine();
			
			while( reader.readLine() != null)
				numLines++;
			
			writer.write(  toParse.getName() + "\t" + getCategory(toParse.getName()) + "\t" + numLines + "\n" );
		}
		
		writer.flush();  writer.close();
	}
	
	private static String getCategory(String name) throws Exception
	{
		if( name.indexOf("241_") == -1 && name.indexOf("242_") == -1 )
			throw new Exception("No");
		
		if( name.indexOf("241_") != -1 && name.indexOf("242_") != -1 )
			return "Cross";
		
		if( name.indexOf("241_") == -1 && name.indexOf("242_") != -1 )
			return "Chs242";
		
		if( name.indexOf("241_") != -1 && name.indexOf("242_") == -1 )
			return "Chs241";

		throw new Exception("Logic error");
	}
}
