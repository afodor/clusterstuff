package chs_snps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.StringTokenizer;

public class WriteSummaryFile
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter("D:\\MelHospital\\snpSummary.txt"));
		
		writer.write("file\tcategory\tnumberOfSnps\taverageSupportingReads\n");
		
		File topDir = new File("D:\\MelHospital\\comparison");
		
		for(String s : topDir.list())
		{
			File toParse = new File(topDir.getAbsolutePath() + File.separator + s);
			
			BufferedReader reader = new BufferedReader(new FileReader(toParse));
			
			int numLines =0;
			
			reader.readLine();
			double sum =0;
			
			for(String s2= reader.readLine() ; s2 != null; s2 = reader.readLine())
			{
				String[] splits = s2.split("\t");
				sum += getSum(splits[1]) + getSum(splits[2]);
				numLines++;
			}
			
			writer.write(  toParse.getName() + "\t" + getCategory(toParse.getName()) + "\t" + numLines + "\t" 
							+ (0.5*sum/numLines) + "\n");
		}
		
		writer.flush();  writer.close();
	}
	
	private static int getSum(String s) throws Exception
	{
		int sum =0;
		
		StringTokenizer sToken = new StringTokenizer(s, "[,]");
		
		for( int x=0; x < 4; x++)
			sum += Integer.parseInt(sToken.nextToken());
		
		if( sToken.hasMoreTokens() )
			throw new Exception("No");
		
		return sum;
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
