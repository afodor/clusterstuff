package chs_snps;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;


public class WriteScriptsForComparePairs
{
	private static List<File> getFiles() throws Exception
	{
		List<File> list = new ArrayList<File>();
		
		for( int x = 1; x <= 8; x++)
		{
			list.add(new File("/projects/afodor_chs/context/chs_241_"+ x + "_context.gz"));
		}
		
		for( int x = 10; x <= 16; x++)
		{
			list.add(new File("/projects/afodor_chs/context/chs_242_"+ x + "_context.gz"));
		}
		
		
		return list;
	}
	
	private static String firstThreeTokensOfName(File f)
	{
		StringTokenizer sToken = new StringTokenizer(f.getName(), "_");
		
		return sToken.nextToken() + "_" + sToken.nextToken() + "_" + sToken.nextToken();
	}
	
	public static void main(String[] args) throws Exception
	{
		BufferedWriter allWriter = new BufferedWriter(new FileWriter(new File("/projects/afodor_chs/runCompare/runAll.sh")));
		
		int index=1;
		List<File> list = getFiles();
		
		for( int x=0; x  < list.size() -1; x++)
		{
			File xFile = list.get(x);
			
			for( int y=x+1 ; y < list.size(); y++)
			{
				File yFile = list.get(y);
				
				File outFile = new File("/projects/afodor_chs/compare/" + firstThreeTokensOfName(xFile) + "_to_" 
						+ firstThreeTokensOfName(yFile) + "_compare.txt");
				
				File outScriptFile = new File("/projects/afodor_chs/runCompare/runC_" + index);
				allWriter.write("qsub -q \"viper\" " + outScriptFile.getName() +  "\n");
				
				BufferedWriter scriptWriter = new BufferedWriter(new FileWriter(outScriptFile));
				scriptWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin -mx2800m " + 
				"coPhylog.WriteSNPFile " + xFile.getAbsolutePath()  + " " + yFile.getAbsolutePath() + " " + 
							outFile.getAbsolutePath());
				
				scriptWriter.flush(); scriptWriter.close();
				
				index++;
			}
		}
		
		allWriter.flush();  allWriter.close();
	}
}
