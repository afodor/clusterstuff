package creOrthologs.kmers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import utils.Avevar;

public class AverageKmers
{
	public static void main(String[] args) throws Exception
	{
		List<Integer> list = new ArrayList<Integer>();
		
		for(String s : GatherDistanceMatrix.GATHERED_DIR.list())
			if( s.startsWith("klebsiella_pneumoniae") && s.endsWith("Key.txt"))
			{
				BufferedReader reader = new BufferedReader(new FileReader(new File(
						GatherDistanceMatrix.GATHERED_DIR.getAbsolutePath() + File.separator + s)));
				
				reader.readLine();
				
				for(String s2 = reader.readLine(); s2 != null; s2 = reader.readLine())
				{
					String[] splits = s2.split("\t");
					if( splits.length != 3)
						throw new Exception("No");
					
					list.add(Integer.parseInt(splits[2]));
				}
				
				reader.close();
			}
		
		Avevar av = new Avevar(list);
		System.out.println(av.getAve() + " " + av.getSD());
	}
}
