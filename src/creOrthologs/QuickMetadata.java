package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class QuickMetadata
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				 "/nobackup/afodor_research/af_broad/strainMeta.txt")));
		
		String[] directories = 
			{
					"/nobackup/afodor_research/af_broad/carolina",
					"/nobackup/afodor_research/af_broad/resistant",
					"/nobackup/afodor_research/af_broad/susceptible"
			};
		
		for(String s : directories)
		{
			File dir = new File(s);
			
			String[] names = dir.list();
			
			for(String name : names)
			{
				if( name.endsWith(".fasta"))
				{
					writer.write(name + "\t" + dir.getName() + "\n");
					writer.flush();
				}
			}
		}
		
		writer.flush();  writer.close();
	}
}
