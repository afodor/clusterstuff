package creOrthologs;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class WriteBroadCategories
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer =new BufferedWriter(new FileWriter(new File(
			"/projects/afodor_research/af_broad/broadCategories.txt")));
		
		writer.write("genome\tcategory\n");
		
		for(String d : RunBlastAll.DIRECTORIES)
		{
			addDirectory(writer, new File( "/projects/afodor_research/af_broad/" + d ), d);
		}
		
		writer.flush();  writer.close();
	}
	
	private static void addDirectory( BufferedWriter writer, File directory , String annotation)
		throws Exception
	{
		String[] list = directory.list();
		
		for(String s : list)
		{
			if (s.endsWith("fasta"))
			{
				writer.write(s.replace(".scaffolds.fasta", "") + "\t" + 
									annotation + "\n");
			}
		}
		
		writer.flush();
	}
}
