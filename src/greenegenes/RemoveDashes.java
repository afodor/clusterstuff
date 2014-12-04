package greenegenes;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;


public class RemoveDashes
{
	public static void main(String[] args) throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				"/projects/afodor/greengenes/isolated_named_strains_gg16S_aligned.fasta")));
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(new File(
				"/projects/afodor/greengenes/isolated_named_strains_gg16S.fasta"
				)));
		
		for(String s= reader.readLine(); s != null; s = reader.readLine())
		{
			writer.write(reader.readLine().replaceAll("-", "") + "\n");
		}
		
		writer.flush();  writer.close();
		reader.close();
	}
}
