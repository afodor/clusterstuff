package kylie_2015;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import parsers.FastQ;

public class CountSeqs
{
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(
				new File("/projects/afodor_research/kylie_2015/counts.txt" )));
		
		writer.write("sample\tcount\n");
		
		for(String s : CreateRDPQSub.FASTQ_DIR.list())
		{
			System.out.println(s);
			File fastQFile = new File( CreateRDPQSub.FASTQ_DIR.getAbsolutePath() + File.separator 
						+ s);
			
			BufferedReader reader =
					new BufferedReader(new InputStreamReader( 
							new GZIPInputStream( new FileInputStream( fastQFile))));
			
			long count =0;
			
			for(FastQ fastq = FastQ.readOneOrNull(reader); fastq != null; 
					fastq = FastQ.readOneOrNull(reader))
			{
				count++;
			}
			
			writer.write(s + "\t" + count + "\n");
		}
		
		writer.flush();  writer.close();
		
	}
}
