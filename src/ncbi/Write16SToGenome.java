package ncbi;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;

import parsers.FastaSequence;

public class Write16SToGenome
{	
	// ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.frn.tar.gz on 11/29/2014
	public static void main(String[] args) throws Exception
	{
		BufferedWriter writer =new BufferedWriter(new FileWriter(new File("/nobackup/afodor_research/ncbi/geneIDToGenome.txt")));
		
		File topDir = new File("/nobackup/afodor_research/ncbi");
		
		String[] list = topDir.list();
		
		for(String s : list)
		{
			File f = new File(topDir.getAbsolutePath() + File.separator + s);
			
			if( f.isDirectory())
			{
				String[] sublist = f.list();
				
				for(String s2 : sublist)
				{
					List<FastaSequence> fastaList = FastaSequence.readFastaFile(f.getAbsolutePath() + File.separator + s2);
					
					for(FastaSequence fs : fastaList)
					{
						if(fs.getHeader().indexOf("16S") != -1)
						{
							writer.write(">" + fs.getHeader()+ "\t");
							writer.write(f.getName() + "\n");
						}
					}
					
					writer.flush();
				}
			}
		}
				
		
	    writer.flush();  writer.close();
	}
}
