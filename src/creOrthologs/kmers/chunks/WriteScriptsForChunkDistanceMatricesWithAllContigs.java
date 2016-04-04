package creOrthologs.kmers.chunks;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;

import creOrthologs.kmers.GatherDistanceMatrix;
import creOrthologs.kmers.WriteSpearmanFromRandom;
import parsers.FastaSequence;

public class WriteScriptsForChunkDistanceMatricesWithAllContigs
{
	private static final File CHUNK_KMER_SCRIPT_DIR = 
			new File( "/nobackup/afodor_research/af_broad/scripts/runKmerChunks");
	
	public static final String GENOME_PATH= 
			"/nobackup/afodor_research/af_broad/carolina/klebsiella_pneumoniae_chs_11.0.scaffolds.fasta";
	

	public static final String LOG_PATH= 
			"/nobackup/afodor_research/af_broad/chunkLog.txt";
	
	
	public static final String CHUNCK_FILE_PATH = 
			"/nobackup/afodor_research/af_broad/pneuOnlyChunksWithContigs_0.85_0.9.txt";
	
	private static void writeOneIfNotThere(BufferedWriter allWriter, String contig,
				int startPos, int endPos, int index, BufferedWriter logWriter, String type) 
		throws Exception
	{
		File shFile = new File(
				CHUNK_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + "run_" + index + ".txt"	);
		

		String outFileBase = GENOME_PATH.substring(GENOME_PATH.lastIndexOf("/")+1)
				.replace(".scaffolds.fasta", "") + "_" + contig + "_" + startPos + "_" + endPos;
		
		File outFile = new File(
				GatherDistanceMatrix.GATHERED_DIR.getAbsolutePath() + File.separator + outFileBase + "_dist.txt");
		
		if( ! outFile.exists() || WriteSpearmanFromRandom.getValsOrNull(outFile, null) == null)
		{
			BufferedWriter aWriter = new BufferedWriter(new FileWriter(shFile));
			
			allWriter.write("qsub -q \"viper\" " +  shFile.getAbsolutePath() + "\n");
			
			aWriter.write("java -cp /users/afodor/gitInstall/clusterstuff/bin "
					+ "creOrthologs.kmers.ConstrainKMersToRegion "  + 
					GENOME_PATH+ " " + contig+ " " +  startPos+ " " 
					+ endPos+ "\n");
			
			allWriter.flush(); 
			
			aWriter.flush();  aWriter.close();
			System.out.println(contig+ " " +  startPos+ " " + endPos);
			
			logWriter.write(outFileBase + "\t" + contig + "\t" + 
					startPos + "\t" + endPos + "\t" + (endPos - startPos) + "\t" + 
						"false" + "\t" + type + "\n");
		}
		else
		{
			System.out.println("found " + outFile.getAbsolutePath() + " skipping ");
			logWriter.write(outFileBase + "\t" + contig + "\t" + 
					startPos + "\t" + endPos + "\t" + (endPos - startPos) + "\t" + 
						"true" + "\t" + type + "\n");
		}
		
		logWriter.flush();
	}
	
	private static class Holder
	{
		private final int start;
		private final int end;
		
		Holder( int start, int end)
		{
			this.start = start;
			this.end = end;
		}
	}
	
	private static HashMap<String, List<Holder>> getAsContigMap() throws Exception
	{
		HashMap<String, List<Holder>>  map = new LinkedHashMap<String,List<Holder>>();
		
		BufferedReader reader = new BufferedReader(new FileReader(new File(
				CHUNCK_FILE_PATH)));
		
		reader.readLine();
		
		for(String s= reader.readLine(); s != null; s= reader.readLine())
		{
			String[] splits = s.split("\t");
			
			if( splits.length != 5)
				throw new Exception("No");
			
			String contig = splits[0].replace("\"", "").replace("Contig_", "");
			
			List<Holder> list = map.get(contig);
			
			if( list == null)
			{
				list = new ArrayList<Holder>();
				map.put(contig, list);
			}
			
			Holder h = new Holder(Integer.parseInt(splits[1]), Integer.parseInt(splits[2]));
			list.add(h);
		}
		
		reader.close();
		
		return map;
	}
	
	public static void main(String[] args) throws Exception
	{
		
		BufferedWriter logWriter = new BufferedWriter(new FileWriter(new File(
				LOG_PATH)));
		
		logWriter.write("basePath\tcontig\tstart\tstop\tsize\tfound\ttype\n");
		
		BufferedWriter allWriter  = new BufferedWriter(new FileWriter(new File(
			CHUNK_KMER_SCRIPT_DIR.getAbsolutePath() + File.separator + 
			"runAll.sh")));
		
		HashMap<String, List<Holder>> contigMap = getAsContigMap();
		
		HashMap<String, FastaSequence> fastaMap = FastaSequence.getFirstTokenSequenceMap(GENOME_PATH);
		
		HashSet<String> includedContigs = new HashSet<String>();
		
		int index = 1;
		for(String contig : contigMap.keySet())
		{
			List<Holder> list = contigMap.get(contig);
			includedContigs.add(contig);
			
			if( list.size() == 1)
			{
				Holder h = list.get(0);
				writeOneIfNotThere(allWriter, contig, h.start, h.end, index,
							logWriter, "singleton");
				index++;
			}
			else
			{
				int listIndex =0;
				int contigLength = fastaMap.get(contig).getSequence().length() -1;
				
				if( list.get(0).start >0)
				{
					writeOneIfNotThere(allWriter, contig, 0, list.get(0).start -1000, index,
							logWriter, "initialBaseline");
					index++;
				}
					
				while(listIndex < list.size() -1)
				{
					
					writeOneIfNotThere(allWriter, contig, 
							list.get(listIndex).start, list.get(listIndex).end, index,
							logWriter, "peak");
					
					index++;
					listIndex++;
					
					if( list.get(listIndex-1).end +5000 < contigLength )
					{
						String type = "endBaseline";
						int endPos = contigLength;
						
						if( listIndex < list.size() -1 )
						{
							endPos = list.get(listIndex).start -1000;
							type = "baseline";
						}
						
						if( endPos - list.get(listIndex-1).end+1000 >= 5000)
						{

							writeOneIfNotThere(allWriter, contig, list.get(listIndex-1).end+1000, 
											endPos, index, logWriter,type);
							
							index++;

						}							
					}
					
				}
			}	
		}
		
		for(String s : fastaMap.keySet())
			if( ! includedContigs.contains(s))
			{
				writeOneIfNotThere(allWriter, s, 0, 
						fastaMap.get(s).getSequence().length()-1, index, logWriter,"singleton");
				
				index++;
			}
		
		logWriter.flush();  logWriter.close();
		allWriter.flush();  allWriter.close();
	}
}
