/**
 * log normalize RDP data
 */
package kw_jobinGA;

import java.io.File;

import parsers.NewRDPParserFileLine;
import parsers.OtuWrapper;
import utils.ConfigReader;

public class logNormalize {
	public static void main(String[] args) throws Exception {
		for( int x=1; x < NewRDPParserFileLine.TAXA_ARRAY.length; x++){
			System.out.println(NewRDPParserFileLine.TAXA_ARRAY[x]);
			OtuWrapper wrapper = new OtuWrapper(ConfigReader.getJobinGAStoolRDPDir() +
					File.separator + "rdp_taxaAsCol_" + NewRDPParserFileLine.TAXA_ARRAY[x] +
					".txt");
			
			wrapper.writeNormalizedLoggedDataToFile(ConfigReader.getJobinGAStoolRDPDir() +
					File.separator + "rdp_taxaAsCol_logNorm_" + NewRDPParserFileLine.TAXA_ARRAY[x] + 
					".txt");
			
			
		}
	}

}
