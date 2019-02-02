package cn.edu.hit.dnmfilterindel.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

public class MergeFiles {

	public static void merge(String outputFile, String[] files) throws IOException {
		PrintWriter writer = new PrintWriter(outputFile);
		writer.write("Father_Allele_Balance" + ","+ "Father_Read_Depth" + ",");
		writer.write("Father_Mean_Mapping_Quality_For_Ref" + "," + "Father_Mean_Mapping_Quality_For_Alt" + "," + "Father_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Father_Mean_Distance_To_Three_Prime_For_Alt" + ",");
		writer.write("Father_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Father_Fraction_Of_MQ0_Reads_For_Alt" + ",");
		writer.write("Father_Fraction_Of_Soft_Clipped_Reads_For_Ref" + "," + "Father_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
		writer.write("Father_Mean_Nearby_Mismatches_For_Ref" + "," + "Father_Mean_Nearby_Mismatches_For_Alt" + ",");
		writer.write("Father_Mean_Nearby_Indels_For_Ref" + "," + "Father_Mean_Nearby_Indels_For_Alt" + ",");
		writer.write("Father_Strand_Direction_For_Ref" + "," + "Father_Strand_Direction_For_Alt" + ",");
		writer.write("Father_Strand_Bias" + ",");
		writer.write("Mother_Allele_Balance" + ","+ "Mother_Read_Depth" + ",");
		writer.write("Mother_Mean_Mapping_Quality_For_Ref" + "," + "Mother_Mean_Mapping_Quality_For_Alt" + "," + "Mother_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Mother_Mean_Distance_To_Three_Prime_For_Alt" + ",");
		writer.write("Mother_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Mother_Fraction_Of_MQ0_Reads_For_Alt" + ",");
		writer.write("Mother_Fraction_Of_Soft_Clipped_Reads_For_Ref" + "," + "Mother_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
		writer.write("Mother_Mean_Nearby_Mismatches_For_Ref" + "," + "Mother_Mean_Nearby_Mismatches_For_Alt" + ",");
		writer.write("Mother_Mean_Nearby_Indels_For_Ref" + "," + "Mother_Mean_Nearby_Indels_For_Alt" + ",");
		writer.write("Mother_Strand_Direction_For_Ref" + "," + "Mother_Strand_Direction_For_Alt" + ",");
		writer.write("Mother_Strand_Bias" + ",");
		writer.write("Offspring_Allele_Balance" + ","+ "Offspring_Read_Depth" + ",");        
		writer.write("Offspring_Mean_Mapping_Quality_For_Ref" + "," + "Offspring_Mean_Mapping_Quality_For_Alt" + "," + "Offspring_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Offspring_Mean_Distance_To_Three_Prime_For_Alt" + ",");
		writer.write("Offspring_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Offspring_Fraction_Of_MQ0_Reads_For_Alt" + ",");
		writer.write("Offspring_Fraction_Of_Soft_Clipped_Reads_For_Ref" + "," + "Offspring_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
		writer.write("Offspring_Mean_Nearby_Mismatches_For_Ref" + "," + "Offspring_Mean_Nearby_Mismatches_For_Alt" + ",");
		writer.write("Offspring_Mean_Nearby_Indels_For_Ref" + "," + "Offspring_Mean_Nearby_Indels_For_Alt" + ",");
		writer.write("Offspring_Strand_Direction_For_Ref" + "," + "Offspring_Strand_Direction_For_Alt" + ",");
		writer.write("Offspring_Strand_Bias" + ",");
		writer.write("PValue_Father_To_Offspring" + "," + "PValue_Mother_To_Offspring" + ",");
		writer.write("Homopolymer_Flag" + "," + "Short_Tandem_Repeat_Flag" + "," + "Class" +"\n");
		for (String file : files) {
			BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(file)));
			String line = null;
			while ((line = bufferedReader.readLine()) != null && line.trim().length() > 0) {
				if (!line.startsWith("Father")) {
					writer.write(line + "\n");
				}
			}
			bufferedReader.close();
		}
		writer.close();
	}
}
