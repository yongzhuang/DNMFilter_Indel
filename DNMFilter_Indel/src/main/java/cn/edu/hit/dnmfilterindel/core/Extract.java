/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cn.edu.hit.dnmfilterindel.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.samtools.SAMFileReader;
import org.apache.log4j.Logger;

import cn.edu.hit.dnmfilterindel.file.ConfigReader;
import cn.edu.hit.dnmfilterindel.file.DNMReader;
import cn.edu.hit.dnmfilterindel.file.DNMRecord;
import cn.edu.hit.dnmfilterindel.file.PEDReader;
import cn.edu.hit.dnmfilterindel.file.Trio;
import cn.edu.hit.dnmfilterindel.sam.MultiPileup;
import cn.edu.hit.dnmfilterindel.sam.MultiSamLocusIterator;

/**
 *
 * @author Yongzhuang Liu
 */
public class Extract {

	private static Logger logger = Logger.getLogger(Extract.class);
	private final static int FATHER_INDEX = 0;
	private final static int MOTHER_INDEX = 1;
	private final static int OFFSPRING_INDEX = 2;
	private Properties properties;
	File referenceSequenceFile;
	private IndexedFastaSequenceFile indexfasta;
	List<Trio> trios;
	Properties bams;
	String outputPath;
	String positivePath;
	String negativePath;

	public Extract(Properties properties) throws IOException {
		this.properties = properties;
		this.referenceSequenceFile = new File(properties.getProperty("reference"));
		this.trios = (new PEDReader(properties.getProperty("pedigree"))).getTrios();
		this.bams = (new ConfigReader(properties.getProperty("bam"))).parse();
		this.positivePath = properties.getProperty("positive");
		this.negativePath = properties.getProperty("negative");
		this.outputPath = properties.getProperty("output");
		this.indexfasta = new IndexedFastaSequenceFile(new File(properties.getProperty("reference")));
	}

	public void run() throws IOException {

		PrintWriter output = new PrintWriter(outputPath);
		output.write("Father_Allele_Balance" + "," + "Father_Read_Depth" + ",");
		output.write("Father_Mean_Mapping_Quality_For_Ref" + "," + "Father_Mean_Mapping_Quality_For_Alt" + ","
				+ "Father_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Father_Mean_Distance_To_Three_Prime_For_Alt"
				+ ",");
		output.write("Father_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Father_Fraction_Of_MQ0_Reads_For_Alt" + ",");
		output.write("Father_Fraction_Of_Soft_Clipped_Reads_For_Ref" + ","
				+ "Father_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
		output.write("Father_Mean_Nearby_Mismatches_For_Ref" + "," + "Father_Mean_Nearby_Mismatches_For_Alt" + ",");
		output.write("Father_Mean_Nearby_Indels_For_Ref" + "," + "Father_Mean_Nearby_Indels_For_Alt" + ",");
		output.write("Father_Strand_Direction_For_Ref" + "," + "Father_Strand_Direction_For_Alt" + ",");
		output.write("Father_Strand_Bias" + ",");
		output.write("Mother_Allele_Balance" + "," + "Mother_Read_Depth" + ",");
		output.write("Mother_Mean_Mapping_Quality_For_Ref" + "," + "Mother_Mean_Mapping_Quality_For_Alt" + ","
				+ "Mother_Mean_Distance_To_Three_Prime_For_Ref" + "," + "Mother_Mean_Distance_To_Three_Prime_For_Alt"
				+ ",");
		output.write("Mother_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Mother_Fraction_Of_MQ0_Reads_For_Alt" + ",");
		output.write("Mother_Fraction_Of_Soft_Clipped_Reads_For_Ref" + ","
				+ "Mother_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
		output.write("Mother_Mean_Nearby_Mismatches_For_Ref" + "," + "Mother_Mean_Nearby_Mismatches_For_Alt" + ",");
		output.write("Mother_Mean_Nearby_Indels_For_Ref" + "," + "Mother_Mean_Nearby_Indels_For_Alt" + ",");
		output.write("Mother_Strand_Direction_For_Ref" + "," + "Mother_Strand_Direction_For_Alt" + ",");
		output.write("Mother_Strand_Bias" + ",");
		output.write("Offspring_Allele_Balance" + "," + "Offspring_Read_Depth" + ",");
		output.write("Offspring_Mean_Mapping_Quality_For_Ref" + "," + "Offspring_Mean_Mapping_Quality_For_Alt" + ","
				+ "Offspring_Mean_Distance_To_Three_Prime_For_Ref" + ","
				+ "Offspring_Mean_Distance_To_Three_Prime_For_Alt" + ",");
		output.write("Offspring_Fraction_Of_MQ0_Reads_For_Ref" + "," + "Offspring_Fraction_Of_MQ0_Reads_For_Alt" + ",");
		output.write("Offspring_Fraction_Of_Soft_Clipped_Reads_For_Ref" + ","
				+ "Offspring_Fraction_Of_Soft_Clipped_Reads_For_Alt" + ",");
		output.write(
				"Offspring_Mean_Nearby_Mismatches_For_Ref" + "," + "Offspring_Mean_Nearby_Mismatches_For_Alt" + ",");
		output.write("Offspring_Mean_Nearby_Indels_For_Ref" + "," + "Offspring_Mean_Nearby_Indels_For_Alt" + ",");
		output.write("Offspring_Strand_Direction_For_Ref" + "," + "Offspring_Strand_Direction_For_Alt" + ",");
		output.write("Offspring_Strand_Bias" + ",");
		output.write("PValue_Father_To_Offspring" + "," + "PValue_Mother_To_Offspring" + ",");
		output.write("Homopolymer_Flag" + "," + "Short_Tandem_Repeat_Flag" + "," + "Class" + "\n");
		logger.info("Extracting features ......");
        int numOfPositive=0;
        int numOfNegative=0;
		if (positivePath != null) {
			numOfPositive = buildTrainingSet(positivePath, "positive", output);
		}
		if (negativePath != null) {
			numOfNegative = buildTrainingSet(negativePath, "negative", output);
		}
		logger.info("Number of true DNMs is " + numOfPositive);
		logger.info("Number of false DNMs is " + numOfNegative);
		logger.info("Total number of DNMs is " + (numOfPositive + numOfNegative));
		output.close();
	}
	
	private int buildTrainingSet(String path, String label, PrintWriter output) throws IOException {
		DNMReader dnmReader = new DNMReader(path);
		HashMap<String, Trio> map = new HashMap();
		for (int i = 0; i < trios.size(); i++) {
			map.put(trios.get(i).getFamilyID(), trios.get(i));
		}
		DNMRecord dnmRecord = dnmReader.getNextRecord();
		int count = 0;
		String familyID = null;
		while (dnmRecord != null) {
			familyID = dnmRecord.getFamilyID();
			Trio trio = map.get(familyID);
			MultiSamLocusIterator trioSamLocusIterator = new MultiSamLocusIterator();
			SAMFileReader[] trioSAMFileReader = new SAMFileReader[3];
			IntervalList[] trioIntervalList = new IntervalList[3];
			SamLocusIterator[] samLocusIterator = new SamLocusIterator[3];
			trioSAMFileReader[0] = new SAMFileReader(new File(bams.getProperty(trio.getFather().getIndividualID())));
			trioSAMFileReader[1] = new SAMFileReader(new File(bams.getProperty(trio.getMother().getIndividualID())));
			trioSAMFileReader[2] = new SAMFileReader(new File(bams.getProperty(trio.getOffspring().getIndividualID())));
			for (int i = 0; i < 3; i++) {
				trioSAMFileReader[i].setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
				trioIntervalList[i] = new IntervalList(trioSAMFileReader[i].getFileHeader());
			}
			Map<String, DNMRecord> dnmMap = new HashMap();
			while (dnmRecord != null && dnmRecord.getFamilyID().equals(familyID)) {
				Interval interval = new Interval(dnmRecord.getChrom(), dnmRecord.getPos(), dnmRecord.getPos());
				for (int i = 0; i < 3; i++) {
					trioIntervalList[i].add(interval);
				}
				dnmMap.put(dnmRecord.getChrom() + "@" + dnmRecord.getPos(), dnmRecord);
				dnmRecord = dnmReader.getNextRecord();
			}
			for (int i = 0; i < 3; i++) {
				samLocusIterator[i] = new SamLocusIterator(trioSAMFileReader[i], trioIntervalList[i], true);
			}
			trioSamLocusIterator.add(samLocusIterator[0], FATHER_INDEX);
			trioSamLocusIterator.add(samLocusIterator[1], MOTHER_INDEX);
			trioSamLocusIterator.add(samLocusIterator[2], OFFSPRING_INDEX);
			ReferenceSequence referenceSequence = null;
			Map<Integer, SamLocusIterator.LocusInfo> tmp = null;
			ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(
					referenceSequenceFile);
			while ((tmp = trioSamLocusIterator.getLocusInfos()) != null) {
				LocusInfo tmplocusInfo = null;
				for (Integer tmpindex : tmp.keySet()) {
					tmplocusInfo = tmp.get(tmpindex);
					break;
				}
				String fordnmmap = tmplocusInfo.getSequenceName() + "@" + tmplocusInfo.getPosition();
				MultiPileup trioPileup = new MultiPileup(tmp, dnmMap.get(fordnmmap).getRef(),
						dnmMap.get(fordnmmap).getVar());
				if (referenceSequence == null || referenceSequence.getContigIndex() != trioPileup.getReferenceIndex()) {
					referenceSequence = referenceSequenceFileWalker.get(trioPileup.getReferenceIndex());
				}
				FeatureSelection featureSelection = new FeatureSelection(referenceSequence, trioPileup, indexfasta);
				output.write(featureSelection.extract() + "," + label + "\n");
				count++;
			}
			for (int i = 0; i < 3; i++) {
				trioSAMFileReader[i].close();
			}
			trioSamLocusIterator.close();
		}
		dnmReader.close();
		return count;
	}
}
