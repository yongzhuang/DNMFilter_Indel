package cn.edu.hit.dnmfilterindel.core;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.Properties;

import cn.edu.hit.dnmfilterindel.file.DNMRecord;
import cn.edu.hit.dnmfilterindel.sam.MultiPileup;
import cn.edu.hit.dnmfilterindel.sam.MultiSamLocusIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.samtools.SAMFileReader;

public class BuildTrainingSet {

	private String referenceSequenceFile;
	private Properties bamFiles;

	public BuildTrainingSet(String referenceSequenceFile, Properties bamFiles) {
		this.referenceSequenceFile = referenceSequenceFile;
		this.bamFiles = bamFiles;
	}

	public void build(String label, String fatherID, String motherID, String offspringID, List<DNMRecord> dnmRecordList,
			String outputPath) throws IOException {
		PrintWriter writer = new PrintWriter(outputPath);
		MultiSamLocusIterator trioSamLocusIterator = new MultiSamLocusIterator();
		SAMFileReader[] trioSAMFileReader = new SAMFileReader[3];
		IntervalList[] trioIntervalList = new IntervalList[3];
		SamLocusIterator[] samLocusIterator = new SamLocusIterator[3];
		trioSAMFileReader[0] = new SAMFileReader(new File(bamFiles.getProperty(fatherID)));
		trioSAMFileReader[1] = new SAMFileReader(new File(bamFiles.getProperty(motherID)));
		trioSAMFileReader[2] = new SAMFileReader(new File(bamFiles.getProperty(offspringID)));
		for (int i = 0; i < 3; i++) {
			trioSAMFileReader[i].setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
			trioIntervalList[i] = new IntervalList(trioSAMFileReader[i].getFileHeader());
		}
		for (DNMRecord dnmRecord : dnmRecordList) {
			Interval interval = new Interval(dnmRecord.getChrom(), dnmRecord.getPos(), dnmRecord.getPos());
			for (int i = 0; i < 3; i++) {
				trioIntervalList[i].add(interval);
			}
		}
		for (int i = 0; i < 3; i++) {
			samLocusIterator[i] = new SamLocusIterator(trioSAMFileReader[i], trioIntervalList[i], true);
		}
		trioSamLocusIterator.add(samLocusIterator[0], 0);
		trioSamLocusIterator.add(samLocusIterator[1], 1);
		trioSamLocusIterator.add(samLocusIterator[2], 2);

		ReferenceSequence referenceSequence = null;
		Map<Integer, SamLocusIterator.LocusInfo> tmp = null;
		ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(
				new File(referenceSequenceFile));
		int i = 0;
		while ((tmp = trioSamLocusIterator.getLocusInfos()) != null) {
			DNMRecord dnmRecord = dnmRecordList.get(i);
			MultiPileup trioPileup = new MultiPileup(tmp, dnmRecord.getRef(), dnmRecord.getVar());
			i++;
			if (referenceSequence == null || referenceSequence.getContigIndex() != trioPileup.getReferenceIndex()) {
				referenceSequence = referenceSequenceFileWalker.get(trioPileup.getReferenceIndex());
			}
			IndexedFastaSequenceFile indexfasta = new IndexedFastaSequenceFile(new File(referenceSequenceFile));
			FeatureSelection featureSelection = new FeatureSelection(referenceSequence, trioPileup, indexfasta);
			writer.write(featureSelection.extract() + "," + label + "\n");
		}
		for (int j = 0; j < 3; j++) {
			trioSAMFileReader[j].close();
		}
		trioSamLocusIterator.close();
		writer.close();
	}
}
