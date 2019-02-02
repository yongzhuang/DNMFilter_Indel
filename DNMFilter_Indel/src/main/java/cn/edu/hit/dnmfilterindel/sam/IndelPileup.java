/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cn.edu.hit.dnmfilterindel.sam;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import cn.edu.hit.dnmfilterindel.core.FisherExact;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.SamLocusIterator.LocusInfo;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.SequenceUtil;

/**
 *
 * @author Yongzhuang Liu
 */
public class IndelPileup {

    private ReferenceSequence referenceSequence;
    private LocusInfo locusInfo;
    private String refAllele;
    private String altAllele;


    private int depth;
    
    private int refAlleleCount;
    private int altAlleleCount;

    public int getForwardRefCount() {
        return forwardRef;
    }

    public int getReverseRefCount() {
        return reverseRef;
    }

    public int getForwardAltCount() {
        return forwardAlt;
    }

    public int getReverseAltCount() {
        return reverseAlt;
    }


    private int forwardRef;
    private int reverseRef;
    private int forwardAlt;
    private int reverseAlt;
    
    private int meanRefMappingQuality;
    private int meanAltMappingQuality;
    
    private int meanRefDistanceToThreePrime;
    private int meanAltDistanceToThreePrime;
    
    private double refFractionOfMQ0Reads;
    private double altFractionOfMQ0Reads;
    
    private double refFractionOfSoftClippedReads;
    private double altFractionOfSoftClippedReads;
    

    
    double meanRefNearbyIndels;
    double meanAltNearbyIndels;
    double meanRefNearbyMismatches;
    double meanAltNearbyMismatches;
   
    public IndelPileup(ReferenceSequence referenceSequence, LocusInfo locusInfo, String refAllele, String altAllele) {
        this.locusInfo = locusInfo;
        this.refAllele = refAllele;
        this.altAllele = altAllele;
        //if (locusInfo != null)
        	//System.out.println(locusInfo.getSequenceName() + "\t" + locusInfo.getSequenceIndex() + "\t" + refAllele + "\t" + altAllele);///////////////
        this.referenceSequence = referenceSequence;
        initialize();
    }

    private void initialize() {
        if (locusInfo == null) {
            return;
        }
        //System.out.println(locusInfo.getSequenceName()+"\t"+locusInfo.getPosition());
        
        /*altAlleleCount = 0;
        refAlleleCount = 0;
        forwardRef = 0;
        forwardAlt = 0;
        reverseAlt = 0;
        reverseRef = 0;
        meanRefMappingQuality = 0;
        meanAltMappingQuality = 0;
        meanRefDistanceToThreePrime = 0;
        meanAltDistanceToThreePrime = 0;
        altFractionOfMQ0Reads = 0;
        refFractionOfMQ0Reads = 0;
        refFractionOfSoftClippedReads = 0;
        altFractionOfSoftClippedReads = 0;
        meanRefNearbyIndels = 0;
        meanAltNearbyIndels = 0;
        meanRefNearbyMismatches = 0;
        meanAltNearbyMismatches = 0;*/
        depth = locusInfo.getRecordAndPositions().size();
        List<RecordAndOffset> list = locusInfo.getRecordAndPositions();
        
        for (int i = 0; i < list.size(); i++) {
            RecordAndOffset recordAndOffset = list.get(i);
            SAMRecord samRecord = recordAndOffset.getRecord();
            int offset = recordAndOffset.getOffset();
            String cigarString = locusInfo.getRecordAndPositions().get(i).getRecord().getCigarString();
            boolean reverse = locusInfo.getRecordAndPositions().get(i).getRecord().getReadNegativeStrandFlag();
            int mappingQuality = samRecord.getMappingQuality();
            int distanceToThreePrime = 0;
            if (!reverse) {
                distanceToThreePrime = samRecord.getReadLength() - offset - 1;
            } else {
                distanceToThreePrime = offset;
            }
            int nearbyIndels = getNearbyIndels(samRecord);
            int nearbyMismatches = getNearbyMismatches(samRecord);

            if (isAltDel(recordAndOffset) || isAltIns(recordAndOffset)) {
                altAlleleCount++;
                if (mappingQuality == 0) {
                    altFractionOfMQ0Reads++;
                }

                if (cigarString.contains("S")) {
                    altFractionOfSoftClippedReads++;
                }
                if (reverse) {
                    reverseAlt++;
                } else {
                    forwardAlt++;
                }
                //System.out.println("AltMappingQuality:"+mappingQuality);
                meanAltMappingQuality += mappingQuality;
                meanAltDistanceToThreePrime += distanceToThreePrime;
                meanAltNearbyIndels += nearbyIndels;
               // System.out.println("altnearbyIndels"+nearbyIndels);
                meanAltNearbyMismatches += nearbyMismatches;

            }
            if (isRefAllele(recordAndOffset)) {
                refAlleleCount++;
                if (mappingQuality == 0) {
                    refFractionOfMQ0Reads++;
                }
                if (cigarString.contains("S")) {
                    refFractionOfSoftClippedReads++;
                }
                if (reverse) {
                    reverseRef++;
                } else {
                    forwardRef++;
                }
                //System.out.println("refmapquality:"+mappingQuality);
                meanRefMappingQuality += mappingQuality;
                meanRefDistanceToThreePrime += distanceToThreePrime;
                meanRefNearbyIndels += nearbyIndels;
                //System.out.println("refnearbyIndels"+nearbyIndels);
                meanRefNearbyMismatches += nearbyMismatches;
            }
        }
        //System.out.println("altAlleleCount:"+altAlleleCount);
        meanAltMappingQuality = (int) divide(meanAltMappingQuality, altAlleleCount);
        //System.out.println("meanAltMappingQuality"+meanAltMappingQuality);
        meanAltDistanceToThreePrime = (int) divide(meanAltDistanceToThreePrime, altAlleleCount);
        meanAltNearbyIndels = divide((int) meanAltNearbyIndels, altAlleleCount);
        meanAltNearbyMismatches = divide((int) meanAltNearbyMismatches, altAlleleCount);
        altFractionOfSoftClippedReads = divide((int) altFractionOfSoftClippedReads, altAlleleCount);
        altFractionOfMQ0Reads = divide((int) altFractionOfMQ0Reads, altAlleleCount);
        
        //System.out.println("refAlleleCount:"+refAlleleCount);
        //System.out.println("meanRefMappingQuality:"+meanRefMappingQuality);
        meanRefMappingQuality = (int) divide(meanRefMappingQuality, refAlleleCount);       
        meanRefDistanceToThreePrime = (int) divide(meanRefDistanceToThreePrime, refAlleleCount);
        meanRefNearbyIndels = divide((int) meanRefNearbyIndels, refAlleleCount);
        meanRefNearbyMismatches = divide((int) meanRefNearbyMismatches, refAlleleCount);
        refFractionOfSoftClippedReads = divide((int) refFractionOfSoftClippedReads, refAlleleCount);
        refFractionOfMQ0Reads = divide((int) refFractionOfMQ0Reads, refAlleleCount);
        //System.out.println(refFractionOfSoftClippedReads + "\t" + refAlleleCount);
        //System.out.println("reverseRef"+reverseRef+"forwardRef"+forwardRef);
    }
    
    public int getDepth() {
        return depth;
    }
    
    public int getMeanRefMappingQuality() {
        return meanRefMappingQuality;
    }

    public int getMeanAltMappingQuality() {
        return meanAltMappingQuality;
    }

    public int getMeanRefDistanceToThreePrime() {
        return meanRefDistanceToThreePrime;
    }

    public int getMeanAltDistanceToThreePrime() {
        return meanAltDistanceToThreePrime;
    }

    public double getRefFractionOfMQ0Reads() {
        return refFractionOfMQ0Reads;
    }

    public double getAltFractionOfMQ0Reads() {
        return altFractionOfMQ0Reads;
    }

    public double getRefFractionOfSoftClippedReads() {
        return refFractionOfSoftClippedReads;
    }

    public double getAltFractionOfSoftClippedReads() {
        return altFractionOfSoftClippedReads;
    }

    public double getMeanRefNearbyIndels() {
        return meanRefNearbyIndels;
    }

    public double getMeanAltNearbyIndels() {
        return meanAltNearbyIndels;
    }

    public double getMeanRefNearbyMismatches() {
        return meanRefNearbyMismatches;
    }

    public double getMeanAltNearbyMismatches() {
        return meanAltNearbyMismatches;
    }
    
    public boolean getRefStrandDirection() {
        if (forwardRef == 0 || reverseRef == 0) {
            return true;
        }
        else
        return false;
    }
    
    
    public boolean getAltStrandDirection() {
        if (forwardAlt == 0 || reverseAlt == 0) {
            return true;
        } else {
            return false;
        }
    }
    
    public int getStrandBias() {
        FisherExact fe = new FisherExact(forwardRef + reverseRef + forwardAlt + reverseAlt);
        double sb = fe.getCumlativeP(forwardRef, reverseRef, forwardAlt, reverseAlt);
        if (sb > 1) {
            sb = 1.0;
        }
        return (int) (Math.log10(sb) * (-10));
    }
    
    private int getNearbyIndels(SAMRecord samRecord) {
        int indelCount = 0;
        Cigar cigar = samRecord.getCigar();
        for (int j = 0; j < cigar.numCigarElements(); j++) {
            CigarOperator operator = cigar.getCigarElement(j).getOperator();
            switch (operator) {
                case M:
                    break;
                case I:
                    indelCount++;
                    break;
                case D:
                    indelCount++;
                    break;
                default:
                    break;
            }
        }
        return indelCount;
    }

    private int getNearbyMismatches(SAMRecord samRecord) {
        byte[] referenceBases = referenceSequence.getBases();
        //System.out.println(new String(referenceBases));
        return SequenceUtil.countMismatches(samRecord, referenceBases);
    }
    
    public boolean isDeletion() {
        if (refAllele.length() > altAllele.length()) {
            return true;
        } else {
            return false;
        }
    }

    public boolean isInsertion() {
        if (refAllele.length() < altAllele.length()) {
            return true;
        } else {
            return false;
        }
    }
    
    public int getRefAlleleCount(){
        return refAlleleCount;  
    }
    
    public int getAltAlleleCount(){
        return altAlleleCount;  
    }

    private boolean isAltIns(RecordAndOffset recordAndOffset) {
        if(!isInsertion()){
            return false;
        }
        SAMRecord samRecord = recordAndOffset.getRecord();
        int offset = recordAndOffset.getOffset();
        Cigar cigar = samRecord.getCigar();
        int length = 0;
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            CigarElement cigarElement = cigar.getCigarElement(i);
            if (length == offset + 1) {
                if (cigarElement.getOperator() == CigarOperator.I && cigarElement.getLength() == (altAllele.length() - refAllele.length())) {
                    return true;
                }
            } else {
                if (cigarElement.getOperator() != CigarOperator.D) {
                    length += cigarElement.getLength();
                }
            }
        }
        return false;
    }

    private boolean isAltDel(RecordAndOffset recordAndOffset) {
        if(!isDeletion()){
            return false;
        }
        SAMRecord samRecord = recordAndOffset.getRecord();
        int offset = recordAndOffset.getOffset();
        Cigar cigar = samRecord.getCigar();
        int length = 0;
        for (int i = 0; i < cigar.numCigarElements(); i++) {
            CigarElement cigarElement = cigar.getCigarElement(i);
            if (length == offset + 1) {
                if (cigarElement.getOperator() == CigarOperator.D) {
                    if (cigarElement.getLength() == (refAllele.length() - altAllele.length())) {
                        return true;
                    }
                }
            } else {
                if (cigarElement.getOperator() != CigarOperator.D) {
                    length += cigarElement.getLength();
                }
            }
        }
        return false;
    }

    private boolean isRefAllele(RecordAndOffset recordAndOffset) {
        SAMRecord samRecord = recordAndOffset.getRecord();
        int offset = recordAndOffset.getOffset();
        Cigar cigar = samRecord.getCigar();
        if (isDeletion()) {
            int delStart = locusInfo.getPosition();
            int delEnd = delStart + refAllele.length() - 1;
            int alignmentStart = samRecord.getAlignmentStart();
            int cigarElementStart = alignmentStart - 1;
            int cigarElementEnd = alignmentStart - 1;
            for (int i = 0; i < cigar.numCigarElements(); i++) {
                CigarElement cigarElement = cigar.getCigarElement(i);
                if (cigarElement.getOperator().consumesReferenceBases()) {
                    cigarElementStart = cigarElementEnd + 1;
                    cigarElementEnd = cigarElementStart + cigarElement.getLength() - 1;
                }
                if (cigarElement.getOperator() == CigarOperator.M) {
                    if (cigarElementStart <= delStart && cigarElementEnd >= delEnd) {
                        String readString = samRecord.getReadString().substring(offset, offset + refAllele.length());
                        if (readString.equals(refAllele)) {
                            return true;
                        } else {
                            return false;
                        }
                    }
                }
            }
            return false;
        }
        if (isInsertion()) {
            int insStart = locusInfo.getPosition();
            int alignmentStart = samRecord.getAlignmentStart();
            int cigarElementStart = alignmentStart - 1;
            int cigarElementEnd = alignmentStart - 1;
            for (int i = 0; i < cigar.numCigarElements(); i++) {
                CigarElement cigarElement = cigar.getCigarElement(i);

                if (cigarElement.getOperator().consumesReferenceBases()) {
                    cigarElementStart = cigarElementEnd + 1;
                    cigarElementEnd = cigarElementStart + cigarElement.getLength() - 1;
                }
                if (cigarElement.getOperator() == CigarOperator.M) {
                    if (cigarElementStart <= insStart && cigarElementEnd > insStart) {
//                        System.out.println("Sample: "+samRecord.getReadGroup().getSample());
//                        System.out.println("Alignment Start: "+samRecord.getAlignmentStart());
//                        System.out.println("POS: "+insStart);
//                        System.out.println("READ: "+ samRecord.getReadString());
//                        System.out.println("REF: "+ refAllele);
//                        System.out.println("ALT: "+ altAllele);
//                        System.out.println("Offset: "+ offset);   
//                        System.out.println("CIGAR: "+ samRecord.getCigarString()); 
                        
                        if ((offset + refAllele.length()) <= samRecord.getReadString().length()) {
                            String readString = samRecord.getReadString().substring(offset, offset + refAllele.length());
                            if (readString.equals(refAllele)) {
                                return true;
                            } else {
                                return false;
                            }
                        }
                    }
                }
            }
        }
        return false;
    }

    public LocusInfo getLocusInfo() {
        return locusInfo;
    }

    private String getReferenceString(int beginIndex, int endIndex) {
        byte[] array = new byte[endIndex - beginIndex];
        byte[] referenceBases = referenceSequence.getBases();
        int j = 0;
        for (int i = beginIndex; i < endIndex; i++) {
            array[j] = referenceBases[i];
            j++;
        }
        return new String(array);
    }
    
    private double divide(int a, int b) {
        if (b == 0) {
            return 0.0;
        } else {
            return (double)a / (double)b;
        }
    }

}
