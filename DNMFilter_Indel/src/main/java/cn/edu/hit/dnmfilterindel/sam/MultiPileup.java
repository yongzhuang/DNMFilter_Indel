/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cn.edu.hit.dnmfilterindel.sam;

import java.util.Map;
import net.sf.picard.util.SamLocusIterator.LocusInfo;

/**
 *
 * @author Yongzhuang Liu
 */
public class MultiPileup {

    private Map<Integer, LocusInfo> map;


    private String refAllele;
    private String altAllele;

    public MultiPileup(Map<Integer, LocusInfo> map, String refAllele, String altAllele) {
        this.map = map;
        this.refAllele=refAllele;
        this.altAllele=altAllele;
    }

    public int getNumberOfPileups() {
        return map.size();
    }

    public int getReferenceIndex() {
        LocusInfo locusInfo = null;
        for (Integer index : map.keySet()) {
            locusInfo = map.get(index);
            break;
        }
        return locusInfo.getSequenceIndex();
    }

    public String getReferenceName() {
        LocusInfo locusInfo = null;
        for (Integer index : map.keySet()) {
            locusInfo = map.get(index);
            break;
        }
        return locusInfo.getSequenceName();
    }

    public int getPosition() {
        LocusInfo locusInfo = null;
        for (Integer index : map.keySet()) {
            locusInfo = map.get(index);
            break;
        }
        return locusInfo.getPosition();
    }

    public int getDepth(int index) {
        if (!map.containsKey(index)) {
            return 0;
        } else {
            return map.get(index).getRecordAndPositions().size();
        }
    }

    public LocusInfo getLocusInfo(int index) {
        if (!map.containsKey(index)) {
            return null;
        } else {
            return map.get(index);
        }
    }

    public Map<Integer, LocusInfo> getLocusInfoMap() {
        return map;
    }
    
    public String getRefAllele() {
        return refAllele;
    }

    public String getAltAllele() {
        return altAllele;
    }
}
