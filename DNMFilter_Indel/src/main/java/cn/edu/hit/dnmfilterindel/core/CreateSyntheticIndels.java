package cn.edu.hit.dnmfilterindel.core;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.Random;

import cn.edu.hit.dnmfilterindel.file.DNMRecord;
import cn.edu.hit.dnmfilterindel.file.PEDReader;
import cn.edu.hit.dnmfilterindel.file.Trio;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class CreateSyntheticIndels {
	
	private String vcfFile;
	private int number;
	private String outputFile;
	private String pedFile;
	private String tag;

	public CreateSyntheticIndels(Properties properties) {
		this.vcfFile = properties.getProperty("vcf");
		this.outputFile = properties.getProperty("output");
		this.number = Integer.parseInt(properties.getProperty("number"));
		this.pedFile = properties.getProperty("pedigree");
		this.tag = properties.getProperty("label");
	}
	
	public void simulateSyntheticIndels() throws IOException {
		List<DNMRecord> variantRecordList = new ArrayList();
		List<Trio> trios = (new PEDReader(pedFile)).getTrios();
		String fatherID = trios.get(0).getFather().getIndividualID();
		String motherID = trios.get(0).getMother().getIndividualID();
		String offspringID = trios.get(0).getOffspring().getIndividualID();
		String familyID = trios.get(0).getFamilyID();
		VCFFileReader vcf = new VCFFileReader(new File(vcfFile));
		CloseableIterator<VariantContext> iterator = vcf.iterator();
		while (iterator.hasNext()) {
			VariantContext variantContext = iterator.next();
			if (variantContext.getChr().equals("chrX") || variantContext.getChr().equals("X")) {
				break;
			}
			if (variantContext.isIndel()) {
				Genotype fatherGenotype = variantContext.getGenotype(fatherID);
				Genotype motherGenotype = variantContext.getGenotype(motherID);
				Genotype offspringGenotype = variantContext.getGenotype(offspringID);
				Genotype varGenotype = null;
				Genotype refGenotype = null;
				if (tag.equals("FATHER")) {
					varGenotype = fatherGenotype;
					refGenotype = motherGenotype;
				} else if (tag.equals("MOTHER")) {
					varGenotype = motherGenotype;
					refGenotype = fatherGenotype;
				} 
				if (varGenotype.isHet() && refGenotype.isHomRef() && offspringGenotype.isHomRef()) {
					String refAllele = variantContext.getReference().getBaseString();
					String genotypeString = varGenotype.getGenotypeString();
					String[] alleles;
					if (varGenotype.isPhased()) {
						alleles = genotypeString.split("|");
					} else {
						alleles = genotypeString.split("/");
					}
					String altAllele = null;
					if (alleles[0].equals(refAllele)) {
						altAllele = alleles[1];
					} else if (alleles[1].equals(refAllele)) {
						altAllele = alleles[0];
					}
					if (altAllele != null) {
						DNMRecord variantRecord = new DNMRecord(familyID, variantContext.getChr(),
								variantContext.getStart(), refAllele, altAllele, null);
						variantRecordList.add(variantRecord);
					}
				}
			}
		}
		vcf.close();
		List<DNMRecord> newDNMList = getRandomDNMList(variantRecordList, number);
		PrintWriter writer = new PrintWriter(outputFile);
		for(DNMRecord dnmRecord:newDNMList) {
			writer.write(dnmRecord.getFamilyID() + "," + dnmRecord.getChrom() + "," + dnmRecord.getPos() + ","
					+ dnmRecord.getRef() + "," + dnmRecord.getVar()+"\n");
		}
		writer.close();
	}
	
	private int[] randomCommon(int size, int n) {
		int[] result = new int[n];
		int count = 0;
		Random random = new Random();
		while (count < n) {
			int num = random.nextInt(size);
			boolean flag = true;
			for (int j = 0; j < n; j++) {
				if (num == result[j]) {
					flag = false;
					break;
				}
			}
			if (flag) {
				result[count] = num;
				count++;
			}
		}
		Arrays.sort(result);
		return result;
	}

	private List<DNMRecord> getRandomDNMList(List<DNMRecord> dnmList, int number) {
		int[] index = randomCommon(dnmList.size(), number);
		List<DNMRecord> newDNMList = new ArrayList();
		for (int i = 0; i < index.length; i++) {
			newDNMList.add(dnmList.get(index[i] - 1));
		}
		return newDNMList;
	}
}
