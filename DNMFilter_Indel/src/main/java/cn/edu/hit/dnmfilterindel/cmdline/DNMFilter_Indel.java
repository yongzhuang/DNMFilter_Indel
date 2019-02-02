/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package cn.edu.hit.dnmfilterindel.cmdline;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Properties;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;

import cn.edu.hit.dnmfilterindel.core.CreateSyntheticIndels;
import cn.edu.hit.dnmfilterindel.core.Extract;
import cn.edu.hit.dnmfilterindel.core.ExtractInheritedIndels;
import cn.edu.hit.dnmfilterindel.core.GBM;

/**
 *
 * @author Yongzhuang Liu
 */
public class DNMFilter_Indel {

	private static Logger logger = Logger.getLogger(DNMFilter_Indel.class);

	public static void main(String[] args) throws IOException {
		String usage = "\nDNMFilter_Indel-0.1.0\n";
		usage = usage + "\nUsage: java -jar DNMFilter_Indel.jar <COMMAND> [OPTIONS]\n\n";
		usage = usage + "COMMANDS:\n" + "\textract\t\tExtract sequence features to build the training set\n"
				+ "\tgbm\t\tUse gradient boosting approach to filter de novo indels\n"
				+ "\tsynthetic\tSimulate synthetic de novo indels\n" + "\tinherited\tExtract inherited indels\n";
		String cmd = null;
		if (args.length > 0) {
			if (args[0].equals("extract") || args[0].equals("gbm") || args[0].equals("synthetic")
					|| args[0].equals("inherited")) {
				cmd = args[0];
			} else {
				logger.error("Command is not recognized!\n" + usage);
				return;
			}
		} else {
			System.out.println(usage);
			return;
		}
		run(cmd, args);
	}

	private static void run(String cmd, String[] args) throws IOException {
		long start = System.currentTimeMillis();
		CommandLineParser parser = new PosixParser();
		CommandLine commandLine = null;
		Options options = createOptions(cmd);
		try {
			if (options != null) {
				commandLine = parser.parse(options, args);
			}
		} catch (ParseException parseException) {
			logger.error("Invalid command line parameters!");
		}
		if (cmd.equals("extract")) {
			if (isValidated(commandLine, "extract")) {
				(new Extract(getProperties(commandLine, "extract"))).run();
			} else {
				printHelp(options, "extract");
				return;
			}
		}
		if (cmd.equals("synthetic")) {
			if (isValidated(commandLine, "synthetic")) {
				(new CreateSyntheticIndels(getProperties(commandLine, "synthetic"))).simulateSyntheticIndels();
			} else {
				printHelp(options, "synthetic");
				return;
			}
		}
		if (cmd.equals("inherited")) {
			if (isValidated(commandLine, "inherited")) {
				(new ExtractInheritedIndels(getProperties(commandLine, "inherited"))).extractInheritedIndels();;
			} else {
				printHelp(options, "inherited");
				return;
			}
		}
		if (cmd.equals("gbm")) {
			if (isValidated(commandLine, "gbm")) {
				if (!(new GBM(getProperties(commandLine, "gbm"))).run()) {
					return;
				}
			} else {
				printHelp(options, "gbm");
				return;
			}
		}
		long end = System.currentTimeMillis();
		logger.info("Total running time is " + (end - start) / 1000 + " seconds");
		logger.info("Done!");
	}

	private static Options createOptions(String cmd) {
		Options options = new Options();
		if (cmd.equals("extract")) {
			options.addOption(OptionBuilder.withLongOpt("reference").withDescription("reference genome file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("bam").withDescription("bam list file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("positive")
					.withDescription("known true positive DNM file (required)").hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("negative")
					.withDescription("known false positive DNM file (required)").hasArg().withArgName("FILE").create());
			return options;
		} else if (cmd.equals("gbm")) {
			options.addOption(OptionBuilder.withLongOpt("reference").withDescription("reference genome file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("bam").withDescription("bam list file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("candidate").withDescription("candidate DNM file (required)")
					.hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("training").withDescription("training set (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("configuration")
					.withDescription("feature configuration file (required)").hasArg().withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("cutoff")
					.withDescription("cutoff to determine a putative DNM (optional, default 0.4)").hasArg()
					.withArgName("DOUBLE").create());
			return options;
		} else if (cmd.equals("synthetic")) {
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("vcf").withDescription("DNM file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("number").withDescription("number of DNMs (required)").hasArg()
					.withArgName("INT").create());
			options.addOption(OptionBuilder.withLongOpt("label").withDescription("exchange 'FATHER' or 'MOTHER' with the offspring (required)").hasArg()
					.withArgName("STR").create());
			return options;
		}else if (cmd.equals("inherited")) {
			options.addOption(OptionBuilder.withLongOpt("pedigree").withDescription("pedigree file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("output").withDescription("output file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("vcf").withDescription("DNM file (required)").hasArg()
					.withArgName("FILE").create());
			options.addOption(OptionBuilder.withLongOpt("number").withDescription("number of DNMs (required)").hasArg()
					.withArgName("INT").create());
			return options;
			}
		else {
			return null;
		}
	}

	private static Properties getProperties(CommandLine line, String cmd) {
		Properties properties = new Properties();
		if (cmd.equals("extract")) {
			properties.put("reference", line.getOptionValue("reference"));
			properties.put("pedigree", line.getOptionValue("pedigree"));
			properties.put("bam", line.getOptionValue("bam"));
			properties.put("output", line.getOptionValue("output"));
			if (line.getOptionValue("positive") != null) {
				System.out.println(line.getOptionValue("positive"));
				properties.put("positive", line.getOptionValue("positive"));
			}
			if (line.getOptionValue("negative") != null) {
				properties.put("negative", line.getOptionValue("negative"));
			}
		}
		if (cmd.equals("gbm")) {
			properties.put("reference", line.getOptionValue("reference"));
			properties.put("pedigree", line.getOptionValue("pedigree"));
			properties.put("bam", line.getOptionValue("bam"));
			properties.put("output", line.getOptionValue("output"));
			properties.put("candidate", line.getOptionValue("candidate"));
			properties.put("training", line.getOptionValue("training"));
			properties.put("configuration", line.getOptionValue("configuration"));
			if (!line.hasOption("cutoff")) {
				properties.put("cutoff", "0.5");
			} else {
				properties.put("cutoff", line.getOptionValue("cutoff"));
			}
		}
		if (cmd.equals("synthetic")) {
			properties.put("pedigree", line.getOptionValue("pedigree"));
			properties.put("output", line.getOptionValue("output"));
			properties.put("vcf", line.getOptionValue("vcf"));
			properties.put("number", line.getOptionValue("number"));
			properties.put("label", line.getOptionValue("label"));
		}
		if (cmd.equals("inherited")) {
			properties.put("pedigree", line.getOptionValue("pedigree"));
			properties.put("output", line.getOptionValue("output"));
			properties.put("vcf", line.getOptionValue("vcf"));
			properties.put("number", line.getOptionValue("number"));
		}
		return properties;
	}

	private static boolean isValidated(CommandLine line, String cmd) {
		boolean tag = true;

		if (!line.hasOption("pedigree") || !(new File(line.getOptionValue("pedigree")).isFile())) {
			logger.error("The pedigree file is not correctly specified!");
			tag = false;
		}

		if (!line.hasOption("output")) {
			logger.error("The output file is not correctly specified!");
			tag = false;
		}
		if (cmd.equals("extract")) {
			if (!line.hasOption("reference") || !(new File(line.getOptionValue("reference")).isFile())) {
				logger.error("The reference genome file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("bam") || !(new File(line.getOptionValue("bam")).isFile())) {
				logger.error("The bam list file is not correctly specified!");
				tag = false;
			}
			if (!(line.hasOption("positive") && (new File(line.getOptionValue("positive")).isFile())
					|| line.hasOption("negative") && (new File(line.getOptionValue("negative")).isFile()))) {
				logger.error("The true or false DNM file is not correctly specified!");
				tag = false;
			}
		}
		if (cmd.equals("gbm")) {
			if (!line.hasOption("reference") || !(new File(line.getOptionValue("reference")).isFile())) {
				logger.error("The reference genome file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("bam") || !(new File(line.getOptionValue("bam")).isFile())) {
				logger.error("The bam list file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("training") || !(new File(line.getOptionValue("training")).isFile())) {
				logger.error("The training set is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("candidate") || !(new File(line.getOptionValue("candidate")).isFile())) {
				logger.error("The candidate DNM file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("configuration") || !(new File(line.getOptionValue("configuration")).isFile())) {
				logger.error("The feature configuration file is not correctly specified!");
				tag = false;
			}
		}
		if (cmd.equals("synthetic")) {
			if (!line.hasOption("vcf") || !(new File(line.getOptionValue("vcf")).isFile())) {
				logger.error("The vcf file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("number")) {
				logger.error("The number of simulated DNMs is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("label")) {
				logger.error("The exchange label is not correctly specified!");
				tag = false;
			}
		}
		if (cmd.equals("inherited")) {
			if (!line.hasOption("vcf") || !(new File(line.getOptionValue("vcf")).isFile())) {
				logger.error("The vcf file is not correctly specified!");
				tag = false;
			}
			if (!line.hasOption("number")) {
				logger.error("The number of simulated DNMs is not correctly specified!");
				tag = false;
			}
		}
		return tag;
	}

	private static void printHelp(Options options, String command) {
		System.out.println();
		String cmdLineSyntax = "java -jar DNMFilter_Indel.jar " + command + " [OPTIONS]\n";
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptionComparator(new DNMFilter_Indel.OptionComarator());
		formatter.printHelp(cmdLineSyntax, options);
	}

	private static class OptionComarator<T extends Option> implements Comparator<T> {
		private static final String OPTS_ORDER = "rpnbktcol";

		public int compare(T o1, T o2) {
			return OPTS_ORDER.indexOf(o1.getLongOpt().charAt(0)) - OPTS_ORDER.indexOf(o2.getLongOpt().charAt(0));
		}
	}
}
