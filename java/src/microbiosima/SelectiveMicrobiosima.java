/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package microbiosima;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import utils.random.MathUtil;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Option.Builder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


public class SelectiveMicrobiosima extends Microbiosima {
	
	private static final String VERSION = "2.0";

	/**
	 * @param args
	 *            the command line arguments
	 * @throws java.io.FileNotFoundException
	 * @throws java.io.UnsupportedEncodingException
	 */
	 
     public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
        int populationSize=5000;//Integer.parseInt(parameters[1]);
        int microSize=1000000000;//Integer.parseInt(parameters[2]);
        int numberOfSpecies=150;//Integer.parseInt(parameters[3]);
        int numberOfGeneration=200000;
        int Ngene=25;
		int numberOfObservation=5000;
		int numberOfReplication=5;
        double Ngenepm=5;
		double pctEnv = 0;
		double pctPool = 0;
		double msCoeff=1;
		double hsCoeff=1;
        boolean HMS_or_TMS=true;
                
		

		Options options = new Options();

		Option help = new Option("h", "help", false, "print this message");
		Option version = new Option("v", "version", false,
				"print the version information and exit");
		options.addOption(help);
		options.addOption(version);
		
		options.addOption(Option.builder("o").longOpt("obs").hasArg()
				.argName("OBS").desc("Number generation for observation [default: 100]")
				.build());
		options.addOption(Option.builder("r").longOpt("rep").hasArg()
				.argName("REP").desc("Number of replication [default: 1]")
				.build());
					
		Builder C = Option.builder("c").longOpt("config")
				.numberOfArgs(6).argName("Pop Micro Spec Gen")
				.desc("Four Parameters in the following orders: "
						+ "(1) population size, (2) microbe size, (3) number of species, (4) number of generation, (5) number of total traits, (6)number of traits per microbe"
						+ " [default: 500 1000 150 10000 10 5]");
		options.addOption(C.build());

		HelpFormatter formatter = new HelpFormatter();
		String syntax = "microbiosima pctEnv pctPool";
		String header = "\nSimulates the evolutionary and ecological dynamics of microbiomes within a population of hosts.\n\n"+
		"required arguments:\n"+"  pctEnv             Percentage of environmental acquisition\n"+
		"  pctPool            Percentage of pooled environmental component\n"+"  msCoeff            Parameter related to microbe selection strength\n"+
		"  hsCoeff            Parameter related to host selection strength\n"+"  HMS_or_TMS         String HMS or TMS to specify host-mediated or trait-mediated microbe selection\n"
		+ "\noptional arguments:\n";
		String footer = "\n";
		
		formatter.setWidth(80);

		
		CommandLineParser parser = new DefaultParser();
		CommandLine cmd = null;

		try {
			cmd = parser.parse(options, args);
			String[] pct_config = cmd.getArgs();
			
			if (cmd.hasOption("h") || args.length == 0) {
				formatter.printHelp(syntax, header, options, footer, true);
				System.exit(0);
			}
			if(cmd.hasOption("v")){
				System.out.println("Microbiosima "+VERSION);
				System.exit(0);
			}
			if (pct_config.length != 5){
				System.out.println("ERROR! Required exactly five argumennts for pct_env, pct_pool, msCoeff, hsCoeff and HMS_or_TMS. It got "
								+ pct_config.length + ": " + Arrays.toString(pct_config));
				formatter.printHelp(syntax, header, options, footer, true);
				System.exit(3);
			}
			else{
				pctEnv = Double.parseDouble(pct_config[0]);
				pctPool = Double.parseDouble(pct_config[1]);
				msCoeff = Double.parseDouble(pct_config[2]);
				hsCoeff = Double.parseDouble(pct_config[3]);
				if (pct_config[4].equals("HMS")) HMS_or_TMS=true;
				if (pct_config[4].equals("TMS")) HMS_or_TMS=false;
				if(pctEnv<0 || pctEnv >1){
					System.out.println(
						"ERROR: pctEnv (Percentage of environmental acquisition) must be between 0 and 1 (pctEnv="
						+ pctEnv + ")! EXIT");
					System.exit(3);
				}
				if(pctPool<0 || pctPool >1){
					System.out.println(
						"ERROR: pctPool (Percentage of pooled environmental component must) must be between 0 and 1 (pctPool="
						+ pctPool + ")! EXIT");
					System.exit(3);
				}
				if(msCoeff<1){
					System.out.println(
						"ERROR: msCoeff (parameter related to microbe selection strength) must be not less than 1 (msCoeff="
						+ msCoeff + ")! EXIT");
					System.exit(3);
				}
				if(hsCoeff<1){
					System.out.println(
						"ERROR: hsCoeff (parameter related to host selection strength) must be not less than 1 (hsCoeff="
						+ hsCoeff + ")! EXIT");
					System.exit(3);
				}
				if (!(pct_config[4].equals("HMS")||pct_config[4].equals("TMS"))){
					System.out.println(
						"ERROR: HMS_or_TMS (parameter specifying host-mediated or trait-mediated selection) must be either 'HMS' or 'TMS' (HMS_or_TMS="
						+ pct_config[4] + ")! EXIT");
					System.exit(3);
				}
				
			}
			if (cmd.hasOption("config")){
				String[] configs = cmd.getOptionValues("config");
				populationSize = Integer.parseInt(configs[0]);
				microSize = Integer.parseInt(configs[1]);
				numberOfSpecies = Integer.parseInt(configs[2]);
				numberOfGeneration = Integer.parseInt(configs[3]);
				Ngene=Integer.parseInt(configs[4]);
				Ngenepm=Double.parseDouble(configs[5]);
				if (Ngenepm>Ngene){
					System.out.println(
						"ERROR: number of traits per microbe must not be greater than number of total traits! EXIT");
					System.exit(3);
				}
			}
			if (cmd.hasOption("obs")){
				numberOfObservation= Integer.parseInt(cmd.getOptionValue("obs"));
			}
			if (cmd.hasOption("rep")){
				numberOfReplication= Integer.parseInt(cmd.getOptionValue("rep"));
			}			
			
			
			
		} catch (ParseException e) {
			e.printStackTrace();
			System.exit(3);
		}
		
		StringBuilder sb = new StringBuilder();
		sb.append("Configuration Summary:")
			.append("\n\tPopulation size: ").append(populationSize)
			.append("\n\tMicrobe size: ").append(microSize)
			.append("\n\tNumber of species: ").append(numberOfSpecies)
			.append("\n\tNumber of generation: ").append(numberOfGeneration)
			.append("\n\tNumber generation for observation: ").append(numberOfObservation)
			.append("\n\tNumber of replication: ").append(numberOfReplication)
			.append("\n\tNumber of total traits: ").append(Ngene)
			.append("\n\tNumber of traits per microbe: ").append(Ngenepm)
			.append("\n");
		System.out.println(sb.toString());
		
		
        double[] environment=new double[numberOfSpecies];
        for (int i=0;i<numberOfSpecies;i++){
            environment[i]=1/(double)numberOfSpecies;
        }
        int[] fitnessToHost=new int[Ngene];
        int[] fitnessToMicrobe=new int[Ngene];
  
        for (int rep=0;rep<numberOfReplication;rep++){
			String prefix = ""+(rep+1)+"_";
                        String sufix;
			if (HMS_or_TMS)
			sufix = "_E" + pctEnv + "_P" + pctPool +"_HS"+hsCoeff+"_HMS"+msCoeff+ ".txt"; 
		    else
			sufix = "_E" + pctEnv + "_P" + pctPool +"_HS"+hsCoeff+"_TMS"+msCoeff+ ".txt";
			System.out.println("Output 5 result files in the format of: "+prefix+"[****]" +sufix);
        try{			
        PrintWriter file1= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"gamma_diversity"+sufix)));
        PrintWriter file2= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"alpha_diversity"+sufix)));
        PrintWriter file3= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"beta_diversity"+sufix)));
        PrintWriter file4= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"sum"+sufix)));
        PrintWriter file5= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"inter_generation_distance"+sufix)));
        PrintWriter file6= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"environment_population_distance"+sufix)));
        PrintWriter file7= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"host_fitness"+sufix)));
        PrintWriter file8= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"cos_theta"+sufix)));
        PrintWriter file9= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"host_fitness_distribution"+sufix)));
        PrintWriter file10= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"microbiome_fitness_distribution"+sufix)));
        PrintWriter file11= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"bacteria_contents"+sufix)));
        PrintWriter file12= new PrintWriter(new BufferedWriter(new FileWriter(prefix+"individual_bacteria_contents"+sufix)));
	MathUtil.setSeed(rep%numberOfReplication);
	for (int i=0;i<Ngene;i++){
          fitnessToMicrobe[i]=MathUtil.getNextInt(2)-1;
          fitnessToHost[i]=MathUtil.getNextInt(2)-1;
        }
        SelectiveSpeciesRegistry ssr=new SelectiveSpeciesRegistry(numberOfSpecies,Ngene, Ngenepm, msCoeff,fitnessToHost,fitnessToMicrobe);
        MathUtil.setSeed();
        SelectivePopulation population=new SelectivePopulation(microSize, environment, populationSize, pctEnv , pctPool,0,0,ssr,hsCoeff,HMS_or_TMS);
        
        while (population.getNumberOfGeneration()<=numberOfGeneration){
            population.sumSpecies();
            if(population.getNumberOfGeneration()%numberOfObservation==0){
                //file1.print(population.gammaDiversity(false));
                //file2.print(population.alphaDiversity(false));
                //file1.print("\t");
                //file2.print("\t");
                file1.println(population.gammaDiversity(true));
                file2.println(population.alphaDiversity(true));                
                //file3.print(population.betaDiversity(true));
                //file3.print("\t");
                file3.println(population.BrayCurtis(true));
                file4.println(population.printOut());
                file5.println(population.interGenerationDistance());
                file6.println(population.environmentPopulationDistance());
                file7.print(population.averageHostFitness());
                file7.print("\t");
                file7.println(population.varianceHostFitness());
                file8.println(population.cosOfMH());
                file9.println(population.printOutHFitness());
                file10.println(population.printOutMFitness());
				file11.println(population.printBacteriaContents());
            }
            population.getNextGen();
        }
        for (SelectiveIndividual host:population.getIndividuals()){
                file12.println(host.printBacteriaContents());
		}
        file1.close();
        file2.close();
        file3.close();
        file4.close();
        file5.close();
        file6.close();
        file7.close();
        file8.close();
        file9.close();
        file10.close();
		file11.close();
		file12.close();
		}catch (IOException e) {
				e.printStackTrace();
			} 
        }
    }
    
}
