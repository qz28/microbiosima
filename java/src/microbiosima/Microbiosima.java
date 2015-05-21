package microbiosima;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

/**
 *
 * @author John
 */
public class Microbiosima {

	/**
	 * @param args
	 *            the command line arguments
	 * @throws java.io.FileNotFoundException
	 * @throws java.io.UnsupportedEncodingException
	 */
	public static void main(String[] args) throws FileNotFoundException,
			UnsupportedEncodingException, IOException {
		int Nh = 500;// Integer.parseInt(parameters[1]);
		int Nm = 1000;// Integer.parseInt(parameters[2]);
		int No = 150;// Integer.parseInt(parameters[3]);
		int Ng = 10000;
		int input = Integer.parseInt(args[0]);
		int x = input / 11;
		int y = input % 11;
		double pooled_or_fixed = 1 - Math.exp(-x);
		double environmental_factor = Math.exp(-y);
		double[] environment = new double[No];
		for (int i = 0; i < No; i++) {
			environment[i] = 1 / (double) No;
		}
		// LogNormalDistribution lgd=new LogNormalDistribution(0,1);
		// environment=lgd.sample(150);
		// double environment_sum=0;
		// for (int i=0;i<No;i++){
		// environment_sum+=environment[i];
		// }
		// for (int i=0;i<No;i++){
		// environment[i]/=environment_sum;
		// }
		for (int rep = 0; rep < 10; rep++) {
			String prefix = args[0] + "_" + Integer.toString(rep);
			PrintWriter file1 = new PrintWriter(new BufferedWriter(
					new FileWriter(prefix + "_gamma_diversity.txt")));
			PrintWriter file2 = new PrintWriter(new BufferedWriter(
					new FileWriter(prefix + "_alpha_diversity.txt")));
			PrintWriter file3 = new PrintWriter(new BufferedWriter(
					new FileWriter(prefix + "_beta_diversity.txt")));
			PrintWriter file4 = new PrintWriter(new BufferedWriter(
					new FileWriter(prefix + "_sum.txt")));
			PrintWriter file5 = new PrintWriter(new BufferedWriter(
					new FileWriter(prefix + "_inter_generation_distance.txt")));
			PrintWriter file6 = new PrintWriter(new BufferedWriter(
					new FileWriter(prefix
							+ "_environment_population_distance.txt")));

			Population population = new Population(Nm, environment, Nh,
					environmental_factor, pooled_or_fixed, 0, 0);

			while (population.getNumberOfGeneration() < Ng) {
				population.sumSpecies();
				if (population.getNumberOfGeneration() % 100 == 0) {
					// file1.print(population.gammaDiversity(false));
					// file2.print(population.alphaDiversity(false));
					// file1.print("\t");
					// file2.print("\t");
					file1.println(population.gammaDiversity(true));
					file2.println(population.alphaDiversity(true));
					population.microbiomeSequenceAlignment();
					file3.print(population.betaDiversity(true));
					file3.print("\t");
					file3.println(population.BrayCurtis(true));
					file4.println(population.printOut());
					file5.println(population.interGenerationDistance());
					file6.println(population.environmentPopulationDistance());
				}
				population.getNextGen();
			}
			file1.close();
			file2.close();
			file3.close();
			file4.close();
			file5.close();
			file6.close();
		}
	}
	// TODO code application logic here
}
