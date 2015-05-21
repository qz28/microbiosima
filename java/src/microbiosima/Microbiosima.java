/*******************************************************************************
 *
 * Copyright (C) 2014, 2015 Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *******************************************************************************/
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package microbiosima;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
/**
 *
 * @author John
 */
public class Microbiosima {

    /**
     * @param args the command line arguments
     * @throws java.io.FileNotFoundException
     * @throws java.io.UnsupportedEncodingException
     */
    public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException, IOException {
        int Nh=50000;//Integer.parseInt(parameters[1]);
        int Nm=100000000;//Integer.parseInt(parameters[2]);
        int No=718;//Integer.parseInt(parameters[3]);
        int Ng=(int) 1e1;
        int input=Integer.parseInt(args[0]);
        int x=input/21;
        int y=input%21;
        double pooled_or_fixed=1-Math.exp(-x);
        double environmental_factor=Math.exp(-y);
        double[] environment=new double[No];
        for (int i=0;i<No;i++){
            environment[i]=1/(double)No;
        }
//		Pre java.nio        
//        PrintWriter out
//        = new PrintWriter(new BufferedWriter(new FileWriter("foo.out")));
//		
//        try (BufferedWriter writer = Files.newBufferedWriter(file, charset)) {
//            writer.write(s, 0, s.length());
//        } catch (IOException x) {
//            System.err.format("IOException: %s%n", x);
//        }
//        PrintWriter file1=new PrintWriter(args[0]+"_"+args[1]+"_"+"gamma_diversity.txt","UTF-8");
//        PrintWriter file2=new PrintWriter(args[0]+"_"+args[1]+"_"+"alpha_diversity.txt","UTF-8");
//        PrintWriter file3=new PrintWriter(args[0]+"_"+args[1]+"_"+"beta_diversity.txt","UTF-8");
//        PrintWriter file4=new PrintWriter(args[0]+"_"+args[1]+"_"+"sum.txt","UTF-8");
//        PrintWriter file5=new PrintWriter(args[0]+"_"+args[1]+"_"+"inter_generation_distance.txt","UTF-8");
//        PrintWriter file6=new PrintWriter(args[0]+"_"+args[1]+"_"+"environment_population_distance.txt","UTF-8");
        Population population=new Population(Nm, environment, Nh, environmental_factor , pooled_or_fixed,args[1],10,5);
        
        for (int i = 0; i < 10; i++) {
			
		
        long timeStart = System.nanoTime();
        while (population.getNumberOfGeneration()<Ng){
            population.sumSpecies();
            //if(population.number_of_generation%1==0){
//                file1.print(population.gammaDiversity(false));
//                file2.print(population.alphaDiversity(false));
//                file1.print("\t");
//                file2.print("\t");
//                file1.println(population.gammaDiversity(true));
//                file2.println(population.alphaDiversity(true));                
//                population.microbiomeSequenceAlignment();
//                file3.println(population.betaDiversity(false));
//                file4.println(population.printOut());
//                file5.println(population.interGenerationDistance());
//                file6.println(population.environmentPopulationDistance());
            //}
            population.getNextGen();
        }
        population.resetGeneration();
        long timeEnd = System.nanoTime();
        System.out.println((timeEnd - timeStart)/1e9);
        
        }
//        System.out.println(Binomial.count1 +"\t"+ (Binomial.time1/Binomial.count1));
//        System.out.println(Binomial.count2 +"\t"+ (Binomial.time2/Binomial.count2));
        /*
        1e6 Generation: 26s (both sumSpecies and getNextGen), 2s without parentalInheritanceAndEnvironmentalAcquisition()
        				RandomSample.randomSample(host_index, numberOfSamples)); take about ~2s
        
28.763721101
27.47352359
27.603628558
28.219436034
27.913342197
28.552429274
28.314822672
28.457028852
28.76946679
28.451131119

vs 
colt
25.300478373
25.063111414
24.287014842
24.714187388
23.888524814
24.180974722
24.776706936
26.975339013
25.035381938
25.004082058

vs static random
22.639185693
22.930234335
23.05024542
22.564727753
22.617750041
22.055159096
22.068620118
22.450818429
22.098670509
22.465345232
        
        
test full dataset
Orig
49.489644297
49.380118422
50.872091345

mod
37.105326458
37.129677668
37.270540841

exclud 0
30.159822112
31.111805308
31.43730276


Sampling_v5 exclude 0
38.520022664
40.292304405
40.886912877
41.49681042


        *
        */
        
        
//        file1.close();
//        file2.close();
//        file3.close();
//        file4.close();
//        file5.close();
//        file6.close();
        }
    
        // TODO code application logic here
    };
