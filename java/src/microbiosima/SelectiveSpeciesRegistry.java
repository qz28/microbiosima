/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package microbiosima;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import utils.RandomSample;
import utils.random.Binomial;
import utils.random.MathUtil;


public class SelectiveSpeciesRegistry extends SpeciesRegistry {
    int numOfTolGenes;
    private double numOfGenesPerMicrobe;
    double msc;
	int[] hfr;
	int[] mfr;
    private ArrayList<Integer> geneCompositionRecords=new ArrayList<>();

    public SelectiveSpeciesRegistry(int initialNumberOfSpecies,int notg, double nogpm, double microbeSelectionCoef,int[] fitnessToHost, int[] fitnessToMicrobe) {
        super(initialNumberOfSpecies);
        numOfTolGenes=notg;
        numOfGenesPerMicrobe=nogpm;
        msc=microbeSelectionCoef;
		hfr=fitnessToHost;
		mfr=fitnessToMicrobe;
        List<Integer> geneIndex = new ArrayList<>();
        for (int i=0;i<numOfTolGenes;i++){
            geneIndex.add(i);
        }
        for (int i=0;i<initialNumberOfSpecies;i++){
            Set<Integer> tempGenes=RandomSample.randomSample(geneIndex, (int)numOfGenesPerMicrobe);
            int geneBinaryComp=0;
            for (int gene : tempGenes){
                geneBinaryComp+=Math.pow(2,gene);
            }
            geneCompositionRecords.add(geneBinaryComp);
        }
        
    }        
    public double getTotalFitness(double[] microbiome, double[] fitnessRecords ){
        double totalFitness=0;
        double abundance=0;
        for (int i=0;i<microbiome.length;i++){
            abundance+=microbiome[i];
            totalFitness+=microbiome[i]*fitnessRecords[i];
        }
        return totalFitness/abundance;
    }
        
    public void getFitness(double[] Fitness, int[] geneFitness){
        for(int i=0;i<Fitness.length;i++){
            int genotype=geneCompositionRecords.get(i);
            double fitness=0;
            int index=0;
                while(genotype>0){
                    if(genotype%2==1)
                        fitness+=geneFitness[index];
                    index++;
                    genotype=genotype>>1;
                }
            Fitness[i]=fitness/numOfGenesPerMicrobe;
        }
    }
    
    public void getFitness(double[] Fitness1,double[]Fitness2, int[] geneFitness1, int[] geneFitness2){
        for(int i=0;i<Fitness1.length;i++){
            int genotype=geneCompositionRecords.get(i);
            double fitness1=0;
            double fitness2=0;
            int index=0;
                while(genotype>0){
                    if(genotype%2==1){
                        fitness1+=geneFitness1[index];
                        fitness2+=geneFitness2[index];}
                    index++;
                    genotype=genotype>>1;
                }
            Fitness1[i]=fitness1/numOfGenesPerMicrobe;
            Fitness2[i]=fitness2/numOfGenesPerMicrobe;
        }
    }
    
    public void getMFitnessSelection(double[] microbiome,double[] fitnessToMicrobeRecords){
        double totalFitness=0;
        for (int i=0;i<microbiome.length;i++){
            microbiome[i]=microbiome[i]*Math.pow(msc,fitnessToMicrobeRecords[i]);
            totalFitness+=microbiome[i];
        }
        for (int i=0;i<microbiome.length;i++){
            microbiome[i]/=totalFitness;
        }
    }
    
    public ArrayList<Integer> getGeneComposition(){
        return geneCompositionRecords;
    }
    
    
    
}
