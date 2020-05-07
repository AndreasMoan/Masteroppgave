package HGSADCwSO;


import HGSADCwSO.protocols.FitnessEvaluationProtocol;

import javax.sound.midi.Soundbank;
import java.sql.SQLOutput;
import java.util.ArrayList;
import java.util.Random;
import java.util.Collections;

public class HGSADCwSOmain {

    private ArrayList<Individual> feasiblePopulation, infeasiblePopulation;
    private Individual bestFeasibleIndividual;
    private long startTime, stopTime;
    private ProblemData problemData;
    private Process process;
    private double mutationChance;

    private double bestCost = Double.POSITIVE_INFINITY;
    private double scheduleCost = Double.POSITIVE_INFINITY;

    private String[] args;

    private int iteration;

    private double repairChance;
    private double educationChance;

    public static int COST_EDUCATIONS = 0;
    public static int INFEASIBLE_EDUCATIONS = 0;

    public static void main(String[] args) {
        HGSADCwSOmain main = new HGSADCwSOmain();
        main.initialize(args);
        main.fullEvolutionaryRun();
    }


    private void initialize(String[] changeParameters) {
        startTime = System.nanoTime();
        //io = new IO(inputFileName); //TODO
        //this.args = changeParameter;
        problemData = HackInitProblemData.hack();
    }

    private void fullEvolutionaryRun(){

        process = new Process(problemData);
        feasiblePopulation = new ArrayList<Individual>();
        infeasiblePopulation = new ArrayList<Individual>();
        bestFeasibleIndividual = null;

        mutationChance = problemData.getHeuristicParameterDouble("Mutation chance");

        iteration = 1;
        problemData.printProblemData(); //TODO
        System.out.println("Creating initial population...");
        createInitialPopulation();

        process.updateIterationsSinceImprovementCounter(true);

        runEvolutionaryLoop();
        terminate();
    }

    private void createInitialPopulation() {
        int initialPopulationSize = problemData.getHeuristicParameterInt("Population size");
        for (int i = 0; i < initialPopulationSize; i++){
            Individual kid = process.createIndividual();
            process.educate(kid);
            if (! kid.isFeasible()) {
                process.repair(kid, repairChance);
            }

            addToSubpopulation(kid);
        }
    }

    private void runEvolutionaryLoop() {
        process.recordRunStatistics(0, feasiblePopulation, infeasiblePopulation, bestFeasibleIndividual);
        while (!stoppingCriterion()) {
            // System.out.println("Iteration " + iteration + "                 Best cost thus far: " + bestCost);
            evolve();
            process.recordRunStatistics(iteration, feasiblePopulation, infeasiblePopulation, bestFeasibleIndividual);
            iteration++;
        }
    }

    private boolean stoppingCriterion() {
        return process.isStoppingIteration();
    }

    private Individual getNextIndividual() {

        Individual kid;

        if (new Random().nextDouble() < mutationChance) {
            kid = process.cloneGoodIndividual(feasiblePopulation,infeasiblePopulation);
            process.mutate(kid);
            process.mutate(kid);
            process.mutate(kid);
        }
        else {
            ArrayList<Individual> parents = process.selectParents(feasiblePopulation, infeasiblePopulation);
            kid = process.mate(parents);
        }
        return kid;
    }

    private void evolve() {

        Individual kid = getNextIndividual();

        process.educate(kid);
        process.repair(kid);
        boolean isImprovingSolution = addToSubpopulation(kid);
        if (kid.getPenalizedCost() < bestCost) {
            bestCost = kid.getPenalizedCost();
            scheduleCost = kid.getScheduleCost();
        }
        System.out.println();
        System.out.println("Iteration " + iteration + ", Schedule cost: " + scheduleCost + ",  Best cost thus far: " + bestCost + ",      Cost this iteration: " + kid.getPenalizedCost() + ",      Chromosome: " + kid.getVesselTourChromosome());
        process.updateIterationsSinceImprovementCounter(isImprovingSolution);
        process.adjustPenaltyParameters(feasiblePopulation, infeasiblePopulation);

        if (process.isDiversifyIteration()) {
            diversify(feasiblePopulation, infeasiblePopulation);
        }
    }

    private void diversify(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {

        System.out.println("Diversifying...");

        genocide(feasiblePopulation, infeasiblePopulation, 2.0/3.0);
        genocide(infeasiblePopulation, feasiblePopulation,  2.0/3.0);

        System.out.println("Breeding new population...");

        createInitialPopulation();

        process.recordDiversification(iteration);
    }

    private void genocide(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation, double proportionToKill){

        Collections.sort(subpopulation, Collections.reverseOrder(Utilities.getBiasedFitnessComparator()));


        /*for (Individual ind : subpopulation) {
            System.out.println(ind.getBiasedFitness() + ", with penalized cost: " + ind.getPenalizedCost());
        }*/


        int numberOfIndividualsToKill = (int) Math.round(subpopulation.size()*proportionToKill);
        ArrayList<Individual> individualsToBeKilled = new ArrayList<>(subpopulation.subList(0,numberOfIndividualsToKill));
        removeFromSubpopulation(subpopulation, otherSubpopulation, individualsToBeKilled);
    }


    public boolean addToSubpopulation(Individual kid) {
        process.updatePenaltyAdjustmentCounter(kid);
        boolean isImprovingSolution = false;
        process.evaluate(kid);

        if (kid.isFeasible()) {
            if ((bestFeasibleIndividual == null) || (kid.getPenalizedCost() < bestFeasibleIndividual.getPenalizedCost())) {
                bestFeasibleIndividual = kid;
                isImprovingSolution = true;
            }
            feasiblePopulation.add(kid);
        }
        else {
            if ((bestFeasibleIndividual == null) || (kid.getPenalizedCost() < bestFeasibleIndividual.getPenalizedCost())) {
                bestFeasibleIndividual = kid;
                isImprovingSolution = true;
            }
            infeasiblePopulation.add(kid);
        }
        process.addDiversityDistance(kid);
        process.updateBiasedFitness(feasiblePopulation, infeasiblePopulation);

        if (kid.isFeasible()) {
            checkSubpopulationSize(feasiblePopulation, infeasiblePopulation);
        }
        else {
            checkSubpopulationSize(infeasiblePopulation, feasiblePopulation);
        }
        return isImprovingSolution;
    }

    private void checkSubpopulationSize(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation) {
        int maxPopulationSize = problemData.getHeuristicParameterInt("Population size")
                + problemData.getHeuristicParameterInt("Number of offspring in a generation");

        if (subpopulation.size() + otherSubpopulation.size() >= maxPopulationSize) {
            genocide(subpopulation, otherSubpopulation, 3.0/4.0);
            genocide(otherSubpopulation, subpopulation, 3.0/4.0);
        }
    }

    private void removeFromSubpopulation(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation, ArrayList<Individual> individualsToKill) {

        for (Individual individual : individualsToKill) {
            HGSADCwSOmain.removeFromSubpopulation(subpopulation, individual, otherSubpopulation, process.getFitnessEvaluationProtocol(), false);
        }
    }

    public static void removeFromSubpopulation(ArrayList<Individual> subpopulation, Individual individual, ArrayList<Individual> otherSubpopulation, FitnessEvaluationProtocol fitnessEvaluationProtocol, boolean updateFitness) {
        subpopulation.remove(individual);
        fitnessEvaluationProtocol.removeDiversityDistance(individual);

        if (updateFitness){ // Updates fitness for all individuals
            fitnessEvaluationProtocol.updateBiasedFitness(Utilities.getAllElements(subpopulation, otherSubpopulation));
        }
    }

    private void terminate() {
        System.out.println("Final population: ");
        printPopulation();
        printRunStatistics();
        printBestSolution();
        stopTime = System.nanoTime();
    }

    private void printPopulation(){

    }

    private void printRunStatistics(){

    }

    private void printBestSolution() {
        Individual winner = bestFeasibleIndividual;

        System.out.println("Vessel tour chromosome" + winner.getVesselTourChromosome());

        //process.getFitnessEvaluationProtocol().getSolutionFromIndividual(winner);
        // process.getTourOfVessel

    }

}
