package main.java.HGSADCwSO.files;


import main.java.HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Collections;

public class HGSADCwSOmain {

    private ArrayList<Individual> feasiblePopulation, infeasiblePopulation;
    private Individual bestFeasibleIndividual;
    private long startTime, stopTime;
    private ProblemData problemData;
    private Process process;

    private double bestCost = Double.POSITIVE_INFINITY;

    private int iteration;
    private double duration;

    private int run_id;

    private IO io;

    public HGSADCwSOmain(int run_id, int scenario_number , String disco_param) throws IOException {
        io = new IO(scenario_number);
        initialize(run_id, scenario_number, disco_param);
        io.writeRunInfo(problemData);
    }

    private void initialize(int run_id, int scenario_number, String disco_param) throws IOException {
        startTime = System.nanoTime();
        problemData = HackInitProblemData.hack(scenario_number);
        problemData.setHeuristicParameter("Number of time periods per hour", disco_param);
        this.run_id = run_id;
    }

    public String fullEvolutionaryRun() throws IOException {

        process = new Process(problemData);
        feasiblePopulation = new ArrayList<Individual>();
        infeasiblePopulation = new ArrayList<Individual>();
        bestFeasibleIndividual = null;
        bestCost = Double.POSITIVE_INFINITY;
        iteration = 1;
        problemData.printProblemData();


        // System.out.println("Creating initial population...");
        createInitialPopulation();
        process.updateIterationsSinceImprovementCounter(true);
        runEvolutionaryLoop();
        terminate();
        return getResultString();
    }

    private String getResultString() {
        String s = "";
        s += run_id + "| ";
        s += duration/1000000000 + "| ";
        s += bestCost;
        return s;
    }

    private void createInitialPopulation() throws IOException {
        int initialPopulationSize = problemData.getHeuristicParameterInt("Population size");
        for (int i = 0; i < initialPopulationSize; i++){
            System.out.println("======================== Iteration: " + iteration + " ========================" );
            Individual kid = process.createIndividual();
            process.educate(kid);
            if (! kid.isFeasible()) {
                double repairChance = problemData.getHeuristicParameterDouble("Repair rate");
                process.repair(kid, repairChance);
            }
            addToSubpopulation(kid);
            iteration++;
        }
    }

    private void runEvolutionaryLoop() throws IOException {
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

        ArrayList<Individual> parents = process.selectParents(feasiblePopulation, infeasiblePopulation);
        kid = process.mate(parents);

        return kid;
    }

    private void evolve() throws IOException {
        System.out.println("======================== Iteration: " + iteration + " ========================" );


        Individual kid = getNextIndividual();

        process.educate(kid);
        process.repair(kid);
        boolean isImprovingSolution = addToSubpopulation(kid);
        System.out.println("FSP size: " + feasiblePopulation.size() + " IFSP size: " + infeasiblePopulation.size());
        System.out.println("Chromosome: " + kid.getVesselTourChromosome() + "  |  this cost: " + kid.getPenalizedCost() + "  ||  Best cost thus far: " + bestCost + "  |  " + bestFeasibleIndividual.getVesselTourChromosome());
        System.out.println("isImpovingSolution: " + isImprovingSolution);
        process.updateIterationsSinceImprovementCounter(isImprovingSolution);
        process.adjustPenaltyParameters(feasiblePopulation, infeasiblePopulation);

        if (process.isDiversifyIteration()) {
            diversify(feasiblePopulation, infeasiblePopulation);
        }
    }

    private void diversify(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) throws IOException {

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


    public boolean addToSubpopulation(Individual kid) throws IOException {
        process.updatePenaltyAdjustmentCounter(kid);
        boolean isImprovingSolution = false;
        process.evaluate(kid);
        if (kid.isFeasible()) {
            if (bestFeasibleIndividual == null || kid.getPenalizedCost() < bestFeasibleIndividual.getPenalizedCost()) {
                bestFeasibleIndividual = kid;
                bestCost = kid.getPenalizedCost();
                // System.out.println("New best solution found!");
                double time = System.nanoTime() - startTime;
                io.writeImprovementIteration(iteration, time, bestCost);
                isImprovingSolution = true;
            }
            feasiblePopulation.add(kid);
        }
        else {
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

        if (subpopulation.size() >= maxPopulationSize) {
            process.survivorSelection(subpopulation, otherSubpopulation);
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

    private void terminate() throws IOException {
        stopTime = System.nanoTime();
        double duration = stopTime - startTime;
        this.duration = duration;
        io.writeFinalSolution(iteration, duration, bestFeasibleIndividual);
        // System.out.println("Final population: ");
        printPopulation();
        printRunStatistics();
        printBestSolution();
    }

    private void printPopulation(){

    }

    private void printRunStatistics(){

    }

    private void printBestSolution() {

        System.out.println("============ BEST SOLUTION: ============");
        System.out.println("(run id: " + run_id + " )");
        Individual winner = bestFeasibleIndividual;

        System.out.println("Vessel tour chromosome: " + winner.getVesselTourChromosome());
        System.out.println("Penalized cost:         " + winner.getPenalizedCost());
        System.out.println();
        process.print_schedule(winner);

        //process.getFitnessEvaluationProtocol().getSolutionFromIndividual(winner);
        // process.getTourOfVessel

    }

}
