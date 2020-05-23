package HGSADCwSO.files;

import HGSADCwSO.implementations.*;
import HGSADCwSO.implementations.DAG.FitnessEvaluationDAG;
import HGSADCwSO.protocols.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class Process {

    private ProblemData problemData;

    private InitialPopulationProtocol initialPopulationProtocol;
    private ParentSelectionProtocol parentSelectionProtocol;
    private EducationProtocol educationProtocol;

    private FitnessEvaluationProtocol fitnessEvaluationProtocol;
    private FitnessEvaluationProtocol heuristicFitnessEvaluationProtocol;

    private ReproductionProtocol reproductionProtocol;
    private PenaltyAdjustmentProtocol penaltyAdjustmentProtocol;
    private DiversificationProtocol diversificationProtocol;
    private SurvivorSelectionProtocol survivorSelectionProtocol;

    private int processIteration;
    private HashMap<Integer, ArrayList<Individual>> feasibleSubPopulationByIteration = new HashMap<Integer, ArrayList<Individual>>();
    private HashMap<Integer, ArrayList<Individual>> infeasibleSubPopulationByIteration = new HashMap<Integer, ArrayList<Individual>>();
    private HashMap<Integer, Individual> bestFeasibleIndividualByIteration = new HashMap<Integer, Individual>();

    public Process(ProblemData problemData) {
        this.problemData = problemData;
        selectProtocols();
    }

    public Individual createIndividual() {
        return initialPopulationProtocol.createIndividual();
    }

    public ArrayList<Individual> selectParents(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
        ArrayList<Individual> entirePopulation = Utilities.getAllElements(feasiblePopulation, infeasiblePopulation);
        return parentSelectionProtocol.selectParents(entirePopulation);
    }

    public Individual mate(ArrayList<Individual> parents) {
        return reproductionProtocol.crossover(parents);
    }

    public void repair(Individual individual, double probability) {
        if (!individual.isFeasible()) {
            double randomDouble = new Random().nextDouble();
            if (randomDouble < probability) {
                int penaltyMultiplier = 10;
                educationProtocol.repairEducate(individual, penaltyMultiplier);
                if (!individual.isFeasible()) {
                    penaltyMultiplier = 100;
                    educationProtocol.repairEducate(individual, penaltyMultiplier);
                }
            }
        }
    }

    public void evaluate(Individual individual) {
        fitnessEvaluationProtocol.evaluate(individual);
        System.out.println("Violations: " + individual.getCapacityViolation() + " " + individual.getDurationViolation() + " " + individual.getDeadlineViolation());
    }

    public void repair(Individual individual) {
        double probability = problemData.getHeuristicParameterDouble("Repair rate");
        repair(individual, probability);
    }

    public void educate(Individual individual) {
        double chance = problemData.getHeuristicParameterDouble("Education rate");
        double luck = new Random().nextDouble();
        if (luck < chance){
            educationProtocol.educate(individual);
        }
    }

    public void mutate(Individual individual) {
        educationProtocol.inter_voyage_improvement(individual);
    }


    public Individual cloneGoodIndividual(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
        ArrayList<Individual> parents = selectParents(feasiblePopulation, infeasiblePopulation);
        if (parents.get(0).getPenalizedCost() < parents.get(1).getPenalizedCost()) {
            return new Individual(Utilities.deepCopyVesselTour(parents.get(0).getVesselTourChromosome()),fitnessEvaluationProtocol);
        }
        else {
            return new Individual(Utilities.deepCopyVesselTour(parents.get(1).getVesselTourChromosome()),fitnessEvaluationProtocol);
        }
    }

    public void survivorSelection(ArrayList<Individual> subpopulation, ArrayList<Individual> otherSubpopulation) {
        survivorSelectionProtocol.selectSurvivors(subpopulation, otherSubpopulation, fitnessEvaluationProtocol);
    }

    public void elite_training(Individual kid) { // TODO: Kuttt - gir ikke markant forbedring
        System.out.println();
        System.out.println("Welcome to the elite training programme!");
        System.out.println();
        int iterations_without_improvement = 0;
        double penalized_cost = kid.getPenalizedCost();
        while (iterations_without_improvement < 10) {
            educate(kid);
            if (kid.getPenalizedCost() < penalized_cost) {
                iterations_without_improvement = 0;
                penalized_cost = kid.getPenalizedCost();
                System.out.println("New penalized cost: " + penalized_cost);
            }
            else {
                iterations_without_improvement++;
            }
        }
        System.out.println();
        System.out.println("elite training over");
        System.out.println();
    }

    private void selectProtocols() {
        selectFitnessEvaluationProtocol();
        selectInitialPopulationProtocol();
        selectParentSelectionProtocol();
        selectPenaltyAdjustmentProtocol();
        selectEducationProtocol();
        selectReproductionProtocol();
        selectDiversificationProtocol();
        selectSurvivorSelectionProtocol();
    }

    private void selectSurvivorSelectionProtocol(){
        survivorSelectionProtocol = new SurvivorSelectionStandard(problemData);
    }

    private void selectDiversificationProtocol() {
        diversificationProtocol = new DiversificationAndStopping(problemData);
    }

    private void selectPenaltyAdjustmentProtocol() {
        penaltyAdjustmentProtocol = new PenaltyAdjustmentProtocol(problemData);
    }

    private void selectInitialPopulationProtocol() {
        switch (problemData.getHeuristicParameters().get("Initial population protocol")) {
            case "standard":
                initialPopulationProtocol = new InitialPopulationStandard(problemData, fitnessEvaluationProtocol);
                break;
            default:
                initialPopulationProtocol = null;
                break;
        }
    }

    private void selectParentSelectionProtocol(){
        switch (problemData.getHeuristicParameters().get("Parents selection protocol")) {
            case "binary tournament":
                parentSelectionProtocol = new ParentsSelectionBinaryTournament();
                break;
            default:
                parentSelectionProtocol = null;
                break;
        }
    }

    private void selectEducationProtocol(){
        switch (problemData.getHeuristicParameters().get("Education protocol")) {
            case "cost":
                educationProtocol = new EducationStandard(problemData, fitnessEvaluationProtocol, penaltyAdjustmentProtocol);
                break;
            default:
                educationProtocol = null;
                break;
        }
    }

    private void selectFitnessEvaluationProtocol() {
        switch (problemData.getHeuristicParameters().get("Fitness evaluation protocol")) {
            case "heuristic":
                fitnessEvaluationProtocol = new FitnessEvaluationHeuristic(problemData);
                heuristicFitnessEvaluationProtocol = new FitnessEvaluationHeuristic(problemData);
                break;
            case "dag":
                fitnessEvaluationProtocol = new FitnessEvaluationDAG(problemData);
                heuristicFitnessEvaluationProtocol = new FitnessEvaluationHeuristic(problemData);
                break;
            default:
                fitnessEvaluationProtocol = null;
                heuristicFitnessEvaluationProtocol = null;
                break;
        }
    }

    private void selectReproductionProtocol() {
        switch (problemData.getHeuristicParameters().get("Reproduction protocol")) {
            case "standard":
                reproductionProtocol = new ReproductionStandard(problemData, fitnessEvaluationProtocol);
                break;
            default:
                reproductionProtocol = null;
                break;
        }
    }

    public void adjustPenaltyParameters(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
        ArrayList<Individual> entirePopulation = Utilities.getAllElements(feasiblePopulation, infeasiblePopulation);
        penaltyAdjustmentProtocol.adjustPenalties(entirePopulation, fitnessEvaluationProtocol);
        penaltyAdjustmentProtocol.adjustPenalties(entirePopulation, heuristicFitnessEvaluationProtocol);
    }

    public boolean isDiversifyIteration() {
        return diversificationProtocol.isDiversifyIteration();
    }

    public boolean isStoppingIteration() {
        return diversificationProtocol.isStoppingIteration();
    }

    public void recordRunStatistics(int iteration, ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation, Individual bestFeasibleIndividual) {
        processIteration = iteration;
        feasibleSubPopulationByIteration.put(iteration, feasiblePopulation);
        infeasibleSubPopulationByIteration.put(iteration, infeasiblePopulation);
        bestFeasibleIndividualByIteration.put(iteration, bestFeasibleIndividual);

        // System.out.println("Iteration: " + iteration + " Feasible pop size: " + feasiblePopulation.size() + " Infeasible pop size: " + infeasiblePopulation.size());
    }

    public void addDiversityDistance(Individual kid) {
        fitnessEvaluationProtocol.addDiversityDistance(kid);
    }

    public void updateBiasedFitness(ArrayList<Individual> feasiblePopulation, ArrayList<Individual> infeasiblePopulation) {
        ArrayList<Individual> entirePopulation = Utilities.getAllElements(feasiblePopulation, infeasiblePopulation);
        fitnessEvaluationProtocol.updateBiasedFitness(entirePopulation);
    }

    public ArrayList<Individual> getClones(ArrayList<Individual> subpopulation){
        return survivorSelectionProtocol.getClones(subpopulation, fitnessEvaluationProtocol);
    }

    public FitnessEvaluationProtocol getFitnessEvaluationProtocol() {
        return fitnessEvaluationProtocol;
    }

    public void updateIterationsSinceImprovementCounter(boolean improvingSolutionFound) {
        diversificationProtocol.updateIterationsSinceImprovementCounter(improvingSolutionFound);
    }

    public void recordDiversification(int iteration) {
        diversificationProtocol.addDiversification(iteration);
        diversificationProtocol.resetDiversificationCounter();
    }

    public void updatePenaltyAdjustmentCounter(Individual individual) {
        penaltyAdjustmentProtocol.countAddedIndividual(individual);
    }

}
