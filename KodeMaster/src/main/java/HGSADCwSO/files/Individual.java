package main.java.HGSADCwSO.files;

import main.java.HGSADCwSO.protocols.FitnessEvaluationProtocol;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class Individual {

    private double penalizedCost;
    private double scheduleCost;
    private double durationViolation;
    private double capacityViolation;
    private double deadlineViolation;
    private double durationViolationCost;
    private double capacityViolationCost;
    private double deadlineViolationCost;
    private double fitness;

    private double heuristicCost;

    private Set<Integer> departingVessels;

    private boolean feasibility;
    private boolean capacityFeasibility = true;
    private boolean durationFeasibility = true;
    private boolean deadlineFeasibility = true;

    private int costRank;
    private int diversityRank;
    private double biasedFitness;
    private double diversityContribution;

    private Genotype genotype;

    private FitnessEvaluationProtocol fitnessEvaluationProtocol;

    public Individual(HashMap<Integer, ArrayList<Integer>>  vesselTourChromosome, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        this.genotype = new Genotype(vesselTourChromosome);
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;
        setDepartingVessels();
    }

    private void setDepartingVessels(){
        Set<Integer> departingSet = new HashSet<>(); //initializing departingVessels
        HashMap<Integer, ArrayList<Integer>> vesselTourChromosome = new HashMap<>(getVesselTourChromosome());
        for (Integer vessel : vesselTourChromosome.keySet()){
            if (vesselTourChromosome.get(vessel).size()>0){
                departingSet.add(vessel);
            }
        }
        this.departingVessels = new HashSet<>(departingSet);
    }

    public Set<Integer> getDepartingVessels() {return departingVessels;}

    public HashMap<Integer, ArrayList<Integer>> getVesselTourChromosome() {
        return genotype.getVesselTourChromosome();
    }

    public void setVesselTourChromosome(HashMap<Integer, ArrayList<Integer>> vesselTour){
        this.genotype.setVesselTourChromosome(vesselTour);
    }

    public Genotype getGenotype() {
        return genotype;
    }

    public double getFitness() {
        return fitness;
    }

    public void setFeasibility(boolean feasibility) {
        this.feasibility = feasibility;
    }

    public void setFitness(double fitness) {
        this.fitness = fitness;
    }

    public boolean isFeasible(){
        return feasibility;
    }

    public double getPenalizedCost() {
        return penalizedCost;
    }

    public void setPenalizedCost() {
        penalizedCost = scheduleCost + durationViolationCost + deadlineViolationCost + capacityViolationCost;
        this.feasibility = (durationFeasibility && deadlineFeasibility && capacityFeasibility);
    }

    /*
    public void updatePenalizedCostForChromosome(FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        fitnessEvaluationProtocol.setPenalizedCostIndividual(this);
    }
    */

    public double getScheduleCost() {return scheduleCost; }

    //public void setScheduleCost() {this.scheduleCost = fitnessEvaluationProtocol.evaluate(this);}

    public double getDurationViolation() {return durationViolation; }

    public void setDurationViolation(double violation, double violationPenalty) {
        this.durationViolation = violation;
        this.durationViolationCost = violationPenalty*violation;
        if (durationViolation > 0) {
            this.durationFeasibility = false;
        }
    }

    public double getCapacityViolation() {return capacityViolation; }

    public void setCapacityViolation(double violation, double violationPenalty) {
        this.capacityViolation = violation;
        this.capacityViolationCost = violationPenalty*violation;
        if (capacityViolation > 0) {
            this.capacityFeasibility = false;
        }
    }

    public double getDeadlineViolation() {
        return deadlineViolation;
    }

    public void setDeadlineViolation(double deadlineViolation, double violationPenalty) {
        this.deadlineViolation = deadlineViolation;
        this.deadlineViolationCost = violationPenalty*deadlineViolation;
        if (capacityViolation > 0) {
            this.capacityFeasibility = false;
        }
    }

    public double getBiasedFitness() {
        return biasedFitness;
    }

    public double getDiversityContribution() {
        return diversityContribution;
    }

    public int getCostRank() {
        return costRank;
    }

    public int getDiversityRank() {
        return diversityRank;
    }

    public void setPenalizedCost(double penalizedCost) {
        this.penalizedCost = penalizedCost;
    }

    public void setScheduleCost(double scheduleCost) {
        this.scheduleCost = scheduleCost;
    }

    public void setBiasedFitness(double biasedFitness) {
        this.biasedFitness = biasedFitness;
    }

    public void setCostRank(int costRank) {
        this.costRank = costRank;
    }

    public void setDiversityContribution(double diversityContribution) {
        this.diversityContribution = diversityContribution;
    }

    public void setDiversityRank(int diversityRank) {
        this.diversityRank = diversityRank;
    }

    public void setHeuristicCost(double heuristicCost) {
        this.heuristicCost = heuristicCost;
    }

    public double getHeuristicCost() {
        return heuristicCost;
    }

    public void resetFeasability() {
        this.durationFeasibility = true;
        this.capacityFeasibility = true;
        this.deadlineFeasibility = true;
    }

    public boolean isCapacityFeasibility() {
        return capacityFeasibility;
    }

    public boolean isDeadlineFeasibility() {
        return deadlineFeasibility;
    }

    public boolean isDurationFeasibility() {
        return durationFeasibility;
    }
}
