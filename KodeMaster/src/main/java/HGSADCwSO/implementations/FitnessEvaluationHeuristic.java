package main.java.HGSADCwSO.implementations;

import main.java.HGSADCwSO.files.*;
import main.java.HGSADCwSO.protocols.SailingLegCalculationsProtocol;

import java.util.ArrayList;
import java.util.HashMap;

public class FitnessEvaluationHeuristic extends FitnessEvaluationBaseline {

    private ProblemData problemData;
    private SailingLegCalculationsProtocol sailingLegCalculationsProtocol;
    private double value;

    public FitnessEvaluationHeuristic(ProblemData problemData){
        super(problemData);
        this.problemData = problemData;
        selectProtocols();
    }

    public void evaluate(Individual individual) {

        double deadline_violation_penalty = super.deadlineViolationPenalty;
        double duration_violation_penalty = super.durationViolationPenalty;
        double capacity_violation_penalty = super.capacityViolationPenalty;

        evaluate(individual, duration_violation_penalty, capacity_violation_penalty, deadline_violation_penalty);
    }

    public void evaluate(Individual individual, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {

        individual.resetFeasability();
        Genotype genotype = individual.getGenotype();
        int nVessels = problemData.getNumberOfVessels();
        double cost = 0;
        double capacityViolation = 0;
        double durationViolation = 0;
        double deadlineViolation = 0;
        for (int i = 0; i < nVessels; i++){
            ArrayList<Integer> route = genotype.getVesselTourChromosome().get(i);
            Vessel vessel = problemData.getVesselByNumber().get(i);
            double[] info = evaluateRoute(route, vessel);
            cost += info[0];
            capacityViolation += info[1];
            deadlineViolation += info[2];
            durationViolation += info[3];
        }

        individual.setScheduleCost(cost);
        individual.setFeasibility(true);
        individual.setCapacityViolation(capacityViolation, capacityViolationPenalty);
        individual.setDeadlineViolation(deadlineViolation, deadlineViolationPenalty);
        individual.setDurationViolation(durationViolation, durationViolationPenalty);
        individual.setPenalizedCost();
    }


    public double[] evaluateRoute(ArrayList<Integer> route, Vessel vessel) {

        double totalConsumption = 0;
        double timeHorizon = 24*vessel.getReturnDay();
        double totalNumberOfHiv = 0;
        double totalDistance = 0;
        double totalCapReq = 0;
        Installation departureInstallation = problemData.getInstallationByNumber().get(0);
        Installation destinationInstallation;
        double fuelPrice = problemData.getProblemInstanceParameterDouble("fuel price");

        for (int i : route){
            totalNumberOfHiv += problemData.getOrdersByNumber().get(i).getDemand();
            destinationInstallation = problemData.getOrdersByNumber().get(i).getInstallation();
            totalDistance += problemData.getDistance(departureInstallation, destinationInstallation);
            departureInstallation = destinationInstallation;
            totalCapReq += problemData.getDemandByOrderNumber(i);
        }
        destinationInstallation = problemData.getInstallationByNumber().get(0);
        totalDistance += problemData.getDistance(departureInstallation, destinationInstallation);
        departureInstallation = destinationInstallation;

        double sailingTime = timeHorizon - totalNumberOfHiv/10;
        double sailingSpeed = Math.min(totalDistance/sailingTime, problemData.getProblemInstanceParameterDouble("Max speed"));


        double timeSpent = 0;
        double totalDeadlineViolatin = 0;

        for (int i = 0; i < route.size(); i++){
            int orderNumber = route.get(i);
            destinationInstallation = problemData.getOrdersByNumber().get(orderNumber).getInstallation();
            int weatherState = problemData.getWeatherStateByHour().get((int)timeSpent);
            double distance = problemData.getDistance(departureInstallation, destinationInstallation);

            sailingLegCalculationsProtocol.calculateSailingLeg(distance, sailingSpeed, timeSpent, problemData.getOrdersByNumber().get(orderNumber).getDemand());

            timeSpent = sailingLegCalculationsProtocol.getArrivalTime();
            totalConsumption += sailingLegCalculationsProtocol.getFuelConsumption();

            totalDeadlineViolatin += Math.max(0, timeSpent - problemData.getOrdersByNumber().get(orderNumber).getDeadline()*24);

            departureInstallation = destinationInstallation;
        }

        double capacityViolation = Math.max(0, totalCapReq - vessel.getCapacity());

        double durationViolation = Math.max(0, timeSpent - vessel.getReturnDay()*24);

        return new double[] {fuelPrice*totalConsumption, capacityViolation, totalDeadlineViolatin, durationViolation};
    }



    private void selectProtocols() {
        selectSailingLegCalculationsProtocol();
    }

    private void selectSailingLegCalculationsProtocol(){
        switch (problemData.getHeuristicParameters().get("Sailing leg calculations protocol")){
            case "quick and dirty":
                sailingLegCalculationsProtocol = new SailingLegCalculationsQuickAndDirty(problemData);
                break;
            default:
                sailingLegCalculationsProtocol = null;
                break;
        }
    }

    @Override
    public void setPenalizedCostIndividual(Individual individual, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {
        evaluate(individual, durationViolationPenalty, capacityViolationPenalty, deadlineViolationPenalty);
    }

    @Override
    public void setPenalizedCostIndividual(Individual individual) {
        evaluate(individual);
    }

    @Override
    public double getPenalizedCostOfVoyage(ArrayList<Integer> orderSequence, int vessel, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {

        HashMap<Integer, ArrayList<Integer>> tempChromosome = new HashMap<>();
        for (int i = 0; i < problemData.getNumberOfVessels(); i ++){
            if (i == vessel) {
                tempChromosome.put(i, orderSequence);
            }
            else {
                tempChromosome.put(i, new ArrayList<>());
            }
        }

        Individual tempIndividual = new Individual(tempChromosome,this);

        evaluate(tempIndividual, durationViolationPenalty, capacityViolationPenalty, deadlineViolationPenalty);

        return tempIndividual.getPenalizedCost();
    }

}
