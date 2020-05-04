package HGSADCwSO.implementations.DAG;

import HGSADCwSO.*;
import HGSADCwSO.implementations.FitnessEvaluationBaseline;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static java.lang.Math.max;


public class FitnessEvaluationDAG extends FitnessEvaluationBaseline { //TODO fix penalty for deadline violation by multiplying cost per hour, not per discrete time period

    private ProblemData problemData;
    
    private DAGSolver solver;
    
    private double nCloseProp;
    protected double nEliteProp;
    private double numberOfOrders;
    private HashMap<Individual, HashMap<Individual, Double>> hammingDistances;

    private int multiplier;

    private HashMap<Integer, HashMap<ArrayList<Integer>, double[]>> cachedVesselTours = new HashMap<>();
    private HashMap<ArrayList<Integer>, Graph> cachedGraphs = new HashMap<>();
    
    // EVT [0] = fitness, [1] = duration violation, [2] = deadline violation, [3] = capacity violation

    public FitnessEvaluationDAG(ProblemData problemData){
        super(problemData);
        this.problemData = problemData;
        for (int i = 0; i < problemData.getNumberOfVessels(); i++){
            cachedVesselTours.put(i, new HashMap<>());
        }
        this.multiplier = (int) problemData.getHeuristicParameterDouble("Number of time periods per hour");
    }

    @Override
    public void evaluate(Individual individual)  {
        double devp = problemData.getHeuristicParameterDouble("Deadline constraint violation penalty");
        double duvp = problemData.getHeuristicParameterDouble("Duration constraint violation penalty");
        double cvp = problemData.getHeuristicParameterDouble("Capacity constraint violation penalty");
        evaluate(individual, devp, cvp, duvp);
    }

    public void evaluate(Individual individual, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {

        Genotype genotype = individual.getGenotype();
        
        boolean feasibility = true;
        double scheduleCost = 0;
        double durationViolation = 0;
        double deadlineViolation = 0;
        double capacityViolation = 0;

        for (Vessel vessel : problemData.getVessels()){
            
            ArrayList<Integer> tour = genotype.getVesselTourChromosome().get(vessel.getNumber());

            if (tour.size() == 0) {
                scheduleCost += 0;
                durationViolation += 0;
                deadlineViolation += 0;
                capacityViolation += 0;
            }
            else {
                if (cachedVesselTours.get(vessel.getNumber()).containsKey(tour)) {
                    scheduleCost += cachedVesselTours.get(vessel.getNumber()).get(tour)[0];
                    durationViolation += cachedVesselTours.get(vessel.getNumber()).get(tour)[1];
                    deadlineViolation += cachedVesselTours.get(vessel.getNumber()).get(tour)[2];
                    capacityViolation += cachedVesselTours.get(vessel.getNumber()).get(tour)[3];
                }
                else {
                    Graph graph = getDAG(tour);
                    doDijkstra(graph, vessel.getReturnDay()*24*multiplier, durationViolationPenalty, deadlineViolationPenalty); //TODO check correct return time
                    double[] vesselTourInfo = getTourInfo(graph, vessel.getReturnDay()*24*multiplier);

                    scheduleCost = vesselTourInfo[0];
                    durationViolation = vesselTourInfo[1];
                    deadlineViolation = vesselTourInfo[2];


                    double capacityReqOfTour = 0;
                    for (int orderNumber : tour) {
                        capacityReqOfTour += problemData.getDemandByOrderNumber(orderNumber);
                    }
                    // System.out.println("Capacuty: " + vessel.getCapacity() + ", tour req: " + capacityReqOfTour);
                    capacityViolation += Math.max(0, capacityReqOfTour - vessel.getCapacity());

                    cachedVesselTours.get(vessel.getNumber()).put(tour, new double[] {scheduleCost, durationViolation, deadlineViolation, capacityViolation});
                }
            }
        }

        double durationViolationCost = durationViolation * durationViolationPenalty;
        double deadlineViolationCost = deadlineViolation * deadlineViolationPenalty;
        double capacityViolationCost = capacityViolation * capacityViolationPenalty;

        individual.setScheduleCost(scheduleCost);
        individual.setDurationViolation(durationViolation, durationViolationCost);
        individual.setDeadlineViolation(deadlineViolation, deadlineViolationCost);
        individual.setCapacityViolation(capacityViolation, capacityViolationCost);
        individual.setPenalizedCost();
    }
    
    public Graph getDAG(ArrayList<Integer> tour) {
        if (cachedGraphs.containsKey(tour)) {
            return cachedGraphs.get(tour);
        }
        else {
            // System.out.println("Tour : " + tour);
            Graph graph = new Graph(tour, problemData);
            cachedGraphs.put(tour, graph);
            return graph;
        }
    }


    //--------------------------------------- DIJKSTRA ----------------------------------------

    private void doDijkstra(Graph graph, int vesselReturnTime, double durationViolationPenalty, double deadlineViolationPenalty){

        // System.out.println("================================== DIJKSTRA ======================================");

        for (int i = 0; i < graph.getSize(); i++) {
            // System.out.println("1 ------");
            for (Map.Entry<Integer, Node> pair : graph.getGraph().get(i).entrySet()){
                // System.out.println("2");
                // System.out.println(pair.getValue().getChildEdges().size());

                expand(pair.getValue(), deadlineViolationPenalty);
            }
        }
        for (Map.Entry<Integer, Node> entry : graph.getGraph().get(graph.getSize()-1).entrySet()){
            double nodeTime = entry.getValue().getTime();
            double nodeCost = entry.getValue().getBestCost();
            if (nodeTime > vesselReturnTime) {
                entry.getValue().setBestCost(nodeCost + (nodeTime - vesselReturnTime)*durationViolationPenalty );
                entry.getValue().setFeasibility(false);
            }
            else {
                entry.getValue().setFeasibility(true);
            }
        }
    }

    private void expand(Node node, double deadlineViolationPenalty) {
        for (Edge childEdge : node.getChildEdges()){
            // System.out.println("3 !!!!");

            Node childNode = childEdge.getChildNode();

            double childNodeDeadlinePenaltyCost = (childNode.getBestTotalDeadlineViolation() + node.getDeadlineViolation()) * deadlineViolationPenalty;

            if (childNode.getBestCost() > node.getBestCost() + childEdge.getCost() + childNodeDeadlinePenaltyCost) {
                childNode.setBestCost(node.getBestCost() + childEdge.getCost() + childNodeDeadlinePenaltyCost);
                childNode.setBestParentEdge(childEdge);
                childNode.setBestTotalDeadlineViolation(childNode.getDeadlineViolation() + node.getBestTotalDeadlineViolation());
                // System.out.println("Parent node cost: " + node.getBestCost() + ", Edge cost: " + childEdge.getCost() + ", Child node cost: " + childEdge.getChildNode().getBestCost());
            }
        }
    }


    //--------------------------------------- SET COST ----------------------------------------



    private double[] getTourInfo(Graph graph, int vesselReturnTime) { //TODO REDO THIS CODE & INCLUDE FEASIBILITY

        double leastCost = Double.POSITIVE_INFINITY;
        double deadlineViolation = 0;
        double durationViolation = 0;

        if (graph.getGraph().get(graph.getSize()-1).containsKey(vesselReturnTime)) {
            leastCost = graph.getGraph().get(graph.getSize()-1).get(vesselReturnTime).getBestCost();
        }

        for ( Map.Entry<Integer, Node> entry : graph.getGraph().get(graph.getSize()-1).entrySet()){
            Node node = entry.getValue();
            if (node.getBestCost() < leastCost) {
                boolean nodeFeasibility =node.getFeasibility();
                leastCost = node.getBestCost();

                deadlineViolation = node.getBestTotalDeadlineViolation();

                durationViolation = node.getTime() - vesselReturnTime;
            }



        }

        // System.out.println("LFC : " + leastFeasibleCost);
        // System.out.println("LIC : " + leastInfeasibleCost);

        return new double[] {leastCost, durationViolation, deadlineViolation};
    }


    // -------------------------------------------- OVERRIDED FUNCTIONS ------------------------------------------------

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
