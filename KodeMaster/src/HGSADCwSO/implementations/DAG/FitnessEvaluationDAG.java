package HGSADCwSO.implementations.DAG;

import HGSADCwSO.*;
import HGSADCwSO.implementations.FitnessEvaluationBaseline;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import static java.lang.Math.max;


public class FitnessEvaluationDAG extends FitnessEvaluationBaseline { //TODO fix penalty for deadline violation by multiplying cost per hour, not per discrete time period

    private ProblemData problemData;
    
    private double nCloseProp;
    protected double nEliteProp;
    private double numberOfOrders;
    private HashMap<Individual, HashMap<Individual, Double>> hammingDistances;

    private int multiplier;
    private int maxSizeCachedVesselTours;
    private int maxSizeCachedGraphs;

    private HashMap<Integer, LinkedHashMap<ArrayList<Integer>, double[]>> cachedVesselTours= new HashMap<>();

    final LinkedHashMap<ArrayList<Integer>, Graph> cachedGraphs = new LinkedHashMap<ArrayList<Integer>, Graph>() {
        @Override
        protected boolean removeEldestEntry(final Map.Entry eldest) {
            return size() > maxSizeCachedGraphs;
        }
    };
    
    // EVT [0] = fitness, [1] = duration violation, [2] = deadline violation, [3] = capacity violation

    public FitnessEvaluationDAG(ProblemData problemData){
        super(problemData);
        this.problemData = problemData;
        this.multiplier = (int) problemData.getHeuristicParameterDouble("Number of time periods per hour");

        maxSizeCachedVesselTours = problemData.getHeuristicParameterInt("Max cached tours per vessel");
        maxSizeCachedGraphs = problemData.getHeuristicParameterInt("Max cached graphs");

        for (int i = 0; i < problemData.getNumberOfVessels(); i++) {

            cachedVesselTours.put(i,
                new LinkedHashMap<ArrayList<Integer>, double[]>() {
                    @Override
                    protected boolean removeEldestEntry(final Map.Entry eldest) {
                        return size() > maxSizeCachedVesselTours;
                    }
                }
            );
        }
    }

    @Override
    public void evaluate(Individual individual)  {
        double deadline_violation_penalty = super.deadlineViolationPenalty;
        double duration_violation_penalty = super.durationViolationPenalty;
        double capacity_violation_penalty = super.capacityViolationPenalty;

        evaluate(individual, duration_violation_penalty, capacity_violation_penalty, deadline_violation_penalty);
    }

    public void evaluate(Individual individual, double durationViolationPenalty, double capacityViolationPenalty, double deadlineViolationPenalty) {

        // System.out.println("dur: " + durationViolationPenalty + "   dead:  " + deadlineViolationPenalty + "   cap: " + capacityViolationPenalty);

        individual.resetFeasability();
        Genotype genotype = individual.getGenotype();

        double scheduleCost = 0;
        double durationViolation = 0;
        double deadlineViolation = 0;
        double capacityViolation = 0;
        double duration_spot_ship_rental = 0;

        for (Vessel vessel : problemData.getVessels()){
            
            ArrayList<Integer> tour = genotype.getVesselTourChromosome().get(vessel.getNumber());

            if (tour.size() == 0) {
                continue;
            }

            int vessel_number = vessel.getNumber();

            if (cachedVesselTours.get(vessel_number).containsKey(tour)) {
                scheduleCost += cachedVesselTours.get(vessel_number).get(tour)[0];
                durationViolation += cachedVesselTours.get(vessel_number).get(tour)[1];
                deadlineViolation += cachedVesselTours.get(vessel_number).get(tour)[2];
                capacityViolation += cachedVesselTours.get(vessel_number).get(tour)[3];
                duration_spot_ship_rental += cachedVesselTours.get(vessel_number).get(tour)[4];
            }
            else {
                Graph graph = getDAG(tour);
                doDijkstra(graph, vessel.getReturnDay()*24*multiplier, durationViolationPenalty, deadlineViolationPenalty); //TODO check correct return time
                double[] vesselTourInfo = getTourInfo(graph, vessel.getReturnDay()*24*multiplier);

                scheduleCost += vesselTourInfo[0];
                durationViolation += vesselTourInfo[1];
                deadlineViolation += vesselTourInfo[2];
                duration_spot_ship_rental += vessel.is_spot_vessel() ? vesselTourInfo[3] : 0;


                double capacityReqOfTour = 0;
                for (int orderNumber : tour) {
                    capacityReqOfTour += problemData.getDemandByOrderNumber(orderNumber);
                }
                // System.out.println("Capacuty: " + vessel.getCapacity() + ", tour req: " + capacityReqOfTour);
                capacityViolation += Math.max(0, capacityReqOfTour - vessel.getCapacity());

                cachedVesselTours.get(vessel.getNumber()).put(tour, new double[] {vesselTourInfo[0], vesselTourInfo[1], vesselTourInfo[2], Math.max(0, capacityReqOfTour - vessel.getCapacity()), vessel.is_spot_vessel() ? vesselTourInfo[3] : 0});

                //System.out.println("Size of vessel tour cache for vessel number " + vessel.getNumber() + " is: " + cachedVesselTours.get(vessel.getNumber()).size());
            }
        }

        double spot_vessel_cost = duration_spot_ship_rental * problemData.getProblemInstanceParameterDouble("Hourly spot rate");


        individual.setScheduleCost(scheduleCost + spot_vessel_cost);
        individual.setDurationViolation(durationViolation, durationViolationPenalty);
        individual.setDeadlineViolation(deadlineViolation, deadlineViolationPenalty);
        individual.setCapacityViolation(capacityViolation, capacityViolationPenalty);
        individual.setPenalizedCost();

        //individual.setSchedule(schedule)

        // System.out.println("Schedule cost: " + scheduleCost + " penalized: " + individual.getPenalizedCost() + " || Violations | cap: " + capacityViolation + ", dead: "+ deadlineViolation + ", dur: "+ durationViolation);
    }
    
    public Graph getDAG(ArrayList<Integer> tour) {
        if (cachedGraphs.containsKey(tour)) {
            return cachedGraphs.get(tour);
        }
        else {
            // System.out.println("Tour : " + tour);
            Graph graph = new Graph(tour, problemData);
            cachedGraphs.put(tour, graph);

            // System.out.println("Size of graph cache: " + cachedGraphs.size());
            return graph;
        }
    }


    //--------------------------------------- DIJKSTRA ----------------------------------------

    private void doDijkstra(Graph graph, int vesselReturnTime, double hourlyDurationViolationPenalty, double hourlyDeadlineViolationPenalty){

        double deadlineViolationPenalty = hourlyDeadlineViolationPenalty/multiplier;
        double durationViolationPenalty = hourlyDurationViolationPenalty/multiplier;

        // System.out.println("================================== DIJKSTRA ======================================");

        for (int i = 0; i < graph.getSize(); i++) {
            for (Map.Entry<Integer, Node> pair : graph.getGraph().get(i).entrySet()){
                expand(pair.getValue(), deadlineViolationPenalty);
            }
        }

        for (Map.Entry<Integer, Node> entry : graph.getGraph().get(graph.getSize()-1).entrySet()){
            Node node = entry.getValue();
            double nodeTime = node.getTime();
            double nodePenalizedCost = node.getBestPenalizedCost();
            if (nodeTime > vesselReturnTime) {
                node.setBestPenalizedCost(nodePenalizedCost + (nodeTime - vesselReturnTime)*durationViolationPenalty );
            }
        }
    }

    private void expand(Node node, double deadlineViolationPenalty) {

        for (Edge childEdge : node.getChildEdges()){
            // System.out.println("3 !!!!");


            Node childNode = childEdge.getChildNode();

            // System.out.println("Expanding: " + node.getOrderNumber() + " " + node.getTime() + " | " + childNode.getOrderNumber() + " " + childNode.getTime() + " | " + node.getBestCost() + " " + node.getBestPenalizedCost());

            double childNodeDeadlinePenaltyCost = childNode.getDeadlineViolation() * deadlineViolationPenalty;

            if (childNode.getBestPenalizedCost() > node.getBestPenalizedCost() + childEdge.getCost() + childNodeDeadlinePenaltyCost) {

                // System.out.println("Pc: " + childNode.getBestPenalizedCost() + " c: " + childNode.getBestCost());
                // System.out.println("Edge: " + childEdge.getCost());
                // System.out.println();



                childNode.setBestPenalizedCost(node.getBestPenalizedCost() + childEdge.getCost() + childNodeDeadlinePenaltyCost);
                childNode.setBestCost(node.getBestCost() + childEdge.getCost());
                childNode.setBestParentEdge(childEdge);
                childNode.setBestTotalDeadlineViolation(childNode.getDeadlineViolation() + node.getBestTotalDeadlineViolation());
                // System.out.println("Parent node cost: " + node.getBestCost() + ", Edge cost: " + childEdge.getCost() + ", Child node cost: " + childEdge.getChildNode().getBestCost());
            }
        }
    }


    //--------------------------------------- SET COST ----------------------------------------



    private double[] getTourInfo(Graph graph, int vesselReturnTime) {

        double leastPenalizedCost = Double.POSITIVE_INFINITY;
        double leastCost = Double.POSITIVE_INFINITY;
        double deadlineViolation = 0;
        double durationViolation = 0;
        double duration_in_hours = 0;

        for ( Map.Entry<Integer, Node> entry : graph.getGraph().get(graph.getSize()-1).entrySet()){

            // System.out.println("Least penalized cost:   " + leastPenalizedCost);
            // System.out.println("Least cost:             " + leastCost);

            Node node = entry.getValue();

            // System.out.println(node.getBestPenalizedCost() + " < " + leastPenalizedCost);
            if (node.getBestPenalizedCost() < leastPenalizedCost) {


                leastPenalizedCost = node.getBestPenalizedCost();

                // System.out.println(leastPenalizedCost + "    !!!!!!!!!!!!!!!! ");
                leastCost = node.getBestCost();
                deadlineViolation = node.getBestTotalDeadlineViolation()/multiplier;
                durationViolation = Math.max(0, (node.getTime() - vesselReturnTime)/multiplier);
                duration_in_hours = ((double) node.getTime() - 8 * multiplier)/multiplier;
            }
            // System.out.println("G - PC: " + node.getBestPenalizedCost());
            // System.out.println("G - C:  " + node.getBestCost());
        }

        // // System.out.println("LFC : " + leastFeasibleCost);
        // // System.out.println("LIC : " + leastInfeasibleCost);

        // System.out.println();
        // System.out.println("FInal: " + leastCost + "   " + deadlineViolation + "   " + durationViolation);
        // System.out.println();

        return new double[] {leastCost, durationViolation, deadlineViolation, duration_in_hours};
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
/*
    public Schedule getSolutionFromIndividual(Individual individual){
        HashMap<Integer, ArrayList<Integer>> chromosome = individual.getVesselTourChromosome();

        HashMap<Integer, ArrayList<SailingLeg>> sol = new HashMap<>();

        for (int i = 0; i < chromosome.size(); i++){
            sol.put(i, new ArrayList<SailingLeg>());
            Graph graph = getDAG(chromosome.get(i));
            doDijkstra(graph, problemData.getVesselByNumber().get(i).getReturnDay()*24*multiplier,);
        }
    }

 */


}
