package HGSADCwSO.implementations.DAG;

import HGSADCwSO.ProblemData;

import java.util.*;

import static java.lang.Math.max;

public class DAG {

    private double shortestPathCost;
    private double shortestFeasiblePathCost;

    private ArrayList<Integer> voyage;
    private ArrayList<Integer> voyageWithDepot;
    private ProblemData problemData;

    private HashMap<Integer, HashMap<Integer, Node>> graph;

    private double timePeriodLength;
    private double maxSpeed;
    private double maxSpeedWS2;
    private double maxSpeedWS3;
    private double minSpeed;
    private int vesselReturnTime;

    public DAG(ArrayList<Integer> voyage, int vesselReturnTime, ProblemData problemData) {
        this.voyage = voyage;
        this.problemData = problemData;
        this.graph = new HashMap<Integer, HashMap<Integer, Node>>();

        createVoyageWithDepot();

        this.timePeriodLength = this.problemData.getHeuristicParameterDouble("Length of time period");
        this.maxSpeed = this.problemData.getProblemInstanceParameterDouble("Max speed");
        this.minSpeed = this.problemData.getProblemInstanceParameterDouble("Min speed");
        this.maxSpeedWS2 = maxSpeed - problemData.getProblemInstanceParameterDouble("Impact on sailing from weather state 2");
        this.maxSpeedWS3 = maxSpeed - problemData.getProblemInstanceParameterDouble("Impact on sailing from weather state 3");
        this.vesselReturnTime = vesselReturnTime;

        buildGraph();
        doDijkstra();
        setCost();
    }

    private void createVoyageWithDepot(){
        voyageWithDepot = new ArrayList<Integer>(voyage);
        voyageWithDepot.add(0,0);
        voyageWithDepot.add(0);
    }

    private void buildGraph(){

        System.out.println("___________________________________________________");

        for (int i = 0; i < voyageWithDepot.size(); i++) {
            graph.put(i, new HashMap<Integer, Node>());
        }

        graph.get(0).put(0, new Node(0,0));

        for (int j = 0; j < voyageWithDepot.size()-1; j++ ) {


            for (int time : graph.get(j).keySet()){
                buildEdgesFromNode(graph.get(j).get(time), voyageWithDepot.get(j + 1));
            }
        }
    }

    private void buildEdgesFromNode(Node node, int destinationOrder) {
        double distance = problemData.getDistanceByIndex(problemData.getInstallationNumberByOrderNumber(node.getOrderNumber()),problemData.getInstallationNumberByOrderNumber(destinationOrder));
        double multiplier = 1/timePeriodLength;
        double startTime = node.getTime()/multiplier;

        double earliestTheoreticalArrivalTime = distance/maxSpeed;
        double latestArrivalTime = distance/minSpeed;


        double arrivalTime = latestArrivalTime;
        while (arrivalTime <= latestArrivalTime) {

            double consumption = 0;

            double[] sailingInfo = evaluateSailing(startTime, distance, arrivalTime);

            double speed_ = sailingInfo[0];
            double timeWS3 = sailingInfo[1];
            double timeWS2 = sailingInfo[2];

            if (speed_ > maxSpeed) {
                break;
            }

            consumption += calculateConsumptionInWeatherState(speed_, timeWS3, 3);
            consumption += calculateConsumptionInWeatherState(speed_, timeWS2, 2);
            consumption += calculateConsumptionInWeatherState(speed_, arrivalTime - startTime - timeWS3 - timeWS2, 0);

            double[] idlingInfo = idlingCalculation(arrivalTime, destinationOrder);

            consumption += idlingInfo[0];

            double[] servicingInfo = servicingCalculations(idlingInfo[1], destinationOrder);

            consumption += servicingInfo[0];
            double endTime = servicingInfo[1];
            if ((endTime*multiplier)%1 != 0) {
                endTime = Math.ceil(endTime*multiplier)/multiplier;
            }

            Node childNode;

            if (!graph.get(voyage.indexOf(node.getOrderNumber())).containsKey(((int)Math.round(node.getTime()*multiplier)))) {
                childNode = new Node((int) Math.round(endTime*multiplier), destinationOrder);
            }
            else {
                childNode = graph.get(voyage.indexOf(node.getOrderNumber())).get(((int)Math.round(node.getTime()*multiplier)));
            }


            node.addChildEdge(new Edge(node, childNode, consumption));

            arrivalTime += timePeriodLength;
        }
    }

    private double calculateConsumptionInWeatherState(double speed, double duration, int weatherState){
        if (weatherState == 3 && speed >= maxSpeedWS3) {
            return duration*cunsumptionFunction(maxSpeedWS3);
        }
        else if (weatherState == 2 && speed >= maxSpeedWS2) {
            return duration*cunsumptionFunction(maxSpeedWS2);
        }
        else return duration*cunsumptionFunction(speed);
    }

    private double cunsumptionFunction(double speed) {
        return (0.8125*speed*speed-13.00*speed+72.75);
    }

    private double[] idlingCalculation(double time, int orderNumber) { //TODO: fix so that this function also accounts for closing time
        double consumptionFactor = problemData.getProblemInstanceParameterDouble("Consumption factor idling");

        if (problemData.getWeatherStateByHour().get((int)time) != 3){
            return new double[] {0, time};
        }
        else {
            double consumption = 0;
            int i = (int)Math.ceil(time);
            while (problemData.getWeatherStateByHour().get(i) == 3) {
                i ++;
                consumption += consumptionFactor;
            }
            return new double[] {consumption, i};
        }
    }

    private double[] servicingCalculations(double time, int order) {
        double servicingTimeLeft = problemData.getOrdersByNumber().get(order).getDemand()*problemData.getHeuristicParameterDouble("Time per Hiv");
        double consumptionFactor = problemData.getProblemInstanceParameterDouble("Consumption factor servicing");

        double consumption = 0;

        while (servicingTimeLeft -1/problemData.getWeatherImpactByHour((int)Math.floor(time)) > 0) {
            if (time - Math.floor(time) > 0) {
                consumption += consumptionFactor * (time - Math.floor(time)) * problemData.getWeatherImpactByHour((int) Math.floor(time));
                servicingTimeLeft -= (1 - (time - Math.floor(time))) / problemData.getWeatherImpactByHour((int) Math.floor(time));
                time = Math.ceil(time);
            } else {
                consumption += consumptionFactor * problemData.getWeatherImpactByHour((int) Math.floor(time));
                servicingTimeLeft -= 1 / problemData.getWeatherImpactByHour((int) Math.floor(time));
                time++;
            }
        }
        if (servicingTimeLeft > 0){
            consumption += consumptionFactor*servicingTimeLeft*problemData.getWeatherImpactByHour((int)Math.floor(time));
            time += servicingTimeLeft;
        }

        return new double[] {consumption, time};
    }

    private double[] evaluateSailing(double startTime, double distance, double arrivalTime){
        double speed = distance/(arrivalTime - startTime);

        // First loop for W = 3

        double durationWS3 = getTimeInWS(startTime, arrivalTime, 3);
        double durationWS2 = getTimeInWS(startTime, arrivalTime, 2);

        if (durationWS2 + durationWS3 == arrivalTime - startTime && speed > (durationWS2*maxSpeedWS2 + durationWS3*maxSpeedWS3)/(durationWS2 + durationWS3)) {
            speed = maxSpeed+1;
        }

        else {
            if (speed > maxSpeedWS3) {
                speed += (durationWS3 * (speed - maxSpeedWS3)) / (arrivalTime - startTime - durationWS3);
            }


            if (speed > maxSpeedWS3) {
                speed += (durationWS2 * (speed - maxSpeedWS2)) / (arrivalTime - startTime - durationWS2 - durationWS3);
            }
        }

        return new double[] {speed, durationWS3, durationWS2};

    }

    private double getTimeInWS(double t1, double t2, int weatherState){
        double time = 0;
        time += isWeatherState(weatherState, (int) t1) ? 1 - t1%1 : 0;
        time += isWeatherState(weatherState, (int) t2) ? t2%1 : 0;
        for (int i = (int) Math.ceil(t1); i < (int) t2; i++ ){
            time += isWeatherState(weatherState, i) ? 1 : 0;
        }
        return time;
    }

    public boolean isWeatherState(int weatherState, int hour) {
        return problemData.getWeatherStateByHour().get(hour) == weatherState;
    }


    private void doDijkstra(){
        for (int i = 0; i < graph.size(); i++) {
            for (Map.Entry<Integer, Node> pair : graph.get(i).entrySet()){
                expand(pair.getValue());
            }
        }
        for (Map.Entry<Integer, Node> entry : graph.get(graph.size()-1).entrySet()){
            double nodeTime = entry.getValue().getTime();
            double nodeCost = entry.getValue().getBestCost();
            if (nodeTime > vesselReturnTime) {
                entry.getValue().setBestCost(nodeCost - (nodeTime - vesselReturnTime)*problemData.getHeuristicParameterDouble("Cost penalty per excessive time period") );
                entry.getValue().setFeasibility(false);
                }
            else {
                entry.getValue().setFeasibility(true);
            }
        }
    }

    private void expand(Node node) {
        for (Edge childEdge : node.getChildEdges()){
            if (childEdge.getChildNode().getBestCost() < node.getBestCost() + childEdge.getCost()) {
                childEdge.getChildNode().setBestCost(node.getBestCost() + childEdge.getCost());
                childEdge.getChildNode().setBestParentEdge(childEdge);
            }
        }
    }

    private void setCost() {
        if (graph.get(graph.size()).containsKey(vesselReturnTime)) {
            shortestPathCost = graph.get(graph.size()-1).get(vesselReturnTime).getBestCost();

        }

        double leastFeasibleCost = Double.POSITIVE_INFINITY;
        double leastInfeasileCost = Double.POSITIVE_INFINITY;

        for ( Map.Entry<Integer, Node> entry : graph.get(graph.size()).entrySet()){

            boolean nodeFeasibility = entry.getValue().getFeasibility();
            double nodeCost = entry.getValue().getBestCost();

            if (nodeCost < leastFeasibleCost || nodeCost < leastInfeasileCost) {
                if (nodeFeasibility) {
                    leastFeasibleCost = nodeCost;
                }
                else {
                    leastInfeasileCost = nodeCost;
                }
            }
        }
        shortestPathCost = max(leastFeasibleCost, leastInfeasileCost);
        shortestFeasiblePathCost = leastFeasibleCost;
    }

    public double getShortestPathCost(){
        return shortestPathCost;
    }

    public double getShortestFeasiblePathCost() {
        return shortestFeasiblePathCost;
    }
}
