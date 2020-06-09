package main.java.HGSADCwSO.implementations.DAG;

import main.java.HGSADCwSO.files.Order;
import main.java.HGSADCwSO.files.ProblemData;

import java.util.ArrayList;
import java.util.HashMap;

import static java.lang.Math.ceil;
import static java.lang.Math.floor;

public class Graph {

    private double scheduleCost;
    private boolean feasibility;

    private ArrayList<Integer> voyage;
    private ArrayList<Integer> voyageWithDepot;
    private ProblemData problemData;

    private HashMap<Integer, HashMap<Integer, Node>> graph;

    private double maxSpeed;
    private double maxSpeedWS2;
    private double maxSpeedWS3;
    private double minSpeed;
    private int vesselReturnTime;
    private double multiplier;
    private double timePerHiv;
    private double speedImpactWS2;
    private double speedImpactWS3;
    private double idlingConsumption;
    private double servicingConsumption;

    private HashMap<Integer, Order> orderByNumber;


    public Graph(ArrayList<Integer> voyage, ProblemData problemData) {
        this.voyage = voyage;
        // // System.out.println(voyage);
        this.problemData = problemData;
        this.graph = new HashMap<Integer, HashMap<Integer, Node>>();

        createVoyageWithDepot();

        this.maxSpeed = this.problemData.getProblemInstanceParameterDouble("Max speed");
        this.minSpeed = this.problemData.getProblemInstanceParameterDouble("Min speed");
        this.speedImpactWS2 = problemData.getProblemInstanceParameterDouble("Impact on sailing from weather state 2");
        this.speedImpactWS3 = problemData.getProblemInstanceParameterDouble("Impact on sailing from weather state 3");
        this.idlingConsumption = problemData.getProblemInstanceParameterDouble("Idling consumption");
        this.servicingConsumption = problemData.getProblemInstanceParameterDouble("Servicing consumption");
        this.maxSpeedWS2 = maxSpeed - speedImpactWS2;
        this.maxSpeedWS3 = maxSpeed - speedImpactWS3;
        this.timePerHiv = problemData.getProblemInstanceParameterDouble("Time per hiv");
        this.multiplier = this.problemData.getHeuristicParameterDouble("Number of time periods per hour");


        this.orderByNumber = problemData.getOrdersByNumber();

        buildGraph();
    }

    private void createVoyageWithDepot(){
        voyageWithDepot = new ArrayList<Integer>(voyage);
        voyageWithDepot.add(0,0);
        voyageWithDepot.add(0);
    }


    //--------------------------------------- BUILD GRAPH ----------------------------------------



    private void buildGraph(){


        for (int i = 0; i < voyageWithDepot.size(); i++) {
            graph.put(i, new HashMap<Integer, Node>());
        }

        graph.get(0).put(0, new Node( 16 * (int) multiplier ,0, 0));
        graph.get(0).get(0).setBestCost(0);
        graph.get(0).get(0).setBestPenalizedCost(0);
        // System.out.println("------- graph: ---------");
        for (int j = 0; j < voyageWithDepot.size() - 1; j++ ) {

            for (int time : graph.get(j).keySet()){


                //// System.out.println();
                //// System.out.println("Building: " + graph.get(j).get(time).getOrderNumber() + " " + voyageWithDepot.get(j + 1));
                //// System.out.println();
                buildEdgesFromNode(graph.get(j).get(time), voyageWithDepot.get(j + 1), j);
            }
        }
    }

    private void buildEdgesFromNode(Node node, int destinationOrderNumber, int legNumber) {

        // System.out.println("Creating edges from time: " + node.getTime());


        double distance = problemData.getDistanceByIndex(problemData.getInstallationNumberByOrderNumber(node.getOrderNumber()),problemData.getInstallationNumberByOrderNumber(destinationOrderNumber));

        int nodeStartTime = node.getTime();
        double real_start_time = convertNodeTimeToRealTime(nodeStartTime);

        // System.out.println("From node Iteration  -  order number: " +  node.getOrderNumber());


        int earliestTheoreticalEndTime = nodeStartTime + (int) ceil((distance/maxSpeed + problemData.getDemandByOrderNumber(destinationOrderNumber)*timePerHiv)*multiplier);
        int latestTheoreticalEndTime = nodeStartTime + (int) ceil((distance/minSpeed + problemData.getDemandByOrderNumber(destinationOrderNumber)*timePerHiv*1.3)*multiplier);
        // System.out.println("a " + earliestTheoreticalEndTime + " " + latestTheoreticalEndTime);

        int finServicingTime = earliestTheoreticalEndTime;

        int counter = 0;
        while (/*continue_ && */ finServicingTime <= latestTheoreticalEndTime ) {

            finServicingTime = getEarliestFeasibleSercivingFinishingTime(finServicingTime , destinationOrderNumber, 0);

            // System.out.println("fst: " + finServicingTime);

            double[] serviceInfo = servicingCalculations(finServicingTime, destinationOrderNumber);

            double servicingCost = serviceInfo[0];
            double real_service_time = serviceInfo[1];

            // System.out.println("sc:  " + servicingCost);
            // System.out.println("rst: " + real_service_time);

            if (!isArrivalPossible(real_start_time, distance, real_service_time)) {
                finServicingTime++;
                continue;
            }

            double[] idlingInfo = idlingCalculations(real_start_time, distance, real_service_time);

            double idlingCost = idlingInfo[0];
            double real_arrival_time = idlingInfo[1];

            // System.out.println("ic:  " + idlingCost);
            // System.out.println("rat: " + real_arrival_time);

            double sailingCost = 0;

            if (distance != 0) {
                double[] timeInAllWeatherStates = getTimeInAllWS(real_start_time, real_arrival_time);
                // System.out.println("tws: " + timeInAllWeatherStates[0] + " " + timeInAllWeatherStates[1] + " " + timeInAllWeatherStates[2] + " " + timeInAllWeatherStates[3]);
                double adjustedAverageSpeed = calculateAdjustedAverageSpeed(timeInAllWeatherStates, distance);
                // System.out.println("aas: " + adjustedAverageSpeed);
                sailingCost = (real_arrival_time == real_start_time) ? 0 : sailingCalculations(timeInAllWeatherStates, adjustedAverageSpeed);
                // System.out.println("sc:  " + sailingCost);
                // System.out.println("d:   " + distance);
            }

            // System.out.println("vessel: " + "?" + " dep order:" + node.getOrderNumber() + " dest order: " + destinationOrderNumber + " start: " + nodeStartTime + " finish: " + finServicingTime + " c sail: " + sailingCost + " c idle: " + idlingCost + " c ser: " + servicingCost);
            Node childNode;

            if (isNodeInGraph(legNumber + 1, finServicingTime)) {
                childNode = graph.get(legNumber + 1).get(finServicingTime);
            }
            else {
                int destination_order_number = voyageWithDepot.get(legNumber + 1);

                int deadlineNodeTime = problemData.getOrderDeadlineByNumber(destination_order_number) * 24 * (int) multiplier;
                int deadlineViolation = Math.max(0, finServicingTime - deadlineNodeTime);

                childNode = new Node(finServicingTime, destinationOrderNumber, deadlineViolation);
                graph.get(legNumber + 1).put(finServicingTime, childNode);
            }

            double total_consumption_tonnes = (sailingCost + idlingCost + servicingCost)/1000;
            double total_fuel_cost = total_consumption_tonnes * problemData.getProblemInstanceParameterDouble("Fuel price");

            double real_end_time = convertNodeTimeToRealTime(finServicingTime);

            Edge currEdge = new Edge(node, childNode, total_fuel_cost, real_start_time, real_arrival_time, real_service_time, real_end_time, distance);

            childNode.addParentEdge(currEdge);
            node.addChildEdge(currEdge);

            counter ++;

            finServicingTime++;
        }
    }

    private double convertNodeTimeToRealTime(int nodeTime) {
        return nodeTime/multiplier;
    }

    private int convertRealTimeToNodeTimeFloor(double time) {
        return (int) floor(time/multiplier);
    }

    private boolean isServicePossible(int finServicingNodeTime, int destinationOrder) {

        // // System.out.println("Service possible check");

        double realTime = convertNodeTimeToRealTime(finServicingNodeTime);

        double servicingTimeLeft = problemData.getDemandByOrderNumber(destinationOrder) * timePerHiv;

        while (servicingTimeLeft > 0) {
            if (problemData.getWeatherStateByHour().get((int) realTime) == 3 || problemData.isInstallationByOrderIndexClosed(destinationOrder, realTime)) {
                return false;
            }
            if (realTime % 1 > 0) {
                if (servicingTimeLeft < (realTime % 1)/ problemData.getWeatherImpactByHour((int) floor(realTime))) {
                    return true;
                }
                else {
                    servicingTimeLeft -= (realTime % 1) / problemData.getWeatherImpactByHour((int) floor(realTime));
                    realTime = floor(realTime);
                }
            }
            else {
                if (servicingTimeLeft < 1 / problemData.getWeatherImpactByHour((int) floor(realTime))){
                    return true;
                }
                else {
                    servicingTimeLeft -= 1 / problemData.getWeatherImpactByHour((int) floor(realTime));
                    realTime--;
                }
            }
        }

        return true;
    }

    private int getEarliestFeasibleSercivingFinishingTime(int finServicingTime, int destinationOrderNumber, int nIterationsIn) {

        int earliestFeasibleSercivingFinishingTime = finServicingTime;
        if (!isServicePossible(finServicingTime, destinationOrderNumber)) {
            earliestFeasibleSercivingFinishingTime = getEarliestFeasibleSercivingFinishingTime(finServicingTime +1, destinationOrderNumber, nIterationsIn +1);
        }
        return earliestFeasibleSercivingFinishingTime;
    } //TODO fix bug: How do we deal with long periods of bad weather? How long waiting do we accept?

    private double[] servicingCalculations(int time, int order) {

        double realTime = convertNodeTimeToRealTime(time);
        double servicingTimeLeft = problemData.getDemandByOrderNumber(order)*timePerHiv;
        double consumption = 0;

        // System.out.println("dem: " + problemData.getDemandByOrderNumber(order));


        // System.out.println(" - " + realTime + " - " + consumption + " - " + servicingTimeLeft);

        while (servicingTimeLeft > 0) {

            if (realTime % 1 > 0) {

                double weather_impact = problemData.getWeatherImpactByHour((int) floor(realTime));
                // System.out.println("wi:  " + weather_impact);

                if (servicingTimeLeft < (realTime % 1)/ weather_impact) {
                    consumption += servicingConsumption * servicingTimeLeft * weather_impact;
                    realTime -= servicingTimeLeft * weather_impact;
                    servicingTimeLeft = -1;
                }
                else {
                    consumption += servicingConsumption * (realTime % 1) * weather_impact;
                    servicingTimeLeft -= (realTime % 1) / weather_impact;
                    realTime = floor(realTime);
                }
            }
            else {

                double weather_impact = problemData.getWeatherImpactByHour((int) floor(realTime -1));
                // System.out.println("wi:  " + weather_impact);

                if (servicingTimeLeft < 1 / weather_impact) {
                    consumption += servicingConsumption * servicingTimeLeft * weather_impact;
                    realTime -= servicingTimeLeft * weather_impact;
                    servicingTimeLeft = - 1;
                }
                else {
                    consumption += servicingConsumption * weather_impact;
                    servicingTimeLeft -= 1 / weather_impact;
                    realTime--;
                }
            }
            // System.out.println(" - " + realTime + " - " + consumption + " - "  + servicingTimeLeft);
        }
        return new double[] {consumption, realTime};
    }

    private boolean isArrivalPossible(double startTime, double distance, double finIdlingRealTime) {

        double[] tiws = getTimeInAllWS(startTime,finIdlingRealTime);

        double maxDistance = tiws[0]*maxSpeed + tiws[1]*maxSpeed + tiws[2]*maxSpeedWS2 + tiws[3]*maxSpeedWS3; //TODO fix get time in all ws and its dependencies

        return maxDistance >= distance;
    }


    private double[] idlingCalculations(double startTime, double distance, double finIdlingRealTime) {

        double longestSailingTime = distance/minSpeed;

        if (longestSailingTime >= finIdlingRealTime - startTime) {
            return new double[] {0 , finIdlingRealTime};
        }
        double idling_duration = finIdlingRealTime - longestSailingTime - startTime;
        double real_arrival_time = finIdlingRealTime - idling_duration;
        double[] time_in_weather_states = getTimeInAllWS(real_arrival_time, finIdlingRealTime);
        double consumption = 0;
        for (int i = 0; i < time_in_weather_states.length; i++) {
            consumption += time_in_weather_states[i] * problemData.getWeatherImpactByState().get(i) * idlingConsumption;
        }
        return new double[] {consumption, real_arrival_time};
    }


    private double[] getTimeInAllWS(double t1, double t2) {
        double tiws3 = getTimeInWS(t1, t2, 3);
        double tiws2 = getTimeInWS(t1, t2, 2);
        double tiws1 = getTimeInWS(t1, t2, 1);
        double tiw0 = t2 - t1 - tiws3 - tiws2 - tiws1;
        return new double[] {tiw0, tiws1, tiws2, tiws3};
    }

    private double getTimeInWS(double t1, double t2, int weatherState){
        double time = 0;
        double _t1 = t1;
        while (_t1 < t2) {
            if (_t1 % 1 > 0) {
                if (t2 - _t1 < 1 - (_t1 % 1)) {
                    time += isWeatherState(weatherState, (int) _t1) ? t2 - _t1 : 0;
                    _t1 = t2 + 1;
                }
                else {
                    time += isWeatherState(weatherState, (int) _t1) ? 1 - (_t1 % 1) : 0;
                    _t1 = ceil(_t1);
                }
            }
            else {
                if (t2 - _t1 <= 1) {
                    time += isWeatherState(weatherState, (int) _t1) ? t2 - _t1 : 0;
                    _t1 = t2 + 1;
                }
                else {
                    time += isWeatherState(weatherState, (int) _t1) ? 1 : 0;
                    _t1 += 1;
                }
            }
        }
        return time;
    }

    private boolean isWeatherState(int weatherState, int hour) {
        return problemData.getWeatherStateByHour().get(hour) == weatherState;
    }

    private double calculateAdjustedAverageSpeed(double[] timeInAllWeatherStates, double distance) {
        double durationWS3 = timeInAllWeatherStates[3];
        double durationWS2 = timeInAllWeatherStates[2];
        double durationWS01 = timeInAllWeatherStates[0] + timeInAllWeatherStates[1];
        double duration = durationWS01 + durationWS2 + durationWS3;

        double speed = distance/duration;

        if (durationWS01 == 0 && speed > (durationWS2*maxSpeedWS2 + durationWS3*maxSpeedWS3)/(durationWS2 + durationWS3)) {
            speed = maxSpeed+1;
        }

        else {
            if (speed > maxSpeedWS3) {
                speed += (durationWS3 * (speed - maxSpeedWS3)) / (durationWS01 + durationWS2);
            }

            if (speed > maxSpeedWS2) {
                speed += (durationWS2 * (speed - maxSpeedWS2)) / (durationWS01);
            }
        }

        return speed;
    }

    private double sailingCalculations(double[] timeInAllWeatherStates, double adjustedAverageSpeed) {
        double cost = 0;
        cost += (timeInAllWeatherStates[0] + timeInAllWeatherStates[1])*consumption_per_hour_at_speed(adjustedAverageSpeed);
        // // System.out.println(cost);
        cost += (adjustedAverageSpeed < maxSpeedWS2) ? timeInAllWeatherStates[2]*consumption_per_hour_at_speed(adjustedAverageSpeed + speedImpactWS2) : timeInAllWeatherStates[2]*consumption_per_hour_at_speed(maxSpeed);
        // // System.out.println(cost);
        cost += (adjustedAverageSpeed < maxSpeedWS3) ? timeInAllWeatherStates[3]*consumption_per_hour_at_speed(adjustedAverageSpeed + speedImpactWS3) : timeInAllWeatherStates[3]*consumption_per_hour_at_speed(maxSpeed);
        // // System.out.println(cost);
        return cost;
    }

    private double consumption_per_hour_at_speed(double speed) {
        return 11.111*speed*speed - 177.78*speed + 1011.1;
    }

    private boolean isNodeInGraph(int columnNumber, int finServicingTime) {
        return graph.get(columnNumber).containsKey(finServicingTime);
    }

    public int getSize() {
        return graph.size();
    }

    public HashMap<Integer, HashMap<Integer, Node>> getGraph() {
        return graph;
    }
}
