package main.java.HGSADCwSO.files;

import main.java.HGSADCwSO.implementations.DAG.Edge;
import main.java.HGSADCwSO.implementations.DAG.Node;
import main.java.HGSADCwSO.protocols.FitnessEvaluationProtocol;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Voyage {


    public static void print_schedule(Individual individual, FitnessEvaluationProtocol fitnessEvaluationProtocol, ProblemData problemData) {
        HashMap<Integer, Node> solution_nodes = fitnessEvaluationProtocol.getSolutionNodes(individual);
        System.out.println("============ Printing schedule: ============");
        System.out.println(individual.getVesselTourChromosome());
        for (Integer vessel_number : solution_nodes.keySet()) {
            printVoyage(vessel_number, solution_nodes.get(vessel_number), individual, problemData);
        }
    }

    public static void printVoyage(int vessel_number, Node best_solution_node, Individual individual, ProblemData problemData) {
        System.out.println();
        System.out.println("------------ vessel number " + vessel_number + " -----------");
        ArrayList<Edge> solution_edges = best_solution_node.getAllBestParentEdges();
        Collections.sort(solution_edges, Utilities.getEdgeTimeComparator());
        for (int i = 0; i < solution_edges.size(); i++) {
            Edge edge = solution_edges.get(i);
            int installation_number = 0;
            if (i < solution_edges.size() -1) {
                installation_number = problemData.getInstallationNumberByOrderNumber(individual.getVesselTourChromosome().get(vessel_number).get(i));
            }
            System.out.println("Leg nr. " + i + " to installation " + installation_number
                    + ", depatures: " + convertTime(edge.getTime_start())
                    + ", arrives: " + convertTime(edge.getTime_arrival())
                    + ", " + "servicing: " + convertTime(edge.getTime_service()));
        }
    }

    public static String convertTime(double time) {
        int day = (int) Math.floor(time / 24);
        int hour = (int) Math.floor(time % 24);
        int minute = (int) Math.floor((time % 1)*60);
        int second = (int) Math.round(((time*60) % 1) * 60);
        String string = day + "/" + timeToString(hour) + ":" + timeToString(minute);
        return string;
    }

    public static String timeToString(int time) {
        String string = "";
        if (time < 10) {
            string += "0" + time;
        }
        else string += time;
        return string;
    }
}
