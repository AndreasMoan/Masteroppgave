package main.java.HGSADCwSO.implementations.DAG;

import java.util.ArrayList;

public class Edge {

    private Node parentNode;
    private Node childNode;

    private double cost;

    private double time_start;
    private double time_arrival;
    private double time_service;
    private double time_end;

    private double aveage_speed;

    public Edge(Node parentNode, Node childNode, double cost, double time_start, double time_arrival, double time_service, double time_end, double distance){
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.cost = cost;
        this.time_start = time_start;
        this.time_arrival = time_arrival;
        this.time_service = time_service;
        this.time_end = time_end;
        this.aveage_speed = distance/(time_arrival - time_start);
    }

    public double getCost() {
        return cost;
    }

    public Node getChildNode() {
        return childNode;
    }

    public Node getParentNode() {
        return parentNode;
    }

    public double getAveage_speed() {
        return aveage_speed;
    }

    public double getTime_arrival() {
        return time_arrival;
    }

    public double getTime_end() {
        return time_end;
    }

    public double getTime_service() {
        return time_service;
    }

    public double getTime_start() {
        return time_start;
    }
}
