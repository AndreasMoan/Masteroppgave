package main.java.HGSADCwSO.implementations.DAG;

import java.util.ArrayList;

public class Node {

    private int time;
    private int orderNumber;
    private ArrayList<Edge> parentEdges;
    private ArrayList<Edge> childEdges;

    private double bestCost;
    private double bestPenalizedCost;
    private Edge bestParentEdge;

    ArrayList<Edge> allBestParentEdges;

    private double bestTotalDeadlineViolation;
    private double deadlineViolation;

    public Node(int time, int orderNumber, double deadlineViolation){
        this.time = time;
        this.orderNumber = orderNumber;
        this.deadlineViolation = deadlineViolation;
        this.childEdges = new ArrayList<Edge>();
        this.parentEdges = new ArrayList<Edge>();
        this.allBestParentEdges = new ArrayList<>();
        this.setBestCost(Double.POSITIVE_INFINITY);
        this.setBestPenalizedCost(Double.POSITIVE_INFINITY);
    }

    public void addParentEdge(Edge edge){
        parentEdges.add(edge);
    }

    public void addChildEdge(Edge edge){
        childEdges.add(edge);
    }

    public ArrayList<Edge> getParentEdges() {
        return parentEdges;
    }

    public ArrayList<Edge> getChildEdges() {
        return childEdges;
    }

    public int getTime() {
        return time;
    }

    public int getOrderNumber() {
        return orderNumber;
    }

    public double getBestCost() {
        return bestCost;
    }

    public void setBestCost(double bestCost) {
        this.bestCost = bestCost;
    }

    public Edge getBestParentEdge() {
        return bestParentEdge;
    }

    public void setBestParentEdge(Edge bestParentEdge) {
        this.bestParentEdge = bestParentEdge;
        if (this.allBestParentEdges.size() > 0) {
            this.allBestParentEdges.clear();
        }
        for (Edge edge: bestParentEdge.getParentNode().getAllBestParentEdges()) {
            this.allBestParentEdges.add(edge);
        }
        this.allBestParentEdges.add(bestParentEdge);
    }

    public double getDeadlineViolation() {
        return deadlineViolation;
    }

    public void setDeadlineViolation(int deadlineViolation) {
        this.deadlineViolation = deadlineViolation;
    }

    public double getBestTotalDeadlineViolation() {
        return bestTotalDeadlineViolation;
    }

    public void setBestTotalDeadlineViolation(double totalDeadlineViolation) {
        this.bestTotalDeadlineViolation = bestTotalDeadlineViolation;
    }

    public void setBestPenalizedCost(double bestPenalizedCost) {
        this.bestPenalizedCost = bestPenalizedCost;
    }

    public double getBestPenalizedCost() {
        return bestPenalizedCost;
    }

    public ArrayList<Edge> getAllBestParentEdges() {
        return allBestParentEdges;
    }
}
