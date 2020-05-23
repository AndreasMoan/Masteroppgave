package main.java.HGSADCwSO.files;

import java.util.ArrayList;
import java.util.HashMap;

public class Schedule {

    private HashMap<Integer, ArrayList<SailingLeg>> schedule;

    private double cost;

    public Schedule(HashMap<Integer, ArrayList<SailingLeg>> schedule, double cost) {
        this.schedule = schedule;
        this.cost = cost;
    }

    public double getCost() {
        return cost;
    }

    public HashMap<Integer, ArrayList<SailingLeg>> getSchedule() {
        return schedule;
    }

    public void setCost(double cost) {
        this.cost = cost;
    }

    public void setSchedule(HashMap<Integer, ArrayList<SailingLeg>> schedule) {
        this.schedule = schedule;
    }
}
