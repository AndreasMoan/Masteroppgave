package main.java.HGSADCwSO.files;

import main.java.HGSADCwSO.implementations.DAG.Edge;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.*;

public class Utilities {

    public static <T> T pickRandomElementFromList(List<T> list) {
        int randomIndex = new Random().nextInt(list.size());
        return list.get(randomIndex);
    }

    public static <T> T pickAndRemoveRandomElementFromList(ArrayList<T> list) {
        T element = pickRandomElementFromList(list);
        list.remove(element);
        return element;
    }

    public static <T> T pickRandomElementFromSet(Set<T> set) {
        int randomIndex = new Random().nextInt(set.size());
        int counter = 0;
        for (T t : set) {
            if (counter == randomIndex){
                return t;
            }
            counter++;
        }
        return null;
    }

    public static <T> T pickAndRemoveRandomElementFromSet(Set<T> set) {
        int randomIndex = new Random().nextInt(set.size());
        int counter = 0;
        for (T t : set) {
            if (counter == randomIndex){
                set.remove(t);
                return t;
            }
            counter++;
        }
        return null;
    }

    public static double parseDouble(String commaSeparatedDouble) {
        DecimalFormat df = new DecimalFormat();
        DecimalFormatSymbols symbols = new DecimalFormatSymbols();
        symbols.setDecimalSeparator(',');
        df.setDecimalFormatSymbols(symbols);
        try {
            return df.parse(commaSeparatedDouble).doubleValue();
        } catch (ParseException e) {
            e.printStackTrace();
            System.out.println("Something went wrong when parsing doubles");
            return -1.0;
        }
    }

    public static <T> ArrayList<T> getAllElements(ArrayList<T> list1, ArrayList<T> list2) {
        ArrayList<T> allElements = new ArrayList<T>();
        allElements.addAll(list1);
        allElements.addAll(list2);
        return allElements;
    }

    /*public static Comparator<Individual> getFitnessComparator() { //sorts in descending order
        return new Comparator<Individual>() {
            @Override
            public int compare(Individual individual1, Individual individual2) {
                if (individual1.getPenalizedCost() < individual2.getPenalizedCost()){
                    return  1;
                }
                else if (individual1.getPenalizedCost() > individual2.getPenalizedCost()) {
                    return  -1;
                }
                else {
                    return 0;
                }
            }
        };
    }*/



    public static <K> Comparator <Map.Entry<K, Double>> getMapEntryWithDoubleComparator() {
        return new Comparator<Map.Entry<K, Double>>() {
            public int compare(Map.Entry<K, Double> o1, Map.Entry<K, Double> o2) {
                return (o1.getValue()).compareTo(o2.getValue());
            }
        };
    }

    public static <T> Set<Set<T>> cartesianProduct(Set<T> set) { //returns cartesian product of itself, i.e. set x set, excluding those pairs where both elements in a pair are equal
        HashSet<Set<T>> cartProduct = new HashSet<>();
        for (T element : set) {
            for (T element2 : set){
                if (element != element2){
                    HashSet<T> pair = new HashSet<T>();
                    pair.add(element);
                    pair.add(element2);
                    cartProduct.add(pair);
                }
            }
        }
        return cartProduct;
    }

    public static HashMap<Integer, ArrayList<Integer>> deepCopyVesselTour(HashMap<Integer, ArrayList<Integer>> vesselTour) {

        HashMap<Integer, ArrayList<Integer>> vesselTourCopy = new HashMap<>();

        for (Integer vessel : vesselTour.keySet()){
            ArrayList<Integer> visitSequence= new ArrayList<Integer>(vesselTour.get(vessel));
            vesselTourCopy.put(vessel, visitSequence);
        }
        return vesselTourCopy;
    }

    public static Comparator<Order> getDeadlineComparator() { //testet! Den gir tidligste deadline sist OBS!!! SNUDDD
        return new Comparator<Order>() {
            @Override
            public int compare(Order order1, Order order2) {
                if (order1.getDeadline() > order2.getDeadline()){
                    return  -1;
                }
                else if (order1.getDeadline() < order2.getDeadline()) {
                    return  1;
                }
                else {
                    return 0;
                }
            }
        };
    }

    public static Comparator<Individual> getBiasedFitnessComparator() { //sorts in ascending order
        return new Comparator<Individual>() {
            public int compare(Individual ind1, Individual ind2) {
                if (ind1.getBiasedFitness() < ind2.getBiasedFitness()) {
                    return -1;
                }
                else if (ind1.getBiasedFitness() > ind2.getBiasedFitness()) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        };
    }

    public static Comparator<Individual> getPenalizedCostComparator() {
        return new Comparator<Individual>() {
            public int compare(Individual ind1, Individual ind2) {
                if (ind1.getPenalizedCost() < ind2.getPenalizedCost()) {
                    return -1;
                }
                else if (ind1.getPenalizedCost() > ind2.getPenalizedCost()) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        };
    }

    public static Comparator<Individual> getDiversityContributionComparator() { //sorts in descending order
        return new Comparator<Individual>() {
            public int compare(Individual ind1, Individual ind2) {
                if (ind1.getDiversityContribution() < ind2.getDiversityContribution()) {
                    return 1;
                }
                else if (ind1.getDiversityContribution() > ind2.getDiversityContribution()) {
                    return -1;
                }
                return 0;
            }
        };
    }

    public static Comparator<Edge> getEdgeTimeComparator() { //sorts in descending order
        return new Comparator<Edge>() {
            public int compare(Edge e1, Edge e2) {
                if (e1.getTime_start() < e2.getTime_start()) {
                    return -1;
                }
                else if (e1.getTime_start() > e2.getTime_start()) {
                    return 1;
                }
                return 0;
            }
        };
    }

    public static ArrayList<Double> normalize(ArrayList<Double> array) {

        double sum = 0;

        for (double value : array) {
            sum += value;
        }

        ArrayList<Double> normalized_array = new ArrayList<>();

        for (int i = 0; i < array.size(); i++) {
            double prob = array.get(i)/sum;
            normalized_array.add(prob);
        }

        return normalized_array;
    }

    public static int getRandomElementFromDistribution(ArrayList<Double> probabilities) {

        double random_var = new Random().nextDouble();

        for (int i = 0; i < probabilities.size(); i++) {
            random_var -= probabilities.get(i);
            if (random_var <= 0) {
                return i;
            }
        }

        return probabilities.size() -1;
    }

    public static int getBestElementFromDistribution(ArrayList<Double> probabilities) {
        double best = 0;
        int index = 0;
        for (int i = 0; i < probabilities.size(); i++) {
            if (probabilities.get(i) > best) {
                best = probabilities.get(i);
                index = i;
            }
        }

        return index;
    }


}


