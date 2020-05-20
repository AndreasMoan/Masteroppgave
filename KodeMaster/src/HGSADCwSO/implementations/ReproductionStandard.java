package HGSADCwSO.implementations;

import HGSADCwSO.Individual;
import HGSADCwSO.ProblemData;
import HGSADCwSO.Utilities;
import HGSADCwSO.protocols.FitnessEvaluationProtocol;
import HGSADCwSO.protocols.ReproductionProtocol;

import java.lang.reflect.Array;
import java.util.*;

public class ReproductionStandard implements ReproductionProtocol {

    private ProblemData problemData;
    private static int NumberOfCrossoverRestarts = 0;
    private FitnessEvaluationProtocol fitnessEvaluationProtocol;

    public ReproductionStandard(ProblemData problemData, FitnessEvaluationProtocol fitnessEvaluationProtocol) {
        this.problemData = problemData;
        this.fitnessEvaluationProtocol = fitnessEvaluationProtocol;
    }


    public Individual crossover(ArrayList<Individual> parents) {

        //lager denne én mer enn
        Set<Integer> allOrders = new HashSet<Integer>(problemData.getOrdersByNumber().keySet());

        // STEP 0: Inheritance rule
        int nVessels = problemData.getVessels().size();
        int nInstallations = problemData.getInstallations().size();

        int randomNumber1 = new Random().nextInt(nVessels);
        int randomNumber2 = new Random().nextInt(nVessels);

        //n1 must be smaller than n2
        int n1 = Math.min(randomNumber1, randomNumber2);
        int n2 = Math.max(randomNumber1, randomNumber2);

        ArrayList<Integer> vessels = new ArrayList<Integer>();
        for (int i = 0; i < nVessels; i++) {
            vessels.add(i);
        }
        Collections.shuffle(vessels);

        Set<Integer> vesselsFromDad = new HashSet<Integer>();
        Set<Integer> vesselsFromMom = new HashSet<Integer>();
        Set<Integer> vesselsFromBoth = new HashSet<Integer>();

        for (int i = 0; i < nVessels; i++) {
            if (i < n1) {
                vesselsFromDad.add(vessels.get(i));
            }
            else if (i < n2) {
                vesselsFromMom.add(vessels.get(i));
            }
            else {
                vesselsFromBoth.add(vessels.get(i));
            }
        }

        HashMap<Integer, ArrayList<Integer>> father = Utilities.deepCopyVesselTour(parents.get(0).getVesselTourChromosome());
        HashMap<Integer, ArrayList<Integer>> mother = Utilities.deepCopyVesselTour(parents.get(1).getVesselTourChromosome());

        HashMap<Integer, ArrayList<Integer>> kid = new HashMap<Integer, ArrayList<Integer>>();
        for (int i = 0; i < problemData.getNumberOfVessels(); i++) {
            kid.put(i, new ArrayList<>());
        }

        // STEP 1: inherit from father

        for (Integer vesselNumber : vesselsFromDad) {
            ArrayList<Integer> ordersToCopy = father.get(vesselNumber);
            kid.put(vesselNumber, ordersToCopy);
            for (int i : father.get(vesselNumber)){
                allOrders.remove(i);
            }
        }

        for (Integer vesselNumber : vesselsFromBoth) {

            ArrayList<Integer> ordersToCopy = new ArrayList<>();
            ArrayList<Integer> fatherOrders = father.get(vesselNumber);

            if (fatherOrders.size() > 0) {

                int alpha = new Random().nextInt(fatherOrders.size());
                int beta = new Random().nextInt(fatherOrders.size());

                if (alpha <= beta) {
                    ordersToCopy = new  ArrayList<>(fatherOrders.subList(alpha, beta));
                    for (int i : fatherOrders.subList(alpha, beta)) {
                        allOrders.remove(i);
                    }
                }
                else {
                    for (int i = 0; i < beta; i++) {
                        ordersToCopy.add(fatherOrders.get(i));
                        allOrders.remove(fatherOrders.get(i));
                    }
                    for (int i = alpha; i < fatherOrders.size(); i++) {
                        ordersToCopy.add(fatherOrders.get(i));
                        allOrders.remove(fatherOrders.get(i));
                    }
                }

            }
            kid.put(vesselNumber, ordersToCopy);
        }

        // STEP 2: inherit from mother

        Set<Integer> vesselsFromMomOrBoth = new HashSet<>(vesselsFromMom);
        vesselsFromMomOrBoth.addAll(vesselsFromBoth);

        while (!vesselsFromMomOrBoth.isEmpty()) {
            int vesselNumber = Utilities.pickAndRemoveRandomElementFromSet(vesselsFromMomOrBoth);

            for (int orderNumber : mother.get(vesselNumber)) {
                for (int i = 0; i < mother.size(); i++) {
                    if (allOrders.contains(orderNumber)) {
                        kid.get(vesselNumber).add(orderNumber);
                        allOrders.remove(orderNumber);
                    }
                }
            }
        }

        //STEP 3:

        Individual offspring = new Individual(kid, fitnessEvaluationProtocol);
        fitnessEvaluationProtocol.evaluate(offspring);
        double unfinished_penalized_cost = offspring.getPenalizedCost();
        double candidate_offspring_penalized_cost = 0;

        for (int orderNumber : allOrders) {

            double least_added_cost =  Double.POSITIVE_INFINITY;
            int bestInsertionVessel = 0;
            int bestInsertionPosition = 0;

            for (int vesselNumber = 0; vesselNumber < kid.size(); vesselNumber++) {

                for (int insertionPosition = 0; insertionPosition <= kid.get(vesselNumber).size(); insertionPosition++){

                    HashMap<Integer, ArrayList<Integer>> copy_of_kid = Utilities.deepCopyVesselTour(kid);
                    copy_of_kid.get(vesselNumber).add(insertionPosition, orderNumber);
                    Individual candidate_offspring = new Individual(copy_of_kid, fitnessEvaluationProtocol);
                    fitnessEvaluationProtocol.evaluate(candidate_offspring);
                    candidate_offspring_penalized_cost = candidate_offspring.getPenalizedCost();

                    double added_cost = Math.max(0, candidate_offspring_penalized_cost - unfinished_penalized_cost);

                    if (added_cost < least_added_cost) {
                        least_added_cost = added_cost;
                        bestInsertionPosition = insertionPosition;
                        bestInsertionVessel = vesselNumber;
                    }

                }
            }
            kid.get(bestInsertionVessel).add(bestInsertionPosition,orderNumber);
            unfinished_penalized_cost = candidate_offspring_penalized_cost;
        }
        return new Individual(kid, fitnessEvaluationProtocol);
    }



    @Override
    public int getNumberOfCrossoverRestarts(){
        return NumberOfCrossoverRestarts;
    }
}
