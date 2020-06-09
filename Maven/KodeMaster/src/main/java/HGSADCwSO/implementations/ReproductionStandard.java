package main.java.HGSADCwSO.implementations;

import main.java.HGSADCwSO.files.Individual;
import main.java.HGSADCwSO.files.ProblemData;
import main.java.HGSADCwSO.files.Utilities;
import main.java.HGSADCwSO.protocols.FitnessEvaluationProtocol;
import main.java.HGSADCwSO.protocols.ReproductionProtocol;

import java.util.*;
import java.util.concurrent.*;

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

        for (int orderNumber : allOrders) {

            ExecutorService executor = Executors.newFixedThreadPool(40); //TODO 40
            List<Callable<double[]>> tasksList = new ArrayList<>();

            for (int vesselNumber = 0; vesselNumber < kid.size(); vesselNumber++) {

                for (int insertionPosition = 0; insertionPosition <= kid.get(vesselNumber).size(); insertionPosition++){


                    int vn = vesselNumber;
                    int ip = insertionPosition;

                    Callable<double[]> task = () -> {

                        HashMap<Integer, ArrayList<Integer>> copy_of_kid = Utilities.deepCopyVesselTour(kid);
                        copy_of_kid.get(vn).add(ip, orderNumber);

                        Individual candidate_offspring = new Individual(copy_of_kid, fitnessEvaluationProtocol);
                        fitnessEvaluationProtocol.evaluate(candidate_offspring);

                        double candidate_offspring_penalized_cost = candidate_offspring.getPenalizedCost();
                        return new double[] {candidate_offspring_penalized_cost, vn, ip};
                    };
                    tasksList.add(task);
                }
            }
            int best_vn = -1;
            int best_ip = -1;
            double best_copc = Double.POSITIVE_INFINITY;

            try {
                List<Future<double[]>> results = executor.invokeAll(tasksList);

                for (Future<double[]> result : results) {
                    if (result.get()[0] < best_copc) {
                        best_copc = result.get()[0];
                        best_vn = (int) result.get()[1];
                        best_ip = (int) result.get()[2];
                    }
                }
            }
            catch (InterruptedException | ExecutionException e1) {
                e1.printStackTrace();
            }
            executor.shutdownNow();

            kid.get(best_vn).add(best_ip, orderNumber);
        }
        return new Individual(kid, fitnessEvaluationProtocol);
    }



    @Override
    public int getNumberOfCrossoverRestarts(){
        return NumberOfCrossoverRestarts;
    }
}
