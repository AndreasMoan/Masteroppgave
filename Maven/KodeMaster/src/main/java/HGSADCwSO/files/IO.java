package main.java.HGSADCwSO.files;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class IO {

    int scenario_number;
    FileWriter file_writer;
    File file;

    public IO(int scenario_number) throws IOException {
        this.scenario_number = scenario_number;
        String filename = "Run_" + scenario_number + ".txt";
        try {
            file = new File(filename);
            if (file.createNewFile()) {
                System.out.println("File created: " + file.getName());
            } else {
                System.out.println("File already exists.");
            }
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
        file_writer = new FileWriter(filename);
    }

    private double convert_time(double nano_time) {
        return nano_time/1000000000.0;
    }

    public void writeRunInfo(ProblemData problemData) throws IOException {
        int min_pop_sz = problemData.getHeuristicParameterInt("Population size");
        int gen_sz = problemData.getHeuristicParameterInt("Number of offspring in a generation");
        double education_rate = problemData.getHeuristicParameterDouble("Education rate");
        double repair_rate = problemData.getHeuristicParameterDouble("Repair rate");
        double move_chance = problemData.getHeuristicParameterDouble("Move chance");
        String line = "\n=================== Changing parameters: ===================" +
                "\nMininum subpopulation size: " + min_pop_sz +
                "\nGeneration size:            " + gen_sz +
                "\nEducation rate:             " + education_rate +
                "\nRepair rate:                " + repair_rate +
                "\nMove chance:                " + move_chance +
                "\n\nIteration - time spent - final objective value\n";
        file_writer.write(line);
        file_writer.flush();
    }

    public void writeImprovementIteration(int iteration, double nano_time, double objective) throws IOException {
        double time = convert_time(nano_time);
        String line = iteration + " - " + time + " - " + objective + "\n";
        file_writer.write(line);
        file_writer.flush();
    }

    public void writeFinalSolution(int iteration, double nano_time, Individual individual) throws IOException {
        double time = convert_time(nano_time);
        HashMap<Integer, ArrayList<Integer>> chromosome = individual.getVesselTourChromosome();
        String line = "\nFinal solution:" +
                "\nIteration:  " + iteration +
                "\nTime:       " + time +
                "\nCost:       " + individual.getPenalizedCost() +
                "\nChromosome: " + individual.getVesselTourChromosome() + "\n\n";
        file_writer.write(line);
        file_writer.flush();
    }
}
