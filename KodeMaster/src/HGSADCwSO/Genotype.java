package HGSADCwSO;

import com.sun.media.jfxmediaimpl.platform.gstreamer.GSTPlatform;

import java.util.ArrayList;
import java.util.HashMap;

public class Genotype {


    private HashMap<Integer, ArrayList<Integer>> vesselTourChromosome;

    public Genotype(HashMap<Integer, ArrayList<Integer>> vesselTourChromosome) {

        this.vesselTourChromosome = vesselTourChromosome;

        System.out.println(vesselTourChromosome);

        //System.out.println("1: " + vesselTourChromosome);

        if(vesselTourChromosome.get(0).size() + vesselTourChromosome.get(1).size() > 9){
            System.out.println("Found error 1");
            System.exit(0);
        }

    }

    public void setVesselTourChromosome(HashMap<Integer, ArrayList<Integer>> vesselTourChromosome) {
        this.vesselTourChromosome = vesselTourChromosome;

        System.out.println(vesselTourChromosome);

        if(vesselTourChromosome.get(0).size() + vesselTourChromosome.get(1).size() > 9){
            System.out.println("Found error 2");
            System.exit(0);
        }

        //System.out.println("2: " + vesselTourChromosome);
    }

    public HashMap<Integer, ArrayList<Integer>> getVesselTourChromosome() {
        //System.out.println("3" + vesselTourChromosome);
        return vesselTourChromosome;

    }
}
