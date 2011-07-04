package ru.spbau.bioinf.shift;

import java.util.ArrayList;

public class Protein {
    private int proteinId;
    private String acids;
    private String simplifiedAcids = null;
    private double[] spectrum = null;

    ArrayList<Double> shifts = new ArrayList<Double>();
    ArrayList<Integer> shiftScores = new ArrayList<Integer>();

    public Protein(int proteinId, String acids) {
        this.proteinId = proteinId;
        this.acids = acids;
    }

    public int getProteinId() {
        return proteinId;
    }

    public String getAcids() {
        return acids;
    }

    public String getSimplifiedAcids() {
        if (simplifiedAcids == null) {
            simplifiedAcids = acids.replaceAll("L", "I").replaceAll("Z", "Q").replaceAll("B", "E").replaceAll("X", "I");
        }
        return simplifiedAcids;
    }

    public double[] getSpectrum() {
        String sa = getSimplifiedAcids();
        if (spectrum == null) {
            spectrum = new double[sa.length()  +1];
            for (int i = 0; i < sa.length(); i++) {
                  spectrum[i+1] = spectrum[i] + Acids.acids.get(sa.charAt(i));
            }
        }
        return spectrum;
    }
}
