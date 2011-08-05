package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SpectrumProteinMatch {
    
    private Spectrum spectrum;
    private Protein protein;
    private double bestScore = 0;
    private double bestShift = 0;
    private Map<Double, Double> shifts = new HashMap<Double, Double>();
    private int firstResidue;
    private double nMass;

    public SpectrumProteinMatch(Spectrum spectrum, Protein protein, ScoringFunction scoringFunction) {
        this.spectrum = spectrum;
        this.protein = protein;
        double[] sd = spectrum.getData();
        double[] pd = protein.getSpectrum();
        Map<String,List<Double>> positions = ProteinFinder.getPositions(spectrum);
        List<Double> shiftsList = ProteinFinder.getShifts(protein, positions);
        firstResidue = 10000;
        for (double shift : shiftsList) {
            double  score = scoringFunction.getScore(spectrum, protein, shift);
            shifts.put(shift, score);
            if (score > bestScore) {
                bestScore = score;
                bestShift = shift;
            }
        }
        if (bestScore > 0) {
            findFirstResidue(sd, pd);
            findFirstResidue(spectrum.getAdditionalSpectrum(), pd);
        }
    }

    private void findFirstResidue(double[] sd, double[] pd) {
        int i = 0;
        int j = 0;
        do {
            double diff = sd[j] - pd[i] + bestShift;
            if (diff < -0.1) {
                j++;
            } else if (diff > 0.1) {
                i++;
            } else {
                break;
            }
        } while (i < pd.length && j < sd.length);

        if (firstResidue > i) {
            firstResidue = i;
            nMass = sd[j];
        }
    }

    public double getBestShift() {
        return bestShift;
    }

    public Element toXml(){
        Element match = new Element("match");
        match.addContent(spectrum.toXml(new HashMap<Integer, List<Break>>()));
        XmlUtil.addElement(match, "first-residue", firstResidue);
        XmlUtil.addElement(match, "n-mass", nMass);
        XmlUtil.addElement(match, "best-score", bestScore);
        XmlUtil.addElement(match, "best-shift", bestShift);
        return match;
    }
}
