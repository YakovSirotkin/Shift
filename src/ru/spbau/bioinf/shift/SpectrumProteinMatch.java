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

    public SpectrumProteinMatch(Spectrum spectrum, Protein protein) {
        this.spectrum = spectrum;
        this.protein = protein;
        double[] sd = spectrum.getData();
        double[] pd = protein.getSpectrum();
        Map<String,List<Double>> positions = ProteinFinder.getPositions(spectrum);
        List<Double> shiftsList = ProteinFinder.getShifts(protein, positions);
        for (double shift : shiftsList) {
            double  score = ProteinFinder.getScore(sd, pd, shift);
            shifts.put(shift, score);
            if (score > bestScore) {
                bestScore = score;
                bestShift = shift;
            }
        }
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
        if (bestScore > 0) {
            firstResidue = i;
            nMass = sd[j];
        }
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
