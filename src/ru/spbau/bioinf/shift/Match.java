package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.util.ArrayList;

public class Match {
    private double score;
    private String protein;
    private int proteinId;
    private double shift;

    ArrayList<Double> shifts = new ArrayList<Double>();
    ArrayList<Double> shiftsScores = new ArrayList<Double>();

    public Match(Protein protein) {
        this.protein = protein.getAcids();
        this.proteinId = protein.getProteinId();
    }

    public Match(Protein protein, double score, double shift) {
        this.score = score;
        this.protein = protein.getAcids();
        this.proteinId = protein.getProteinId();
        this.shift = shift;
    }

    public void update(double score, double shift) {
        if (score > this.score) {
            this.score = score;
            this.shift = shift;
        }
    }

    public void addShift(double v, double score) {
        shifts.add(v);
        shiftsScores.add(score);
    }
    public int getProteinId() {
        return proteinId;
    }


    public double getScore() {
        return score;
    }

    public double getShift() {
        return shift;
    }

    public Element toXml() {
        Element match = new Element("match");
        XmlUtil.addElement(match, "protein", protein);
        XmlUtil.addElement(match, "protein-id", proteinId);
        Element sh = new Element("shifts");
        match.addContent(sh);
        for (int i = 0; i < shifts.size(); i++) {
            Element shift = new Element("shift");
            XmlUtil.addElement(shift, "value", shifts.get(i));
            XmlUtil.addElement(shift, "score", shiftsScores.get(i));
            sh.addContent(shift);
        }
        return match;
    }
}
