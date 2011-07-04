package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.util.ArrayList;

public class Match {
    private int score;
    private String protein;
    private int proteinId;

    ArrayList<Double> shifts = new ArrayList<Double>();
    ArrayList<Integer> shiftsScores = new ArrayList<Integer>();

    public Match(Protein protein) {
        this.protein = protein.getAcids();
        this.proteinId = protein.getProteinId();
    }

    public Match(Protein protein, int score) {
        this.score = score;
        this.protein = protein.getAcids();
        this.proteinId = protein.getProteinId();
    }

    public void addShift(double v, int score) {
        shifts.add(v);
        shiftsScores.add(score);
    }
    public int getProteinId() {
        return proteinId;
    }


    public int getScore() {
        return score;
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
