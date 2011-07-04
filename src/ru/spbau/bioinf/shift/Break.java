package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.XmlUtil;

public class Break {
    private char type;

    private int peakId;
    private int charge;
    private int matchId;

    private float recalibrationMass;
    private float fragmentMass;
    private float massShift;
    private int ionPosition;
    private int residueId;

    public Break(char type) {
        this.type = type;
    }

    public char getType() {
        return type;
    }

    public int getPeakId() {
        return peakId;
    }

    public void setPeakId(int peakId) {
        this.peakId = peakId;
    }

    public int getCharge() {
        return charge;
    }

    public void setCharge(int charge) {
        this.charge = charge;
    }

    public int getMatchId() {
        return matchId;
    }

    public void setMatchId(int matchId) {
        this.matchId = matchId;
    }

    public boolean isReverse() {
        switch (type) {
            case 'Y' : return true;
            case 'Z' : return true;
        }
        return false;
    }

    public float getRecalibrationMass() {
        return recalibrationMass;
    }

    public void setRecalibrationMass(float recalibrationMass) {
        this.recalibrationMass = recalibrationMass;
    }

    public float getFragmentMass() {
        return fragmentMass;
    }

    public void setFragmentMass(float fragmentMass) {
        this.fragmentMass = fragmentMass;
    }

    public float getMassError() {
        return recalibrationMass - fragmentMass;
    }

    public float getPPM() {
        return getMassError() * 1000000 / recalibrationMass;
    }

    public float getMassShift() {
        return massShift;
    }

    public void setMassShift(float massShift) {
        this.massShift = massShift;
    }

    public int getIonPosition() {
        return ionPosition;
    }

    public void setIonPosition(int ionPosition) {
        this.ionPosition = ionPosition;
    }

    public int getResidueId() {
        return residueId;
    }

    public void setResidueId(int residueId) {
        this.residueId = residueId;
    }

    public void addToXml(Element parent) {
        Element match = new Element("break");
        parent.addContent(match);
        XmlUtil.addElement(match, "type", "" + getType());
        XmlUtil.addElement(match, "mass-shift", getMassShift());
        XmlUtil.addElement(match, "recalibration-mass", getRecalibrationMass());
        XmlUtil.addElement(match, "match-id", getMatchId());
        XmlUtil.addElement(match, "ion-position", getIonPosition());
        XmlUtil.addElement(match, "residue-id", getResidueId());
        XmlUtil.addElement(match, "mass-error", getMassError());
        XmlUtil.addElement(match, "ppm", getPPM());
    }
}
