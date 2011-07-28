package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.util.List;

public class Peak {
    private double monoisotopicMass;
    private double intensity;
    private int charge;

    public Peak(double monoisotopicMass, double intensity, int charge) {
        this.monoisotopicMass = monoisotopicMass;
        this.intensity = intensity;
        this.charge = charge;
    }

    public double getMonoisotopicMass() {
        return monoisotopicMass;
    }

    public double getIntensity() {
        return intensity;
    }

    public int getCharge() {
        return charge;
    }

    public void addToXml(Element peak, List<Break> breaks) {
        XmlUtil.addElement(peak, "monoisotopic-mass", getMonoisotopicMass());
        XmlUtil.addElement(peak, "intensity", getIntensity());
        XmlUtil.addElement(peak, "charge", getCharge());
        Element matchesTag = new Element("breaks");
        peak.addContent(matchesTag);
        if (breaks != null) {
            for (Break aBreak : breaks) {
                aBreak.addToXml(matchesTag);
            }
        }
    }
}
