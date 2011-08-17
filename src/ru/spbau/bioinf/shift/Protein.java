package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.util.ArrayList;
import java.util.List;

public class Protein {
    private int proteinId;
    private String acids;
    private String simplifiedAcids = null;
    private double[] spectrum = null;
    private List<SpectrumProteinMatch> matches = new ArrayList<SpectrumProteinMatch>();
    private String name;

    public Protein(int proteinId, String acids, String name) {
        this.proteinId = proteinId;
        this.acids = acids;
        this.name = name;
    }

    public int getProteinId() {
        return proteinId;
    }

    public String getAcids() {
        return acids;
    }

    public String getName() {
        return name;
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

    public void addSpectrumMatch(SpectrumProteinMatch match) {
          matches.add(match);
    }

    public Element toXml() {
        Element p = new Element("protein");
        XmlUtil.addElement(p, "acids", acids);
        XmlUtil.addElement(p, "protein-id", proteinId);
        XmlUtil.addElement(p, "name", name);
        Element spectrums = new Element("matches");
        for (SpectrumProteinMatch match : matches) {
            spectrums.addContent(match.toXmlWithoutProtein());
        }
        p.addContent(spectrums);
        return p;
    }
}
