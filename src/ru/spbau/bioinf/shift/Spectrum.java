package ru.spbau.bioinf.shift;

import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Properties;

public class Spectrum {

    private int id;
    private List<Peak> peaks = new ArrayList<Peak>();
    private String scans;
    private double precursorMz;
    private int precursorCharge;
    private double precursorMass;

    private double[] data = null;
    private double[] additionalSpectrum = null;

    public Spectrum(Properties prop, BufferedReader input) throws IOException {
        id = ReaderUtil.getIntValue(prop, "ID");
        scans = ReaderUtil.getValue(prop, "SCANS");
        precursorMz = ReaderUtil.getDoubleValue(prop, "PRECURSOR_MZ");
        precursorCharge = ReaderUtil.getIntValue(prop, "PRECURSOR_CHARGE");
        precursorMass = ReaderUtil.getDoubleValue(prop, "PRECURSOR_MASS");
        List<String[]> datas = ReaderUtil.readDataUntil(input, "END IONS");
        for (String[] data : datas) {
            double mass = Double.parseDouble(data[0]);
            if (mass < precursorMass - 50) {
                peaks.add(new Peak(mass, Double.parseDouble(data[1]), Integer.parseInt(data[2])));
            }
        }
    }

    public Spectrum(List<Peak> peaks, double precursorMass) {
        this.precursorMass = precursorMass;
        this.peaks = peaks;
    }

    public int getId() {
        return id;
    }

    public double getPrecursorMass() {
        return precursorMass;
    }

    public List<Peak> getPeaks() {
        return peaks;
    }

    public int getPeakCharge(int peakId) {
        return peaks.get(peakId).getCharge();
    }

    public double[] getData() {
        if (data == null) {
            List<Double> d = new ArrayList<Double>();
            d.add(0d);
            d.add(precursorMass - Consts.WATER);
            data = generateMasses(d, Modifications.CORE);
        }
        return data;
    }

    private double[] generateMasses(List<Double> d, Modifications mod) {
        for (Peak peak : peaks) {
            double p = peak.getMonoisotopicMass();
            List<Double> values = mod.getModifications(p, precursorMass);
            for (Double value : values) {
                if (value > 0 && value < precursorMass) {
                    d.add(value);
                }
            }
        }
        Collections.sort(d);
        d = ProteinFinder.merge(d);
        double[] a = new double[d.size()];
        for (int i = 0; i < a.length; i++) {
            a[i] = d.get(i);
        }
        return a;
    }

    public double[] getAdditionalSpectrum() {
        if (additionalSpectrum == null) {
           List<Double> d = new ArrayList<Double>();
            additionalSpectrum = generateMasses(d, Modifications.ADDITIONAL);
        }
        return additionalSpectrum;
    }

    public Spectrum getLeft(double[] p, double shift) {
        List<Peak> left = new ArrayList<Peak>();
        for (Peak peak : peaks) {
            List<Double> mods = Modifications.CORE.getModifications(peak.getMonoisotopicMass(), precursorMass);
            boolean found = false;
            for (double mod : mods) {
                for (double v : p) {
                    if (Math.abs(mod - v + shift) < 0.1) {
                        found = true;
                        break;
                    }
                }
                if (found) {
                    break;
                }
            }
            if (!found) {
                left.add(peak);
            }
        }
        return new Spectrum(left, precursorMass);
    }

    public void clearData() {
        data = null;
        additionalSpectrum = null;
    }

    public Element toXml(Map<Integer, List<Break>> breaks) {
        Element spectrum = new Element("spectrum");
        XmlUtil.addElement(spectrum, "scans", scans);
        XmlUtil.addElement(spectrum, "spectrum-id", id);
        XmlUtil.addElement(spectrum, "precursor-mz", precursorMz);
        XmlUtil.addElement(spectrum, "precursor-charge", precursorCharge);
        XmlUtil.addElement(spectrum, "precursor-mass", precursorMass);
        Element peaksTag = new Element("peaks");
        spectrum.addContent(peaksTag);
        for (int i = 0; i < peaks.size(); i++) {
            Peak peak = peaks.get(i);
            Element peakTag = new Element("peak");
            peaksTag.addContent(peakTag);
            XmlUtil.addElement(peakTag, "peak-id", i);
            peak.addToXml(peakTag, breaks.get(i));
        }
        return spectrum;
    }
}
